#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "quac.h"
#include "operators.h"
#include "error_correction.h"
#include "solver.h"
#include "dm_utilities.h"
#include "quantum_gates.h"
#include "quantum_circuits.h"
#include "petsc.h"
#include "qasm_parser.h"
#include "qsystem.h"


PetscScalar omega(PetscReal time,void *ctx);
PetscScalar delta(PetscReal time,void *ctx);

PetscErrorCode ts_monitor(TS,PetscInt,PetscReal,Vec,void*);
PetscErrorCode ts_monitor_std(TS,PetscInt,PetscReal,Vec,void*);
qvec dm_dummy,dm32;
operator op_list[6];
vec_op *atoms;
operator *atomsstd;
FILE *data_fp = NULL;

typedef struct {
  PetscScalar omega,delta;
  PetscReal deltat,stime;
} PulseParams;

int main(int argc,char *args[]){


  PetscInt n_atoms,i,n_levels,n_seqgroups,max_seqgroupsize,val_init,pos1,pos2,dmpos,dmstdpos;
  PetscScalar tmp_scalar = 1.0,b_field,b_dr,b_0r,b_1r,gamma_r,valpar,diagsum;
  PetscInt steps_max;
  qvec dm,dmstd;
  qsystem qsys,qsysstd,qsys32;
  PetscReal dt,time_max;
  PetscReal single_qubit_gate_time=0.1,two_qubit_gate_time=0.1,var,fidelity;
  circuit circ;
  //PulseParams pulse_params[2]; //moved to line 61
  //State identifiers
  
  dmstdpos=atoi(args[1]);
	dmpos=atoi(args[2]);
	printf("dmstd=%d, dm=%d\n",dmpos,dmstdpos);
	
  enum STATE {zero=0,d,one,r};

  /* Initialize QuaC */
  QuaC_initialize(argc,args);

  //Initialize the qsystem
  initialize_system(&qsys);
  initialize_system(&qsysstd);
  initialize_system(&qsys32);

  n_atoms = 5;
  n_levels = 4;
  
  //Parameters for doing sequential operation
  
  //number of sequential operations in time domain
  n_seqgroups = 4;
  //number of qubits in each sequential operation
  PetscInt seqgroupsize[4] = {2,2,2,2};
  //maximum element in n_seqgroupsize
  max_seqgroupsize=2;
  //indices of atoms in each sequential operation
  PetscInt seqgroup[n_seqgroups][max_seqgroupsize];  
  seqgroup[0][0]=0, seqgroup[0][1]=1;
  seqgroup[1][0]=0, seqgroup[1][1]=2; 
  seqgroup[2][0]=0, seqgroup[2][1]=3;
  seqgroup[3][0]=0, seqgroup[3][1]=4; 
  //parameters for each pulse in the sequence
  PulseParams pulse_params[n_seqgroups];
  for(i=0;i<n_seqgroups;i++){

    pulse_params[i].stime = i* 2;
    
    pulse_params[i].omega = 17.0*(2*3.1416); //MHz 2pi?
  	pulse_params[i].deltat = 0.2; //us
  	pulse_params[i].delta = -0.50*17.0*(2*3.1416); //MHz 2pi?
  }
  //branching ratios
  //Careful to use 1.0 instead of 1 because of integer division
  b_1r = 1.0/16.0;
  b_0r = 1.0/16.0;
  b_dr = 7.0/8.0;

  //decay rate
  gamma_r = 1.0/(540.0);

  //magnetic field
  b_field = 600*2*PETSC_PI; //2pi?

  //Create the operators for the atoms
  atoms = malloc(n_atoms*sizeof(vec_op));
  atomsstd = malloc(n_atoms*sizeof(operator));
  
  for(i=0;i<n_atoms;i++){ 
    //Create an n_level system which is stored in atoms[i]
    create_vec_op_sys(qsys,n_levels,&(atoms[i]));
    create_op_sys(qsysstd,2,&(atomsstd[i]));
    create_op_sys(qsys32,2,&(atomsstd[i]));
  }

  //data_fp = fopen("neutral_atom_3atom_seq111.dat","w");
  // PetscFPrintf(PETSC_COMM_WORLD,data_fp,"#Step_num time omega delta |10><10| |r0><r0| |d0><d0|\n");
  
  //Add hamiltonian terms

  
  for(int i=0;i<n_seqgroups;i++){
  	for(int j=0;j<seqgroupsize[i];j++){
  		tmp_scalar = 1.0;
    	//1.0 * omega(t) * |r><1|
    	add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params[i],omega,2,atoms[seqgroup[i][j]][r],atoms[seqgroup[i][j]][one]);
    	//1.0 * omega(t) * |1><r|
    	add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params[i],omega,2,atoms[seqgroup[i][j]][one],atoms[seqgroup[i][j]][r]);
   		//1.0 * delta(t) * |r><r|
    	add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params[i],delta,2,atoms[seqgroup[i][j]][r],atoms[seqgroup[i][j]][r]);
	  }
  }
  /* 
  for(i=0;i<n_atoms;i++){
    tmp_scalar = 1.0;
    //1.0 * omega(t) * |r><1|
    add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params,omega,2,atoms[i][r],atoms[i][one]);
    //1.0 * omega(t) * |1><r|
    add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params,omega,2,atoms[i][one],atoms[i][r]);

    //1.0 * delta(t) * |r><r|
    add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params,delta,2,atoms[i][r],atoms[i][r]);
  }
  */
  //Coupling term
  //b_field * (|r_0> <r_0|) (|r_1><r_1|) = (|r_0 r_1><r_0 r_1|)
  add_ham_term(qsys,b_field,4,atoms[0][r],atoms[0][r],atoms[1][r],atoms[1][r]);
  add_ham_term(qsys,b_field,4,atoms[0][r],atoms[0][r],atoms[2][r],atoms[2][r]);
  add_ham_term(qsys,b_field,4,atoms[0][r],atoms[0][r],atoms[3][r],atoms[3][r]);
  add_ham_term(qsys,b_field,4,atoms[0][r],atoms[0][r],atoms[4][r],atoms[4][r]);

  //Add lindblad terms
  for(i=0;i<n_atoms;i++){
    tmp_scalar = b_0r*gamma_r;
    //tmp_scalar * L(|0><r|) - no sqrt needed because lin term wants squared term
    add_lin_term(qsys,tmp_scalar,2,atoms[i][zero],atoms[i][r]);

    tmp_scalar = b_1r*gamma_r;
    //L(|1><r|)
    add_lin_term(qsys,tmp_scalar,2,atoms[i][one],atoms[i][r]);

    tmp_scalar = b_dr*gamma_r;
    //L(|d><r|)
    add_lin_term(qsys,tmp_scalar,2,atoms[i][d],atoms[i][r]);

  }
  //Now that we've added all the terms, we construct the matrix
  construct_matrix(qsys);

  add_ham_term(qsysstd,0.0,1,atomsstd[0]->n);
  add_lin_term(qsysstd,0.0,1,atomsstd[0]->n);  
  //
  construct_matrix(qsysstd);
  add_ham_term(qsys32,0.0,1,atomsstd[0]->n);
  //add_lin_term(qsysstd,0.0,1,atomsstd[0]->n);
  construct_matrix(qsys32);
  /* print_mat_sparse(qsys->mat_A); */
  //Create a density matrix object
  
  create_dm_sys(qsysstd,&(dmstd));
  //create_qvec_sys()
  //crete_wf_sys()
  create_dm_sys(qsysstd,&(dm32));
  create_qvec_sys(qsys,&(dm));
  create_qvec_sys(qsys,&(dm_dummy));
  /*
  val_init=0;
  for(i=0;i<4;i++){
  	int r=rand()%2;
  	val_init+=r*pow(2,2*i+1);
  	printf("%d %d %d\n",4-i,r,val_init);
	}
	printf("%d\n",val_init);
	*/

  add_to_qvec(dm,1.0,dmpos,dmpos); //start in the |111><11| state
  add_to_qvec(dmstd,1.0,dmstdpos,dmstdpos); //start in the |111><11| state
  assemble_qvec(dm); 
  assemble_qvec(dmstd);
  time_max  = 3;
  dt        = 0.00001;
  steps_max = 100000;

  /* Set the ts_monitor to print results at each time step */
  set_ts_monitor_sys(qsys,ts_monitor,&pulse_params);
  
  //trace: sum the diagonal
  //population in r and d states
  
  //Do the timestepping
  //  print_qvec(dm);
  printf("---------------------\n");
//  print_qvec(dm);

  time_step_sys(qsys,dm,0.0,time_max/4,dt,steps_max);
  time_step_sys(qsys,dm,time_max/4,time_max/2,dt,steps_max);
  time_step_sys(qsys,dm,time_max/2,3*time_max/4,dt,steps_max);
  time_step_sys(qsys,dm,3*time_max/4,time_max,dt,steps_max);

  int enum_list[2];
  enum_list[0] = zero;
  enum_list[1] = one;
  valpar=0.0;
  diagsum=0.0;
  for(int l_00=0;l_00<2;l_00++){ //atom 0, 0th element  
    for(int l_01=0;l_01<2;l_01++){ //atom 0, 1st element
      for(int l_10=0;l_10<2;l_10++){ //atom 1, 0th element
        for(int l_11=0;l_11<2;l_11++){ //atom 1, 1st element
          for(int l_20=0;l_20<2;l_20++){ //atom 2, 0th element
            for(int l_21=0;l_21<2;l_21++){ //atom 2, 1st element
                for(int l_30=0;l_30<2;l_30++){ //atom 3, 0th element
            			for(int l_31=0;l_31<2;l_31++){ //atom 3, 1st element
            			  for(int l_40=0;l_40<2;l_40++){ //atom 4, 0th element
            					for(int l_41=0;l_41<2;l_41++){ //atom 4, 1st element
              					op_list[0] = atoms[0][enum_list[l_00]]; //atom 0 0th element
              					op_list[1] = atoms[0][enum_list[l_01]]; //atom 0 1st element
             						op_list[2] = atoms[1][enum_list[l_10]]; //atom 1 0th element
              					op_list[3] = atoms[1][enum_list[l_11]]; //atom 1 1st element
              					op_list[4] = atoms[2][enum_list[l_20]]; //atom 2 0th element
              					op_list[5] = atoms[2][enum_list[l_21]]; //atom 2 1st element
              					op_list[6] = atoms[3][enum_list[l_30]]; //atom 2 0th element
              					op_list[7] = atoms[3][enum_list[l_31]]; //atom 2 1st element
              					op_list[8] = atoms[4][enum_list[l_40]]; //atom 2 0th element
              					op_list[9] = atoms[4][enum_list[l_41]]; //atom 2 1st element
              					get_expectation_value_qvec_list(dm,&valpar,10,op_list);
              					pos1=l_00*16+l_10*8+l_20*4+l_30*2+l_40*1;
												pos2=l_01*16+l_11*8+l_21*4+l_31*2+l_41*1;
												if (pos1==pos2){
													printf("%d %d %f\n",pos1,pos2,valpar);
													diagsum=diagsum+valpar;
												}
              					add_to_qvec(dm32,valpar,pos1,pos2);
              				}
              		}
               }
             }
            }
          }
        }
      }
    }
  }
  assemble_qvec(dm32);
  printf("dm32 constructed\n");
  //print_qvec(dm32);  
//----------------------------------------------------------------------------------------------
  set_ts_monitor_sys(qsysstd,ts_monitor_std,&pulse_params);
  create_circuit(&circ,5);
  //Add some gates
  add_gate_to_circuit_sys(&circ,0.0,CZ_ARP,0,1);
  add_gate_to_circuit_sys(&circ,0.5,CZ_ARP,0,2);
  add_gate_to_circuit_sys(&circ,1.0,CZ_ARP,0,3);
  add_gate_to_circuit_sys(&circ,1.5,CZ_ARP,0,4);
  single_qubit_gate_time = 0.1;
  two_qubit_gate_time = 0.2;
  for (i=0;i<circ.num_gates;i++){
    if (circ.gate_list[i].num_qubits==1){
      circ.gate_list[i].run_time = single_qubit_gate_time;
      time_max += single_qubit_gate_time;
    } else if (circ.gate_list[i].num_qubits==2){
      circ.gate_list[i].run_time = two_qubit_gate_time;
      time_max += two_qubit_gate_time;
    }
  }
  schedule_circuit_layers(qsysstd,&circ);
  //Start out circuit at time 0.0, first gate will be at 0 
  //printf("Before\n");
  //print_qvec(dmstd);  
  apply_circuit_to_qvec(qsysstd,circ,dmstd);
  //printf("dmstd\n");
 // print_qvec(dmstd);
//-----------------------------------------------------------------------------------------------------------
  get_fidelity_qvec(dm32,dmstd,&fidelity,&var);
  printf("fidelity between seq(32*32) and std is %lf\n",fidelity);
  printf("sum of the diag is %f\n",diagsum);
  //print_qvec(dm);
  //clean up memory
  data_fp = fopen("seq_sp.dat","a");
  
  PetscFPrintf(PETSC_COMM_WORLD,data_fp,"%d %d %f %f\n",dmpos,dmstdpos,fidelity,diagsum);
  
  destroy_system(&qsys);
  destroy_qvec(&(dm));
  //destroy_qvec(&(dm_dummy));
  destroy_system(&qsysstd);
  destroy_qvec(&(dmstd));
  for(i=0;i<n_atoms;i++){
    destroy_vec_op_sys(&atoms[i]);
    destroy_op_sys(&atomsstd[i]);
  }

  return;
}  

PetscErrorCode ts_monitor(TS ts,PetscInt step,PetscReal time,Vec rho_data,void *ctx){
  PulseParams *pulse_params = (PulseParams*) ctx;   /* user-defined struct */
  PetscScalar trace_val,trace_val2,trace_val3,trace_val4,trace_val5,trace_val6,trace_val7,trace_val8;
  Vec tmp_data;
  enum STATE {zero=0,d,one,r};
  //Print out things at each time step, if desired
  //tmp data necessary for technical reasons
  tmp_data = dm_dummy->data;
  dm_dummy->data = rho_data;
  /*
  //ev( (|1><1|)(|0><0|) )= ev(|10><10|)
  //get_expectation_value_qvec_list
  get_expectation_value_qvec(dm_dummy,&trace_val,6,atoms[0][one],atoms[0][one],atoms[1][one],atoms[1][one],atoms[2][one],atoms[2][one]);
  get_expectation_value_qvec(dm_dummy,&trace_val2,6,atoms[0][one],atoms[0][one],atoms[1][one],atoms[1][one],atoms[2][r],atoms[2][r]);
  get_expectation_value_qvec(dm_dummy,&trace_val3,6,atoms[0][one],atoms[0][one],atoms[1][r],atoms[1][r],atoms[2][one],atoms[2][one]);
  get_expectation_value_qvec(dm_dummy,&trace_val4,6,atoms[0][r],atoms[0][r],atoms[1][one],atoms[1][one],atoms[2][one],atoms[2][one]);
  get_expectation_value_qvec(dm_dummy,&trace_val5,6,atoms[0][one],atoms[0][one],atoms[1][r],atoms[1][r],atoms[2][r],atoms[2][r]);
  get_expectation_value_qvec(dm_dummy,&trace_val6,6,atoms[0][r],atoms[0][r],atoms[1][one],atoms[1][one],atoms[2][r],atoms[2][r]);
  get_expectation_value_qvec(dm_dummy,&trace_val7,6,atoms[0][r],atoms[0][r],atoms[1][r],atoms[1][r],atoms[2][one],atoms[2][one]);
  get_expectation_value_qvec(dm_dummy,&trace_val8,6,atoms[0][r],atoms[0][r],atoms[1][r],atoms[1][r],atoms[2][r],atoms[2][r]);

  
  PetscFPrintf(PETSC_COMM_WORLD,data_fp,"%d %f %f %f %f %f %f %f %f %f\n",step,time,PetscRealPart(trace_val),PetscRealPart(trace_val2),PetscRealPart(trace_val3),PetscRealPart(trace_val4),PetscRealPart(trace_val5),PetscRealPart(trace_val6),PetscRealPart(trace_val7),PetscRealPart(trace_val8));
  dm_dummy->data = tmp_data;
  /* print_qvec(rho); */
  PetscFunctionReturn(0);
}
PetscErrorCode ts_monitor_std(TS ts,PetscInt step,PetscReal time,Vec rho_data,void *ctx){
	  PetscFunctionReturn(0);
}

//Define time dependent pulses

PetscScalar omega(PetscReal time,void *ctx){
  PetscScalar pulse_value;
  PetscReal tau,dt,p,a,ts;
  PulseParams *pulse_params = (PulseParams*) ctx;   /* user-defined struct */


  //I don't know what the true form is, so I used two gaussians,
  //one centered at l/4, one centered at 3l/4
  dt = pulse_params->deltat;
  ts = pulse_params->stime;

  p=exp(-pow((time-ts-5*dt),2)/pow(dt,2));

  pulse_value = pulse_params->omega/2.0*p;

  return pulse_value;
}


PetscScalar delta(PetscReal time,void *ctx){
  PetscScalar pulse_value;
  PulseParams *pulse_params = (PulseParams*) ctx;   /* user-defined struct */

  pulse_value = pulse_params->delta;
  return pulse_value;
}


