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
  PetscInt steps_max,n_ens=-1;
  qvec dm,dmstd;
  qsystem qsys,qsysstd;
  PetscReal dt,time_max;
  PetscReal single_qubit_gate_time=0.1,two_qubit_gate_time=0.1,var,fidelity;
  char bitstr[PETSC_MAX_PATH_LEN] = "111"; //Default bitstr to start with
  circuit circ;
  int length;
  //PulseParams pulse_params[2]; //moved to line 61
  //State identifiers

  enum STATE {zero=0,one,r};

  /* Initialize QuaC */ 
  QuaC_initialize(argc,args);
  //Get the bitstring we want to simulate
  PetscOptionsGetString(NULL,NULL,"-bitstr",bitstr,PETSC_MAX_PATH_LEN,NULL);
  PetscOptionsGetInt(NULL,NULL,"-n_ens",&n_ens,NULL);

  dmpos = 0;
  dmstdpos = 0;
  length = strlen(bitstr);
  if(length!=3){
    PetscPrintf(PETSC_COMM_WORLD,"ERROR: bitstr must be of length 3!\n");
    exit(8);
  }
  //Convert from the bitstr to the dmpos and dmstdpos
  for(i=0;i<length;i++){
    //We use length-1-i to go through the list in reverse, because we want 00001 to be dmpos=2
    if(bitstr[length-1-i]=='0'){ //Must use single apostrophe for character equality
      dmstdpos += 0*pow(2,i);
      dmpos += 0*pow(3,i);
    } else if(bitstr[length-1-i]=='1') {
      dmstdpos += 1*pow(2,i);
      dmpos += 1*pow(3,i);
    } else {
      PetscPrintf(PETSC_COMM_WORLD,"ERROR: Must be 0 or 1\n");
      exit(0);
    }
  }
  PetscPrintf(PETSC_COMM_WORLD,"Simulating bitstr %s with dmpos=%d and dmstdpos=%d\n",bitstr,dmpos,dmstdpos);

  //Initialize the qsystem
  initialize_system(&qsys);
  initialize_system(&qsysstd);

  n_atoms = 3;
  n_levels = 3;

  //Parameters for doing sequential operation

  //number of sequential operations in time domain
  n_seqgroups = 2;
  //number of qubits in each sequential operation
  PetscInt seqgroupsize[2] = {2,2};
  //maximum element in n_seqgroupsize
  max_seqgroupsize=2;
  //indices of atoms in each sequential operation
  PetscInt seqgroup[n_seqgroups][max_seqgroupsize];
  seqgroup[0][0]=0, seqgroup[0][1]=1;
  seqgroup[1][0]=0, seqgroup[1][1]=2;

  //parameters for each pulse in the sequence
  PulseParams pulse_params[n_seqgroups];
  for(i=0;i<n_seqgroups;i++){

    pulse_params[i].stime = i* 0.0;

    pulse_params[i].omega = 17.0*(2*3.1416); //MHz 2pi?
  	pulse_params[i].deltat = 0.2; //us
  	pulse_params[i].delta = 0.50*17.0*(2*3.1416); //MHz 2pi?
  }
  //branching ratios
  //Careful to use 1.0 instead of 1 because of integer division
  b_1r = 0.0/16.0;
  b_0r = 0.0/16.0;
  b_dr = 0.0*7.0/8.0; 

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
  }

  data_fp = fopen("na_1117_3a.txt","w");
  
  //Add hamiltonian terms

  for(i=0;i<n_atoms;i++){
    tmp_scalar = 1.0;
    //1.0 * omega(t) * |r><1|
    add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params,omega,2,atoms[i][r],atoms[i][one]);
    //1.0 * omega(t) * |1><r|
    add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params,omega,2,atoms[i][one],atoms[i][r]);

    //1.0 * delta(t) * |r><r|
    add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params,delta,2,atoms[i][r],atoms[i][r]);
  }
  /*
  for(int i=0;i<n_seqgroups;i++){
  	for(int j=0;j<seqgroupsize[i];j++){
  		tmp_scalar = 1.0;
    	//1.0 * omega(t) * |r><1|
    	add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params[i],omega,2,atoms[seqgroup[i][j]][r],atoms[seqgroup[i][j]][one]);
    	//1.0 * omega(t) * |1><r|
    	add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params[i],omega,2,atoms[seqgroup[i][j]][one],atoms[seqgroup[i][j]][r]);
   		//1.0 * delta(t) * |r><r|
      tmp_scalar = pulse_params[i].delta;
      //MJO: Delta has no time dependence, so we use a normal ham term
      add_ham_term(qsys,tmp_scalar,2,atoms[seqgroup[i][j]][r],atoms[seqgroup[i][j]][r]);
	  }
  }
  */
	float ddfac=0;
	printf("ddfac is %f\n",ddfac);
  //Coupling term
  //b_field * (|r_0> <r_0|) (|r_1><r_1|) = (|r_0 r_1><r_0 r_1|)
  add_ham_term(qsys,b_field,4,atoms[0][r],atoms[0][r],atoms[1][r],atoms[1][r]); 
  add_ham_term(qsys,b_field,4,atoms[0][r],atoms[0][r],atoms[2][r],atoms[2][r]);
  add_ham_term(qsys,ddfac*b_field,4,atoms[1][r],atoms[1][r],atoms[2][r],atoms[2][r]);

  //Add lindblad terms
  for(i=0;i<n_atoms;i++){
    tmp_scalar = b_0r*gamma_r;
    //tmp_scalar * L(|0><r|) - no sqrt needed because lin term wants squared term
    add_lin_term(qsys,tmp_scalar,2,atoms[i][zero],atoms[i][r]);

    tmp_scalar = b_1r*gamma_r;
    //L(|1><r|)
    add_lin_term(qsys,tmp_scalar,2,atoms[i][one],atoms[i][r]);

    //L(|0><1|) for testing purposes
    //tmp_scalar = 10*gamma_r;
    //add_lin_term(qsys,tmp_scalar,2,atoms[i][zero],atoms[i][one]);

    //tmp_scalar = b_dr*gamma_r;
    //L(|d><r|)
    /* add_lin_term(qsys,tmp_scalar,2,atoms[i][d],atoms[i][r]); */
    // Instead of a lindblad term, we do an imaginary hamiltonian to allow for leakage
    // need to have a different tmp_scalar based on the repumping simulations
    tmp_scalar = b_dr*gamma_r*PETSC_i;
    add_ham_term(qsys,tmp_scalar,2,atoms[i][r],atoms[i][r]);

  }
  //Now that we've added all the terms, we construct the matrix
  if(n_ens>0){
    use_mcwf_solver(qsys,n_ens,NULL);
  }
  construct_matrix(qsys);

  add_ham_term(qsysstd,0.0,1,atomsstd[0]->n);
  add_lin_term(qsysstd,0.0,1,atomsstd[0]->n);
  //

  construct_matrix(qsysstd);

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
  time_max  = 4;
  dt        = 0.0001;
  steps_max = 100000;

  /* Set the ts_monitor to print results at each time step */
  set_ts_monitor_sys(qsys,ts_monitor,&pulse_params); 


  //Do the timestepping
  PetscPrintf(PETSC_COMM_WORLD,"---------------------\n");

  time_step_sys(qsys,dm,0.0,time_max,dt,steps_max);
  PetscPrintf(PETSC_COMM_WORLD,"Timestep 1 done\n");

  /* time_step_sys(qsys,dm,time_max/4,time_max/2,dt,steps_max); */
  /* PetscPrintf(PETSC_COMM_WORLD,"Timestep 2 done\n"); */

  /* time_step_sys(qsys,dm,time_max/2,3*time_max/4,dt,steps_max); */
  /* PetscPrintf(PETSC_COMM_WORLD,"Timestep 3 done\n"); */

  /* time_step_sys(qsys,dm,3*time_max/4,time_max,dt,steps_max); */
  /* PetscPrintf(PETSC_COMM_WORLD,"Timestep 4 done\n"); */

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
    					op_list[0] = atoms[0][enum_list[l_00]]; //atom 0 0th element
    					op_list[1] = atoms[0][enum_list[l_01]]; //atom 0 1st element
   						op_list[2] = atoms[1][enum_list[l_10]]; //atom 1 0th element
    					op_list[3] = atoms[1][enum_list[l_11]]; //atom 1 1st element
    					op_list[4] = atoms[2][enum_list[l_20]]; //atom 2 0th element
    					op_list[5] = atoms[2][enum_list[l_21]]; //atom 2 1st element
    					get_expectation_value_qvec_list(dm,&valpar,6,op_list);
    					pos1=l_00*4+l_10*2+l_20*1;
							pos2=l_01*4+l_11*2+l_21*1;
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

  assemble_qvec(dm32);
  PetscPrintf(PETSC_COMM_WORLD,"dm32 constructed\n");
  /* print_qvec(dm32); */
//----------------------------------------------------------------------------------------------
  create_circuit(&circ,3);
  //Add some gates
  add_gate_to_circuit_sys(&circ,0.0,CZ_ARP,0,1);
  add_gate_to_circuit_sys(&circ,0.5,CZ_ARP,0,2);

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
  apply_circuit_to_qvec(qsysstd,circ,dmstd);

//-----------------------------------------------------------------------------------------------------------
  get_fidelity_qvec(dm32,dmstd,&fidelity,&var);
  printf("fidelity between seq(32*32) and std is %lf\n",fidelity);
  printf("sum of the diag is %f\n",diagsum);

  //print_qvec(dm);
  //clean up memory
  /* data_fp = fopen("seq_sp_3lvl.dat","a"); */

  /* PetscFPrintf(PETSC_COMM_WORLD,data_fp,"%d %d %f %f\n",dmpos,dmstdpos,fidelity,diagsum); */
  /* fclose(data_fp); */

  destroy_qvec(&(dm));
  destroy_qvec(&(dm_dummy));
  destroy_qvec(&(dmstd));
  destroy_qvec(&(dm32));

  for(i=0;i<n_atoms;i++){
    destroy_vec_op_sys(&atoms[i]);
    destroy_op_sys(&atomsstd[i]);
  }
  destroy_system(&qsys);
  destroy_system(&qsysstd);

  QuaC_finalize();
  return;
}

PetscErrorCode ts_monitor(TS ts,PetscInt step,PetscReal time,Vec rho_data,void *ctx){
  PulseParams *pulse_params = (PulseParams*) ctx;   /* user-defined struct */
  PetscScalar trace_val,trace_val2,trace_val3,trace_val4,trace_val5,trace_val6,trace_val7,trace_val8;
  Vec tmp_data;
  enum STATE {zero=0,one,r};
  //Print out things at each time step, if desired
  //tmp data necessary for technical reasons
  tmp_data = dm_dummy->data;
  dm_dummy->data = rho_data;

  //ev( (|1><1|)(|0><0|) )= ev(|10><10|)
  //get_expectation_value_qvec_list
  get_expectation_value_qvec(dm_dummy,&trace_val,6,atoms[0][one],atoms[0][one],atoms[1][one],atoms[1][one],atoms[2][one],atoms[2][one]);
  get_expectation_value_qvec(dm_dummy,&trace_val2,6,atoms[0][r],atoms[0][r],atoms[1][one],atoms[1][one],atoms[2][one],atoms[2][one]);
  get_expectation_value_qvec(dm_dummy,&trace_val3,6,atoms[0][one],atoms[0][one],atoms[1][r],atoms[1][r],atoms[2][one],atoms[2][one]);
  get_expectation_value_qvec(dm_dummy,&trace_val4,6,atoms[0][one],atoms[0][one],atoms[1][one],atoms[1][one],atoms[2][r],atoms[2][r]);
  get_expectation_value_qvec(dm_dummy,&trace_val5,6,atoms[0][r],atoms[0][r],atoms[1][r],atoms[1][r],atoms[2][one],atoms[2][one]);
  get_expectation_value_qvec(dm_dummy,&trace_val6,6,atoms[0][one],atoms[0][one],atoms[1][r],atoms[1][r],atoms[2][r],atoms[2][r]);
  get_expectation_value_qvec(dm_dummy,&trace_val7,6,atoms[0][r],atoms[0][r],atoms[1][one],atoms[1][one],atoms[2][r],atoms[2][r]);
  get_expectation_value_qvec(dm_dummy,&trace_val8,6,atoms[0][r],atoms[0][r],atoms[1][r],atoms[1][r],atoms[2][r],atoms[2][r]);
  
  //if(step%1000==0){
  //	PetscPrintf(PETSC_COMM_WORLD,"%d %f %f %f %f\n",step,time,PetscRealPart(trace_val),PetscRealPart(trace_val2),PetscRealPart(trace_val3));	
	//}
	
  PetscFPrintf(PETSC_COMM_WORLD,data_fp,"%d %f %f %f %f %f %f %f %f %f\n",step,time,PetscRealPart(trace_val),PetscRealPart(trace_val2),PetscRealPart(trace_val3),PetscRealPart(trace_val4),PetscRealPart(trace_val5),PetscRealPart(trace_val6),PetscRealPart(trace_val7),PetscRealPart(trace_val8));

  /* print_qvec(rho); */

  dm_dummy->data = tmp_data;
  PetscFunctionReturn(0);
}

//Define time dependent pulses

PetscScalar omega(PetscReal time,void *ctx){
  PetscScalar pulse_value;
  PetscReal tau,dt,p,a,ts;
  PulseParams *pulse_params = (PulseParams*) ctx;   /* user-defined struct */

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
