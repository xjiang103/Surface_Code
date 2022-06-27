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
qvec dm_dummy;
vec_op *atoms;
FILE *data_fp = NULL;
operator op_list[4];
operator *atoms_arp;

typedef struct {
  PetscScalar omega,delta;
  PetscReal deltat;
} PulseParams;

int main(int argc,char **args){
  PetscInt n_atoms,i,n_levels,pos1,pos2;
  PetscScalar tmp_scalar = 1.0,b_field,b_dr,b_0r,b_1r,gamma_r,tmp_val,val_arp;
  PetscInt steps_max;
  qvec dm,dm_std,dm_arp;
  qsystem qsys,qsys_arp;
  PetscReal dt,time_max,val,fidelity;
  PulseParams pulse_params;
  circuit circ;
  //State identifiers
  enum STATE {zero=0,d,one,r};

  /* Initialize QuaC */
  QuaC_initialize(argc,args);

  //Initialize the qsystem
  initialize_system(&qsys);
  initialize_system(&qsys_arp);


  n_atoms = 2;
  n_levels = 4;

  //branching ratios
  //Careful to use 1.0 instead of 1 because of integer division
  b_1r = 1.0/16.0;
  b_0r = 1.0/16.0;
  b_dr = 7.0/8.0;

  //decay rate
  gamma_r = 1.0/(130.0);

  //magnetic field
  b_field = 60*(2*3.1416); //2pi?

  //Create the operators for the atoms
  atoms = malloc(n_atoms*sizeof(vec_op));
  atoms_arp = malloc(n_atoms*sizeof(operator));

  for(i=0;i<n_atoms;i++){
    //Create an n_level system which is stored in atoms[i]
    create_vec_op_sys(qsys,n_levels,&(atoms[i]));
    create_op_sys(qsys_arp,2,&(atoms_arp[i]));
  }

  data_fp = fopen("na1_sp.dat","w");
  PetscFPrintf(PETSC_COMM_WORLD,data_fp,"#Step_num time omega delta |10><10| |r0><r0| |d0><d0|\n");
  //Add hamiltonian terms
  pulse_params.omega = 17.0*(2*3.1416); //MHz 2pi?
  pulse_params.deltat = 0.2; //us
  pulse_params.delta = -0.5*17.0*(2*3.1416); //MHz 2pi?


  for(i=0;i<n_atoms;i++){
    tmp_scalar = 1.0;
    //1.0 * omega(t) * |r><1|
    add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params,omega,2,atoms[i][r],atoms[i][one]);
    //1.0 * omega(t) * |1><r|
    add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params,omega,2,atoms[i][one],atoms[i][r]);

    //1.0 * delta(t) * |r><r|
    add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params,delta,2,atoms[i][r],atoms[i][r]);
  }

  //Coupling term
  //b_field * (|r_0> <r_0|) (|r_1><r_1|) = (|r_0 r_1><r_0 r_1|)
  add_ham_term(qsys,b_field,4,atoms[0][r],atoms[0][r],atoms[1][r],atoms[1][r]);

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
  add_ham_term(qsys_arp,0.0,1,atoms_arp[0]->n);
  add_lin_term(qsys_arp,0.0,1,atoms_arp[0]->n);
  construct_matrix(qsys_arp);
  /* print_mat_sparse(qsys->mat_A); */
  //Create a density matrix object
  create_qvec_sys(qsys,&(dm));
  create_dm_sys(qsys_arp,&(dm_std)); 
  create_dm_sys(qsys_arp,&(dm_arp));
  create_qvec_sys(qsys,&(dm_dummy));
/*
	FILE *fp_read=NULL;
	fp_read=fopen("init_11.txt","r");
	PetscInt init_enum[4]={0,2,8,10};
	PetscReal init_real[4]={0};
	PetscReal init_imag[4]={0};
	PetscScalar init_state[4];
	PetscScalar init_state_conj[4];

	for(int i=0;i<4;i++){
		fscanf(fp_read,"%lf %lf",&init_real[i],&init_imag[i]);
		//printf("%d\n",i);
		//printf("%lf + %lf j\n", init_real[i], init_imag[i]);
		init_state[i]=init_real[i] + init_imag[i]*PETSC_i;
		init_state_conj[i]=init_real[i] - init_imag[i]*PETSC_i;
	}

  for(int i=0;i<4;i++){
  	for (int j=0;j<4;j++){
      tmp_val = init_state[i]*init_state_conj[j];
      printf("Added dm[%d][%d]: %f %f\n",init_enum[i],init_enum[j],tmp_val);
  		add_to_qvec(dm,init_state[i]*init_state_conj[j],init_enum[i],init_enum[j]);

//  	  printf("%lf + %lf j, ",(double)PetscRealPart(init_state[i]*init_state_conj[j]),(double)PetscImaginaryPart(init_state[i]*init_state_conj[j]));
		}
//		printf("\n");
	}
 //start in the |10><10| state
 */
 	
  add_to_qvec(dm,1.0,10,10); //start in the |111><11| state
  assemble_qvec(dm);

	printf("Intial dm is:\n");
	print_qvec(dm);
  time_max  = 3;
  dt        = 0.00001;
  steps_max = 100000;

  /* Set the ts_monitor to print results at each time step */
  set_ts_monitor_sys(qsys,ts_monitor,&pulse_params);

  //Do the timestepping
  time_step_sys(qsys,dm,0.0,time_max,dt,steps_max);

  //print_qvec(dm);
  //clean up memory

  int enum_list[2];
  enum_list[0] = zero;
  enum_list[1] = one;
  val_arp=0.0;
  for(int l_00=0;l_00<2;l_00++){ //atom 0, 0th element
    for(int l_01=0;l_01<2;l_01++){ //atom 0, 1st element
      for(int l_10=0;l_10<2;l_10++){ //atom 1, 0th element
        for(int l_11=0;l_11<2;l_11++){ //atom 1, 1st element
              op_list[0] = atoms[0][enum_list[l_00]]; //atom 0 0th element
              op_list[1] = atoms[0][enum_list[l_01]]; //atom 0 1st element
              op_list[2] = atoms[1][enum_list[l_10]]; //atom 1 0th element
              op_list[3] = atoms[1][enum_list[l_11]]; //atom 1 1st element

              get_expectation_value_qvec_list(dm,&val_arp,4,op_list);
              pos1=l_00*2+l_10*1;
							pos2=l_01*2+l_11*1;
              add_to_qvec(dm_arp,val_arp,pos1,pos2);
            }
          }
        }
      }

  assemble_qvec(dm_arp);
	printf("After CZ_ARP, dm_arp is:\n");  
	print_qvec(dm_arp);
	
  create_circuit(&circ,3);
  add_gate_to_circuit_sys(&circ,0.0,HADAMARD,0);
  schedule_circuit_layers(qsys_arp,&circ);
  apply_circuit_to_qvec(qsys_arp,circ,dm_arp);

	printf("After 2nd Hadamard, dm_arp is:\n");
	print_qvec(dm_arp);
	
  fp_read=fopen("final_11.txt","r");

	PetscReal final_real[16]={0};
	PetscReal final_imag[16]={0};
	PetscScalar final_state[16];

	for(int i=0;i<16;i++){
		fscanf(fp_read,"%lf %lf",&final_real[i],&final_imag[i]);
		//printf("%d\n",i);
		//printf("%lf + %lf j\n", final_real[i], final_imag[i]);
		final_state[i]=final_real[i] + final_imag[i]*PETSC_i;
	}

  for(int i=0;i<4;i++){
  	for (int j=0;j<4;j++){
  		int pos=i*4+j*1;
      tmp_val = final_state[pos];
      //printf("Added dm[%d][%d]: %f %f\n",i,j,tmp_val);
  		add_to_qvec(dm_std,tmp_val,i,j);

//  	  printf("%lf + %lf j, ",(double)PetscRealPart(init_state[i]*init_state_conj[j]),(double)PetscImaginaryPart(init_state[i]*init_state_conj[j]));
		}
//		printf("\n");
	}
	assemble_qvec(dm_std);
	printf("dm_std from qutip is:\n");
	print_qvec(dm_std);

	
	val=1;
	get_fidelity_qvec(dm_arp,dm_std,&fidelity,&val);
  printf("fidelity between dm_arp and dm_std is %lf\n",fidelity);
  
  destroy_system(&qsys);
  destroy_qvec(&(dm));
  destroy_qvec(&(dm_std));
  destroy_qvec(&(dm_arp));
  destroy_qvec(&(dm_dummy));
  for(i=0;i<n_atoms;i++){
    destroy_vec_op_sys(&atoms[i]);
  }

  return;
}

PetscErrorCode ts_monitor(TS ts,PetscInt step,PetscReal time,Vec rho_data,void *ctx){
  PulseParams *pulse_params = (PulseParams*) ctx;   /* user-defined struct */
  PetscScalar trace_val,trace_val2,trace_val3,trace_val4,trace_val5;
  Vec tmp_data;
  enum STATE {zero=0,d,one,r};
  //Print out things at each time step, if desired
  //tmp data necessary for technical reasons
  tmp_data = dm_dummy->data;
  dm_dummy->data = rho_data;

  //ev( (|1><1|)(|0><0|) )= ev(|10><10|)
  get_expectation_value_qvec(dm_dummy,&trace_val,4,atoms[0][one],atoms[0][one],atoms[1][one],atoms[1][one]);
  get_expectation_value_qvec(dm_dummy,&trace_val2,4,atoms[0][one],atoms[0][one],atoms[1][r],atoms[1][r]);
  get_expectation_value_qvec(dm_dummy,&trace_val3,4,atoms[0][r],atoms[0][r],atoms[1][one],atoms[1][one]);
  get_expectation_value_qvec(dm_dummy,&trace_val4,4,atoms[0][r],atoms[0][r],atoms[1][r],atoms[1][r]);
  get_expectation_value_qvec(dm_dummy,&trace_val5,4,atoms[0][d],atoms[0][d],atoms[1][one],atoms[1][one]);
  get_expectation_value_qvec(dm_dummy,&trace_val5,4,atoms[0][d],atoms[0][d],atoms[1][one],atoms[1][one]);

  if(step>0){
     PetscFPrintf(PETSC_COMM_WORLD,data_fp,"%d %f %f %f %f %f %f %f %f\n",step,time,PetscRealPart(omega(time,ctx)),PetscRealPart(delta(time,ctx)),PetscRealPart(trace_val),PetscRealPart(trace_val2),PetscRealPart(trace_val3),PetscRealPart(trace_val4),PetscRealPart(trace_val5));
    }
  dm_dummy->data = tmp_data;
  /* print_qvec(rho); */
  PetscFunctionReturn(0);
}


//Define time dependent pulses
PetscScalar omega(PetscReal time,void *ctx){
  PetscScalar pulse_value;
  PetscReal tau,dt,p,a;
  PulseParams *pulse_params = (PulseParams*) ctx;   /* user-defined struct */


  //I don't know what the true form is, so I used two gaussians,
  //one centered at l/4, one centered at 3l/4
  dt = pulse_params->deltat;

  p=exp(-pow((time-3*dt),2)/pow(dt,2));

  pulse_value = pulse_params->omega/2.0*p;

  return pulse_value;
}


PetscScalar delta(PetscReal time,void *ctx){
  PetscScalar pulse_value;
  PulseParams *pulse_params = (PulseParams*) ctx;   /* user-defined struct */

  pulse_value = pulse_params->delta;
  return pulse_value;
}
