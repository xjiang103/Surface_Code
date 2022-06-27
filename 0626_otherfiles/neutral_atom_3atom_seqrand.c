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
//#include "sprng_f.h"
#include "qasm_parser.h"
#include "qsystem.h"


PetscScalar omega(PetscReal time,void *ctx);
PetscScalar delta(PetscReal time,void *ctx);

PetscErrorCode ts_monitor(TS,PetscInt,PetscReal,Vec,void*);
qvec dm_dummy;
vec_op *atoms;
FILE *data_fp = NULL;

typedef struct {
  PetscScalar omega,delta;
  PetscReal length,stime;
} PulseParams;

int main(int argc,char **args){
  PetscInt n_atoms,i,n_levels,n_seqgroups,max_seqgroupsize;
  PetscScalar tmp_scalar = 1.0,b_field,b_dr,b_0r,b_1r,gamma_r;
  PetscInt steps_max;
  qvec dm;
  qsystem qsys;
  PetscReal dt,time_max,rn;
  //PulseParams pulse_params[2]; //moved to line 61
  //State identifiers
  enum STATE {zero=0,d,one,r};

  /* Initialize QuaC */
  QuaC_initialize(argc,args);

  //Initialize the qsystem
  initialize_system(&qsys);

  n_atoms = 3;
  n_levels = 4;
  
  //Parameters for doing sequential operation
  
  //number of sequential operations in time domain
  n_seqgroups = 2;
  //number of qubits in each sequential operation
  PetscInt seqgroupsize[2] = {2,2};
  //maximum element in n_seqgroupsize
  max_seqgroupsize=2;
  //indices of atoms in each sequential operation
  PetscInt seqgroup[n_seqgroups][max_seqgroupsize];  
  seqgroup[0][0]=0, seqgroup[0][1]=2;
  seqgroup[1][0]=1, seqgroup[1][1]=2; 
  //parameters for each pulse in the sequence
  PulseParams pulse_params[n_seqgroups];
  for(i=0;i<n_seqgroups;i++){
  	pulse_params[i].omega = 17.0*2*PETSC_PI; //MHz 2pi?
    pulse_params[i].length = 0.54; //us
    pulse_params[i].delta = 23.0*2*PETSC_PI; //MHz 2pi?
    pulse_params[i].stime = i* 0.54;
  }
  //branching ratios
  //Careful to use 1.0 instead of 1 because of integer division
  b_1r = 1.0/16.0;
  b_0r = 1.0/16.0;
  b_dr = 7.0/8.0;

  //decay rate
  gamma_r = 1.0/(540.0);

  //magnetic field
  b_field = 100*2*PETSC_PI; //2pi?

  //Create the operators for the atoms
  atoms = malloc(n_atoms*sizeof(vec_op));
  
  for(i=0;i<n_atoms;i++){
    //Create an n_level system which is stored in atoms[i]
    create_vec_op_sys(qsys,n_levels,&(atoms[i]));
  }

  data_fp = fopen("neutral_atom_3atom_seqrand.dat","w");
  PetscFPrintf(PETSC_COMM_WORLD,data_fp,"#Step_num time omega delta |10><10| |r0><r0| |d0><d0|\n");
  
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
  add_ham_term(qsys,b_field*0,4,atoms[0][r],atoms[0][r],atoms[1][r],atoms[1][r]);
  add_ham_term(qsys,b_field,4,atoms[0][r],atoms[0][r],atoms[2][r],atoms[2][r]);
  add_ham_term(qsys,b_field,4,atoms[1][r],atoms[1][r],atoms[2][r],atoms[2][r]);
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
  /* print_mat_sparse(qsys->mat_A); */
  //Create a density matrix object
  create_qvec_sys(qsys,&(dm));
  create_qvec_sys(qsys,&(dm_dummy));

  /* Initialize with random state
	0: |000>
	2: |001>
	8: |010>
	10: |011>
	32: |100>
	34: |101>
	40: |110>
	42: |111>
	*/
	
	FILE *fp_read=NULL;
	fp_read=fopen("randket8.txt","r");
	PetscInt init_enum[8]={0,2,8,10,32,34,40,42};
	PetscReal init_real[8]={0};
	PetscReal init_imag[8]={0};
	PetscScalar init_state[8];
	PetscScalar init_state_conj[8];
	
	for(int i=0;i<8;i++){
		fscanf(fp_read,"%lf %lf",&init_real[i],&init_imag[i]);
//		printf("%lf + %lf j\n", init_real[i], init_imag[i]);
		init_state[i]=init_real[i] + init_imag[i]*PETSC_i;
		init_state_conj[i]=init_real[i] - init_imag[i]*PETSC_i;
//		printf("%lf + %lf j\n", (double)PetscRealPart(init_state[i]*init_state_conj[i]), (double)PetscImaginaryPart(init_state[i]*init_state_conj[i]));
	}

  for(int i=0;i<8;i++){
  	for (int j=0;j<8;j++){
  		add_to_qvec(dm,init_state[i]*init_state_conj[j],init_enum[i],init_enum[j]);
  	  printf("%lf + %lf j, ",(double)PetscRealPart(init_state[i]*init_state_conj[j]),(double)PetscImaginaryPart(init_state[i]*init_state_conj[j]));
		}
		printf("\n");
	}
	
	//add_to_qvec(dm,1.0,42,42); //start in the |111><11| state
  assemble_qvec(dm);
  

  time_max  = 1.08;
  dt        = 0.05;
  steps_max = 1000;

  /* Set the ts_monitor to print results at each time step */
  set_ts_monitor_sys(qsys,ts_monitor,&pulse_params);
  //Do the timestepping
  time_step_sys(qsys,dm,0.0,time_max,dt,steps_max);
  time_step_sys(qsys,dm,time_max, time_max,dt,steps_max);
  print_qvec(dm);
  //clean up memory
  destroy_system(&qsys);
  destroy_qvec(&(dm));
  destroy_qvec(&(dm_dummy));
  for(i=0;i<n_atoms;i++){
    destroy_vec_op_sys(&atoms[i]);
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


//Define time dependent pulses
PetscScalar omega(PetscReal time,void *ctx){
  PetscScalar pulse_value;
  PetscReal tau,t1,t2,ts,p,a;
  PulseParams *pulse_params = (PulseParams*) ctx;   /* user-defined struct */


  //I don't know what the true form is, so I used two gaussians,
  //one centered at l/4, one centered at 3l/4
  t1 = pulse_params->length/4.0;
  t2 = 3*pulse_params->length/4.0;
  ts = pulse_params->stime;
  tau = 0.175*pulse_params->length;
  a = exp(-pow(t1,4)/pow(tau,4));
  if(time-ts>0 && time-ts<pulse_params->length/2){
    p = (exp(-pow((time-ts-t1),4)/pow(tau,4)) - a)/(1-a);
  } else if(time-ts>pulse_params->length/2 && time-ts<pulse_params->length){
    p = (exp(-pow((time-ts-t2),4)/pow(tau,4)) - a)/(1-a);
  } else{
    p = 0.0;
  }

  pulse_value = pulse_params->omega/2.0*p;

  return pulse_value;
}


PetscScalar delta(PetscReal time,void *ctx){
  PetscScalar pulse_value;
  PulseParams *pulse_params = (PulseParams*) ctx;   /* user-defined struct */
  PetscReal sigma,p,ts;
  
  ts = pulse_params->stime;
  
  //I don't know what the true form is, so I used two gaussians,
  //one centered at l/2, one centered at t, the first stops at t/2

  if(time-ts>0 && time-ts<pulse_params->length/2){
    p = sin(2*PETSC_PI*(time-ts + 3*pulse_params->length/4)/pulse_params->length);
  } else if(time-ts>pulse_params->length/2 && time-ts<pulse_params->length){
    p = sin(2*PETSC_PI*(time-ts + pulse_params->length/4)/pulse_params->length);
  } else {
    p=0.0;
  }
  pulse_value = pulse_params->delta * p;
  return pulse_value;
}


