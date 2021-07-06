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

PetscErrorCode ts_monitor_par(TS,PetscInt,PetscReal,Vec,void*);
PetscErrorCode ts_monitor_seq(TS,PetscInt,PetscReal,Vec,void*);
qvec dmpar_dummy,dmseq_dummy;
vec_op *atomspar,*atomsseq;
FILE *data_fppar = NULL;
FILE *data_fpseq = NULL;

typedef struct {
  PetscScalar omega,delta;
  PetscReal length,stime;
} PulseParams;

int main(int argc,char **args){
  PetscInt n_atoms,i,n_levels,n_seqgroups,max_seqgroupsize;
  PetscScalar tmp_scalar = 1.0,b_field,b_dr,b_0r,b_1r,gamma_r;
  PetscInt steps_max;
  qvec dmpar,dmseq;
  qsystem qsyspar,qsysseq;
  PetscReal dt,time_max,fidelity,var;
  //PulseParams pulse_params[2]; //moved to line 61
  //State identifiers
  enum STATE {zero=0,d,one,r};

  /* Initialize QuaC */
  QuaC_initialize(argc,args);

  //Initialize the qsystem
  initialize_system(&qsyspar);
  initialize_system(&qsysseq);
  
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
  printf("test 1\n");
  //magnetic field
  b_field = 3000*2*PETSC_PI; //2pi?

  //Create the operators for the atoms
  atomspar = malloc(n_atoms*sizeof(vec_op));
  atomsseq = malloc(n_atoms*sizeof(vec_op));
	  
  for(i=0;i<n_atoms;i++){
    //Create an n_level system which is stored in atoms[i]
    create_vec_op_sys(qsyspar,n_levels,&(atomspar[i]));
    create_vec_op_sys(qsysseq,n_levels,&(atomsseq[i]));
  }

  data_fppar = fopen("na_3_comp_par.dat","w");
  data_fpseq = fopen("na_3_comp_seq.dat","w"); 
  
  
  //Add hamiltonian terms  
  for(int i=0;i<n_seqgroups;i++){
  	for(int j=0;j<seqgroupsize[i];j++){
  		tmp_scalar = 1.0;
    	//1.0 * omega(t) * |r><1|
    	add_ham_term_time_dep(qsysseq,tmp_scalar,&pulse_params[i],omega,2,atomsseq[seqgroup[i][j]][r],atomsseq[seqgroup[i][j]][one]);
    	//1.0 * omega(t) * |1><r|
    	add_ham_term_time_dep(qsysseq,tmp_scalar,&pulse_params[i],omega,2,atomsseq[seqgroup[i][j]][one],atomsseq[seqgroup[i][j]][r]);
   		//1.0 * delta(t) * |r><r|
    	add_ham_term_time_dep(qsysseq,tmp_scalar,&pulse_params[i],delta,2,atomsseq[seqgroup[i][j]][r],atomsseq[seqgroup[i][j]][r]);
	  }
  }
  
  for(i=0;i<n_atoms;i++){
    tmp_scalar = 1.0;
    //1.0 * omega(t) * |r><1|
    add_ham_term_time_dep(qsyspar,tmp_scalar,&pulse_params[0],omega,2,atomspar[i][r],atomspar[i][one]);
    //1.0 * omega(t) * |1><r|
    add_ham_term_time_dep(qsyspar,tmp_scalar,&pulse_params[0],omega,2,atomspar[i][one],atomspar[i][r]);

    //1.0 * delta(t) * |r><r|
    add_ham_term_time_dep(qsyspar,tmp_scalar,&pulse_params[0],delta,2,atomspar[i][r],atomspar[i][r]);
  }
  
  //Coupling term
  //b_field * (|r_0> <r_0|) (|r_1><r_1|) = (|r_0 r_1><r_0 r_1|)
  add_ham_term(qsyspar,b_field*0,4,atomspar[0][r],atomspar[0][r],atomspar[1][r],atomspar[1][r]);
  add_ham_term(qsyspar,b_field,4,atomspar[0][r],atomspar[0][r],atomspar[2][r],atomspar[2][r]);
  add_ham_term(qsyspar,b_field,4,atomspar[1][r],atomspar[1][r],atomspar[2][r],atomspar[2][r]);
  
  add_ham_term(qsysseq,b_field*0,4,atomsseq[0][r],atomsseq[0][r],atomsseq[1][r],atomsseq[1][r]);
  add_ham_term(qsysseq,b_field,4,atomsseq[0][r],atomsseq[0][r],atomsseq[2][r],atomsseq[2][r]);
  add_ham_term(qsysseq,b_field,4,atomsseq[1][r],atomsseq[1][r],atomsseq[2][r],atomsseq[2][r]);
  
  //Add lindblad terms
  for(i=0;i<n_atoms;i++){
    tmp_scalar = b_0r*gamma_r;
    //tmp_scalar * L(|0><r|) - no sqrt needed because lin term wants squared term
    add_lin_term(qsyspar,tmp_scalar,2,atomspar[i][zero],atomspar[i][r]);
    add_lin_term(qsysseq,tmp_scalar,2,atomsseq[i][zero],atomsseq[i][r]);
    tmp_scalar = b_1r*gamma_r;
    //L(|1><r|)
    add_lin_term(qsyspar,tmp_scalar,2,atomspar[i][one],atomspar[i][r]);
    add_lin_term(qsysseq,tmp_scalar,2,atomsseq[i][one],atomsseq[i][r]);
    tmp_scalar = b_dr*gamma_r;
    //L(|d><r|)
    add_lin_term(qsyspar,tmp_scalar,2,atomspar[i][d],atomspar[i][r]);
    add_lin_term(qsysseq,tmp_scalar,2,atomsseq[i][d],atomsseq[i][r]);
  }

  //Now that we've added all the terms, we construct the matrix
  construct_matrix(qsyspar);
  construct_matrix(qsysseq);
  /* print_mat_sparse(qsys->mat_A); */
  //Create a density matrix object
  create_qvec_sys(qsyspar,&(dmpar));
  create_qvec_sys(qsyspar,&(dmpar_dummy));
  create_qvec_sys(qsysseq,&(dmseq));
  create_qvec_sys(qsysseq,&(dmseq_dummy));
/* 
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
*/
	add_to_qvec(dmpar,1.0,42,42); //start in the |111><11| state
	add_to_qvec(dmseq,1.0,42,42); //start in the |111><11| state
  assemble_qvec(dmpar);
  assemble_qvec(dmseq);

  time_max  = 1.08;
  dt        = 0.02;
  steps_max = 1000;

  /* Set the ts_monitor to print results at each time step */
  set_ts_monitor_sys(qsyspar,ts_monitor_par,&pulse_params);
  time_step_sys(qsyspar,dmpar,0.0,time_max,dt,steps_max);
  
  set_ts_monitor_sys(qsysseq,ts_monitor_seq,&pulse_params);
  time_step_sys(qsysseq,dmseq,0.0,time_max,dt,steps_max);
  //clean up memory
  
  var=1;
  get_fidelity_qvec(dmpar,dmseq,&fidelity,&var);
  printf("fidelity between par and seq is %lf\n",fidelity);
  
  destroy_system(&qsyspar);
  destroy_qvec(&(dmpar));
  destroy_qvec(&(dmpar_dummy));
  for(i=0;i<n_atoms;i++){
    destroy_vec_op_sys(&atomspar[i]);
  }
  destroy_system(&qsysseq);
  destroy_qvec(&(dmseq));
  destroy_qvec(&(dmseq_dummy));
  for(i=0;i<n_atoms;i++){
    destroy_vec_op_sys(&atomsseq[i]);
  }
  return;
}

PetscErrorCode ts_monitor_par(TS ts,PetscInt step,PetscReal time,Vec rho_data,void *ctx){
  PulseParams *pulse_params = (PulseParams*) ctx;   /* user-defined struct */
  PetscScalar trace_val,trace_val2,trace_val3,trace_val4,trace_val5,trace_val6,trace_val7,trace_val8;
  Vec tmp_data;
  enum STATE {zero=0,d,one,r};
  //Print out things at each time step, if desired
  //tmp data necessary for technical reasons
  tmp_data = dmpar_dummy->data;
  dmpar_dummy->data = rho_data;

  //ev( (|1><1|)(|0><0|) )= ev(|10><10|)
  //get_expectation_value_qvec_list
  get_expectation_value_qvec(dmpar_dummy,&trace_val,6,atomspar[0][one],atomspar[0][one],atomspar[1][one],atomspar[1][one],atomspar[2][one],atomspar[2][one]);
  get_expectation_value_qvec(dmpar_dummy,&trace_val2,6,atomspar[0][one],atomspar[0][one],atomspar[1][one],atomspar[1][one],atomspar[2][r],atomspar[2][r]);
  get_expectation_value_qvec(dmpar_dummy,&trace_val3,6,atomspar[0][one],atomspar[0][one],atomspar[1][r],atomspar[1][r],atomspar[2][one],atomspar[2][one]);
  get_expectation_value_qvec(dmpar_dummy,&trace_val4,6,atomspar[0][r],atomspar[0][r],atomspar[1][one],atomspar[1][one],atomspar[2][one],atomspar[2][one]);
  get_expectation_value_qvec(dmpar_dummy,&trace_val5,6,atomspar[0][one],atomspar[0][one],atomspar[1][r],atomspar[1][r],atomspar[2][r],atomspar[2][r]);
  get_expectation_value_qvec(dmpar_dummy,&trace_val6,6,atomspar[0][r],atomspar[0][r],atomspar[1][one],atomspar[1][one],atomspar[2][r],atomspar[2][r]);
  get_expectation_value_qvec(dmpar_dummy,&trace_val7,6,atomspar[0][r],atomspar[0][r],atomspar[1][r],atomspar[1][r],atomspar[2][one],atomspar[2][one]);
  get_expectation_value_qvec(dmpar_dummy,&trace_val8,6,atomspar[0][r],atomspar[0][r],atomspar[1][r],atomspar[1][r],atomspar[2][r],atomspar[2][r]);

  
  PetscFPrintf(PETSC_COMM_WORLD,data_fppar,"%d %f %f %f %f %f %f %f %f %f\n",step,time,PetscRealPart(trace_val),PetscRealPart(trace_val2),PetscRealPart(trace_val3),PetscRealPart(trace_val4),PetscRealPart(trace_val5),PetscRealPart(trace_val6),PetscRealPart(trace_val7),PetscRealPart(trace_val8));
  dmpar_dummy->data = tmp_data;
  /* print_qvec(rho); */
  PetscFunctionReturn(0);
}
PetscErrorCode ts_monitor_seq(TS ts,PetscInt step,PetscReal time,Vec rho_data,void *ctx){
  PulseParams *pulse_params = (PulseParams*) ctx;   /* user-defined struct */
  PetscScalar trace_val,trace_val2,trace_val3,trace_val4,trace_val5,trace_val6,trace_val7,trace_val8;
  Vec tmp_data;
  enum STATE {zero=0,d,one,r};
  //Print out things at each time step, if desired
  //tmp data necessary for technical reasons
  tmp_data = dmseq_dummy->data;
  dmseq_dummy->data = rho_data;

  //ev( (|1><1|)(|0><0|) )= ev(|10><10|)
  //get_expectation_value_qvec_list
  get_expectation_value_qvec(dmseq_dummy,&trace_val,6,atomsseq[0][one],atomsseq[0][one],atomsseq[1][one],atomsseq[1][one],atomsseq[2][one],atomsseq[2][one]);
  get_expectation_value_qvec(dmseq_dummy,&trace_val2,6,atomsseq[0][one],atomsseq[0][one],atomsseq[1][one],atomsseq[1][one],atomsseq[2][r],atomsseq[2][r]);
  get_expectation_value_qvec(dmseq_dummy,&trace_val3,6,atomsseq[0][one],atomsseq[0][one],atomsseq[1][r],atomsseq[1][r],atomsseq[2][one],atomsseq[2][one]);
  get_expectation_value_qvec(dmseq_dummy,&trace_val4,6,atomsseq[0][r],atomsseq[0][r],atomsseq[1][one],atomsseq[1][one],atomsseq[2][one],atomsseq[2][one]);
  get_expectation_value_qvec(dmseq_dummy,&trace_val5,6,atomsseq[0][one],atomsseq[0][one],atomsseq[1][r],atomsseq[1][r],atomsseq[2][r],atomsseq[2][r]);
  get_expectation_value_qvec(dmseq_dummy,&trace_val6,6,atomsseq[0][r],atomsseq[0][r],atomsseq[1][one],atomsseq[1][one],atomsseq[2][r],atomsseq[2][r]);
  get_expectation_value_qvec(dmseq_dummy,&trace_val7,6,atomsseq[0][r],atomsseq[0][r],atomsseq[1][r],atomsseq[1][r],atomsseq[2][one],atomsseq[2][one]);
  get_expectation_value_qvec(dmseq_dummy,&trace_val8,6,atomsseq[0][r],atomsseq[0][r],atomsseq[1][r],atomsseq[1][r],atomsseq[2][r],atomsseq[2][r]);

  
  PetscFPrintf(PETSC_COMM_WORLD,data_fpseq,"%d %f %f %f %f %f %f %f %f %f\n",step,time,PetscRealPart(trace_val),PetscRealPart(trace_val2),PetscRealPart(trace_val3),PetscRealPart(trace_val4),PetscRealPart(trace_val5),PetscRealPart(trace_val6),PetscRealPart(trace_val7),PetscRealPart(trace_val8));
  dmseq_dummy->data = tmp_data;
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


