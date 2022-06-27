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

qvec dm;

vec_op *atom;
FILE *data_fp = NULL;



int main(int argc,char *args[]){


  PetscInt n_atoms,i,n_levels;
  PetscScalar tmp_scalar = 1.0,b_field,b_dr,b_0r,b_1r,gamma_r,valpar,diagsum;
  PetscInt steps_max;
  qsystem qsys;
  PetscReal dt,time_max;
  PetscReal var,fidelity;
  //PulseParams pulse_params[2]; //moved to line 61
  //State identifiers
	
  //enum STATE {zero=0,d,one,r};

  /* Initialize QuaC */
  QuaC_initialize(argc,args);

  //Initialize the qsystem
  initialize_system(&qsys);

  n_atoms = 1;
  n_levels = 16;
  atom = malloc(n_atoms*sizeof(vec_op));
  

  create_vec_op_sys(qsys,n_levels,&(atoms[i]));

/*
  //data_fp = fopen("neutral_atom_3atom_seq111.dat","w");
  // PetscFPrintf(PETSC_COMM_WORLD,data_fp,"#Step_num time omega delta |10><10| |r0><r0| |d0><d0|\n");
  
  //Add hamiltonian terms

/* 
  for(int i=0;i<16;i++){
  	for(int j=0;j<16;j++){
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
  /*
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
  */
  construct_matrix(qsys);

  create_qvec_sys(qsys,&(dm));
  
	dmpos1=3;
	dmpos2=11;
	leakage=0.001
	initval=(1.0-leakage)/2.0;
	
  add_to_qvec(dm,initval,dmpos1,dmpos1); 
  add_to_qvec(dm,initval,dmpos2,dmpos2); 
  for(i=0;i<16;i++){
  	if((i!=4)&&(i!=12)):
  		add_to_qvec(dm,leakage/14.0,i,i); 
	}
  assemble_qvec(dm); 

  time_max  = 8;
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

  time_step_sys(qsys,dm,0.0,time_max,dt,steps_max);

	get_expectation_value_qvec_list(dm,&val4,2,atom[3],atom[3])
	get_expectation_value_qvec_list(dm,&val12,2,atom[11],atom[11])
	diagsum=val4+val12
  printf("sum of the diag is %f\n",diagsum);
  //print_qvec(dm);
  //clean up memory
  data_fp = fopen("s16.dat","a");
  
  PetscFPrintf(PETSC_COMM_WORLD,data_fp,"%f %f\n",leakage,diagsum);
  
  destroy_system(&qsys);
  destroy_qvec(&(dm));
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


