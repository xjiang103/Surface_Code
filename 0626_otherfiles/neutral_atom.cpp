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

typedef struct {
  PetscScalar omega,delta;
  PetscReal length;
} PulseParams;

int main(int argc,char **args){
  PetscInt n_atoms,i,n_levels;
  PetscScalar tmp_scalar = 1.0,b_field,b_dr,b_0r,b_1r,gamma_r;
  PetscInt steps_max;
  qvec dm;
  qsystem qsys;
  PetscReal dt,time_max;
  PulseParams pulse_params;
  //State identifiers
  enum STATE {zero=0,d,one,r};

  /* Initialize QuaC */
  QuaC_initialize(argc,args);

  //Initialize the qsystem
  initialize_system(&qsys);

  n_atoms = 2;
  n_levels = 4;

  //branching ratios
  //Careful to use 1.0 instead of 1 because of integer division
  b_1r = 1.0/16.0;
  b_0r = 1.0/16.0;
  b_dr = 7.0/8.0;

  //decay rate
  gamma_r = 1.0/(540.0);

  //magnetic field
  b_field = 100; //2pi?

  //Create the operators for the atoms
  atoms = malloc(n_atoms*sizeof(vec_op));

  for(i=0;i<n_atoms;i++){
    //Create an n_level system which is stored in atoms[i]
    create_vec_op_sys(qsys,n_levels,&(atoms[i]));
  }

  data_fp = fopen("neutral_atom.dat","w");
  PetscFPrintf(PETSC_COMM_WORLD,data_fp,"#Step_num time omega delta |10><10| |r0><r0| |d0><d0|\n");
  //Add hamiltonian terms
  pulse_params.omega = 17.; //MHz 2pi?
  pulse_params.length = 0.54; //us
  pulse_params.delta = 23.; //MHz 2pi?


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
  /* print_mat_sparse(qsys->mat_A); */
  //Create a density matrix object
  create_qvec_sys(qsys,&(dm));
  create_qvec_sys(qsys,&(dm_dummy));

  add_to_qvec(dm,1.0,2,2); //start in the |01><01| state
  assemble_qvec(dm);

  time_max  = 1;
  dt        = 0.05;
  steps_max = 1000;

  /* Set the ts_monitor to print results at each time step */
  set_ts_monitor_sys(qsys,ts_monitor,&pulse_params);

  //Do the timestepping
  time_step_sys(qsys,dm,0.0,time_max,dt,steps_max);

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
  PetscScalar trace_val,trace_val2,trace_val3;
  Vec tmp_data;
  enum STATE {zero=0,d,one,r};
  //Print out things at each time step, if desired
  //tmp data necessary for technical reasons
  tmp_data = dm_dummy->data;
  dm_dummy->data = rho_data;

  //ev( (|1><1|)(|0><0|) )= ev(|10><10|)
  get_expectation_value_qvec(dm_dummy,&trace_val,4,atoms[0][one],atoms[0][one],atoms[1][zero],atoms[1][zero]);
  get_expectation_value_qvec(dm_dummy,&trace_val2,4,atoms[0][r],atoms[0][r],atoms[1][zero],atoms[1][zero]);
  get_expectation_value_qvec(dm_dummy,&trace_val3,4,atoms[0][d],atoms[0][d],atoms[1][zero],atoms[1][zero]);

  get_expectation_value_qvec(dm_dummy,&trace_val,4,atoms[1][one],atoms[1][one],atoms[0][zero],atoms[0][zero]);
  get_expectation_value_qvec(dm_dummy,&trace_val2,4,atoms[1][r],atoms[1][r],atoms[0][zero],atoms[0][zero]);
  get_expectation_value_qvec(dm_dummy,&trace_val3,4,atoms[1][d],atoms[1][d],atoms[0][zero],atoms[0][zero]);
  if(step>0){
     PetscFPrintf(PETSC_COMM_WORLD,data_fp,"%d %f %f %f %f %f %f\n",step,time,PetscRealPart(omega(time,ctx)),PetscRealPart(delta(time,ctx)),PetscRealPart(trace_val),PetscRealPart(trace_val2),PetscRealPart(trace_val3));
   }
  dm_dummy->data = tmp_data;
  /* print_qvec(rho); */
  PetscFunctionReturn(0);
}


//Define time dependent pulses
PetscScalar omega(PetscReal time,void *ctx){
  PetscScalar pulse_value;
  PetscReal tau,t1,t2,p,a;
  PulseParams *pulse_params = (PulseParams*) ctx;   /* user-defined struct */


  //I don't know what the true form is, so I used two gaussians,
  //one centered at l/4, one centered at 3l/4
  t1 = pulse_params->length/4.0;
  t2 = 3*pulse_params->length/4.0;
  tau = 0.175*pulse_params->length;
  a = exp(-pow(t1,4)/pow(tau,4));
  if(time>0 && time<pulse_params->length/2){
    p = (exp(-pow((time-t1),4)/pow(tau,4)) - a)/(1-a);
  } else if(time>pulse_params->length/2 && time<pulse_params->length){
    p = (exp(-pow((time-t2),4)/pow(tau,4)) - a)/(1-a);
  } else{
    p = 0.0;
  }

  pulse_value = pulse_params->omega/2.0*p;

  return pulse_value;
}


PetscScalar delta(PetscReal time,void *ctx){
  PetscScalar pulse_value;
  PulseParams *pulse_params = (PulseParams*) ctx;   /* user-defined struct */
  PetscReal sigma,p;
  //I don't know what the true form is, so I used two gaussians,
  //one centered at l/2, one centered at t, the first stops at t/2

  if(time>0 && time<pulse_params->length/2){
    p = sin(2*PETSC_PI*(time + 3*pulse_params->length/4)/pulse_params->length);
  } else if(time>pulse_params->length/2 && time<pulse_params->length){
    p = sin(2*PETSC_PI*(time + pulse_params->length/4)/pulse_params->length);
  } else {
    p=0.0;
  }
  pulse_value = pulse_params->delta * p;
  return pulse_value;
}
