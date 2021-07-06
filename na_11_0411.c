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
FILE *data_fp = NULL;
typedef struct {
  PetscScalar omega,delta;
  PetscReal length;
  qvec dm,dm_dummy;
  vec_op *atoms;
  qsystem qsys;
} PulseParams;

int main(int argc,char **args){
  PetscInt n_atoms,i,n_levels;
  PetscScalar tmp_scalar = 1.0,b_field,b_dr,b_0r,b_1r,gamma_r;
  PetscInt steps_max;
  PetscReal dt,time_max;
  PulseParams pulse_params;
  //State identifiers
  enum STATE {zero=0,d,one,r};

  /* Initialize QuaC */
  QuaC_initialize(argc,args);

  //Initialize the qsystem
  initialize_system(&pulse_params.qsys);

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
  pulse_params.atoms = malloc(n_atoms*sizeof(vec_op));

  for(i=0;i<n_atoms;i++){
    //Create an n_level system which is stored in atoms[i]
    create_vec_op_sys(pulse_params.qsys,n_levels,&(pulse_params.atoms[i]));
  }


  //Add hamiltonian terms
  pulse_params.omega = 17.; //MHz 2pi?
  pulse_params.length = 0.54; //us
  pulse_params.delta = 23.; //MHz 2pi?


  for(i=0;i<n_atoms;i++){
    tmp_scalar = 1.0;
    //1.0 * omega(t) * |r><1|
    add_ham_term_time_dep(pulse_params.qsys,tmp_scalar,&pulse_params,omega,2,pulse_params.atoms[i][r],pulse_params.atoms[i][one]);
    //1.0 * omega(t) * |1><r|
    add_ham_term_time_dep(pulse_params.qsys,tmp_scalar,&pulse_params,omega,2,pulse_params.atoms[i][one],pulse_params.atoms[i][r]);

    //1.0 * delta(t) * |r><r|
    add_ham_term_time_dep(pulse_params.qsys,tmp_scalar,&pulse_params,delta,2,pulse_params.atoms[i][r],pulse_params.atoms[i][r]);
  }

  //Coupling term
  //b_field * (|r_0> <r_0|) (|r_1><r_1|) = (|r_0 r_1><r_0 r_1|)
  add_ham_term(pulse_params.qsys,b_field,4,pulse_params.atoms[0][r],pulse_params.atoms[0][r],pulse_params.atoms[1][r],pulse_params.atoms[1][r]);

  //Add lindblad terms
  for(i=0;i<n_atoms;i++){
    tmp_scalar = b_0r*gamma_r;
    //tmp_scalar * L(|0><r|) - no sqrt needed because lin term wants squared term
    add_lin_term(pulse_params.qsys,tmp_scalar,2,pulse_params.atoms[i][zero],pulse_params.atoms[i][r]);

    tmp_scalar = b_1r*gamma_r;
    //L(|1><r|)
    add_lin_term(pulse_params.qsys,tmp_scalar,2,pulse_params.atoms[i][one],pulse_params.atoms[i][r]);

    tmp_scalar = b_dr*gamma_r;
    //L(|d><r|)
    add_lin_term(pulse_params.qsys,tmp_scalar,2,pulse_params.atoms[i][d],pulse_params.atoms[i][r]);

  }

  //Now that we've added all the terms, we construct the matrix
  construct_matrix(pulse_params.qsys);
  /* print_mat_sparse(qsys->mat_A); */
  //Create a density matrix object
  create_qvec_sys(pulse_params.qsys,&(pulse_params.dm));
  create_qvec_sys(pulse_params.qsys,&(pulse_params.dm_dummy));

  add_to_qvec(pulse_params.dm,1.0,10,10); //start in the |10><10| state
  assemble_qvec(pulse_params.dm);

  time_max  = 1;
  dt        = 0.05;
  steps_max = 1000;
  
  data_fp = fopen("na_11_0411.dat","w");

  /* Set the ts_monitor to print results at each time step */
  set_ts_monitor_sys(pulse_params.qsys,ts_monitor,&pulse_params);

  //Do the timestepping
  time_step_sys(pulse_params.qsys,pulse_params.dm,0.0,time_max,dt,steps_max);

  print_qvec(pulse_params.dm);
  //clean up memory
  destroy_system(&pulse_params.qsys);
  destroy_qvec(&(pulse_params.dm));
  destroy_qvec(&(pulse_params.dm_dummy));
  for(i=0;i<n_atoms;i++){
    destroy_vec_op_sys(&pulse_params.atoms[i]);
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
  tmp_data = pulse_params->dm_dummy->data;
  pulse_params->dm_dummy->data = rho_data;

  //ev( (|1><1|)(|0><0|) )= ev(|10><10|)
  get_expectation_value_qvec(pulse_params->dm_dummy,&trace_val,4,pulse_params->atoms[0][one],pulse_params->atoms[0][one],pulse_params->atoms[1][one],pulse_params->atoms[1][one]);
  get_expectation_value_qvec(pulse_params->dm_dummy,&trace_val2,4,pulse_params->atoms[0][one],pulse_params->atoms[0][one],pulse_params->atoms[1][r],pulse_params->atoms[1][r]);
  get_expectation_value_qvec(pulse_params->dm_dummy,&trace_val3,4,pulse_params->atoms[0][r],pulse_params->atoms[0][r],pulse_params->atoms[1][one],pulse_params->atoms[1][one]);
  get_expectation_value_qvec(pulse_params->dm_dummy,&trace_val4,4,pulse_params->atoms[0][r],pulse_params->atoms[0][r],pulse_params->atoms[1][r],pulse_params->atoms[1][r]);
  get_expectation_value_qvec(pulse_params->dm_dummy,&trace_val5,4,pulse_params->atoms[0][d],pulse_params->atoms[0][d],pulse_params->atoms[1][one],pulse_params->atoms[1][one]);
  PetscFPrintf(PETSC_COMM_WORLD,"%d %f %f %f %f %f %f %f %f\n",step,time,PetscRealPart(omega(time,ctx)),PetscRealPart(delta(time,ctx)),PetscRealPart(trace_val),PetscRealPart(trace_val2),PetscRealPart(trace_val3));
  pulse_params->dm_dummy->data = tmp_data;
  /* print_qvec(rho); */
  PetscFunctionReturn(0);
}


//Define time dependent pulses
PetscScalar omega(PetscReal time,void *ctx){
  PetscScalar pulse_value;
  PetscReal sigma,t1,t2,g1,g2;
  PulseParams *pulse_params = (PulseParams*) ctx;   /* user-defined struct */


  //I don't know what the true form is, so I used two gaussians,
  //one centered at l/4, one centered at 3l/4
  t1 = pulse_params->length/4.0;
  t2 = 3*pulse_params->length/4.0;

  //Sigma was approximated by eye
  sigma = pulse_params->length/(5.4*2*sqrt(2*log(2)));

  g1 = exp(-(time-t1)*(time-t1)/(2*sigma*sigma));
  g2 = exp(-(time-t2)*(time-t2)/(2*sigma*sigma));

  pulse_value = pulse_params->omega*(g1+g2);

  return pulse_value;
}


PetscScalar delta(PetscReal time,void *ctx){
  PetscScalar pulse_value;
  PulseParams *pulse_params = (PulseParams*) ctx;   /* user-defined struct */
  PetscReal sigma,g1,t1;
  //I don't know what the true form is, so I used two gaussians,
  //one centered at l/2, one centered at t, the first stops at t/2

  if(time>pulse_params->length){
    pulse_value = 0.0;
  } else {
    if(time<pulse_params->length/2){
      t1 = pulse_params->length/2;
      sigma = pulse_params->length/(2*2*sqrt(2*log(2)));
      g1 = exp(-(time-t1)*(time-t1)/(2*sigma*sigma));
    } else {
      t1 = pulse_params->length;
      sigma = pulse_params->length/(2*2*sqrt(2*log(2)));
      g1 = exp(-(time-t1)*(time-t1)/(2*sigma*sigma));
    }
    pulse_value = -pulse_params->delta + 2*pulse_params->delta*g1;
  }
  return pulse_value;
}
