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

PetscErrorCode ts_monitor(TS,PetscInt,PetscReal,Vec,void*);
qvec dm_dummy;
vec_op *atoms;
FILE *data_fp = NULL;

int main(int argc,char **args){
  PetscInt n_atoms,i,n_levels;

  PetscInt steps_max;
  qvec dm;
  qsystem qsys;
  PetscReal dt,time_max;
  
  //State identifiers
  enum STATE {zero=0,d,one,r};

  /* Initialize QuaC */
  QuaC_initialize(argc,args);

  //Initialize the qsystem
  initialize_system(&qsys);

  n_atoms = 3;
  n_levels = 4;

  //Create the operators for the atoms
  atoms = malloc(n_atoms*sizeof(vec_op));

  for(i=0;i<n_atoms;i++){
    //Create an n_level system which is stored in atoms[i]
    create_vec_op_sys(qsys,n_levels,&(atoms[i]));
  }

  data_fp = fopen("neutral_atom_3atom.dat","w");
  PetscFPrintf(PETSC_COMM_WORLD,data_fp,"#Step_num time omega delta |10><10| |r0><r0| |d0><d0|\n");

  //Now that we've added all the terms, we construct the matrix
  construct_matrix(qsys);

  /* print_mat_sparse(qsys->mat_A); */
  //Create a density matrix object
  create_qvec_sys(qsys,&(dm));
  create_qvec_sys(qsys,&(dm_dummy));

  add_to_qvec(dm,1.0,2,2); //start in the |100><100| state
  assemble_qvec(dm);

  /* Set the ts_monitor to print results at each time step */
  set_ts_monitor_sys(qsys,ts_monitor,&pulse_params);
  
  create_circuit(&circ,3);
  //Add some gates
  add_gate_to_circuit_sys(&circ,0.3,CNOT,2,0);
  add_gate_to_circuit_sys(&circ,0.6,CNOT,2,1);

  time_max  = 1;
  dt        = 0.05;
  steps_max = 1000;
	 
  time_max = 0;
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
  schedule_circuit_layers(system,&circ);

  //Start out circuit at time 0.0, first gate will be at 0
  apply_circuit_to_sys(system,&circ,0.0);

  //Do the timestepping
  time_step_sys(qsys,dm,0.0,time_max,dt,steps_max);
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
  PetscScalar trace_val,trace_val2,trace_val3,trace_val4,trace_val5;
  Vec tmp_data;
  enum STATE {zero=0,d,one,r};
  //Print out things at each time step, if desired
  //tmp data necessary for technical reasons
  tmp_data = dm_dummy->data;
  dm_dummy->data = rho_data;

  //ev( (|1><1|)(|0><0|) )= ev(|10><10|)
  //get_expectation_value_qvec_list
  get_expectation_value_qvec(dm_dummy,&trace_val,6,atoms[0][zero],atoms[0][zero],atoms[1][zero],atoms[1][zero],atoms[2][one],atoms[2][one]);
  get_expectation_value_qvec(dm_dummy,&trace_val2,6,atoms[0][zero],atoms[0][zero],atoms[1][zero],atoms[1][zero],atoms[2][r],atoms[2][r]);
  get_expectation_value_qvec(dm_dummy,&trace_val3,6,atoms[0][zero],atoms[0][zero],atoms[1][zero],atoms[1][zero],atoms[2][d],atoms[2][d]);
  PetscFPrintf(PETSC_COMM_WORLD,data_fp,"%d %f %f %f %f %f %f\n",step,time,PetscRealPart(omega(time,ctx)),PetscRealPart(delta(time,ctx)),PetscRealPart(trace_val),PetscRealPart(trace_val2),PetscRealPart(trace_val3));
  dm_dummy->data = tmp_data;
  /* print_qvec(rho); */
  PetscFunctionReturn(0);
}



