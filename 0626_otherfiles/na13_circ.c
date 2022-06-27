#include "neutral_atom.h"

operator *atoms;
FILE *data_fp = NULL;

int main(int argc,char *args[]){
  QuaC_initialize(argc,args);
  PetscInt n_atoms;
  circuit circ;
  qvec dm;
  qsystem qsys;
  initialize_system(&qsys);
  n_atoms = 13;
  atoms = malloc(n_atoms*sizeof(operator));
  
  for(int i=0;i<n_atoms;i++){
    create_op_sys(qsys,2,&(atoms[i]));
  }
  
  add_ham_term(qsys,1.0,1,atoms[0]->n);
  
  construct_matrix(qsys);
  
  create_wf_sys(qsys,&(dm));
  
  //print_qvec(dmstd);
  
  add_to_qvec(dm,1.0,0,0); 
  
  assemble_qvec(dm);

  create_circuit(&circ,52);
   
  //initialize 9,10,11,12 to |0>
  add_gate_to_circuit_sys(&circ,0.1,HADAMARD,10);
  add_gate_to_circuit_sys(&circ,0.1,HADAMARD,11);
  
  add_gate_to_circuit_sys(&circ,0.2,CNOT,0,9);
  add_gate_to_circuit_sys(&circ,0.2,CNOT,10,1);
  add_gate_to_circuit_sys(&circ,0.2,CNOT,11,6);
  add_gate_to_circuit_sys(&circ,0.2,CNOT,5,12);
  
  add_gate_to_circuit_sys(&circ,0.3,CNOT,3,9);
  add_gate_to_circuit_sys(&circ,0.3,CNOT,10,2);
  add_gate_to_circuit_sys(&circ,0.3,CNOT,11,7);
  add_gate_to_circuit_sys(&circ,0.3,CNOT,8,12); 
  
  add_gate_to_circuit_sys(&circ,0.4,HADAMARD,10);
  add_gate_to_circuit_sys(&circ,0.4,HADAMARD,11);
  
  //measurement on 9,10,11,12
  //initialize 9,10,11,12 to |0>
  
  add_gate_to_circuit_sys(&circ,0.6,HADAMARD,9);
  
  add_gate_to_circuit_sys(&circ,0.7,CNOT,9,0);
  add_gate_to_circuit_sys(&circ,0.7,CNOT,1,10);  
  
  add_gate_to_circuit_sys(&circ,0.8,CNOT,9,3);
  add_gate_to_circuit_sys(&circ,0.8,CNOT,4,10);  
  
  add_gate_to_circuit_sys(&circ,0.9,CNOT,9,1);
  add_gate_to_circuit_sys(&circ,0.9,CNOT,2,10);  
  
  add_gate_to_circuit_sys(&circ,1.0,CNOT,9,4);
  add_gate_to_circuit_sys(&circ,1.0,CNOT,5,10); 
  
  add_gate_to_circuit_sys(&circ,1.1,HADAMARD,9);   

  //measurement on 9,10
        
  add_gate_to_circuit_sys(&circ,1.3,HADAMARD,12);
  
  add_gate_to_circuit_sys(&circ,1.4,CNOT,4,11);
  add_gate_to_circuit_sys(&circ,1.4,CNOT,5,12); 
  
  add_gate_to_circuit_sys(&circ,1.5,CNOT,7,11);
  add_gate_to_circuit_sys(&circ,1.5,CNOT,8,12); 
  
  add_gate_to_circuit_sys(&circ,1.6,CNOT,3,11);
  add_gate_to_circuit_sys(&circ,1.6,CNOT,4,12); 
  
  add_gate_to_circuit_sys(&circ,1.7,CNOT,6,11);
  add_gate_to_circuit_sys(&circ,1.7,CNOT,7,12); 

  add_gate_to_circuit_sys(&circ,1.8,HADAMARD,12);

  schedule_circuit_layers(qsys,&circ);
  printf("test1\n");
  apply_circuit_to_qvec(qsys,circ,dm);
  
  printf("test2\n");  
  destroy_qvec(&(dm));

  for(int i=0;i<n_atoms;i++){
    destroy_op_sys(&atoms[i]);
    
  }
  destroy_system(&qsys);

  QuaC_finalize();
  return;
}
