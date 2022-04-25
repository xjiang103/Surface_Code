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
#include "neutral_atom.h"

PetscErrorCode ts_monitor(TS,PetscInt,PetscReal,Vec,void*);
qvec dm_dummy,dm32;
operator op_list[10];
vec_op *atoms;
operator *atomsstd;
FILE *data_fp = NULL;

int main(int argc,char *args[]){
  PetscInt n_atoms,i,j,n_levels,n_seqgroups,max_seqgroupsize,val_init,pos1,pos2,dmpos,dmstdpos;
  PetscScalar tmp_scalar = 1.0,b_term=600,b_dr,b_0r,b_1r,gamma_r,valpar,diagsum,meas_val,gamma_t1,gamma_t2s;
  PetscInt steps_max,n_ens=0,seed=12;
  qvec dm,dmstd;
  qsystem qsys,qsysstd;
  PetscReal dt,time_max,measurement_time,phase_qb0=0,phase_qb1=0,phase_qb2=0,phase_qb3=0,phase_qb4=0,dd_fac=0;
  PetscReal single_qubit_gate_time=0.1,two_qubit_gate_time=0.1,var,fidelity;
  char bitstr[PETSC_MAX_PATH_LEN] = "11111"; //Default bitstr to start with
  char pulse_type[PETSC_MAX_PATH_LEN] = "ARP",filename[PETSC_MAX_PATH_LEN]="dm.dat";
  circuit circ;
  int length;
	PetscScalar omega,delta;
  PetscReal pulse_length,deltat;
  //State identifiers

  enum STATE {zero=0,one,r};

  /* Initialize QuaC */
  QuaC_initialize(argc,args);
  //rydberg coupling
  b_term = 600;
  n_atoms = 5;

  //Get the bitstring we want to simulate
  PetscOptionsGetString(NULL,NULL,"-bitstr",bitstr,PETSC_MAX_PATH_LEN,NULL);
  PetscOptionsGetInt(NULL,NULL,"-n_ens",&n_ens,NULL);
  PetscOptionsGetInt(NULL,NULL,"-seed",&seed,NULL);
  PetscOptionsGetReal(NULL,NULL,"-b_term",&b_term,NULL);
  PetscOptionsGetString(NULL,NULL,"-pulse_type",pulse_type,PETSC_MAX_PATH_LEN,NULL);
  PetscOptionsGetString(NULL,NULL,"-file",filename,PETSC_MAX_PATH_LEN,NULL);

  //Set default parameters for the two pulse_types
  if(strcmp("ARP",pulse_type)==0){
    pulse_length = 0.54;
    omega = 17.0;
    delta = 23.0;
  } else if(strcmp("SP",pulse_type)==0){
    pulse_length = 2.0;
    omega = 17.0;
    delta = -0.50;
    deltat = 0.2;
  } else {
    PetscPrintf(PETSC_COMM_WORLD,"Pulse type not understood!\n");
  }

  PetscOptionsGetReal(NULL,NULL,"-omega",&omega,NULL);
  PetscOptionsGetReal(NULL,NULL,"-delta",&delta,NULL);
  PetscOptionsGetReal(NULL,NULL,"-deltat",&deltat,NULL);
  PetscOptionsGetReal(NULL,NULL,"-length",&length,NULL);
  PetscOptionsGetReal(NULL,NULL,"-pulse_length",&pulse_length,NULL);
  PetscOptionsGetReal(NULL,NULL,"-phase_qb0",&phase_qb0,NULL);
  PetscOptionsGetReal(NULL,NULL,"-phase_qb1",&phase_qb1,NULL);
  PetscOptionsGetReal(NULL,NULL,"-phase_qb2",&phase_qb2,NULL);
  PetscOptionsGetReal(NULL,NULL,"-phase_qb3",&phase_qb3,NULL);
  PetscOptionsGetReal(NULL,NULL,"-phase_qb4",&phase_qb4,NULL);
  PetscOptionsGetReal(NULL,NULL,"-dd_fac",&dd_fac,NULL);




  b_term = b_term*2*PETSC_PI;
  dmpos = 0;
  dmstdpos = 0;
  length = strlen(bitstr);
  if(length!=5){
    PetscPrintf(PETSC_COMM_WORLD,"ERROR: bitstr must be of length 5!\n");
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

  measurement_time = 1;//1000us = 1ms
  n_levels = 3;

  //Parameters for doing sequential operation

  //number of sequential operations in time domain
  n_seqgroups = 1;
  //number of qubits in each sequential operation
  PetscInt seqgroupsize[1] = {5};
  //maximum element in n_seqgroupsize
  max_seqgroupsize=5;
  //indices of atoms in each sequential operation
  PetscInt seqgroup[n_seqgroups][max_seqgroupsize];
  seqgroup[0][0]=0, seqgroup[0][1]=1, seqgroup[0][2]=2, seqgroup[0][3]=3, seqgroup[0][4]=4;
  //parameters for each pulse in the sequence
  PulseParams pulse_params[n_seqgroups];

  time_max = 0;
  for(i=0;i<n_seqgroups;i++){

    if(strcmp("ARP",pulse_type)==0){

      pulse_params[i].stime = i*pulse_length;
      pulse_params[i].omega = omega*2*PETSC_PI; //MHz 2pi?
      pulse_params[i].length = pulse_length; //us
      pulse_params[i].delta = delta*2*PETSC_PI; //MHz 2pi?
      time_max = time_max + pulse_length;

    } else if(strcmp("SP",pulse_type)==0){

      pulse_params[i].stime = i*pulse_length;
      pulse_params[i].omega = omega*(2*PETSC_PI); //MHz 2pi?
      pulse_params[i].deltat = deltat; //us
      pulse_params[i].delta = delta*omega*(2*PETSC_PI); //MHz 2pi?
      time_max = time_max + pulse_length;

    }
  }

  //branching ratios
  //Careful to use 1.0 instead of 1 because of integer division
  b_1r = 1.0/16.0;
  b_0r = 1.0/16.0;
  b_dr = 7.0/8.0;

  //decay rate
  gamma_r = 1.0/(540.0);//us
  gamma_t1 = 1.0/(10000000);//T1 = 10s
  gamma_t2s = 1.0/(1000000);//T2s = 1s


  //Create the operators for the atoms
  atoms = malloc(n_atoms*sizeof(vec_op));
  atomsstd = malloc(n_atoms*sizeof(operator));

  for(i=0;i<n_atoms;i++){
    //Create an n_level system which is stored in atoms[i]
    create_vec_op_sys(qsys,n_levels,&(atoms[i]));
  }

  for(i=0;i<n_atoms;i++){
    create_op_sys(qsysstd,2,&(atomsstd[i]));
  }
  data_fp = fopen("neutral_atom_3atom_seq111.dat","w");
  PetscFPrintf(PETSC_COMM_WORLD,data_fp,"#Step_num time omega delta |10><10| |r0><r0| |d0><d0|\n");
  
  //Add hamiltonian terms


  for(int i=0;i<n_seqgroups;i++){
  	for(int j=0;j<seqgroupsize[i];j++){

      if(strcmp("ARP",pulse_type)==0){
        tmp_scalar = 1.0;
        //1.0 * omega(t) * |r><1|
        add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params[i],omega_arp,2,atoms[seqgroup[i][j]][r],atoms[seqgroup[i][j]][one]);
        //1.0 * omega(t) * |1><r|
        add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params[i],omega_arp,2,atoms[seqgroup[i][j]][one],atoms[seqgroup[i][j]][r]);

        //1.0 * delta(t) * |r><r|
        //delta is time dependent in this case
        add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params[i],delta_arp,2,atoms[seqgroup[i][j]][r],atoms[seqgroup[i][j]][r]);

      } else if(strcmp("SP",pulse_type)==0){
        tmp_scalar = 1.0;
        //1.0 * omega(t) * |r><1|
        add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params[i],omega_sp,2,atoms[seqgroup[i][j]][r],atoms[seqgroup[i][j]][one]);
        //1.0 * omega(t) * |1><r|
        add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params[i],omega_sp,2,atoms[seqgroup[i][j]][one],atoms[seqgroup[i][j]][r]);

        tmp_scalar = pulse_params[i].delta;
        if(seqgroup[i][j]!=0){
          //MJO: Delta has no time dependence, so we use a normal ham term
          add_ham_term(qsys,tmp_scalar,2,atoms[seqgroup[i][j]][r],atoms[seqgroup[i][j]][r]);
        }
      }
	  }
  }
  if(strcmp("SP",pulse_type)==0){
    tmp_scalar = pulse_params[0].delta;
    //MJO: Delta has no time dependence, so we use a normal ham term
    add_ham_term(qsys,tmp_scalar,2,atoms[0][r],atoms[0][r]);
  }
  /* } */
  //Coupling term
  //b * (|r_0> <r_0|) (|r_1><r_1|) = (|r_0 r_1><r_0 r_1|)
  add_ham_term(qsys,b_term,4,atoms[0][r],atoms[0][r],atoms[1][r],atoms[1][r]);
  add_ham_term(qsys,b_term,4,atoms[0][r],atoms[0][r],atoms[2][r],atoms[2][r]);
  add_ham_term(qsys,b_term,4,atoms[0][r],atoms[0][r],atoms[3][r],atoms[3][r]);
  add_ham_term(qsys,b_term,4,atoms[0][r],atoms[0][r],atoms[4][r],atoms[4][r]);  
  add_ham_term(qsys,dd_fac*b_term,4,atoms[1][r],atoms[1][r],atoms[2][r],atoms[2][r]);
  add_ham_term(qsys,dd_fac*b_term,4,atoms[1][r],atoms[1][r],atoms[3][r],atoms[3][r]);
  add_ham_term(qsys,dd_fac*b_term,4,atoms[1][r],atoms[1][r],atoms[4][r],atoms[4][r]);
  add_ham_term(qsys,dd_fac*b_term,4,atoms[2][r],atoms[2][r],atoms[3][r],atoms[3][r]);
  add_ham_term(qsys,dd_fac*b_term,4,atoms[2][r],atoms[2][r],atoms[4][r],atoms[4][r]);
  add_ham_term(qsys,dd_fac*b_term,4,atoms[3][r],atoms[3][r],atoms[4][r],atoms[4][r]);
 

  //Add lindblad terms
  if(n_ens!=0){
    for(i=0;i<n_atoms;i++){
      tmp_scalar = b_0r*gamma_r;
      //tmp_scalar * L(|0><r|) - no sqrt needed because lin term wants squared term
      add_lin_term(qsys,tmp_scalar,2,atoms[i][zero],atoms[i][r]);

      tmp_scalar = b_1r*gamma_r;
      //L(|1><r|)
      add_lin_term(qsys,tmp_scalar,2,atoms[i][one],atoms[i][r]);

      /* T1 time - L(|0><1|) */
      tmp_scalar = gamma_t1;
      add_lin_term(qsys,tmp_scalar,2,atoms[i][zero],atoms[i][one]);

      /* T2s time - L(|1><1|) */
      tmp_scalar = gamma_t2s;
      add_lin_term(qsys,tmp_scalar,2,atoms[i][one],atoms[i][one]);

      tmp_scalar = b_dr*gamma_r;
      //L(|d><r|)
      /* add_lin_term(qsys,tmp_scalar,2,atoms[i][d],atoms[i][r]); */
      // Instead of a lindblad term, we do an imaginary hamiltonian to allow for leakage
      // need to have a different tmp_scalar based on the repumping simulations
      if(n_ens==-1){
        //density matrix - two is needed
        tmp_scalar = -1*b_dr*gamma_r*PETSC_i/2;
      } else {
        //wf ensemble, do not need 2 !! Check this
        tmp_scalar = -1*b_dr*gamma_r*PETSC_i;
      }
      add_ham_term(qsys,tmp_scalar,2,atoms[i][r],atoms[i][r]);

      //Add terms for t1 and t2 time

    }
  }
  //Now that we've added all the terms, we construct the matrix
  if(n_ens>0){
    use_mcwf_solver(qsys,n_ens,seed);
  }
  printf("test1\n");
  construct_matrix(qsys);
  printf("test2\n");
  add_ham_term(qsysstd,1.0,1,atomsstd[0]->n);
  add_lin_term(qsysstd,1.0,1,atomsstd[0]->n);
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

  apply_1q_na_gate_to_qvec(dm,HADAMARD,atoms[0][zero]);
  apply_1q_na_gate_to_qvec(dm,HADAMARD,atoms[1][zero]);
  apply_1q_na_gate_to_qvec(dm,HADAMARD,atoms[2][zero]);
  apply_1q_na_gate_to_qvec(dm,HADAMARD,atoms[3][zero]);
  apply_1q_na_gate_to_qvec(dm,HADAMARD,atoms[4][zero]);

  //time_max  = 5;
  time_max = time_max + measurement_time; //100us
  dt        = 0.01;
  steps_max = 10000000;

  /* Set the ts_monitor to print results at each time step */
  /* set_ts_monitor_sys(qsys,ts_monitor,&pulse_params); */


  //Do the timestepping
  PetscPrintf(PETSC_COMM_WORLD,"---------------------\n");

  time_step_sys(qsys,dm,0.0,time_max,dt,steps_max);
  PetscPrintf(PETSC_COMM_WORLD,"Timestep 1 done\n");

  op_list[0] = atoms[0][zero]->sig_z;
  op_list[1] = atoms[1][zero]->sig_z;
  op_list[2] = atoms[2][zero]->sig_z;
  op_list[3] = atoms[3][zero]->sig_z;
  op_list[4] = atoms[4][zero]->sig_z;

  //  apply_projective_measurement_tensor_list(dm,&meas_val,5,op_list);
  //I don't think we need to split it up anymore
  PetscPrintf(PETSC_COMM_WORLD,"Full quantum state: \n");
  //  print_qvec(dm);

  PetscPrintf(PETSC_COMM_WORLD,"logical state populations: \n");

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
              					pos1=l_00*16+l_10*8+l_20*4+l_30*2+l_40;
												pos2=l_01*16+l_11*8+l_21*4+l_31*2+l_41;
												if (pos1==pos2){
													PetscPrintf(PETSC_COMM_WORLD,"%d %d %f\n",pos1,pos2,valpar);
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
  PetscPrintf(PETSC_COMM_WORLD,"dm32 constructed\n");
  print_qvec_file(dm32,filename);
  /* print_qvec(dm32); */
//----------------------------------------------------------------------------------------------
  create_circuit(&circ,25);
//  //Add some gates
//  add_gate_to_circuit_sys(&circ,0.09,HADAMARD,0);
//  add_gate_to_circuit_sys(&circ,0.1,HADAMARD,1);
//  add_gate_to_circuit_sys(&circ,0.15,HADAMARD,2);
//  add_gate_to_circuit_sys(&circ,0.16,HADAMARD,3);
//  add_gate_to_circuit_sys(&circ,0.17,HADAMARD,4);
//  add_gate_to_circuit_sys(&circ,0.16,RZ,0,phase_qb0+PETSC_PI);
//  add_gate_to_circuit_sys(&circ,0.17,RZ,1,phase_qb1+PETSC_PI);
//  add_gate_to_circuit_sys(&circ,0.2,CZ_ARP,0,1);
//  add_gate_to_circuit_sys(&circ,0.21,RZ,0,phase_qb0+PETSC_PI);
//  add_gate_to_circuit_sys(&circ,0.22,RZ,2,phase_qb1+PETSC_PI);
//  add_gate_to_circuit_sys(&circ,0.5,CZ_ARP,0,2);
//  add_gate_to_circuit_sys(&circ,0.51,RZ,0,phase_qb0+PETSC_PI);
//  add_gate_to_circuit_sys(&circ,0.52,RZ,3,phase_qb1+PETSC_PI);
//  add_gate_to_circuit_sys(&circ,1.0,CZ_ARP,0,3);
//  add_gate_to_circuit_sys(&circ,1.01,RZ,0,phase_qb0+PETSC_PI);
//  add_gate_to_circuit_sys(&circ,1.02,RZ,4,phase_qb1+PETSC_PI);
//  add_gate_to_circuit_sys(&circ,1.5,CZ_ARP,0,4);
//  single_qubit_gate_time = 0.1;
//  two_qubit_gate_time = 0.2;
//  for (i=0;i<circ.num_gates;i++){
//    if (circ.gate_list[i].num_qubits==1){
//      circ.gate_list[i].run_time = single_qubit_gate_time;
//      time_max += single_qubit_gate_time;
//    } else if (circ.gate_list[i].num_qubits==2){
//      circ.gate_list[i].run_time = two_qubit_gate_time;
//      time_max += two_qubit_gate_time;
//    }
//  }
//  schedule_circuit_layers(qsysstd,&circ);
//  //Start out circuit at time 0.0, first gate will be at 0
//  apply_circuit_to_qvec(qsysstd,circ,dmstd);

//-----------------------------------------------------------------------------------------------------------
  get_fidelity_qvec(dm32,dmstd,&fidelity,&var);
  PetscPrintf(PETSC_COMM_WORLD,"fidelity between seq(32*32) and std is %.20e\n",fidelity);
  PetscPrintf(PETSC_COMM_WORLD,"sum of the diag is %f\n",diagsum);
  PetscPrintf(PETSC_COMM_WORLD,"%.20e\n",fidelity);
  //print_qvec(dm);
  //clean up memory
  /* data_fp = fopen("seq_sp_3lvl.dat","a"); */

  /* for(i=0;i<32;i++){ */
  /*   get_dm_element_qvec(dmstd,i,i,&tmp_scalar); */
  /*   printf("%d %d %f %f\n",i,i,tmp_scalar); */
  /* } */
  /* PetscFPrintf(PETSC_COMM_WORLD,data_fp,"%d %d %f %f\n",dmpos,dmstdpos,fidelity,diagsum); */
  /* fclose(data_fp); */
  /* printf("dm32\n"); */
  /* print_qvec(dm32); */
  /* printf("dmstd\n"); */
  /* print_qvec(dmstd); */

  destroy_qvec(&(dm));
  destroy_qvec(&(dm_dummy));
  destroy_qvec(&(dmstd));
  destroy_qvec(&(dm32));

  for(i=0;i<n_atoms;i++){
    destroy_vec_op_sys(&atoms[i]);
    if(i<5){
      destroy_op_sys(&atomsstd[i]);
    }
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
  get_expectation_value_qvec(dm_dummy,&trace_val,6,atoms[2][one],atoms[2][one],atoms[3][one],atoms[3][one],atoms[4][one],atoms[4][one]);
  get_expectation_value_qvec(dm_dummy,&trace_val2,6,atoms[0][one],atoms[0][one],atoms[1][one],atoms[1][one],atoms[2][r],atoms[2][r]);
  get_expectation_value_qvec(dm_dummy,&trace_val3,6,atoms[0][one],atoms[0][one],atoms[1][r],atoms[1][r],atoms[2][one],atoms[2][one]);
  get_expectation_value_qvec(dm_dummy,&trace_val4,6,atoms[0][r],atoms[0][r],atoms[1][one],atoms[1][one],atoms[2][one],atoms[2][one]);
  get_expectation_value_qvec(dm_dummy,&trace_val5,6,atoms[0][one],atoms[0][one],atoms[1][r],atoms[1][r],atoms[2][r],atoms[2][r]);
  get_expectation_value_qvec(dm_dummy,&trace_val6,6,atoms[0][r],atoms[0][r],atoms[1][one],atoms[1][one],atoms[2][r],atoms[2][r]);
  get_expectation_value_qvec(dm_dummy,&trace_val7,6,atoms[0][r],atoms[0][r],atoms[1][r],atoms[1][r],atoms[2][one],atoms[2][one]);
  get_expectation_value_qvec(dm_dummy,&trace_val8,6,atoms[0][r],atoms[0][r],atoms[1][r],atoms[1][r],atoms[2][r],atoms[2][r]);


  PetscFPrintf(PETSC_COMM_WORLD,data_fp,"%d %f %f %f %f %f %f %f %f %f\n",step,time,PetscRealPart(trace_val),PetscRealPart(trace_val2),PetscRealPart(trace_val3),PetscRealPart(trace_val4),PetscRealPart(trace_val5),PetscRealPart(trace_val6),PetscRealPart(trace_val7),PetscRealPart(trace_val8));

  /* print_qvec(rho); */

  dm_dummy->data = tmp_data;
  PetscFunctionReturn(0);
}
