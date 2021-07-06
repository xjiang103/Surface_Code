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
void fillSeqGroup(const PetscInt n_seqgroups,const PetscInt max_seqgroupsize,PetscInt seqgroup[n_seqgroups][max_seqgroupsize]);
PetscErrorCode ts_monitor(TS,PetscInt,PetscReal,Vec,void*);
qvec dm_dummy,dm_logical;
operator op_list[10];
vec_op *atoms;
operator *atomsstd;
FILE *data_fp = NULL;

typedef struct {
  PetscScalar omega,delta;
  PetscReal length,stime;
} PulseParams;

int main(int argc,char *args[]){

  const PetscInt max_seqgroupsize=8,n_seqgroups=10;
  PetscInt seqgroupsize[10] = {8,8,4,4,4,4,4,4,4,4};
  PetscInt seqgroup[n_seqgroups][max_seqgroupsize];
  PetscInt n_atoms,i,n_levels,val_init,pos1,pos2,dmpos,dmstdpos;
  PetscScalar tmp_scalar = 1.0,rydberg_coupling,b_dr,b_0r,b_1r,gamma_r,valpar,diagsum;
  PetscInt steps_max,n_ens=0,seed=12;
  qvec dm,qvec_std;
  qsystem qsys,qsysstd;
  PetscReal dt,time_max,tm=0.; //tm-> measurement time
  PetscReal single_qubit_gate_time=0.1,two_qubit_gate_time=0.1,var,fidelity;
  char bitstr[PETSC_MAX_PATH_LEN] = "01111"; //Default bitstr to start with
  circuit circ;
  int length;
  PulseParams pulse_params[n_seqgroups];
  //State identifiers
  enum STATE {zero=0,one,r};

  //Below effectively implements a coupling map
  PetscInt aqubit[4]={9,10,11,12};
  PetscInt dqubit[4][4]={
    {0,1,3,4},
    {1,2,4,5},
    {3,4,6,7},
    {4,5,7,8}
  };

  /* Initialize QuaC */
  QuaC_initialize(argc,args);
  //Get the bitstring we want to simulate
  PetscOptionsGetString(NULL,NULL,"-bitstr",bitstr,PETSC_MAX_PATH_LEN,NULL);

  //If n_ens > 0, use MCWF with that many walkers
  //if n_ens < 0, use DM
  //If n_ens = 0, just do the coherent part (good for fast debugging of logical behavior)
  PetscOptionsGetInt(NULL,NULL,"-n_ens",&n_ens,NULL);
  PetscOptionsGetInt(NULL,NULL,"-seed",&seed,NULL);

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

  n_atoms = 13;
  n_levels = 3;

  fillSeqGroup(n_seqgroups,max_seqgroupsize,seqgroup);
  for(i=0;i<n_seqgroups;i++){

    pulse_params[i].stime = i* 0.54;
    if(1<i<6){
			pulse_params[i].stime += tm;
		}
		if(i>=6){
			pulse_params[i].stime += tm;
		}

    pulse_params[i].omega = 17.0*2*PETSC_PI; //MHz 2pi?
 		pulse_params[i].length = 0.54; //us
  	pulse_params[i].delta = 23.0*2*PETSC_PI; //MHz 2pi?
  }

  //branching ratios
  //Careful to use 1.0 instead of 1 because of integer division
  b_1r = 1.0/16.0;
  b_0r = 1.0/16.0;
  b_dr = 7.0/8.0;

  //decay rate
  gamma_r = 1.0/(540.0);

  //rydberg coupling
  rydberg_coupling = 600*2*PETSC_PI; //2pi?

  //Create the operators for the atoms
  atoms = malloc(n_atoms*sizeof(vec_op));
  atomsstd = malloc(n_atoms*sizeof(operator));

  for(i=0;i<n_atoms;i++){
    //Create an n_level system which is stored in atoms[i]
    create_vec_op_sys(qsys,n_levels,&(atoms[i]));
  }

  for(i=0;i<5;i++){
    create_op_sys(qsysstd,2,&(atomsstd[i]));
  }
  data_fp = fopen("neutral_atom_3atom_seq111.dat","w");
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
      //delta is time dependent in this case
      add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params[i],delta,2,atoms[seqgroup[i][j]][r],atoms[seqgroup[i][j]][r]);

      //tmp_scalar = pulse_params[i].delta;
      //MJO: Delta has no time dependence, so we use a normal ham term
      //add_ham_term(qsys,tmp_scalar,2,atoms[seqgroup[i][j]][r],atoms[seqgroup[i][j]][r]);
    }
  }


  //Coupling term
  //rydberg_coupling * (|r_0> <r_0|) (|r_1><r_1|) = (|r_0 r_1><r_0 r_1|)
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      add_ham_term(qsys,rydberg_coupling,4,atoms[aqubit[i]][r],atoms[aqubit[i]][r],atoms[dqubit[i][j]][r],atoms[dqubit[i][j]][r]);
    }
  }


  if(n_ens!=0){
    //We only add linblad terms if n_ens = 0
    //If n_ens > 0, use MCWF
    //if n_ens < 0, use DM
    //If n_ens = 0, just do the coherent part (good for fast debugging of logical behavior)
    //Add lindblad terms
    for(i=0;i<n_atoms;i++){
      tmp_scalar = b_0r*gamma_r;
      //tmp_scalar * L(|0><r|) - no sqrt needed because lin term wants squared term
      add_lin_term(qsys,tmp_scalar,2,atoms[i][zero],atoms[i][r]);

      tmp_scalar = b_1r*gamma_r;
      //L(|1><r|)
      add_lin_term(qsys,tmp_scalar,2,atoms[i][one],atoms[i][r]);

      //L(|0><1|) for testing purposes
      tmp_scalar = 10*gamma_r;
      add_lin_term(qsys,tmp_scalar,2,atoms[i][zero],atoms[i][one]);

      tmp_scalar = b_dr*gamma_r;
      //L(|d><r|)
      /* add_lin_term(qsys,tmp_scalar,2,atoms[i][d],atoms[i][r]); */
      // Instead of a lindblad term, we do an imaginary hamiltonian to allow for leakage
      // need to have a different tmp_scalar based on the repumping simulations
      // Ignore for now
      /* tmp_scalar = b_dr*gamma_r*PETSC_i; */
      /* add_ham_term(qsys,tmp_scalar,2,atoms[i][r],atoms[i][r]); */

    }
  }
  //Now that we've added all the terms, we construct the matrix
  if(n_ens>0){
    use_mcwf_solver(qsys,n_ens,seed);
  }
  construct_matrix(qsys);

  add_ham_term(qsysstd,0.0,1,atomsstd[0]->n);

  //Do we need the std to be a dm?
  if(n_ens!=0){
    add_lin_term(qsysstd,0.0,1,atomsstd[0]->n);
  }
  //

  construct_matrix(qsysstd);

  //Create a density matrix object
  create_qvec_sys(qsysstd,&(qvec_std));
  //create_qvec_sys()
  //crete_wf_sys()
  create_dm_sys(qsysstd,&(dm_logical)); // need a dm because we are effectively doing some partial trace like operation?

  //These are for the real system
  create_qvec_sys(qsys,&(dm));
  create_qvec_sys(qsys,&(dm_dummy));

  add_to_qvec(dm,1.0,dmpos,dmpos); //start in the state defined by the input bitstring
  add_to_qvec(qvec_std,1.0,dmstdpos,dmstdpos); //start in the corresponding state for the qubit

  assemble_qvec(dm);
  assemble_qvec(qvec_std);
  time_max  = 8;
  dt        = 0.01;
  steps_max = 100000;

  /* Set the ts_monitor to print results at each time step */
  /* set_ts_monitor_sys(qsys,ts_monitor,&pulse_params); */


  //Do the timestepping
  PetscPrintf(PETSC_COMM_WORLD,"---------------------\n");
  PetscPrintf(PETSC_COMM_WORLD,"Timestep starting\n");
  time_step_sys(qsys,dm,0.0,time_max,dt,steps_max);
  PetscPrintf(PETSC_COMM_WORLD,"Timestep done\n");
  QuaC_finalize();
  return;
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
              					pos1=l_00*16+l_10*8+l_20*4+l_30*2+l_40*1;
												pos2=l_01*16+l_11*8+l_21*4+l_31*2+l_41*1;
												if (pos1==pos2){
													printf("%d %d %f\n",pos1,pos2,valpar);
													diagsum=diagsum+valpar;
												}
              					add_to_qvec(dm_logical,valpar,pos1,pos2);
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
  assemble_qvec(dm_logical);
  PetscPrintf(PETSC_COMM_WORLD,"dm_logical constructed\n");
  /* print_qvec(dm_logical); */
//----------------------------------------------------------------------------------------------
  create_circuit(&circ,5);
  //Add some gates
  add_gate_to_circuit_sys(&circ,0.0,CZ_ARP,0,1);
  add_gate_to_circuit_sys(&circ,0.5,CZ_ARP,0,2);
  add_gate_to_circuit_sys(&circ,1.0,CZ_ARP,0,3);
  add_gate_to_circuit_sys(&circ,1.5,CZ_ARP,0,4);
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
  apply_circuit_to_qvec(qsysstd,circ,qvec_std);

//-----------------------------------------------------------------------------------------------------------
  get_fidelity_qvec(dm_logical,qvec_std,&fidelity,&var);
  printf("fidelity between seq(32*32) and std is %lf\n",fidelity);
  printf("sum of the diag is %f\n",diagsum);
  //print_qvec(dm);
  //clean up memory
  /* data_fp = fopen("seq_sp_3lvl.dat","a"); */

  /* PetscFPrintf(PETSC_COMM_WORLD,data_fp,"%d %d %f %f\n",dmpos,dmstdpos,fidelity,diagsum); */
  /* fclose(data_fp); */

  destroy_qvec(&(dm));
  destroy_qvec(&(dm_dummy));
  destroy_qvec(&(qvec_std));
  destroy_qvec(&(dm_logical));

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

void fillSeqGroup(const PetscInt n_seqgroups,const PetscInt max_seqgroupsize,PetscInt seqgroup[n_seqgroups][max_seqgroupsize]){
  seqgroup[0][0]=0, seqgroup[0][1]=9, seqgroup[0][2]=10, seqgroup[0][3]=1, seqgroup[0][4]=11, seqgroup[0][5]=6, seqgroup[0][6]=5, seqgroup[0][7]=12;
  seqgroup[1][0]=3, seqgroup[1][1]=9, seqgroup[1][2]=10, seqgroup[1][3]=2, seqgroup[1][4]=11, seqgroup[1][5]=7, seqgroup[1][6]=8, seqgroup[1][7]=12;

  seqgroup[2][0]=9, seqgroup[2][1]=0, seqgroup[2][2]=1, seqgroup[2][3]=10;
  seqgroup[3][0]=9, seqgroup[3][1]=3, seqgroup[3][2]=4, seqgroup[3][3]=10;
  seqgroup[4][0]=9, seqgroup[4][1]=1, seqgroup[4][2]=2, seqgroup[4][3]=10;
  seqgroup[5][0]=9, seqgroup[5][1]=4, seqgroup[5][2]=5, seqgroup[5][3]=10;

  seqgroup[6][0]=4, seqgroup[6][1]=11, seqgroup[6][2]=12, seqgroup[6][3]=5;
  seqgroup[7][0]=7, seqgroup[7][1]=11, seqgroup[7][2]=12, seqgroup[7][3]=8;
  seqgroup[8][0]=3, seqgroup[8][1]=11, seqgroup[8][2]=12, seqgroup[8][3]=4;
  seqgroup[9][0]=6, seqgroup[9][1]=11, seqgroup[9][2]=12, seqgroup[9][3]=7;
  return;
}
