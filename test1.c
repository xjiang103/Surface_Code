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

qvec dm_dummy,dm32;

operator *atomsstd;
FILE *data_fp = NULL;

typedef struct {
  PetscScalar omega,delta;
  PetscReal stime,length;
} PulseParams;

int main(int argc,char **args){
  PetscInt n_atoms,i,n_levels;
  PetscScalar tmp_scalar = 1.0,b_field,b_dr,b_0r,b_1r,gamma_r,valpar,diagsum;
  PetscInt steps_max;
  qvec dm,dmstd;
  qsystem qsys;
  PetscReal dt,time_max;

  //PulseParams pulse_params[2]; //moved to line 61
  //State identifiers
  /* Initialize QuaC */
  QuaC_initialize(argc,args);
  //Initialize the qsystem
  initialize_system(&qsys);

  n_atoms = 1;
  n_levels = 32;
  
  atoms = malloc(n_atoms*sizeof(vec_op));
    
  PetscInt lvlinfo[32][2];  
  //6s1/2 levels(f,mf):
  lvlinfo[0][0]=3, lvlinfo[0][1]=-3;
  lvlinfo[1][0]=3, lvlinfo[1][1]=-2; 
  lvlinfo[2][0]=3, lvlinfo[2][1]=-1;
  lvlinfo[3][0]=3, lvlinfo[3][1]=0;
  lvlinfo[4][0]=3, lvlinfo[3][1]=1; 
  lvlinfo[5][0]=3, lvlinfo[3][1]=2; 
  lvlinfo[6][0]=3, lvlinfo[3][1]=3; 
  lvlinfo[7][0]=4, lvlinfo[3][1]=-4;  
  lvlinfo[8][0]=4, lvlinfo[3][1]=-3;
  lvlinfo[9][0]=4, lvlinfo[3][1]=-2;
  lvlinfo[10][0]=4, lvlinfo[3][1]=-1;
  lvlinfo[11][0]=4, lvlinfo[3][1]=0;
  lvlinfo[12][0]=4, lvlinfo[3][1]=1;
  lvlinfo[13][0]=4, lvlinfo[3][1]=2;
  lvlinfo[14][0]=4, lvlinfo[3][1]=3;
  lvlinfo[15][0]=4, lvlinfo[3][1]=4;
  
  //7p1/2 levels(f,mf):
  lvlinfo[16][0]=3, lvlinfo[0][1]=-3;
  lvlinfo[17][0]=3, lvlinfo[1][1]=-2; 
  lvlinfo[18][0]=3, lvlinfo[2][1]=-1;
  lvlinfo[19][0]=3, lvlinfo[3][1]=0;
  lvlinfo[20][0]=3, lvlinfo[3][1]=1; 
  lvlinfo[21][0]=3, lvlinfo[3][1]=2; 
  lvlinfo[22][0]=3, lvlinfo[3][1]=3; 
  lvlinfo[23][0]=4, lvlinfo[3][1]=-4;  
  lvlinfo[24][0]=4, lvlinfo[3][1]=-3;
  lvlinfo[25][0]=4, lvlinfo[3][1]=-2;
  lvlinfo[26][0]=4, lvlinfo[3][1]=-1;
  lvlinfo[27][0]=4, lvlinfo[3][1]=0;
  lvlinfo[28][0]=4, lvlinfo[3][1]=1;
  lvlinfo[29][0]=4, lvlinfo[3][1]=2;
  lvlinfo[30][0]=4, lvlinfo[3][1]=3;
  lvlinfo[31][0]=4, lvlinfo[3][1]=4;
  //parameters for each pulse in the sequence
  
  PulseParams pulse_params[n_seqgroups];
  for(i=0;i<n_seqgroups;i++){

    pulse_params[i].stime = stime;
  	pulse_params[i].length = pulse_length; //
    pulse_params[i].omega = pulse_omega; //
  	pulse_params[i].delta = delta; //need to be calculated
  }
  
  for(int i=0;i<n_seqgroups;i++){
  	for(int j=0;j<seqgroupsize[i];j++){
  		tmp_scalar = 1.0;

    	add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params[i],omega,2,atoms[0][r],atoms[0][one]);

    	add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params[i],omega,2,atoms[0][one],atoms[0][r]);

    	add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params[i],delta,2,atoms[0][r],atoms[0][r]);
  //branching ratios
  //Careful to use 1.0 instead of 1 because of integer division
  b_1r = 1.0/16.0;
  b_0r = 1.0/16.0;
  b_dr = 7.0/8.0;

  //decay rate
  gamma_r = 1.0/(0.0304);

  //magnetic field
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

  
  //create_qvec_sys()
  //crete_wf_sys()

  create_qvec_sys(qsys,&(dm));

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
  time_max  = n_seqgroups*10*deltat;
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

  time_step_sys(qsys,dm,0.0,time_max/4,dt,steps_max);
  time_step_sys(qsys,dm,time_max/4,time_max/2,dt,steps_max);
  time_step_sys(qsys,dm,time_max/2,3*time_max/4,dt,steps_max);
  time_step_sys(qsys,dm,3*time_max/4,time_max,dt,steps_max);

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
  printf("dm32 constructed\n");
  //print_qvec(dm32);  
//----------------------------------------------------------------------------------------------
  set_ts_monitor_sys(qsysstd,ts_monitor_std,&pulse_params);
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
  //printf("Before\n");
  //print_qvec(dmstd);  
  apply_circuit_to_qvec(qsysstd,circ,dmstd);
  //printf("dmstd\n");
 // print_qvec(dmstd);
//-----------------------------------------------------------------------------------------------------------
  get_fidelity_qvec(dm32,dmstd,&fidelity,&var);
  data_fp = fopen("seq_sp.dat","a");
  
  PetscFPrintf(PETSC_COMM_WORLD,data_fp,"%d %d %f %f\n",dmpos,dmstdpos,fidelity,diagsum);
   
	printf("sum of the diag is %f\n",diagsum);
  printf("fidelity between seq(32*32) and std is %lf\n",fidelity);

  //print_qvec(dm);
  //clean up memory

  destroy_system(&qsys);
  destroy_qvec(&(dm));
  //destroy_qvec(&(dm_dummy));
  destroy_system(&qsysstd);
  destroy_qvec(&(dmstd));
  for(i=0;i<n_atoms;i++){
    destroy_vec_op_sys(&atoms[i]);
    destroy_op_sys(&atomsstd[i]);
  }

  return;
}  
