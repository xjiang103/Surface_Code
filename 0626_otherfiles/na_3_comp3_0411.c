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

PetscErrorCode ts_monitor_std(TS,PetscInt,PetscReal,Vec,void*);
PetscErrorCode ts_monitor_par(TS,PetscInt,PetscReal,Vec,void*);

qvec dmpar_dummy,dmpar8;
vec_op *atomspar;
operator op_listpar[6];
operator *atomsstd;
FILE *data_fppar = NULL;

typedef struct {
  PetscScalar omega,delta;
  PetscReal length,stime;
} PulseParams;

int main(int argc,char **args){
  PetscInt n_atoms,i,n_levels,n_seqgroups,max_seqgroupsize,pos1,pos2;
  PetscScalar tmp_scalar = 1.0,b_field,b_dr,b_0r,b_1r,gamma_r,valpar,valseq,tmp_val;
  PetscInt steps_max;
  qvec dmpar,dmstd;
  qsystem qsyspar,qsysstd,qsyspar8;
  PetscReal dt,time_max,fidelity,var;
  PetscReal single_qubit_gate_time=0.1,two_qubit_gate_time=0.1;
  circuit circ;
  //PulseParams pulse_params[2]; //moved to line 61
  //State identifiers
  enum STATE {zero=0,d,one,r};

  /* Initialize QuaC */
  QuaC_initialize(argc,args);

  //Initialize the qsystem
  initialize_system(&qsyspar);
  initialize_system(&qsyspar8);
  initialize_system(&qsysstd);

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
  gamma_r = 0.0/(540.0);

  //magnetic field
  b_field = 600*2*PETSC_PI; //2pi?

  //Create the operators for the atoms
  atomspar = malloc(n_atoms*sizeof(vec_op));
	atomsstd = malloc(n_atoms*sizeof(operator));

  for(i=0;i<n_atoms;i++){
    //Create an n_level system which is stored in atoms[i]
    create_vec_op_sys(qsyspar,n_levels,&(atomspar[i]));
    create_op_sys(qsysstd,2,&(atomsstd[i]));
    create_op_sys(qsyspar8,2,&(atomsstd[i]));
  }

  data_fppar = fopen("na_3_comp3_par.dat","w");


  //Add hamiltonian terms

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
  add_ham_term(qsyspar,b_field,4,atomspar[0][r],atomspar[0][r],atomspar[1][r],atomspar[1][r]);
  add_ham_term(qsyspar,b_field,4,atomspar[0][r],atomspar[0][r],atomspar[2][r],atomspar[2][r]);
  add_ham_term(qsyspar,b_field,4,atomspar[1][r],atomspar[1][r],atomspar[2][r],atomspar[2][r]);

  //add_ham_term(qsysstd,0.0,1,atomsstd[0][0]);

  //Add lindblad terms
  for(i=0;i<n_atoms;i++){
    tmp_scalar = b_0r*gamma_r;
    //tmp_scalar * L(|0><r|) - no sqrt needed because lin term wants squared term
    add_lin_term(qsyspar,tmp_scalar,2,atomspar[i][zero],atomspar[i][r]);

    tmp_scalar = b_1r*gamma_r;
    //L(|1><r|)
    add_lin_term(qsyspar,tmp_scalar,2,atomspar[i][one],atomspar[i][r]);

    tmp_scalar = b_dr*gamma_r;
    //L(|d><r|)
    add_lin_term(qsyspar,tmp_scalar,2,atomspar[i][d],atomspar[i][r]);

  }

  //Now that we've added all the terms, we construct the matrix
  construct_matrix(qsyspar);
  add_ham_term(qsysstd,0.0,1,atomsstd[0]->n);
  add_lin_term(qsysstd,0.0,1,atomsstd[0]->n);
  construct_matrix(qsysstd);
  add_ham_term(qsyspar8,0.0,1,atomsstd[0]->n);
  construct_matrix(qsyspar8);
  /* print_mat_sparse(qsys->mat_A); */
  //Create a density matrix object
  create_qvec_sys(qsyspar,&(dmpar));
  create_qvec_sys(qsyspar,&(dmpar_dummy));


  create_dm_sys(qsysstd,&(dmstd));
  create_dm_sys(qsysstd,&(dmpar8));

	FILE *fp_read=NULL;
	fp_read=fopen("init_0111.txt","r");
	PetscInt init_enum[8]={0,2,8,10,32,34,40,42};
	PetscReal init_real[8]={0};
	PetscReal init_imag[8]={0};
	PetscScalar init_state[8];
	PetscScalar init_state_conj[8];

	for(int i=0;i<8;i++){
		fscanf(fp_read,"%lf %lf",&init_real[i],&init_imag[i]);
		//printf("%d\n",i);
		printf("%lf + %lf j\n", init_real[i], init_imag[i]);
		init_state[i]=init_real[i] + init_imag[i]*PETSC_i;
		init_state_conj[i]=init_real[i] - init_imag[i]*PETSC_i;
	}

  for(int i=0;i<8;i++){
  	for (int j=0;j<8;j++){
      tmp_val = init_state[i]*init_state_conj[j];
      printf("Added dm[%d][%d]: %f %f\n",i,j,tmp_val);
  		add_to_qvec(dmpar,init_state[i]*init_state_conj[j],init_enum[i],init_enum[j]);
//  	  printf("%lf + %lf j, ",(double)PetscRealPart(init_state[i]*init_state_conj[j]),(double)PetscImaginaryPart(init_state[i]*init_state_conj[j]));
		}
//		printf("\n");

	fp_read=fopen("final_0111.txt","r");
	}
	for(int i=0;i<8;i++){
		fscanf(fp_read,"%lf %lf",&init_real[i],&init_imag[i]);
		//printf("%d\n",i);
		printf("%lf + %lf j\n", init_real[i], init_imag[i]);
		init_state[i]=init_real[i] + init_imag[i]*PETSC_i;
		init_state_conj[i]=init_real[i] - init_imag[i]*PETSC_i;
	}

  for(int i=0;i<8;i++){
  	for (int j=0;j<8;j++){
      tmp_val = init_state[i]*init_state_conj[j];
      printf("Added dm[%d][%d]: %f %f\n",i,j,tmp_val);
  		add_to_qvec(dmstd,init_state[i]*init_state_conj[j],i,j);
//  	  printf("%lf + %lf j, ",(double)PetscRealPart(init_state[i]*init_state_conj[j]),(double)PetscImaginaryPart(init_state[i]*init_state_conj[j]));
		}
//		printf("\n");
	}
//	add_to_qvec(dmpar,1.0,34,34); //start in the |101><101| state
//	add_to_qvec(dmseq,1.0,34,34); //start in the |101><101| state
//	add_to_qvec(dmstd,1.0,5,5); //start in the |101><101| state
  assemble_qvec(dmpar);
  assemble_qvec(dmstd);

  printf("Initialization\n");
  print_qvec(dmstd);
  //   exit(9);

  time_max  = 1.08;
  dt        = 0.00000001;
  steps_max = 1000000;

  /* Set the ts_monitor to print results at each time step */
  set_ts_monitor_sys(qsyspar,ts_monitor_par,&pulse_params);
  time_step_sys(qsyspar,dmpar,0.0,time_max,dt,steps_max);

  int enum_list[2];
  enum_list[0] = zero;
  enum_list[1] = one;
  valpar=0.0;
	valseq=0.0;
  for(int l_00=0;l_00<2;l_00++){ //atom 0, 0th element
    for(int l_01=0;l_01<2;l_01++){ //atom 0, 1st element
      for(int l_10=0;l_10<2;l_10++){ //atom 1, 0th element
        for(int l_11=0;l_11<2;l_11++){ //atom 1, 1st element
          for(int l_20=0;l_20<2;l_20++){ //atom 2, 0th element
            for(int l_21=0;l_21<2;l_21++){ //atom 2, 1st element
              op_listpar[0] = atomspar[0][enum_list[l_00]]; //atom 0 0th element
              op_listpar[1] = atomspar[0][enum_list[l_01]]; //atom 0 1st element
              op_listpar[2] = atomspar[1][enum_list[l_10]]; //atom 1 0th element
              op_listpar[3] = atomspar[1][enum_list[l_11]]; //atom 1 1st element
              op_listpar[4] = atomspar[2][enum_list[l_20]]; //atom 2 0th element
              op_listpar[5] = atomspar[2][enum_list[l_21]]; //atom 2 1st element
              get_expectation_value_qvec_list(dmpar,&valpar,6,op_listpar);
              pos1=l_00*4+l_10*2+l_20*1;
							pos2=l_01*4+l_11*2+l_21*1;
              add_to_qvec(dmpar8,valpar,pos1,pos2);
            }
          }
        }
      }
    }
  }
  assemble_qvec(dmpar8);

  printf("dmpar8\n");
  print_qvec(dmpar8);

  //Do the timestepping
  //time_step_sys(qsysstd,dmstd,0.0,time_max,dt,steps_max);
  var=1;
  get_fidelity_qvec(dmpar8,dmstd,&fidelity,&var);
  printf("fidelity between par(8*8) and std is %lf\n",fidelity);

  destroy_system(&qsyspar);
  destroy_qvec(&(dmpar));
  destroy_qvec(&(dmpar_dummy));
  for(i=0;i<n_atoms;i++){
    destroy_vec_op_sys(&atomspar[i]);
  }
  destroy_system(&qsysstd);
  destroy_qvec(&(dmstd));
  for(i=0;i<n_atoms;i++){
    destroy_op_sys(&atomsstd[i]);
  }
  destroy_system(&qsyspar8);
  destroy_qvec(&(dmpar8));
  return;
}

PetscErrorCode ts_monitor_std(TS ts,PetscInt step,PetscReal time,Vec rho_data,void *ctx){
	  PetscFunctionReturn(0);
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
