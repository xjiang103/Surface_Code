tm=0.0 //measurement time

//Parameters for doing sequential operation

//number of sequential operations in time domain
n_seqgroups = 10;
//number of qubits in each sequential operation
PetscInt seqgroupsize[2] = {8,8,4,4,4,4,4,4,4,4};
//maximum element in n_seqgroupsize
max_seqgroupsize=8;
//indices of atoms in each sequential operation
PetscInt seqgroup[n_seqgroups][max_seqgroupsize];
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

//parameters for each pulse in the sequence
PulseParams pulse_params[n_seqgroups];

  for(i=0;i<n_seqgroups;i++){

    pulse_params[i].stime = i* 2.0;
    if(1<i<6){
			pulse_params[i].stime += tm;
		}
		if(i>=6){
			pulse_params.stime += tm;
		}

    pulse_params[i].omega = 17.0*(2*3.1416); //MHz 2pi?
  	pulse_params[i].deltat = 0.2; //us
  	pulse_params[i].delta = -0.50*17.0*(2*3.1416); //MHz 2pi?
  }
  //branching ratios
  //Careful to use 1.0 instead of 1 because of integer division
  b_1r = 1.0/16.0;
  b_0r = 1.0/16.0;
  b_dr = 1.0*7.0/8.0;

  //decay rate
  gamma_r = 1.0/(540.0);

  //magnetic field
  b_field = 600*2*PETSC_PI; //2pi?

  //Create the operators for the atoms
  atoms = malloc(n_atoms*sizeof(vec_op));
  atomsstd = malloc(n_atoms*sizeof(operator));

  for(i=0;i<n_atoms;i++){
    //Create an n_level system which is stored in atoms[i]
    create_vec_op_sys(qsys,n_levels,&(atoms[i]));
    create_op_sys(qsysstd,2,&(atomsstd[i]));
  }

  data_fp = fopen("na13_3lvl.dat","w");
  
for(int i=0;i<n_seqgroups;i++){
	for(int j=0;j<seqgroupsize[i];j++){
		tmp_scalar = 1.0;
  	//1.0 * omega(t) * |r><1|
  	add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params[i],omega,2,atoms[seqgroup[i][j]][r],atoms[seqgroup[i][j]][one]);
  	//1.0 * omega(t) * |1><r|
  	add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params[i],omega,2,atoms[seqgroup[i][j]][one],atoms[seqgroup[i][j]][r]);
 		//1.0 * delta(t) * |r><r|
  	add_ham_term_time_dep(qsys,tmp_scalar,&pulse_params[i],delta,2,atoms[seqgroup[i][j]][r],atoms[seqgroup[i][j]][r]);
  }
}

//Coupling term
//b_field * (|r_0> <r_0|) (|r_1><r_1|) = (|r_0 r_1><r_0 r_1|)
PetscInt aqubit[4]={9,10,11,12}
PetscInt dqubit[4][4]={
	{0,1,3,4},
	{1,2,4,5},
	{3,4,6,7},
	{4,5,7,8}
};
for(int i=0;i<4;i++){
	for(int j=0;j<4;j++){
		add_ham_term(qsys,b_field,4,atoms[aqubit[i]][r],atoms[aqubit[i]][r],atoms[dqubit[i][j]][r],atoms[dqubit[i][j]][r]);
  }
}
//Add lindblad terms
for(i=0;i<n_atoms;i++){
  tmp_scalar = b_0r*gamma_r;
  //tmp_scalar * L(|0><r|) - no sqrt needed because lin term wants squared term
  add_lin_term(qsys,tmp_scalar,2,atoms[i][zero],atoms[i][r]);

  tmp_scalar = b_1r*gamma_r;
  //L(|1><r|)
  add_lin_term(qsys,tmp_scalar,2,atoms[i][one],atoms[i][r]);
 
  //L(|0><1|) for testing purposes
  //tmp_scalar = 10*gamma_r;
  //add_lin_term(qsys,tmp_scalar,2,atoms[i][zero],atoms[i][one]);

  //tmp_scalar = b_dr*gamma_r;
  //L(|d><r|)
  /* add_lin_term(qsys,tmp_scalar,2,atoms[i][d],atoms[i][r]); */
  // Instead of a lindblad term, we do an imaginary hamiltonian to allow for leakage
  // need to have a different tmp_scalar based on the repumping simulations
  tmp_scalar = b_dr*gamma_r*PETSC_i;
  add_ham_term(qsys,tmp_scalar,2,atoms[i][r],atoms[i][r]);

}
