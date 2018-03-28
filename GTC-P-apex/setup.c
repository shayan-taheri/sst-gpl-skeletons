#include "RngStream.h"
#include "bench_gtc.h"

void usage (const char* progname) {

  fprintf(stderr, "usage: %s <filename> <micell> <ntoroidal>\n", progname);
  fprintf(stderr,
	  "\n"
	  "GTC benchmark.\n"
	  "  <filename>\n"
	  "    Input file.\n"
	  "  <micell>\n"
	  "    Number of particles per cell.\n"
	  "  <ntoroidal>\n"
	  "    1D domain decomposition (< mzetamax).\n"
	  "\n");
}

int read_input_file (char *filename, gtc_global_params_t *global_params, 
		     gtc_particle_decomp_t *parallel_decomp, 
		     gtc_radial_decomp_t *radial_decomp) {

  FILE *fp;
  char real_scan_str[100], int_scan_str[100], buf[100];
  int mype;

  gtc_global_params_t* gp = global_params;

  mype = parallel_decomp->mype;

#if SINGLE_PRECISION 
  strcpy(real_scan_str, "%*s %f");    
#else
  strcpy(real_scan_str, "%*s %lf");
#endif
  strcpy(int_scan_str, "%*s %d");

  if (mype == 0) {
    fp = fopen(filename, "r");
    while (fgets(buf, 100, fp) != NULL)  {
  
      /* comment symbols */
      if (buf[0] == '!' || buf[0] == '#' || buf[0] == '/')
	continue;

      if (strncmp(buf,"irun",4) == 0)
	sscanf(buf, int_scan_str, &(global_params->irun));

      if (strncmp(buf,"mstep",5) == 0)
	sscanf(buf, int_scan_str, &(global_params->mstep));

      if (strncmp(buf,"mpsi",4) == 0) 
	sscanf(buf, int_scan_str, &(global_params->mpsi));
    
      if (strncmp(buf,"mthetamax",9) == 0)
	sscanf(buf, int_scan_str, &(global_params->mthetamax));

      if (strncmp(buf,"mzetamax",8) == 0) 
	sscanf(buf, int_scan_str, &(global_params->mzetamax));

      if (strncmp(buf,"npe_radiald",11) == 0)
	sscanf(buf, int_scan_str, &(radial_decomp->npe_radiald));

      if (strncmp(buf,"r0",2) == 0)
	sscanf(buf, real_scan_str, &(global_params->r0));

      if (strncmp(buf, "nbound", 6) == 0)
	sscanf(buf, int_scan_str, &(global_params->nbound));

    }
    fclose(fp);
    fprintf(stderr, "Finished reading input file.\n");
  }

#if USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  /* mype=0 broadcasts the values read from input file */
  MPI_Bcast(&(global_params->irun), 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&(global_params->mstep), 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&(global_params->mpsi), 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&(global_params->mthetamax), 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&(global_params->mzetamax), 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&(radial_decomp->npe_radiald),1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&(global_params->nbound), 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&(global_params->r0), 1, MPI_MYREAL, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
    
  return 0;
}

int gtc_mem_free(gtc_bench_data_t* gtc_input) {

  gtc_field_data_t        *field_data;
  gtc_particle_data_t     *particle_data;
  gtc_aux_particle_data_t *aux_particle_data;
  gtc_radial_decomp_t     *radial_decomp;
  gtc_particle_decomp_t   *parallel_decomp;
  gtc_diagnosis_data_t    *diagnosis_data;

  field_data        = &(gtc_input->field_data);
  aux_particle_data = &(gtc_input->aux_particle_data);
  particle_data     = &(gtc_input->particle_data);
  parallel_decomp   = &(gtc_input->parallel_decomp);
  radial_decomp = &(gtc_input->radial_decomp);
  diagnosis_data = &(gtc_input->diagnosis_data);

  _mm_free(field_data->igrid);
  _mm_free(field_data->qtinv); 
  _mm_free(field_data->mtheta); 
  _mm_free(field_data->pgyro); 
  _mm_free(field_data->tgyro); 
  _mm_free(field_data->densityi); 
  _mm_free(field_data->densityi_local); 
  _mm_free(field_data->delt); 
  _mm_free(field_data->den00);
  _mm_free(field_data->difft);
  _mm_free(field_data->rtemi);
  _mm_free(field_data->temp); 
  _mm_free(field_data->dtemp);
  _mm_free(field_data->evector); 
  _mm_free(field_data->heatflux);
  _mm_free(field_data->itran);
  _mm_free(field_data->rden);
  _mm_free(field_data->dnitmp);
  _mm_free(field_data->jtp1);
  _mm_free(field_data->jtp2);
  _mm_free(field_data->wtp1); 
  _mm_free(field_data->wtp2);
  _mm_free(field_data->markeri);
  _mm_free(field_data->phitmps);

  _mm_free(field_data->hfluxpsi);
  _mm_free(field_data->pfluxpsi);
  _mm_free(field_data->vdrtmp);
  _mm_free(field_data->zonali); 
  _mm_free(field_data->sendl);
  _mm_free(field_data->recvr); 
  _mm_free(field_data->dtemper);
  _mm_free(field_data->phi);
  _mm_free(field_data->pmarki);
  _mm_free(field_data->phi00);
  _mm_free(field_data->phip00);
  _mm_free(field_data->gradt);
  _mm_free(field_data->deltat);

  _mm_free(field_data->indexp);
  _mm_free(field_data->nindex);
  _mm_free(field_data->ring);

  _mm_free(field_data->sendlf);
  _mm_free(field_data->sendrf);
  _mm_free(field_data->recvlf);
  _mm_free(field_data->recvrf);
  _mm_free(field_data->sendrsf);
  _mm_free(field_data->recvlsf);
   
  _mm_free(field_data->perr);
  _mm_free(field_data->ptilde);
  _mm_free(field_data->phitmp);
  _mm_free(field_data->dentmp);

  _mm_free(field_data->drdpa);
  _mm_free(field_data->diffta);
  _mm_free(field_data->idx1a);
  _mm_free(field_data->idx2a);
  _mm_free(field_data->recvl_index);
  _mm_free(field_data->recvr_index);

  _mm_free(particle_data->z0); 
  _mm_free(particle_data->z1); 
  _mm_free(particle_data->z2); 
  _mm_free(particle_data->z3); 
  _mm_free(particle_data->z4); 
  _mm_free(particle_data->z5);

  _mm_free(particle_data->z00); 
  _mm_free(particle_data->z01); 
  _mm_free(particle_data->z02); 
  _mm_free(particle_data->z03); 
  _mm_free(particle_data->z04); 
  _mm_free(particle_data->z05);
  _mm_free(particle_data->ztmp);
  _mm_free(particle_data->ztmp2);
  _mm_free(particle_data->psi_count);
  _mm_free(particle_data->psi_offsets);
  _mm_free(aux_particle_data->kzion); 
  _mm_free(aux_particle_data->jtion0); 
  _mm_free(aux_particle_data->jtion1); 
  _mm_free(aux_particle_data->wzion); 
  _mm_free(aux_particle_data->wpion); 
  _mm_free(aux_particle_data->wtion0); 
  _mm_free(aux_particle_data->wtion1);
  _mm_free(aux_particle_data->kzi);
#if USE_MPI
  _mm_free(parallel_decomp->recvbuf);
  _mm_free(parallel_decomp->sendbuf);
#endif
  _mm_free(radial_decomp->ri_pe);
  _mm_free(radial_decomp->ri_pe2);
  _mm_free(radial_decomp->ghost_comm_list);
  _mm_free(radial_decomp->ghost_start);
  _mm_free(radial_decomp->ghost_end);
  _mm_free(radial_decomp->ghost_sendrecvbuf);

  _mm_free(radial_decomp->nghost_comm_list);
  _mm_free(radial_decomp->nghost_start);
  _mm_free(radial_decomp->nghost_end);
  _mm_free(radial_decomp->nghost_sendrecvbuf);
    
  _mm_free(diagnosis_data->scalar_data);
  _mm_free(diagnosis_data->eflux);
  _mm_free(diagnosis_data->rmarker);
  _mm_free(diagnosis_data->dmark);
  _mm_free(diagnosis_data->dden);
  _mm_free(diagnosis_data->rdtemi);
  _mm_free(diagnosis_data->rdteme);
  _mm_free(diagnosis_data->flux_data);
  _mm_free(diagnosis_data->amp_mode);
  _mm_free(diagnosis_data->eigenmode);
  _mm_free(diagnosis_data->nmode);
  _mm_free(diagnosis_data->mmode);

  return 0;
}

int setup(gtc_bench_data_t* gtc_input) {

  gtc_global_params_t     *global_params;
  gtc_field_data_t        *field_data;
  gtc_aux_particle_data_t *aux_particle_data;
  gtc_particle_data_t     *particle_data;
  gtc_particle_decomp_t   *parallel_decomp;
  gtc_radial_decomp_t     *radial_decomp;
  gtc_diagnosis_data_t    *diagnosis_data;

  /* parameters that were read from the command-line/input file */
  int irun, mstep, numberpe, mype, ntoroidal, micell, mpsi, mthetamax, mzetamax, npe_radiald;
  real r0;

  /* parameters that are initialized in this routine */
  /* temporary variables */
#if USE_MPI
  MPI_Comm *toroidal_comm, *partd_comm, *radiald_comm, *radial_partd_comm;
  //MPI_Comm *comm2d;
  //int rank2d; //YKT on BGQ
#endif

  int mtdiag;
  int msnap;
  real nonlinear;
  int mode00;
  int mi, mgrid, mzeta;
  real zetamin, zetamax,  smu_inv;
  real delr, delz, pi2_inv;
  real a, a0, a1, q0, q1, q2, rc, rw, 
    aion, qion, b0, temperature, umax;
  int ndiag;


  int *kzion, *jtion0, *jtion1;
  real *wzion, *wpion, *wtion0, *wtion1;
  int *kzi;

  int *igrid, *mtheta;
  real *delt, *qtinv, *pgyro, *tgyro;
  wreal *densityi;
  wreal *densityi_local;

  real *evector; real *rtemi;
  real *temp, *dtemp;

  real *den00;
  real *difft;
  real *rden, *pfluxpsi, *vdrtmp, *hfluxpsi, *zonali, *wtp1, *wtp2;
  real *phi, *dtemper, *heatflux;
  int *jtp1, *jtp2;

  real *drdpa, *diffta;
  int *idx1a, *idx2a;
  int *recvl_index, *recvr_index;

  int *psi_count;
  int *psi_offsets;

  int *indexp, *nindex;
  real *ring;
  wreal *recvr; wreal *sendl;
  wreal *dnitmp;

  real *perr, *ptilde, *phitmp, *dentmp;

  real *pmarki, *phi00, *phip00, *gradt;
  real *markeri; 
  real *phitmps;

  real *sendrf, *recvrf, *sendlf, *recvlf, *sendrsf, *recvlsf;

  real *z0, *z1, *z2, *z3, *z4, *z5;
  real *z00, *z01, *z02, *z03, *z04, *z05;
  real *ztmp, *ztmp2;

  // variables for diagnosis
  real *scalar_data;
  real *eflux, *rmarker, *dmark, *dden, *rdtemi, *rdteme;
  real *flux_data;
  real *amp_mode, *eigenmode;
  int *nmode, *mmode;

  // variables for toroidal decomposition
  int nradial_dom, npartdom, particle_domain_location, toroidal_domain_location;
  int nproc_partd, myrank_partd, nproc_toroidal, myrank_toroidal, left_pe, right_pe; 

  // variables for radial decomposition 
  int n_rad_buf;
  int nproc_radiald, myrank_radiald, left_radial_pe, right_radial_pe;
  int nproc_radial_partd, myrank_radial_partd;
  int radial_domain_location, radial_part_domain_location;
   
  int ipsi_nover_in, ipsi_nover_out, ipsi_in, ipsi_out,
    ipsi_valid_in, ipsi_valid_out, igrid_in, igrid_out,
    igrid_nover_in, igrid_nover_out, igrid_valid_in,
    igrid_valid_out, nloc_nover, nloc_over;
  int ipsi_nover_in_radiald, ipsi_nover_out_radiald, igrid_nover_in_radiald, igrid_nover_out_radiald;

  real a_in, a_out, a_nover_in, a_nover_out, a_valid_in, a_valid_out, rho_max;               
  int *ri_pe, *ri_pe2, *tempi, *arr, *arr2;
  real *adum;
    
  // variables for radial domain grid ghost cell communication
  int *ghost_comm_list, *ghost_start, *ghost_end, *nghost_comm_list, *nghost_start, *nghost_end;
  real *ghost_sendrecvbuf,*nghost_sendrecvbuf;
  int ghost_comm_num, nghost_comm_num, ghost_bufsize, nghost_bufsize;
  int i_loc;

  // temparory variables
  int *ghost_start_local, *ghost_end_local, *ghost_start_global, *ghost_end_global;

  int nthreads;
  /* int ntracer; */
  int mflux;

  RngStream rs;
  double *rng_seed;

  real deltaz, deltar, pi, ulength, gyroradius, rad, zdum, tdum;
  int mthetatmp, mimax;
  int *itran;
  real *deltat;
  int i, j, k, ip, indp, indt, ij, jt, m, mt;
  uint64_t mi_local;
  real b, q, dtheta_dx, w_initial, zt1, zt2, vthi, z3tmp, wt;
  real c0, c1, c2, d1, d2, d3, rmi;
  real kappati, tstep;
  real r, rmax, rmin;
#if USE_MPI
  real *recvbuf, *sendbuf;
  int recvbuf_size, sendbuf_size;
#endif
  /* 
   * ************************************************************************
   */

  global_params     = &(gtc_input->global_params);
  field_data        = &(gtc_input->field_data);
  aux_particle_data = &(gtc_input->aux_particle_data);
  particle_data     = &(gtc_input->particle_data);
  parallel_decomp   = &(gtc_input->parallel_decomp);
  radial_decomp     = &(gtc_input->radial_decomp);
  diagnosis_data    = &(gtc_input->diagnosis_data);
    
  /* load the values that were previously read from input file */
  mype       = parallel_decomp->mype;
  ntoroidal  = parallel_decomp->ntoroidal;
  numberpe   = parallel_decomp->numberpe;

  irun       = global_params->irun;
  micell     = global_params->micell;
  mstep      = global_params->mstep;
  mpsi       = global_params->mpsi;
  mthetamax  = global_params->mthetamax; 
  mzetamax   = global_params->mzetamax;

  gtc_global_params_t* gp = global_params;

  /* initialize other global parameters */
  //BMA: set diagnotic period to run only at end of job 
  //ndiag = global_params->ndiag = mstep+1;
  ndiag = global_params->ndiag = mstep;
  msnap = global_params->msnap = 1;
  nonlinear = global_params->nonlinear = 1.0;
  global_params->paranl = 1.0;
  mode00 = global_params->mode00 = 1;
  tstep = global_params->tstep = 0.15;
  global_params->ncycle = 0;
  a = global_params->a = 0.358; 
  a0 = global_params->a0  = 0.1; 
  a1 = global_params->a1 = 0.9; 
  q0 = global_params->q0 = 0.854;
  q1 = global_params->q1 = 0.0;
  q2 = global_params->q2 = 2.184; 
  rc = global_params->rc = 0.5;
  rw = global_params->rw = 0.35; 
  aion = global_params->aion = 1.0; 
  qion = global_params->qion = 1.0; 
  kappati = global_params->kappati = 6.9;
  global_params->aelectron = 5.443658e-4;
  global_params->qelectron = -1.0;
  global_params->kappan = 2.2;
  global_params->tite = 1.0;
  global_params->flow0 = 0.0;
  global_params->flow1 = 0.0;
  global_params->flow2 = 0.0;
  r0 = global_params->r0; 
  b0 = global_params->b0 = 19100.0; 
  temperature = global_params->temperature = 2500.0; 
  umax = global_params->umax = 4.0;
  mflux = global_params->mflux = 5;
  global_params->mmomentsoutput = 5;
  npe_radiald = radial_decomp->npe_radiald;
  if (irun==0) global_params->mstepall = 0;

#if USE_MPI
  toroidal_comm = &(parallel_decomp->toroidal_comm);
  partd_comm = &(parallel_decomp->partd_comm);
  radiald_comm = &(radial_decomp->radiald_comm);
  radial_partd_comm = &(radial_decomp->radial_partd_comm);
#endif

#pragma omp parallel
  {
#ifdef _OPENMP
#pragma omp master
    nthreads = omp_get_num_threads();
#else
    nthreads = 1;
#endif
#pragma omp barrier
    //assert(nthreads>=NUM_THREADS);
  }

  /* Change the units of a0 and a1 from 
   * units of "a" to units of "R_0" */
  a0 = a0*a;
  a1 = a1*a;

  /* Numerical constant */
  pi = 4.0 * atan(1.0);
    
  if (mstep < 2)
    mstep = 2;
    
  if (mstep/ndiag < msnap)
    msnap = mstep/ndiag;

  if (nonlinear < 0.5) {
    global_params->paranl = 0.0;
    mode00 = 0;
  }

  rc = rc*(a0+a1);
  rw = 1.0/(rw*(a1-a0));

  /* Particle decomposition */
  /* ntoroidal * npartdom * nradial_dom = numberpe */
  /* ntoroidal and nradial_dom are specified in the input file */
  if (numberpe == 1)
    npartdom = ntoroidal = nradial_dom = 1;
  else
    npartdom = numberpe/ntoroidal;
      
  if (mype == 0) {
    fprintf(stderr, "Using npartdom %d, ntoroidal %d, nthreads %d\n", 
	    npartdom, ntoroidal, nthreads);
  }
  assert(npartdom >= 1);

  /* Make sure that mzetamax is a multiple of ntoroidal */
  mzetamax = ntoroidal * (mzetamax/ntoroidal);
  if (mzetamax == 0)
    mzetamax = ntoroidal;

  /* Make sure that mpsi is even */ 
  mpsi = 2*(mpsi/2);
    
#if USE_MPI
#ifdef GTC_COMM2D
  // YKT:optionally use a 2d communicator to set domains and split comms
  if (mype==0) printf("use customized MPI mapping on BGQ\n");
  get_2d_communicator(comm2d);
#else
  if (mype==0) printf("NOT use customized MPI mapping on BGQ\n");
  //*comm2d = MPI_COMM_WORLD;
#endif
  //MPI_Comm_rank(*comm2d, &rank2d);
#endif

  //particle_domain_location = rank2d % npartdom;
  //toroidal_domain_location = rank2d / npartdom;
  /* assign each PE two domain numbers, particle and toroidal */
  particle_domain_location = mype % npartdom;
  toroidal_domain_location = mype / npartdom;

  radial_part_domain_location = particle_domain_location % npe_radiald;
  radial_domain_location = particle_domain_location / npe_radiald;

  /* domain decomposition in toroidal direction */
  mzeta = mzetamax/ntoroidal;
  zetamin = 2.0 * pi * ((real) toroidal_domain_location)/ntoroidal;
  zetamax = 2.0 * pi * ((real) (toroidal_domain_location+1))/ntoroidal;
  deltaz  = (zetamax - zetamin)/mzeta;

#if USE_MPI
  MPI_Comm_split(MPI_COMM_WORLD, toroidal_domain_location,
		 particle_domain_location, partd_comm);
  MPI_Comm_split(MPI_COMM_WORLD, particle_domain_location,
		 toroidal_domain_location, toroidal_comm);
  MPI_Comm_size(*partd_comm, &nproc_partd);
  MPI_Comm_rank(*partd_comm, &myrank_partd);
  MPI_Comm_size(*toroidal_comm, &nproc_toroidal);
  MPI_Comm_rank(*toroidal_comm, &myrank_toroidal);

  MPI_Comm_split(*partd_comm, radial_domain_location,
		 radial_part_domain_location, radial_partd_comm);
  MPI_Comm_split(*partd_comm, radial_part_domain_location,
		 radial_domain_location, radiald_comm);
  MPI_Comm_size(*radial_partd_comm, &nproc_radial_partd);
  MPI_Comm_rank(*radial_partd_comm, &myrank_radial_partd);
  MPI_Comm_size(*radiald_comm, &nproc_radiald);
  MPI_Comm_rank(*radiald_comm, &myrank_radiald);
  nradial_dom = nproc_radiald;
#else
  nproc_partd = npartdom;
  nproc_toroidal = 1;
  myrank_partd = 0;
  myrank_toroidal = 0;
  nradial_dom = 1;
  nproc_radiald = 1;
  nproc_radial_partd = 1;
  myrank_radiald = 0;
  myrank_radial_partd = 0;
#endif

  if (mype == 0) {
    fprintf(stderr, "number of radial parition=%d number of particle copies=%d\n",
	    nproc_radiald, nproc_radial_partd);
  }

  /* find toroidal neighbors */
  left_pe  = (myrank_toroidal-1+ntoroidal) % ntoroidal;
  right_pe = (myrank_toroidal+1) % ntoroidal;

  left_radial_pe = myrank_partd - nproc_radial_partd;
  if (left_radial_pe < 0) left_radial_pe = left_radial_pe + nproc_partd;
  right_radial_pe = myrank_partd + nproc_radial_partd;
  if (right_radial_pe >= nproc_partd) right_radial_pe = right_radial_pe - nproc_partd;
    
  /* Equilibrium unit: length (unit=cm) and time (unit=second) unit */
  ulength = r0;
  gyroradius = 102.0 * sqrt(aion*temperature)/(fabs(qion)*b0)/ulength;
  tstep = tstep * aion/(fabs(qion)*gyroradius*kappati);

  qtinv    = (real *) _mm_malloc((mpsi+1)*sizeof(real), IDEAL_ALIGNMENT);
  itran    = (int  *) _mm_malloc((mpsi+1)*sizeof(int),  IDEAL_ALIGNMENT);
  mtheta   = (int  *) _mm_malloc((mpsi+1)*sizeof(int),  IDEAL_ALIGNMENT);
  deltat   = (real *) _mm_malloc((mpsi+1)*sizeof(real), IDEAL_ALIGNMENT);
  rtemi    = (real *) _mm_malloc((mpsi+1)*sizeof(real), IDEAL_ALIGNMENT);
  rden     = (real *) _mm_malloc((mpsi+1)*sizeof(real), IDEAL_ALIGNMENT);
  igrid    = (int  *) _mm_malloc((mpsi+1)*sizeof(int),  IDEAL_ALIGNMENT);
  pmarki   = (real *) _mm_malloc((mpsi+1)*sizeof(real),  IDEAL_ALIGNMENT);
  phi00    = (real *) _mm_malloc((mpsi+1)*sizeof(real),  IDEAL_ALIGNMENT);
  den00    = (real *) _mm_malloc((mpsi+1)*sizeof(real),  IDEAL_ALIGNMENT);
  phip00   = (real *) _mm_malloc((mpsi+1)*sizeof(real),  IDEAL_ALIGNMENT);
  hfluxpsi = (real *) _mm_malloc((mpsi+1)*sizeof(real), IDEAL_ALIGNMENT);
  zonali   = (real *) _mm_malloc((mpsi+1)*sizeof(real), IDEAL_ALIGNMENT);
  gradt    = (real *) _mm_malloc((mpsi+1)*sizeof(real),  IDEAL_ALIGNMENT);
  delt     = (real *) _mm_malloc((mpsi+1)*sizeof(real), IDEAL_ALIGNMENT);
  pfluxpsi = (real *) _mm_malloc(mflux*sizeof(real), IDEAL_ALIGNMENT);
  vdrtmp   = (real *) _mm_malloc(mflux*sizeof(real), IDEAL_ALIGNMENT);

  temp   = (real *) _mm_malloc((mpsi+1)*sizeof(real), IDEAL_ALIGNMENT);
  dtemp  = (real *) _mm_malloc((mpsi+1)*sizeof(real), IDEAL_ALIGNMENT);
  difft  = (real *) _mm_malloc((mpsi+1)*sizeof(real), IDEAL_ALIGNMENT);

  /* Define poloidal grid */

  /* grid spacing */
  deltar = (a1-a0)/mpsi;
   
  /* grid shift associated with fieldline following coordinates */
  tdum = 2.0*pi*a1/mthetamax;
#pragma omp parallel for private(i,rad,mthetatmp,mt,q)
  for (i=0; i<mpsi+1; i++) {
    rad = a0 + deltar*i;
    mthetatmp = 2 * (int) (pi*rad/tdum+0.5);
    mt = (mthetatmp < mthetamax) ? mthetatmp : mthetamax;
    mtheta[i] = (mt > 2) ? mt : 2;
    deltat[i] = 2.0*pi/mtheta[i];
    difft[i] = 0.5/deltat[i];
    q = q0 + q1*rad/a + q2*rad*rad/(a*a);
    itran[i] = (int) (mtheta[i]/q + 0.5);
    qtinv[i] = ((real) mtheta[i])/itran[i];
    qtinv[i] = 1.0/qtinv[i];
    itran[i] = itran[i] - mtheta[i] * (itran[i]/mtheta[i]);

    rtemi[i] = 1.0;
    zonali[i] = 0.0;
    rden[i]  = 1.0;
    temp[i]  = 0.0;
    dtemp[i] = 0.0;
    phip00[i] = 0.0;
    den00[i] = 0.0;
    pmarki[i] = 0.0;
  }

  mtdiag = (mthetamax/mzetamax)*mzetamax;

  for (i=0; i<mflux; i++) {
    pfluxpsi[i] = 0.0;
    vdrtmp[i] = 0.0;
  }

  igrid[0] = 1;
  for (i=1; i<mpsi+1; i++) {
    igrid[i] = igrid[i-1] + mtheta[i-1] + 1;
  }

  /* Number of grid on a poloidal plane */
  mgrid = 0;
  for (i=0; i<mpsi+1; i++) {
    mgrid += mtheta[i] + 1;
  }

  //mi_local = micell * (mgrid-mpsi) * mzeta;
  mi_local  = micell;
  mi_local *= (mgrid-mpsi);
  mi_local *= mzeta;
  mi = mi_local/npartdom;
  mimax = 2*mi;
  if( particle_domain_location == npartdom-1 ) mi += mi_local % npartdom;

  delr = deltar/gyroradius;
  /* fprintf(stderr, "gyro radius: %lf, deltar: %lf\n", 
     gyroradius, deltar);
  */
  if (mype == 0) 
    fprintf(stderr, "mi: %d, mi_total: %ld, mgrid: %d, mzetamax: %d\n", mi,
	    ((long) mi)*((long) npartdom*ntoroidal), mgrid, mzetamax);

  a_nover_in = sqrt((real)myrank_radiald/(real)nproc_radiald*(a1*a1-a0*a0)+a0*a0);
  a_nover_out = sqrt((real)(myrank_radiald+1)/(real)nproc_radiald*(a1*a1-a0*a0)+a0*a0);
    
  /* pick inner/outer valid (working) radiuss, used to determine valid region */
  a_valid_in = -1.0;                                                                              
  a_valid_out = -1.0;
  for (i=0; i<mpsi+1; i++){
    rad = a0 + deltar*i;
    /* support of local geometric partition */
    if (rad - a_nover_in > -1.0e-15 && a_valid_in == -1.0) {
      a_valid_in = rad;
      ipsi_in = i; // start of domain
      if (i != 0 ) ipsi_in = i - 1;
      if (i != 0 ) a_valid_in = rad - deltar;
    }

    /* find end of local valid region (in mesh) */
    if (a_valid_out == -1.0 && rad - a_nover_out > -1.0e-15){
      a_valid_out = rad;
      ipsi_out = i; // end of working (valid) domain
    }
  }

  if (myrank_radiald==nproc_radiald-1&& ipsi_out != mpsi){
    fprintf(stderr, "ipsi_out %d /= mpsi %d\n", ipsi_out, mpsi);
#if USE_MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#else
    exit(1);
#endif
  }

  ipsi_valid_in = ipsi_in;
  ipsi_valid_out = ipsi_out;

  /* create non-overlapping partition with (beginning of ) hear of domain */
  ri_pe = (int *) _mm_malloc((mpsi+1)*sizeof(int), IDEAL_ALIGNMENT);
  tempi = (int *) _mm_malloc((mpsi+1)*sizeof(int), IDEAL_ALIGNMENT);
  for (i=0; i<mpsi+1; i++){
    tempi[i] = nproc_partd + 1;
    if (i>=ipsi_valid_in&&i<=ipsi_valid_out)
      tempi[i] = myrank_partd;
  }

#if USE_MPI
  MPI_Allreduce(tempi, ri_pe, mpsi+1, MPI_INT, MPI_MIN, *partd_comm);
#else
  for (i=0; i<mpsi+1; i++){
    ri_pe[i] = tempi[i];
  }
#endif

  /* set non-overlap region : ipsi_nover_in, ipsi_nover_out */
  ipsi_nover_in = -1;
  ipsi_nover_out = -1;
  for (i=ipsi_valid_in; i<ipsi_valid_out+1; i++){
    if (ri_pe[i] == myrank_partd && ipsi_nover_in == -1) ipsi_nover_in = i;
    if (ri_pe[i] != myrank_partd && ipsi_nover_in != -1 && ipsi_nover_out == -1)
      ipsi_nover_out = i - 1;
  }

  if (ipsi_nover_out == -1) ipsi_nover_out = ipsi_valid_out;
  if (ipsi_nover_in == -1) {
    ipsi_nover_out = 0;
    ipsi_nover_in = 1;
  }
  igrid_nover_in = igrid[ipsi_nover_in];
  igrid_nover_out = igrid[ipsi_nover_out] + mtheta[ipsi_nover_out];
    
  /* check for correct */
  arr = (int *) _mm_malloc((mpsi+1)*sizeof(int), IDEAL_ALIGNMENT);
  arr2 = (int *) _mm_malloc((mpsi+1)*sizeof(int), IDEAL_ALIGNMENT);
  for (i=0; i<mpsi+1; i++){
    arr[i] = 0;
    if (i>=ipsi_nover_in&&i<=ipsi_nover_out)
      arr[i] = 1;
  }
#if USE_MPI
  MPI_Allreduce(arr, arr2, mpsi+1, MPI_INT, MPI_SUM, *partd_comm);
#else
  for (i=0; i<mpsi+1; i++){
    arr2[i] = arr[i];
  }
#endif
  for (i=0; i<mpsi+1; i++){
    if (arr2[i]!=1){
      fprintf(stderr, "failed radial partition: use more radial grid cells (or less procs)");
#if USE_MPI
      MPI_Abort(MPI_COMM_WORLD,1);
#else
      exit(1);
#endif
    }
  }
    
  /* pick inner/outer radius for ghost cells */
  n_rad_buf = 8.0*gyroradius/deltar + 0.5;
  if (n_rad_buf<3) n_rad_buf = 3; // need buffer for (large) stencil in grad op
  rho_max = n_rad_buf*deltar;
  if (mype==0) printf("Using n_rad_buf=%d buffer cells: gyroradius/deltar=%f\n",
		      n_rad_buf, gyroradius/deltar);
  a_in = a_valid_in - n_rad_buf*deltar;
  ipsi_in = ipsi_valid_in - n_rad_buf;
  if (a_in <= a0) {
    a_in = a0;
    ipsi_in = 0;
  }
  a_out = a_valid_out + n_rad_buf*deltar;
  ipsi_out = ipsi_valid_out + n_rad_buf;
  if (a_out >= a1) {
    a_out = a1;
    ipsi_out = mpsi;
  }

  if (ipsi_in < 0 || ipsi_out > mpsi){
    fprintf(stderr, "mype=%d: in set_radial_decomp:ipsi_in < 0 or ipsi_out > mpsi\n",
	    mype);
#if USE_MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#else
    exit(1);
#endif
  }

  igrid_valid_in = igrid[ipsi_valid_in];
  igrid_valid_out = igrid[ipsi_valid_out] + mtheta[ipsi_valid_out];
    
  /* cache bounds of local grid array */
  igrid_in = igrid[ipsi_in];
  igrid_out = igrid[ipsi_out] + mtheta[ipsi_out];

  nloc_over = igrid_out-igrid_in+1;
  nloc_nover = igrid_nover_out-igrid_nover_in+1;

  /* create ghost points for updating ghost cell information in grid based subroutine */
  /* create non-overlapping partition with (beginning of ) hear of domain */
  ri_pe2 = (int *) _mm_malloc((mpsi+1)*sizeof(int), IDEAL_ALIGNMENT);
  for (i=0; i<mpsi+1; i++){
    tempi[i] = nproc_radiald + 1;
    if (i>=ipsi_valid_in&&i<=ipsi_valid_out)
      tempi[i] = myrank_radiald;
  }

#if USE_MPI
  MPI_Allreduce(tempi, ri_pe2, mpsi+1, MPI_INT, MPI_MIN, *radiald_comm);
#else
  for (i=0; i<mpsi+1; i++){
    ri_pe2[i] = tempi[i];
  }
#endif

  ghost_start_local = (int *) _mm_malloc(nproc_radiald*sizeof(int), IDEAL_ALIGNMENT);
  ghost_end_local = (int *) _mm_malloc(nproc_radiald*sizeof(int), IDEAL_ALIGNMENT);
  ghost_start_global = (int *) _mm_malloc(nproc_radiald*nproc_radiald*sizeof(int), IDEAL_ALIGNMENT);
  ghost_end_global = (int *) _mm_malloc(nproc_radiald*nproc_radiald*sizeof(int), IDEAL_ALIGNMENT); 

  for (i=0; i<nproc_radiald; i++){
    ghost_start_local[i] = 0;
    ghost_end_local[i] = 0;
  }

  for (i=ipsi_in; i<ipsi_out+1; i++){
    int proc_num = ri_pe2[i];
    if (proc_num!=myrank_radiald&&ghost_end_local[proc_num]!=0){
      ghost_end_local[proc_num] = i;
    }
    if (proc_num!=myrank_radiald &&
	ghost_start_local[proc_num]==0 &&
	ghost_end_local[proc_num]==0){
      ghost_start_local[proc_num] = i;	
      ghost_end_local[proc_num] = i;
    }
  }

#if USE_MPI
  MPI_Allgather(ghost_start_local, nproc_radiald, MPI_INT,
		ghost_start_global, nproc_radiald, MPI_INT, *radiald_comm);
  MPI_Allgather(ghost_end_local, nproc_radiald, MPI_INT,
		ghost_end_global, nproc_radiald, MPI_INT, *radiald_comm);
#else
  for (i=0; i<nproc_radiald; i++){
    ghost_start_global[i] = ghost_start_local[i];
  }
#endif    

  nghost_comm_num = 0;
  for (i=0; i<nproc_radiald; i++){
    if (i!=myrank_radiald){
      if (ghost_start_global[i*nproc_radiald+myrank_radiald]!=0){
	nghost_comm_num += 1;
      }
    }
  }
  /* nghost_comm_list: ghost cells from other PEs to non ghost cells of my PE*/
  nghost_comm_list = (int *) _mm_malloc(nghost_comm_num*sizeof(int), IDEAL_ALIGNMENT);
  nghost_start = (int *) _mm_malloc(nghost_comm_num*sizeof(int), IDEAL_ALIGNMENT);
  nghost_end = (int *) _mm_malloc(nghost_comm_num*sizeof(int), IDEAL_ALIGNMENT);
    
  i_loc = 0;
  nghost_bufsize = 0;
  for (i=0; i<nproc_radiald; i++){
    if (i!=myrank_radiald){
      if (ghost_start_global[i*nproc_radiald+myrank_radiald]!=0){
	nghost_comm_list[i_loc] = i;
	nghost_start[i_loc] = ghost_start_global[i*nproc_radiald+myrank_radiald];
	nghost_end[i_loc] = ghost_end_global[i*nproc_radiald+myrank_radiald];
	nghost_bufsize += igrid[nghost_end[i_loc]]+mtheta[nghost_end[i_loc]]-igrid[nghost_start[i_loc]]+1;
	i_loc++;
      }
    }
  }
  assert(nghost_comm_num==i_loc);
  nghost_sendrecvbuf  = (real *) _mm_malloc(nghost_bufsize*(mzeta+1)*3*sizeof(real), IDEAL_ALIGNMENT);

  ghost_comm_num = 0;
  for (i=0; i<nproc_radiald; i++){
    if (i!=myrank_radiald){
      if (ghost_start_global[myrank_radiald*nproc_radiald+i]!=0){
	ghost_comm_num++;
      }
    }
  }
  /* ghost_comm_list: ghost cells of my PEs to non ghost cells of other PEs */
  ghost_comm_list = (int *) _mm_malloc(ghost_comm_num*sizeof(int), IDEAL_ALIGNMENT);
  ghost_start = (int *) _mm_malloc(ghost_comm_num*sizeof(int), IDEAL_ALIGNMENT);
  ghost_end = (int *) _mm_malloc(ghost_comm_num*sizeof(int), IDEAL_ALIGNMENT);

  i_loc = 0;
  ghost_bufsize = 0;
  for(i=0; i<nproc_radiald; i++){
    if (i!=myrank_radiald){
      if (ghost_start_global[myrank_radiald*nproc_radiald+i]!=0){
	ghost_comm_list[i_loc] = i;
	ghost_start[i_loc] = ghost_start_global[myrank_radiald*nproc_radiald+i];
	ghost_end[i_loc] = ghost_end_global[myrank_radiald*nproc_radiald+i];
	ghost_bufsize += igrid[ghost_end[i_loc]]+mtheta[ghost_end[i_loc]]-igrid[ghost_start[i_loc]]+1;
	i_loc++;
      }
    }
  }
  assert(ghost_comm_num==i_loc);
  ghost_sendrecvbuf = (real *) _mm_malloc(ghost_bufsize*(mzeta+1)*3*sizeof(real), IDEAL_ALIGNMENT);

  ipsi_nover_in_radiald = -1;
  ipsi_nover_out_radiald = -1;
  for (i=ipsi_valid_in; i<ipsi_valid_out+1; i++){
    if (ri_pe2[i] == myrank_radiald && ipsi_nover_in_radiald == -1) ipsi_nover_in_radiald = i;
    if (ri_pe2[i] != myrank_radiald && ipsi_nover_in_radiald != -1 && ipsi_nover_out_radiald == -1)
      ipsi_nover_out_radiald = i - 1;
  }

  if (ipsi_nover_out_radiald == -1) ipsi_nover_out_radiald = ipsi_valid_out;
  if (ipsi_nover_in_radiald == -1) {
    ipsi_nover_out_radiald = 0;
    ipsi_nover_in_radiald = 1;
  }
  igrid_nover_in_radiald = igrid[ipsi_nover_in_radiald];
  igrid_nover_out_radiald = igrid[ipsi_nover_out_radiald] + mtheta[ipsi_nover_out_radiald];
    
  pgyro    = (real *) _mm_malloc(4*nloc_over*sizeof(real), IDEAL_ALIGNMENT);
  tgyro    = (real *) _mm_malloc(4*nloc_over*sizeof(real), IDEAL_ALIGNMENT);
  markeri  = (real *) _mm_malloc(nloc_over * mzeta * sizeof(real), IDEAL_ALIGNMENT);

  densityi = (wreal *) _mm_malloc((mzeta+1)*nloc_over*sizeof(wreal), IDEAL_ALIGNMENT);
  densityi_local = (wreal *) _mm_malloc((mzeta+1)*nloc_over*nthreads*sizeof(wreal), IDEAL_ALIGNMENT);
    
  dnitmp   = (wreal *) _mm_malloc((mzeta+1)*nloc_over*sizeof(wreal), IDEAL_ALIGNMENT);
  recvr    = (wreal *) _mm_malloc(nloc_over * sizeof(wreal), IDEAL_ALIGNMENT);
  sendl    = (wreal *) _mm_malloc(nloc_over * sizeof(wreal), IDEAL_ALIGNMENT);
   
  int mumax2 = 8;
  real delvperp = umax/(real)mumax2;
  int m2pi=16;

  int mpsi_loc = ipsi_valid_out-ipsi_valid_in+1;
  assert(mpsi_loc>0);

  psi_count = (int *) _mm_malloc(mpsi_loc*nthreads*sizeof(int), IDEAL_ALIGNMENT);
  psi_offsets = (int *) _mm_malloc(mpsi_loc*nthreads*sizeof(int), IDEAL_ALIGNMENT);

  phi      = (real *) _mm_malloc((mzeta+1)*nloc_over*sizeof(real), IDEAL_ALIGNMENT);
  evector  = (real *) _mm_malloc(3*(mzeta+1)*nloc_over*sizeof(real), IDEAL_ALIGNMENT);

  jtp1     = (int *) _mm_malloc(2 * nloc_over * mzeta * sizeof(int), IDEAL_ALIGNMENT);
  jtp2     = (int *) _mm_malloc(2 * nloc_over * mzeta * sizeof(int), IDEAL_ALIGNMENT);
  wtp1     = (real *) _mm_malloc(2 * nloc_over * mzeta * sizeof(real), IDEAL_ALIGNMENT);
  wtp2     = (real *) _mm_malloc(2 * nloc_over * mzeta * sizeof(real), IDEAL_ALIGNMENT);
  dtemper  = (real *) _mm_malloc(nloc_over * mzeta * sizeof(real), IDEAL_ALIGNMENT);
  heatflux = (real *) _mm_malloc(nloc_over * mzeta * sizeof(real), IDEAL_ALIGNMENT);

  drdpa  = (real *) _mm_malloc(nloc_over*sizeof(real), IDEAL_ALIGNMENT);
  diffta = (real *) _mm_malloc(nloc_over*sizeof(real), IDEAL_ALIGNMENT);
  idx1a  = (int *)  _mm_malloc(nloc_over*sizeof(int), IDEAL_ALIGNMENT);
  idx2a  = (int *)  _mm_malloc(nloc_over*sizeof(int), IDEAL_ALIGNMENT);
  recvl_index = (int *)  _mm_malloc(nloc_over*sizeof(int), IDEAL_ALIGNMENT);
  recvr_index= (int *)  _mm_malloc(nloc_over*sizeof(int), IDEAL_ALIGNMENT);

  sendrf     = (real *) _mm_malloc(nloc_over*sizeof(real), IDEAL_ALIGNMENT);
  sendlf     = (real *) _mm_malloc(nloc_over*sizeof(real), IDEAL_ALIGNMENT);
  recvrf     = (real *) _mm_malloc(nloc_over*sizeof(real), IDEAL_ALIGNMENT);
  recvlf     = (real *) _mm_malloc(nloc_over*sizeof(real), IDEAL_ALIGNMENT);
  sendrsf    = (real *) _mm_malloc(3*nloc_over*sizeof(real), IDEAL_ALIGNMENT);
  recvlsf    = (real *) _mm_malloc(3*nloc_over*sizeof(real), IDEAL_ALIGNMENT);
    
  perr       = (real *) _mm_malloc(nloc_over * sizeof(real), IDEAL_ALIGNMENT);
  ptilde     = (real *) _mm_malloc(nloc_over * sizeof(real), IDEAL_ALIGNMENT);
  dentmp     = (real *) _mm_malloc(nloc_over * sizeof(real), IDEAL_ALIGNMENT);
  phitmp     = (real *) _mm_malloc(nloc_over * sizeof(real), IDEAL_ALIGNMENT);
  phitmps    = (real *) _mm_malloc(nloc_over * (mzeta+1) * sizeof(real), IDEAL_ALIGNMENT);
   
  z0       = (real *) _mm_malloc(mimax*sizeof(real), IDEAL_ALIGNMENT);
  z1       = (real *) _mm_malloc(mimax*sizeof(real), IDEAL_ALIGNMENT);
  z2       = (real *) _mm_malloc(mimax*sizeof(real), IDEAL_ALIGNMENT);
  z3       = (real *) _mm_malloc(mimax*sizeof(real), IDEAL_ALIGNMENT);
  z4       = (real *) _mm_malloc(mimax*sizeof(real), IDEAL_ALIGNMENT);
  z5       = (real *) _mm_malloc(mimax*sizeof(real), IDEAL_ALIGNMENT);
  z00      = (real *) _mm_malloc(mimax*sizeof(real), IDEAL_ALIGNMENT);
  z01      = (real *) _mm_malloc(mimax*sizeof(real), IDEAL_ALIGNMENT);
  z02      = (real *) _mm_malloc(mimax*sizeof(real), IDEAL_ALIGNMENT);
  z03      = (real *) _mm_malloc(mimax*sizeof(real), IDEAL_ALIGNMENT);
  z04      = (real *) _mm_malloc(mimax*sizeof(real), IDEAL_ALIGNMENT);
  z05      = (real *) _mm_malloc(mimax*sizeof(real), IDEAL_ALIGNMENT);
  ztmp     = (real *) _mm_malloc(mimax*sizeof(real), IDEAL_ALIGNMENT);
  ztmp2    = (real *) _mm_malloc(mimax*sizeof(real), IDEAL_ALIGNMENT);
  wzion    = (real *) _mm_malloc(mimax*sizeof(real), IDEAL_ALIGNMENT);
  wpion    = (real *) _mm_malloc(4*mimax*sizeof(real), IDEAL_ALIGNMENT);
  wtion0   = (real *) _mm_malloc(4*mimax*sizeof(real), IDEAL_ALIGNMENT);
  wtion1   = (real *) _mm_malloc(4*mimax*sizeof(real), IDEAL_ALIGNMENT);
  kzion    = (int  *) _mm_malloc(mimax*sizeof(int),  IDEAL_ALIGNMENT);
  jtion0   = (int  *) _mm_malloc(4*mimax*sizeof(int), IDEAL_ALIGNMENT);
  jtion1   = (int  *) _mm_malloc(4*mimax*sizeof(int), IDEAL_ALIGNMENT);
  kzi      = (int  *) _mm_malloc(mimax*sizeof(int), IDEAL_ALIGNMENT);

  assert(pgyro != NULL); assert(tgyro != NULL); 
  assert(markeri != NULL);
  assert(densityi != NULL);
  assert(densityi_local != NULL);

  assert(phi != NULL); assert(evector != NULL);
  assert(jtp1 != NULL); assert(jtp2 != NULL); 
  assert(wtp1 != NULL); assert(wtp2 != NULL);
  assert(dtemper != NULL); assert(heatflux != NULL);

  assert(z0 != NULL); assert(z1 != NULL); assert(z2 != NULL); 
  assert(z3 != NULL); assert(z4 != NULL); assert(z5 != NULL);
  assert(z00 != NULL); assert(z01 != NULL); assert(z02 != NULL);
  assert(z03 != NULL); assert(z04 != NULL); assert(z05 != NULL);
  assert(ztmp != NULL); assert(ztmp2 != NULL); assert(kzi != NULL);

  assert(wzion != NULL); assert(wpion != NULL); 
  assert(wtion0 != NULL); assert(wtion1 != NULL);
  assert(kzion != NULL); assert(jtion0 != NULL); assert(jtion1 != NULL);
    
  assert(drdpa != NULL); assert(diffta != NULL);
  assert(idx1a != NULL); assert(idx2a != NULL);
  assert(recvl_index != NULL); assert(recvr_index != NULL);

  assert(sendrsf != NULL); assert(recvlsf != NULL);
  assert(sendrf != NULL); assert(sendlf != NULL);
  assert(recvrf != NULL); assert(recvlf != NULL);
  assert(perr != NULL); assert(ptilde != NULL);
  assert(phitmp != NULL); assert(dentmp != NULL);
  assert(phitmps != NULL);
 
#pragma omp parallel for private(i) 
  for (i=0; i<mi; i++) {
    z0[i] = 0.0; z1[i] = 0.0; z2[i] = 0.0; z3[i] = 0.0; 
    z4[i] = 0.0; z5[i] = 0.0;
    z00[i] = 0.0; z01[i] = 0.0; z02[i] = 0.0; z03[i] = 0.0; 
    z04[i] = 0.0; z05[i] = 0.0; ztmp[i] = 0.0; ztmp2[i] = 0.0;
    wpion[4*i] = 0.0; wpion[4*i+1] = 0.0; 
    wpion[4*i+2] = 0.0; wpion[4*i+3] = 0.0;
    wtion0[4*i] = 0.0; wtion0[4*i+1] = 0.0; 
    wtion0[4*i+2] = 0.0; wtion0[4*i+3] = 0.0;
    wtion1[4*i] = 0.0; wtion1[4*i+1] = 0.0; 
    wtion1[4*i+2] = 0.0; wtion1[4*i+3] = 0.0;
    kzion[i] = 0; wzion[i] = 0.0;
    jtion0[4*i] = 0; jtion0[4*i+1] = 0; 
    jtion0[4*i+2] = 0; jtion0[4*i+3] = 0;
    jtion1[4*i] = 0; jtion1[4*i+1] = 0; 
    jtion1[4*i+2] = 0; jtion1[4*i+3] = 0;
  }

#pragma omp parallel for private(i)
  for (i=0; i<mimax; i++) {
    kzi[i] = 0;
  }

#pragma omp parallel for private(i)
  for (i=0; i<(mzeta+1)*nloc_over; i++) {
    phi[i] = 0.0;
    densityi[i] = 0.0;
    dnitmp[i] = 0.0;
  }

#pragma omp parallel for private(i)
  for (i=0; i<nthreads*(mzeta+1)*nloc_over; i++) {
    densityi_local[i] = 0.0;
  }

  for (i=0; i<nthreads*mpsi_loc; i++){
    psi_count[i] = 0;
    psi_offsets[i] = 0;
  }

#pragma omp parallel for private(i)
  for (i=0; i<2*mzeta*nloc_over; i++) {
    jtp1[i] = 0;
    jtp2[i] = 0;
    wtp1[i] = 0.0;
    wtp2[i] = 0.0;
  }

#pragma omp parallel for private(i)
  for (i=0; i<mzeta*nloc_over; i++) {
    dtemper[i] = 0.0;
    heatflux[i] = 0.0;
    markeri[i] = 0.0;
  }

#pragma omp parallel for private(i)
  for (i=0; i<nloc_over; i++) {
    pgyro[4*i]   = 0.0;
    pgyro[4*i+1] = 0.0;
    pgyro[4*i+2] = 0.0;
    pgyro[4*i+3] = 0.0;
    tgyro[4*i]   = 0.0;
    tgyro[4*i+1] = 0.0;
    tgyro[4*i+2] = 0.0;
    tgyro[4*i+3] = 0.0;
    sendrf[i]    = 0.0;
    sendlf[i]    = 0.0;
    recvrf[i]    = 0.0;
    recvlf[i]    = 0.0;
    perr[i]      = 0.0;
    ptilde[i]    = 0.0;
    phitmp[i]    = 0.0;
    dentmp[i]    = 0.0;

    drdpa[i]     = 0.0;
    diffta[i]    = 0.0;
    idx1a[i]     = 0;
    idx2a[i]     = 0;
  }

#pragma omp parallel for private(i)
  for (i=0; i<3*(mzeta + 1)*nloc_over; i++) {
    evector[i] = 0.0;
  }

#pragma omp parallel for private(i)
  for (i=0; i<(mzeta + 1)*nloc_over; i++) {
    phitmps[i] = 0.0;
  }

#pragma omp parallel for private(i)
  for (i=0; i<3*nloc_over; i++) {
    sendrsf[i] = 0.0;
    recvlsf[i] = 0.0;
  }

  /* # of marker per grid, Jacobian */
  for (i=ipsi_in; i<ipsi_out+1; i++) {
    r = a0+deltar*i;
    for (j=1; j<mtheta[i]+1; j++) {
      ij = igrid[i]+j;
      for (k=1; k<mzeta+1; k++) {
	zdum = zetamin + k*deltaz;
	tdum = j*deltat[i]+zdum*qtinv[i];
	markeri[(ij-igrid_in)*mzeta+k-1] = pow((1.0 + r*cos(tdum)), 2);
	if (i>=ipsi_nover_in&&i<=ipsi_nover_out)
	  pmarki[i] = pmarki[i]+markeri[(ij-igrid_in)*mzeta+k-1];
      }
    }
  }
    
  adum = (real *)_mm_malloc((mpsi+1)*sizeof(real), IDEAL_ALIGNMENT);
  assert(adum != NULL);
#if USE_MPI
  MPI_Allreduce(pmarki, adum, mpsi+1, MPI_MYREAL, MPI_SUM, *partd_comm);
#else
  for (i=0; i<mpsi+1; i++){
    adum[i] = pmarki[i];
  }
#endif
  for (i=0; i<mpsi+1; i++){
    pmarki[i] = adum[i];
  }

  for (i=ipsi_in; i<ipsi_out+1; i++){
    r = a0 + deltar*i;
    rmax = min(a1,r+0.5*deltar);
    rmin = max(a0,r-0.5*deltar);
    tdum = mi*npartdom*(rmax*rmax-rmin*rmin)/(a1*a1-a0*a0);
      
    for (j=1; j<mtheta[i]+1; j++) {
      ij = igrid[i]+j;
      for (k=1; k<mzeta+1; k++) {
	markeri[(ij-igrid_in)*mzeta+k-1] = tdum*markeri[(ij-igrid_in)*mzeta+k-1]/pmarki[i];
	markeri[(ij-igrid_in)*mzeta+k-1] = 1.0/markeri[(ij-igrid_in)*mzeta+k-1]; // avoid divide operation
      }
    }

    for (k=1; k<mzeta+1; k++) {
      markeri[(igrid[i]-igrid_in)*mzeta+k-1] = markeri[(igrid[i]+mtheta[i]-igrid_in)*mzeta+k-1];
    }
  }

  for (i=0; i<mpsi+1; i++){
    r = a0 + deltar*i;
    rmax = min(a1,r+0.5*deltar);
    rmin = max(a0,r-0.5*deltar);
    tdum = mi*npartdom*(rmax*rmax-rmin*rmin)/(a1*a1-a0*a0);  

    pmarki[i] = 1.0/(tdum*ntoroidal); 
  }
    
  for (i=ipsi_in; i<ipsi_out+1; i++) {
    real rhoi;
    rad = a0 + deltar * i;
    for (j=0; j<mtheta[i]+1; j++) {
      ij = igrid[i] + j - igrid_in;
      tdum = deltat[i]*j;
      q = q0 + q1*rad/a + q2*rad*rad/(a*a);
      b = 1.0/(1.0+rad*cos(tdum));
      dtheta_dx = 1.0/rad;

      /* first two points perpendicular to field line on poloidal surface
       */
      rhoi = sqrt(2.0/b)*gyroradius;
      pgyro[0+4*ij] = -rhoi;
      pgyro[1+4*ij] = rhoi;
        
      /* Non-orthogonality between psi and theta */
      tgyro[0+4*ij] = 0;
      tgyro[1+4*ij] = 0;
            
      /* the other two points tangential to field line */
      tgyro[2+4*ij] = -rhoi * dtheta_dx;
      tgyro[3+4*ij] =  rhoi * dtheta_dx;
            
      pgyro[2+4*ij] = rhoi*0.5*rhoi/rad;
      pgyro[3+4*ij] = rhoi*0.5*rhoi/rad;
    }
  }

  for (i=ipsi_nover_in_radiald; i<ipsi_nover_out_radiald+1; i++){
    real r = a0 + deltar*((real) i);
    real drdp = 1.0/r;
    for (int j=0; j<mtheta[i]+1; j++) {
      int ij = igrid[i]+j-igrid_in;
      //assert(ij>=0);
      //assert(ij<nloc_over);
      int jt = j+1-mtheta[i]*(j/mtheta[i]);
      drdpa[ij] = drdp;
      diffta[ij] = difft[i];
      idx1a[ij] = igrid[i]+jt-igrid_in;
      idx2a[ij] = igrid[i]+j-1-igrid_in;
    }
  }  

  for (int i=ipsi_in; i<ipsi_out+1; i++){
    int ii = igrid[i];
    int jt = mtheta[i];
    int iti = itran[i];

    if (myrank_toroidal == 0){
      for (int j=0; j<jt+1; j++){
	recvr_index[ii+j-igrid_in] = ii+j-igrid_in;
      }
      for (int j=ii; j<ii+iti+1; j++){
	recvl_index[j-igrid_in] = j+jt-iti-igrid_in;
      }
      for (int j=ii+iti+1; j<ii+jt+1; j++){
	recvl_index[j-igrid_in] = j-iti-igrid_in;
      }
    } else if (myrank_toroidal == (ntoroidal - 1)) {
      for (int j=0; j<jt+1; j++){
	recvl_index[ii+j-igrid_in] = ii+j-igrid_in;
      }
      for (int j=ii; j<ii+jt-iti+1; j++){
	recvr_index[j-igrid_in] = j+iti-igrid_in;
      }
      for (int j=ii+jt-iti+1; j<ii+jt+1; j++){
	recvr_index[j-igrid_in] = j+iti-jt-igrid_in;
      }
    } else {
      for (int j=0; j<jt+1; j++){
	recvl_index[ii+j-igrid_in] = ii+j-igrid_in;
	recvr_index[ii+j-igrid_in] = ii+j-igrid_in;
      }
    }
  }

  /* initiate radial interpolation for grid */
  for (k=1; k<mzeta+1; k++) {
    zdum = zetamin + deltaz*k;
    for (i=ipsi_in; i<ipsi_out+1; i++) {
      for (ip=1; ip<3; ip++) {
	indp = min(ipsi_out, i+ip);
	indt = max(ipsi_in, i-ip);
	for (j=0; j<mtheta[i]+1; j++) {
	  ij = igrid[i] + j;
                    
	  /* upward */
	  tdum = (j * deltat[i] 
		  + zdum*(qtinv[i]-qtinv[indp]))/deltat[indp];
	  jt = ((int) floor(tdum));
	  wt = tdum - ((real) jt);
	  jt = (jt+mtheta[indp]) % mtheta[indp];
	  if (ip == 1) {
	    wtp1[(k-1)*2*nloc_over + (ij-igrid_in)*2] = wt;
	    jtp1[(k-1)*2*nloc_over + (ij-igrid_in)*2] = igrid[indp] + jt;
	  } else {
	    wtp2[(k-1)*2*nloc_over + (ij-igrid_in)*2] = wt;
	    jtp2[(k-1)*2*nloc_over + (ij-igrid_in)*2] = igrid[indp] + jt;
	  }

	  /* downward */
	  tdum = (j * deltat[i] + zdum*(qtinv[i]-qtinv[indt]))/deltat[indt];
	  jt = ((int) floor(tdum));
	  wt = tdum - ((real) jt);
	  jt = (jt+mtheta[indt]) % mtheta[indt];
	  if (ip == 1) {
	    wtp1[(k-1)*2*nloc_over + (ij-igrid_in)*2 + 1] = wt;
	    jtp1[(k-1)*2*nloc_over + (ij-igrid_in)*2 + 1] = igrid[indt] + jt;
	  } else {
	    wtp2[(k-1)*2*nloc_over + (ij-igrid_in)*2 + 1] = wt;
	    jtp2[(k-1)*2*nloc_over + (ij-igrid_in)*2 + 1] = igrid[indt] + jt;
	  }

	}

      }
    }
  }

  /* initialize diagnostic data */
  diagnosis_data->ptracer[0] = diagnosis_data->ptracer[1] = diagnosis_data->ptracer[2] = diagnosis_data->ptracer[3] = 0.0;
  diagnosis_data->eflux_average = 0.0;
  int num_mode = 8; int m_poloidal = 9;

  scalar_data = (real *) _mm_malloc(15*sizeof(real), IDEAL_ALIGNMENT);
  eflux = (real *) _mm_malloc(mflux*sizeof(real), IDEAL_ALIGNMENT);
  rmarker = (real *) _mm_malloc(mflux*sizeof(real), IDEAL_ALIGNMENT);
  dmark = (real *) _mm_malloc(mflux*sizeof(real), IDEAL_ALIGNMENT);
  dden = (real *) _mm_malloc(mflux*sizeof(real), IDEAL_ALIGNMENT);
  rdtemi = (real *) _mm_malloc(mflux*sizeof(real), IDEAL_ALIGNMENT);
  rdteme = (real *) _mm_malloc(mflux*sizeof(real), IDEAL_ALIGNMENT);
  flux_data = (real *) _mm_malloc(6*mflux*sizeof(real), IDEAL_ALIGNMENT);
  amp_mode = (real *) _mm_malloc(4*num_mode*sizeof(real), IDEAL_ALIGNMENT);
  eigenmode = (real *) _mm_malloc(num_mode*m_poloidal*mpsi*sizeof(real), IDEAL_ALIGNMENT);
  nmode = (int *) _mm_malloc(num_mode*sizeof(int), IDEAL_ALIGNMENT);
  mmode = (int *) _mm_malloc(num_mode*sizeof(int), IDEAL_ALIGNMENT);

  assert(scalar_data!=NULL); assert(flux_data!=NULL);
  assert(eflux!=NULL); assert(rmarker!=NULL); assert(dmark!=NULL);assert(dden!=NULL);assert(rdtemi!=NULL); assert(rdteme!=NULL);
  assert(amp_mode!=NULL); assert(eigenmode!=NULL); assert(nmode!=NULL); assert(mmode!=NULL);
  for (int i=0; i<15; i++){
    scalar_data[i] = 0.0;
  }

  for (int i=0; i<mflux; i++)
    {
      eflux[i] = 0.0;
      rmarker[i] = 0.0;
      dmark[i] = 0.0;
      dden[i] = 0.0;
      rdtemi[i] = 0.0;
      rdteme[i] = 0.0;
    }

  for (int i=0; i<6*mflux; i++){
    flux_data[i] = 0.0;
  }

  for (int i=0; i<num_mode; i++)
    {
      for (int j=0; j<4; j++)
	{
	  amp_mode[i*4+j] = 0.0;
	}

      for (int j=0; j<m_poloidal*mpsi; j++)
	{
	  eigenmode[i*m_poloidal*mpsi+j] = 0.0;
	}
    }

  nmode[0]=5;nmode[1]=7;nmode[2]=9;nmode[3]=11;nmode[4]=13;
  nmode[5]=15;nmode[6]=18;nmode[7]=20;
  mmode[0]=7;mmode[1]=10;mmode[2]=13;mmode[3]=15;mmode[4]=18;
  mmode[5]=21;mmode[6]=25;mmode[7]=28;

  /* load function */
  c0 = 2.515517;
  c1 = 0.802853;
  c2 = 0.010328;
  d1 = 1.432788;
  d2 = 0.189269;
  d3 = 0.001308;
#if UNIFORM_PARTICLE_LOADING
  rmi = 1.0/mi;
#else
  rmi = 1.0/(mi*npartdom);
#endif
  pi2_inv = 0.5/pi;
  delr =  1.0/deltar;
  w_initial = 1.0e-3;
   
  /* radial: uniformly distributed in r^2, later transform to psi */
  /*
    for (m=0; m<mi; m++) {
    #if UNIFORM_PARTICLE_LOADING
    z0[m] = sqrt(a0*a0+(((real) m+0.5)*(a1*a1-a0*a0)*rmi));
    #else
    z0[m] = sqrt(a0*a0+((real) (m+1+myrank_partd*mi)-0.5)*(a1*a1-a0*a0)*rmi);
    // this doesn't result in correct physics, fix it!
    //z0[m] = sqrt(a0*a0+((real) (myrank_radial_partd+m+1+myrank_radiald*nproc_radial_partd*mi)-0.5)*(a1*a1-a0*a0)*rmi);
    #endif
    }
  */

  /* Initialize RNG */
  rng_seed = (double *) malloc(6 * sizeof(double));
  assert(rng_seed != NULL);
  rng_seed[0] = 12345; rng_seed[1] = 12345; rng_seed[2] = 12345;
  rng_seed[3] = 12345; rng_seed[4] = 12345; rng_seed[5] = 12345;
  //rng_seed[0] = 54321; rng_seed[1] = 54321; rng_seed[2] = 54321;
  //rng_seed[3] = 54321; rng_seed[4] = 54321; rng_seed[5] = 54321;
  //rng_seed[0] = 21354; rng_seed[1] = 21354; rng_seed[2] = 21354;
  //rng_seed[3] = 21354; rng_seed[4] = 21354; rng_seed[5] = 21354;

    //BMA: change RNG initialization so particle generation is independent of concurrency
    //BMA RngStream_ParInit(mype, numberpe, 0, 1, 1, rng_seed);
    //BMA rs = RngStream_CreateStream("", rng_seed);
    RngStream_ParInit( 0, numberpe, 0, 1, 1, rng_seed );
    rs = RngStream_CreateStream("", rng_seed );
    long rng_offset;
    rng_offset = ( mi_local * toroidal_domain_location +
		   mi_local / npartdom * particle_domain_location ) * 6;
    //RngStream_AdvanceState( rs, 0, rng_offset );
    RngStream_FastAdvanceState( rs, 0, rng_offset );

  for (m=0; m<mi; m++) {
    z0[m] = (real) RngStream_RandU01(rs);  
    z1[m] = (real) RngStream_RandU01(rs);
    z2[m] = (real) RngStream_RandU01(rs);
    z3[m] = (real) RngStream_RandU01(rs);
    z4[m] = (real) RngStream_RandU01(rs);
    z5[m] = (real) RngStream_RandU01(rs);
  }

  //if( mype==0 ) printf("nproc_radiald: %d \n", nproc_radiald );
  //real volumn_frac = (a1*a1-a0*a0)/nproc_radiald;
  //real x2_l = a0*a0 +  myrank_radiald  * volumn_frac;
  real volumn_frac = (a1*a1-a0*a0);
  real x2_l = a0*a0;
  // randomly uniform distribute in r^2
  for (m=0; m<mi; m++) {
    z0[m] = sqrt(x2_l + z0[m]*volumn_frac);
  }
    
  for (m=0; m<mi; m++) {
    z1[m] = 2.0*pi*(z1[m]-0.5);
    z01[m] = z1[m];
  }

  for (i=0; i<10; i++) {
    for (m=0; m<mi; m++) {
      z1[m] = z01[m]-2.0*z0[m]*sin(z1[m]);
    }
  }

  for (m=0; m<mi; m++) {
    z1[m] = z1[m]*pi2_inv+10.0;
    z1[m] = 2.0 * pi * (z1[m] - (int) z1[m]);
  }

  for (m=0; m<mi; m++) {
    z3tmp  = z3[m];
    z3[m]  = z3[m]-0.5;
    z03[m] = ((z3[m] > 0) ? 1.0 : -1.0);
    z3[m]  = sqrt(max(1e-20,log(1.0/max(1e-20, pow(z3[m],2)))));
    z3[m]  = z3[m]-(c0+c1*z3[m]+c2*pow(z3[m],2))/
      (1.0+d1*z3[m]+d2*pow(z3[m],2)+d3*pow(z3[m],3));
    if (z3[m] > umax)
      z3[m] = z3tmp;
  }

  for (m=0; m<mi; m++) {
    z2[m] = zetamin + (zetamax-zetamin)*z2[m];
    z3[m] = z03[m]*min(umax,z3[m]);
    z4[m] =
      2.0*w_initial*(z4[m]-0.5)*(1.0+cos(z1[m]));
    zt1 = max(1.0e-20,z5[m]);
    /* (z5[m] > 1.0e-20) ? z5[m] : 1.0e-20; */
    zt2 = min(umax*umax,-log(zt1));
    /* (umax*umax < -log(zt1)) ? umax*umax : -log(zt1); */
    z5[m] = max(1.0e-20,zt2);
  }

  vthi = gyroradius*fabs(qion)/aion;
  for (m=0; m<mi; m++) {
    z00[m] = 1.0/(1.0+z0[m]*cos(z1[m]));
#if !SQRT_PRECOMPUTED
    z0[m] = 0.5 * z0[m]*z0[m];
#endif
    z3[m] = vthi*z3[m]*aion/(qion*z00[m]);
    z5[m] = sqrt(aion*vthi*vthi*z5[m]/z00[m]);
  }

  for (i=0; i<mi; i++) {
    z05[i] = 1.0;
  }

  delr=1.0/deltar; 
 
  delz=1.0/deltaz;
  smu_inv=sqrt(aion)/(fabs(qion)*gyroradius);
  pi2_inv=0.5/pi;

  for (i=0; i<mpsi+1; i++) {
    delt[i] = 2.0*pi/deltat[i];
  }

  for (i=0; i<nloc_over*(mzeta+1); i++) {
    densityi[i] = 0.0;
  }

  /*
    ntracer = 0;
    if (mype == 0) {
    z0[ntracer] = 0.5*pow((0.5*(a0+a1)),2);
    z1[ntracer] = 0.0;
    z2[ntracer] = 0.5*(zetamin+zetamax);
    z3[ntracer] = 0.5*vthi*aion/qion;
    z4[ntracer] = 0.0;
    z5[ntracer] = sqrt(aion*vthi*vthi);
    }
  */

  indexp = (int *)  _mm_malloc(65 * nloc_over * mzeta * sizeof(int),
			       IDEAL_ALIGNMENT);
  ring   = (real *) _mm_malloc(65 * nloc_over * mzeta * sizeof(real),
			       IDEAL_ALIGNMENT);
  nindex = (int *)  _mm_malloc(nloc_over  * mzeta * sizeof(int),
			       IDEAL_ALIGNMENT);
  assert(indexp != NULL); assert(ring != NULL); assert(nindex != NULL);

#pragma omp parallel for private(i)
  for (i=0; i<nloc_over*mzeta; i++) {
    nindex[i] = 0;
  }

#pragma omp parallel for private(i)
  for (i=0; i<nloc_over*mzeta*65; i++) {
    ring[i] = 0.0;
    indexp[i] = 0;
  }



#if USE_MPI
  /* Maximum number of particles sent/received in an iteration */
  recvbuf_size = (mi/5);
  recvbuf = (real *) _mm_malloc(recvbuf_size * 12 * sizeof(real), IDEAL_ALIGNMENT);
  assert(recvbuf != NULL);
  sendbuf_size = (mi/5);
  sendbuf = (real *) _mm_malloc(sendbuf_size * 12 * sizeof(real), IDEAL_ALIGNMENT);
  assert(sendbuf != NULL);
#endif

  /* 
*****************************************************************
*/

  global_params->micell = micell;
  global_params->mi = mi;
  global_params->mimax = mimax;
  global_params->mpsi = mpsi;
  global_params->mgrid = mgrid;
  global_params->zetamax = zetamax;
  global_params->zetamin = zetamin;
  global_params->mthetamax = mthetamax;
  global_params->mzetamax = mzetamax;
  global_params->mzeta = mzeta;
  global_params->smu_inv = smu_inv;
  global_params->mpsi_loc = mpsi_loc;
  global_params->m2pi = m2pi;
  global_params->mumax2 = mumax2;
  global_params->delvperp = delvperp;
  global_params->a = a;
  global_params->a0 = a0;
  global_params->a1 = a1;
  global_params->aion = aion;
  global_params->delr = delr;
  global_params->delz = delz;
  global_params->q0 = q0;
  global_params->q1 = q1;
  global_params->q2 = q2;
  global_params->qion = qion;
  global_params->rc = rc;
  global_params->r0 = r0;
  global_params->rw = rw;
  global_params->b0 = b0;

  global_params->temperature = temperature;
  global_params->umax = umax;
  global_params->pi2_inv = pi2_inv;
  global_params->mtdiag = mtdiag;
  global_params->deltar = deltar;
  global_params->mode00 = mode00;
  global_params->gyroradius = gyroradius;
  global_params->deltaz = deltaz;
  global_params->holecount = 0;
  global_params->miinit = mi;
  global_params->tstep = tstep;
  global_params->pi = pi;

  /* currently unused */
  global_params->hole_remove_freq = 1;

  /* The binning frequency can be modified */
  /* Binning every time step works best on most systems */
  global_params->radial_bin_freq = 1;

  if (mype == 0) {
    fprintf(stderr, "Binning and removing holes every %d steps\n", 
	    global_params->radial_bin_freq);
  }

  aux_particle_data->kzion = kzion;
  aux_particle_data->wzion = wzion;
  aux_particle_data->wpion = wpion;
  aux_particle_data->jtion0 = jtion0;
  aux_particle_data->jtion1 = jtion1;
  aux_particle_data->wtion0 = wtion0;
  aux_particle_data->wtion1 = wtion1;
  aux_particle_data->kzi = kzi;

  diagnosis_data->scalar_data = scalar_data;
  diagnosis_data->eflux = eflux;
  diagnosis_data->rmarker = rmarker;
  diagnosis_data->dmark = dmark;
  diagnosis_data->dden = dden;
  diagnosis_data->rdtemi = rdtemi;
  diagnosis_data->rdteme = rdteme;
  diagnosis_data->flux_data = flux_data;
  diagnosis_data->amp_mode = amp_mode;
  diagnosis_data->eigenmode = eigenmode;
  diagnosis_data->nmode = nmode;
  diagnosis_data->mmode = mmode;

  field_data->igrid = igrid;
  field_data->delt = delt;
  field_data->qtinv = qtinv; 
  field_data->mtheta = mtheta;
  field_data->pgyro = pgyro;
  field_data->tgyro = tgyro; 
  field_data->densityi = densityi; 
  field_data->densityi_local = densityi_local;
  field_data->evector = evector; 
  field_data->rtemi = rtemi;
  field_data->temp = temp;
  field_data->dtemp = dtemp;
  field_data->mmpsi = 22;
  field_data->itran = itran;
  field_data->jtp1 = jtp1;
  field_data->jtp2 = jtp2;
  field_data->rden = rden;
  field_data->zonali = zonali;
  field_data->dtemper = dtemper;
  field_data->wtp1 = wtp1;
  field_data->wtp2 = wtp2;
  field_data->indexp = indexp;
  field_data->nindex = nindex;
  field_data->ring = ring;
  field_data->recvr = recvr;
  field_data->sendl = sendl;
  field_data->dnitmp = dnitmp;
  field_data->heatflux = heatflux;
  field_data->hfluxpsi = hfluxpsi;
  field_data->pfluxpsi = pfluxpsi;
  field_data->vdrtmp = vdrtmp;
  field_data->difft = difft;
  field_data->phi = phi;
  field_data->pmarki = pmarki;
  field_data->den00 = den00;
  field_data->phi00 = phi00;
  field_data->phip00 = phip00;
  field_data->gradt = gradt;
  field_data->markeri = markeri;
  field_data->deltat = deltat;   
 
  field_data->drdpa = drdpa;
  field_data->diffta = diffta;
  field_data->idx1a = idx1a;
  field_data->idx2a = idx2a;
  field_data->recvl_index = recvl_index;
  field_data->recvr_index = recvr_index;

  field_data->sendrf = sendrf;
  field_data->recvrf = recvrf;
  field_data->sendlf = sendlf;
  field_data->recvlf = recvlf;
  field_data->sendrsf = sendrsf;
  field_data->recvlsf = recvlsf;

  field_data->perr = perr;
  field_data->ptilde = ptilde;
  field_data->phitmp = phitmp;
  field_data->phitmps = phitmps;
  field_data->dentmp = dentmp;

  particle_data->z0 = z0;
  particle_data->z1 = z1;
  particle_data->z2 = z2;
  particle_data->z3 = z3;
  particle_data->z4 = z4;
  particle_data->z5 = z5;
  particle_data->z00 = z00;
  particle_data->z01 = z01;
  particle_data->z02 = z02;
  particle_data->z03 = z03;
  particle_data->z04 = z04;
  particle_data->z05 = z05;
  particle_data->ztmp = ztmp;
  particle_data->ztmp2 = ztmp2;
  particle_data->psi_count = psi_count;
  particle_data->psi_offsets = psi_offsets;

  parallel_decomp->nthreads = nthreads;
  parallel_decomp->ntoroidal = ntoroidal;
  parallel_decomp->npartdom  = npartdom;
  parallel_decomp->nproc_partd = nproc_partd;
  parallel_decomp->myrank_partd = myrank_partd;
  parallel_decomp->nproc_toroidal = nproc_toroidal;
  parallel_decomp->myrank_toroidal = myrank_toroidal;
  parallel_decomp->left_pe = left_pe;
  parallel_decomp->right_pe = right_pe;
  parallel_decomp->toroidal_domain_location = toroidal_domain_location;
  parallel_decomp->particle_domain_location = particle_domain_location;

#if USE_MPI
  parallel_decomp->recvbuf = recvbuf;
  parallel_decomp->recvbuf_size = recvbuf_size;
  parallel_decomp->sendbuf = sendbuf;
  parallel_decomp->sendbuf_size = sendbuf_size;
#endif

  radial_decomp->ipsi_nover_in = ipsi_nover_in;
  radial_decomp->ipsi_nover_out = ipsi_nover_out;
  radial_decomp->ipsi_in = ipsi_in;
  radial_decomp->ipsi_out = ipsi_out;
  radial_decomp->ipsi_valid_in = ipsi_valid_in;
  radial_decomp->ipsi_valid_out = ipsi_valid_out;
  radial_decomp->igrid_in = igrid_in;
  radial_decomp->igrid_out = igrid_out;
  radial_decomp->igrid_nover_in = igrid_nover_in;
  radial_decomp->igrid_nover_out = igrid_nover_out;
  radial_decomp->nloc_nover = nloc_nover;
  radial_decomp->nloc_over = nloc_over;
  radial_decomp->ipsi_nover_in_radiald = ipsi_nover_in_radiald;
  radial_decomp->ipsi_nover_out_radiald = ipsi_nover_out_radiald;
  radial_decomp->igrid_nover_in_radiald = igrid_nover_in_radiald;
  radial_decomp->igrid_nover_out_radiald = igrid_nover_out_radiald;
  radial_decomp->a_nover_in = a_nover_in;
  radial_decomp->a_nover_out = a_nover_out;
  radial_decomp->a_valid_in = a_in;
  radial_decomp->a_valid_out = a_out;
  radial_decomp->rho_max = rho_max;
  radial_decomp->ri_pe = ri_pe;
  radial_decomp->ri_pe2 = ri_pe2;
  radial_decomp->npe_radiald = npe_radiald;
  radial_decomp->nradial_dom = nradial_dom;
  radial_decomp->myrank_radiald = myrank_radiald;
  radial_decomp->nproc_radiald = nproc_radiald;
  radial_decomp->myrank_radial_partd = myrank_radial_partd;
  radial_decomp->nproc_radial_partd = nproc_radial_partd;
  radial_decomp->left_radial_pe = left_radial_pe;
  radial_decomp->right_radial_pe = right_radial_pe;
  radial_decomp->radial_domain_location = radial_domain_location;
  radial_decomp->radial_part_domain_location = radial_part_domain_location;
  radial_decomp->ghost_comm_list = ghost_comm_list;
  radial_decomp->ghost_start = ghost_start;
  radial_decomp->ghost_end = ghost_end;
  radial_decomp->ghost_comm_num = ghost_comm_num;
  radial_decomp->ghost_sendrecvbuf = ghost_sendrecvbuf;
  radial_decomp->ghost_bufsize = ghost_bufsize;
  radial_decomp->nghost_comm_list = nghost_comm_list;
  radial_decomp->nghost_start = nghost_start;
  radial_decomp->nghost_end = nghost_end;
  radial_decomp->nghost_comm_num = nghost_comm_num;
  radial_decomp->nghost_sendrecvbuf = nghost_sendrecvbuf;
  radial_decomp->nghost_bufsize = nghost_bufsize;

  RngStream_DeleteStream(rs);
  free(rng_seed);
  _mm_free(adum);
  _mm_free(tempi);
  _mm_free(arr);
  _mm_free(arr2);
  _mm_free(ghost_start_local);
  _mm_free(ghost_end_local);
  _mm_free(ghost_start_global);
  _mm_free(ghost_end_global);
  return 0;
}

