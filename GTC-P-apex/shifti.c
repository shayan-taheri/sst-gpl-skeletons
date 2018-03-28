#include "bench_gtc.h"

int radial_bin_particles(gtc_bench_data_t *gtc_input) {
 
    gtc_global_params_t     *params;
    gtc_field_data_t        *field_data;
    gtc_particle_data_t     *particle_data;
    gtc_aux_particle_data_t *aux_particle_data;
    gtc_particle_decomp_t   *parallel_decomp;
    gtc_radial_decomp_t     *radial_decomp;

    real *z0, *z1, *z2, *z3, *z4, *z5;
    real *z00, *z01, *z02, *z03, *z04, *z05;
   
    params            = &(gtc_input->global_params);
    field_data        = &(gtc_input->field_data);
    particle_data     = &(gtc_input->particle_data);
    aux_particle_data = &(gtc_input->aux_particle_data);
    parallel_decomp   = &(gtc_input->parallel_decomp);
    radial_decomp     = &(gtc_input->radial_decomp);

    if (((istep % params->radial_bin_freq) != 0) || (irk == 1)) {
        return 0;
    }
 
    z0 = particle_data->z0; z1 = particle_data->z1; z2 = particle_data->z2;
    z3 = particle_data->z3; z4 = particle_data->z4; z5 = particle_data->z5;
    z00 = particle_data->z00; z01 = particle_data->z01; 
    z02 = particle_data->z02; z03 = particle_data->z03;
    z04 = particle_data->z04; z05 = particle_data->z05;
    int *psi_count = particle_data->psi_count;
    int *psi_offsets = particle_data->psi_offsets;

    int *ipval = aux_particle_data->kzi;
    real *ztmp  = particle_data->ztmp;
    real *ztmp2 = particle_data->ztmp2;
    int mpsi = params->mpsi;
    real a0 = params->a0;
    real delr = params->delr;

    int mstep = params->mstep;
    int mi = params->mi;
    int mi_new = mi;
    int ipsi_valid_in = radial_decomp->ipsi_valid_in;
    int ipsi_valid_out = radial_decomp->ipsi_valid_out;
    int mpsi_loc = params->mpsi_loc;

#pragma omp parallel
{
    int tid, nthreads;
    int *psi_count_l, *psi_offsets_l;
    int i, j, m;

#ifdef _OPENMP
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
#else
    tid = 0;
    nthreads = 1;
#endif

    psi_count_l   = psi_count + mpsi_loc*tid;
    psi_offsets_l = psi_offsets + mpsi_loc*tid;

    for (i=0; i<mpsi_loc; i++){
        psi_count_l[i] = 0;
        psi_offsets_l[i] = 0;
    }

#pragma omp barrier

#pragma omp for
    for (m = 0; m < mi; m++) {

        real zetatmp = z2[m];
        if (zetatmp != HOLEVAL) {
            real psitmp  = z0[m];
#if SQRT_PRECOMPUTED
	    real r = psitmp;
#else
            real r  = sqrt(2.0*psitmp);
#endif
            int iptmp   = (int) ((r-a0)*delr+0.5);
            int ip      = abs_min_int(mpsi, iptmp);

	    ip = ip - ipsi_valid_in;
            psi_count_l[ip]++;
	    ipval[m] = ip;
            
        } else {
            ipval[m] = -1;
        }
    }

    if (tid == 0) {
      for (i=0; i<mpsi_loc-1; i++){
            int full_bin_count = 0; 
            for (j=0; j<nthreads; j++) {
	      full_bin_count += psi_count[j*mpsi_loc+i];
            }
            psi_offsets[i+1] = psi_offsets[i] + full_bin_count;
        }

      for (i=0; i<mpsi_loc; i++){
            for (j=1; j<nthreads; j++) {
	      psi_offsets[j*mpsi_loc+i] = psi_offsets[(j-1)*mpsi_loc+i] + psi_count[(j-1)*mpsi_loc+i];
            }
        }

        if (psi_offsets[(nthreads-1)*mpsi_loc+mpsi_loc-1]+
            psi_count[(nthreads-1)*mpsi_loc+mpsi_loc-1] !=
            (mi-(params->holecount))) {
	  fprintf(stderr, "psi offset sum %d, mi %d\n",
		  psi_offsets[(nthreads-1)*mpsi_loc+mpsi_loc-1]+psi_count[(nthreads-1)*mpsi_loc+mpsi_loc-1], mi);
	  exit(1);
        }
        mi_new =  mi - params->holecount;
        params->mi = mi_new;
        params->holecount = 0;
    }

#pragma omp barrier

    for (i=0; i<mpsi_loc; i++) {
        psi_count_l[i] = 0;
    }

#pragma omp barrier

#pragma omp for
    for (m = 0; m < mi; m++) {

        int ip = ipval[m];
	
        if (ip != -1) {
	  int new_pos   = psi_offsets_l[ip] + psi_count_l[ip]++;
	  z00[new_pos]  = z0[m];
	  z01[new_pos]  = z1[m];
	  z02[new_pos]  = z2[m];
	  z03[new_pos]  = z3[m];
	  z04[new_pos]  = z4[m];
	  ztmp[new_pos] = z5[m];
	  ztmp2[new_pos]= z05[m];
        }
    }

#pragma omp for
    for (m = 0; m < mi_new; m++) {
        z0[m]  = z00[m];
        z1[m]  = z01[m];
        z2[m]  = z02[m];
        z3[m]  = z03[m];
        z4[m]  = z04[m];
    }

}      /* end of parallel region */

    real *tmp_ptr, *tmp_ptr2;
    /* swap ztmp and z5 */
    tmp_ptr = particle_data->z5;
    particle_data->z5 = particle_data->ztmp;
    particle_data->ztmp = tmp_ptr;

    /* swap ztmp2 and z05 */
    tmp_ptr2 = particle_data->z05;
    particle_data->z05 = particle_data->ztmp2;
    particle_data->ztmp2 = tmp_ptr2;

    return 0;
}


#if USE_MPI
int shifti_toroidal(gtc_bench_data_t *gtc_input) {
    gtc_global_params_t     *params;
    gtc_particle_data_t     *particle_data;
    gtc_aux_particle_data_t *aux_particle_data;
    gtc_particle_decomp_t   *parallel_decomp;
    
    int mi, mimax;
    int ntoroidal;
    real zetamax, zetamin;

    int mi_end;

    real *z0, *z1, *z2, *z3, *z4, *z5;
    real *z00, *z01, *z02, *z03, *z04, *z05, *ztmp;
    int *kzi;

    int *mshift;

    int nparam;
    real pi, pi2, pi_inv;
    
    int m0, iteration;
    int max_threads;

    int mstep;

    /*******/

    params            = &(gtc_input->global_params);
    particle_data     = &(gtc_input->particle_data);
    aux_particle_data = &(gtc_input->aux_particle_data);
    parallel_decomp   = &(gtc_input->parallel_decomp);
    
    mi    = params->mi; mimax = params->mimax;
    pi = params->pi;
    zetamax = params->zetamax; zetamin = params->zetamin;
    mstep = params->mstep;

    z0 = particle_data->z0; z1 = particle_data->z1; 
    z2 = particle_data->z2; z3 = particle_data->z3;
    z4 = particle_data->z4; z5 = particle_data->z5;

    z00 = particle_data->z00; z01 = particle_data->z01; 
    z02 = particle_data->z02; z03 = particle_data->z03;
    z04 = particle_data->z04; z05 = particle_data->z05;
    ztmp = particle_data->ztmp;

    kzi = aux_particle_data->kzi;

    ntoroidal = parallel_decomp->ntoroidal;

    nparam = 12;
    pi_inv = 1.0/pi;
    pi2 = 2.0*pi;

    /********/

    if (ntoroidal == 1)
        return 0;

    m0 = 0;
    mi_end = mi;
    iteration = 0;
    
    //max_threads = 64;
#ifdef _OPENMP
    max_threads = omp_get_max_threads();
#else
    max_threads = 1;
#endif
    mshift = (int *) _mm_malloc(max_threads * 8 * sizeof(int),
            IDEAL_ALIGNMENT);
    assert(mshift != NULL);

    while (iteration <= ntoroidal) {
{

        real *sendright, *sendleft;
        real *recvleft, *recvright;
        int msend, mrecv, msendright, msendleft, pos;
        int mrecvleft, mrecvright;
        int isendtag, irecvtag;
        MPI_Status istatus1, istatus2;

        msend = msendright = msendleft = 0;

#if FINE_TIMER
        double start_t = timer();
#endif

#pragma omp parallel 
{    
        int mlstack_max, mrstack_max; 
        int *lstack, *rstack, *stack_tmp;
        int lpos, rpos;

        int i;        
        int tid, nthreads;
        real zetaright, zetaleft;
        real z2val;

#ifdef _OPENMP
        tid = omp_get_thread_num();
        nthreads = omp_get_num_threads();
#else
        tid = 0;
        nthreads = 1;
#endif
        
        assert(nthreads == parallel_decomp->nthreads);
        assert(nthreads <= max_threads);

        mshift[8*tid] = mshift[8*tid+1] 
            = mshift[8*tid+2] = mshift[8*tid+3] 
            = mshift[8*tid+4] = 0;

        mlstack_max = mrstack_max = 32768;
        lstack = (int *) malloc(mlstack_max * sizeof(int));
        rstack = (int *) malloc(mrstack_max * sizeof(int));
        assert(lstack != NULL); assert(rstack != NULL);

#pragma omp barrier

#pragma omp for 
        for (i=m0; i<mi_end; i++) {

            z2val = z2[i];
            if (z2val == HOLEVAL) {
                continue;
            }

            zetaright = min(pi2, z2val) - zetamax;
            zetaleft  = z2val - zetamin;

            if (zetaright * zetaleft > 0.0) {

                zetaright = zetaright*0.5*pi_inv;
                zetaright = zetaright-floor(zetaright);
                
                mshift[8*tid]++;
                
                /* particle to move right */
                if (zetaright < 0.5) {

                    rpos = mshift[8*tid+1]++;
                    if (rpos == mrstack_max) {
		      //		      printf("mype=%d mrstack_max=%d\n",parallel_decomp->mype,mrstack_max);
                        /* double the size */
                        stack_tmp = (int *)
                            malloc(2*mrstack_max*sizeof(int));
			assert(stack_tmp != NULL);
                        memcpy(stack_tmp, rstack, mrstack_max*sizeof(int));
                        free(rstack);
                        rstack = stack_tmp;
                        mrstack_max = 2 * mrstack_max;
                    }
                    rstack[rpos] = i;

                /* particle to move left */        
                } else {
               
                    lpos = mshift[8*tid+2]++;
                    if (lpos == mlstack_max) {
                      //printf("mype=%d mlstack_max=%d\n",parallel_decomp->mype,mlstack_max);
                        /* double the size */
                        stack_tmp = (int *)
                            malloc(2*mlstack_max*sizeof(int));
			assert(stack_tmp != NULL);
                        memcpy(stack_tmp, lstack, mlstack_max*sizeof(int));
                        free(lstack);
                        lstack = stack_tmp;
                        mlstack_max = 2 * mlstack_max;
                    }
                    lstack[lpos] = i;
                }

            }
        }

        /* Merge partial arrays */
        if (tid == 0) {
            mshift[8*0+3] = 0;
            mshift[8*0+4] = 0;
            for (i=1; i<nthreads; i++) {
                mshift[8*i+3] = mshift[8*(i-1)+3] + mshift[8*(i-1)+1];
                mshift[8*i+4] = mshift[8*(i-1)+4] + mshift[8*(i-1)+2];
            }
            msendright  = mshift[8*(nthreads-1)+3] + mshift[8*(nthreads-1)+1];
            msendleft   = mshift[8*(nthreads-1)+4] + mshift[8*(nthreads-1)+2];
            msend = msendright + msendleft;
	    
            if (msend >= mi) {
                fprintf(stderr, "Error! mype %d, msend %d, left %d, right"
                        " %d, mi %d\n", 
                        parallel_decomp->mype, msend, msendleft, msendright, mi);
                exit(1);
            }
        } 

#pragma omp barrier
        
        memcpy(kzi+mshift[8*tid+3], rstack, mshift[8*tid+1]*sizeof(int));
        memcpy(kzi+msendright+mshift[8*tid+4], lstack,
                mshift[8*tid+2]*sizeof(int));
        free(lstack);
        free(rstack);
}

#if FINE_TIMER
 double end_t = timer();
 shift_t_comp += (end_t - start_t);
#endif
        iteration++;
        if (iteration > 1) {
            mrecv = 0;
            MPI_Allreduce(&msend, &mrecv, 1, MPI_INT, MPI_SUM,
                    parallel_decomp->toroidal_comm );
            if (mrecv == 0) {
                break;
            }
        }

        /* allocate mem for incoming particle data */
        if (msendleft + msendright >= parallel_decomp->sendbuf_size) {
            fprintf(stderr, "Toroidal Error! PE %d, msendleft %d, msendright %d, "
                    "sendbuf_size %d\n", parallel_decomp->mype, 
                    msendleft, msendright,
                    parallel_decomp->sendbuf_size);
	    MPI_Abort(MPI_COMM_WORLD, 1);
        } else {
            sendleft = parallel_decomp->sendbuf;
            sendright = parallel_decomp->sendbuf + nparam*msendleft;
        }

#if 0
        /* Allocate space for sendright and sendleft */
        sendright = (real *) malloc(nparam*msendright*sizeof(real));
        assert(sendright != NULL);

        sendleft = (real *) malloc(nparam*msendleft*sizeof(real));
        assert(sendleft != NULL);
#endif

#if FINE_TIMER
        start_t = timer();
#endif

        /* pack particle data */
#pragma omp parallel for private(pos)
        for (int i=0; i<msendright; i++) {
            pos = kzi[i];
            sendright[nparam*i+0] = z0[pos];
            sendright[nparam*i+1] = z1[pos];
            sendright[nparam*i+2] = z2[pos];
            z2[pos] = HOLEVAL;
            sendright[nparam*i+3] = z3[pos];
            sendright[nparam*i+4] = z4[pos];
            sendright[nparam*i+5] = z5[pos];
            sendright[nparam*i+6] = z00[pos];
            sendright[nparam*i+7] = z01[pos];
            sendright[nparam*i+8] = z02[pos];
            z02[pos] = HOLEVAL;
            sendright[nparam*i+9] = z03[pos];
            sendright[nparam*i+10] = z04[pos];
            sendright[nparam*i+11] = z05[pos];
        }

#pragma omp parallel for private(pos)
        for (int i=0; i<msendleft; i++) {
            pos = kzi[i+msendright];
            sendleft[nparam*i+0] = z0[pos];
            sendleft[nparam*i+1] = z1[pos];
            sendleft[nparam*i+2] = z2[pos];
            z2[pos] = HOLEVAL;
            sendleft[nparam*i+3] = z3[pos];
            sendleft[nparam*i+4] = z4[pos];
            sendleft[nparam*i+5] = z5[pos];
            sendleft[nparam*i+6] = z00[pos];
            sendleft[nparam*i+7] = z01[pos];
            sendleft[nparam*i+8] = z02[pos];
            z02[pos] = HOLEVAL;
            sendleft[nparam*i+9] = z03[pos];
            sendleft[nparam*i+10] = z04[pos];
            sendleft[nparam*i+11] = z05[pos];
        }

#if FINE_TIMER
        end_t = timer();
        shift_t_comm1 += (end_t - start_t);
        start_t = timer();
#endif

        /* send # of particles to move right, get # to recv from left */
        isendtag = parallel_decomp->myrank_toroidal;
        irecvtag = parallel_decomp->left_pe;
        mrecvleft = 0;
        MPI_Sendrecv(&msendright, 1, MPI_INT, parallel_decomp->right_pe, 
                     isendtag, &mrecvleft, 1, MPI_INT, parallel_decomp->left_pe,
                     irecvtag, parallel_decomp->toroidal_comm, &istatus1);
       
        /* send # of particles to move left, get # to recv from right */
        irecvtag = parallel_decomp->right_pe;
        mrecvright = 0;
        MPI_Sendrecv(&msendleft, 1, MPI_INT, parallel_decomp->left_pe, 
               isendtag, &mrecvright, 1, MPI_INT, parallel_decomp->right_pe, 
               irecvtag, parallel_decomp->toroidal_comm, &istatus2);
#if FINE_TIMER
        end_t = timer();
        shift_t_comm2 += (end_t - start_t);
#endif

        /* allocate mem for incoming particle data */
        if (mrecvleft + mrecvright >= parallel_decomp->recvbuf_size) {
            fprintf(stderr, "Toroidal Error! PE %d, mrecvleft %d, mrecvright %d, "
                    "recvbuf_size %d\n", parallel_decomp->mype, 
                    mrecvleft, mrecvright,
                    parallel_decomp->recvbuf_size);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        recvleft = parallel_decomp->recvbuf;
        recvright = parallel_decomp->recvbuf + nparam*mrecvleft;

#if 0
        /* Allocate space for sendright and sendleft */
        recvright = (real *) malloc(nparam*mrecvright*sizeof(real));
        assert(recvright != NULL);

        recvleft = (real *) malloc(nparam*mrecvleft*sizeof(real));
        assert(recvleft != NULL);
#endif
   
#if FINE_TIMER
        start_t = timer();
#endif
        /* send particles to right neighbor, recv from left */
        irecvtag = parallel_decomp->left_pe;
        MPI_Sendrecv(sendright, msendright*nparam, MPI_MYREAL, 
                parallel_decomp->right_pe,
                isendtag, recvleft, mrecvleft*nparam, MPI_MYREAL, 
                parallel_decomp->left_pe,
                irecvtag, parallel_decomp->toroidal_comm, &istatus1);

        /* send particles to left neighbor, recv from right */
        irecvtag = parallel_decomp->right_pe;
        MPI_Sendrecv(sendleft, msendleft*nparam, MPI_MYREAL, 
                parallel_decomp->left_pe, 
                isendtag, recvright, mrecvright*nparam, MPI_MYREAL, 
                parallel_decomp->right_pe,
                irecvtag, parallel_decomp->toroidal_comm, &istatus2);
        /* copy received data to particle arrays */
	assert(mi_end + mrecvleft + mrecvright < mimax);
#if FINE_TIMER
        end_t = timer();
        shift_t_comm2 += (end_t - start_t);
        start_t = timer();
#endif

#pragma omp parallel for private(pos)
        for (int i=0; i<mrecvleft; i++) {
            pos = mi_end + i;
            z0[pos]  = recvleft[nparam*i+0];
            z1[pos]  = recvleft[nparam*i+1];
            z2[pos]  = recvleft[nparam*i+2];
            z3[pos]  = recvleft[nparam*i+3];
            z4[pos]  = recvleft[nparam*i+4];
            z5[pos]  = recvleft[nparam*i+5];
            z00[pos] = recvleft[nparam*i+6];
            z01[pos] = recvleft[nparam*i+7];
            z02[pos] = recvleft[nparam*i+8];
            z03[pos] = recvleft[nparam*i+9];
            z04[pos] = recvleft[nparam*i+10];
            z05[pos] = recvleft[nparam*i+11];
        }

#pragma omp parallel for private(pos)
        for (int i=0; i<mrecvright; i++) {
            pos = mi_end + mrecvleft + i;
            z0[pos]  = recvright[nparam*i+0];
            z1[pos]  = recvright[nparam*i+1];
            z2[pos]  = recvright[nparam*i+2];
            z3[pos]  = recvright[nparam*i+3];
            z4[pos]  = recvright[nparam*i+4];
            z5[pos]  = recvright[nparam*i+5];
            z00[pos] = recvright[nparam*i+6];
            z01[pos] = recvright[nparam*i+7];
            z02[pos] = recvright[nparam*i+8];
            z03[pos] = recvright[nparam*i+9];
            z04[pos] = recvright[nparam*i+10];
            z05[pos] = recvright[nparam*i+11];
        }
#if FINE_TIMER
        end_t = timer()
        shift_t_comm1 += (end_t - start_t);
#endif

        /* update m0, m_end, and other global counts */
        m0 = mi_end;
        params->mi += mrecvleft + mrecvright;
        mi_end = params->mi;
        params->holecount += msend;

#if 0
        free(sendleft);
	free(sendright);
	free(recvleft);
	free(recvright);
#endif


} /* while section */

   }


   if (iteration > (ntoroidal+1)) {
      fprintf(stderr, "Error! Endless particle sorting loop at PE %d\n",
	  parallel_decomp->mype);
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   params->mi = mi_end;

   _mm_free(mshift);

    return 0;
}

int shifti_radial(gtc_bench_data_t *gtc_input){

  gtc_global_params_t     *params;
  gtc_particle_data_t     *particle_data;
  gtc_aux_particle_data_t *aux_particle_data;
  gtc_particle_decomp_t   *parallel_decomp;
  gtc_radial_decomp_t     *radial_decomp;

  int mi, mimax;
  int nradial_dom, myrank_radiald;
  real a_nover_in, a_nover_out;  

  int mi_end;
  real *z0, *z1, *z2, *z3, *z4, *z5;
  real *z00, *z01, *z02, *z03, *z04, *z05, *ztmp;
  int *kzi;

  int *mshift;
  
  int nparam;
  real pi, pi2, pi_inv;

  int m0, iteration;
  int max_threads;

  /************/
  
  params            = &(gtc_input->global_params);
  particle_data     = &(gtc_input->particle_data);
  aux_particle_data = &(gtc_input->aux_particle_data);
  parallel_decomp   = &(gtc_input->parallel_decomp);
  radial_decomp     = &(gtc_input->radial_decomp);

  mi = params->mi; mimax = params->mimax;
  pi = params->pi;
  
  z0 = particle_data->z0; z1 = particle_data->z1;
  z2 = particle_data->z2; z3 = particle_data->z3;
  z4 = particle_data->z4; z5 = particle_data->z5;

  z00 = particle_data->z00; z01 = particle_data->z01;  
  z02 = particle_data->z02; z03 = particle_data->z03;
  z04 = particle_data->z04; z05 = particle_data->z05;
  ztmp = particle_data->ztmp;

  kzi = aux_particle_data->kzi;

  nradial_dom = radial_decomp->nradial_dom;  
  myrank_radiald = radial_decomp->myrank_radiald;

  a_nover_in = radial_decomp->a_nover_in;
  a_nover_out = radial_decomp->a_nover_out;

  nparam = 12;
  pi_inv = 1.0/pi;
  pi2    = 2.0*pi;

  if (nradial_dom == 1) return 0;

  m0 = 0;
  mi_end = mi;
  iteration = 0;

  //max_threads = 64;
#ifdef _OPENMP
  max_threads = omp_get_max_threads();
#else
  max_threads = 1;
#endif

  mshift = (int *) _mm_malloc(max_threads * 8 * sizeof(int), IDEAL_ALIGNMENT);
  assert(mshift != NULL);

  while (iteration <= nradial_dom) {
    {
      real *sendright, *sendleft;
      real *recvleft, *recvright;
      int msend, mrecv, msendright, msendleft, pos;
      int mrecvleft, mrecvright;
      int isendtag, irecvtag;
      MPI_Status istatus1, istatus2;

      msend = msendright = msendleft = 0;

#pragma omp parallel
      {
        int mlstack_max, mrstack_max;
        int *lstack, *rstack,*stack_tmp;
        int lpos, rpos;

        int i;
        int tid, nthreads;
        real radialright, radialleft;
        real z2val, psitmp, r;

#ifdef _OPENMP
        tid = omp_get_thread_num();
        nthreads = omp_get_num_threads();
#else
        tid = 0;
        nthreads = 1;
#endif

        assert(nthreads == parallel_decomp->nthreads);
        assert(nthreads <= max_threads);

        mshift[8*tid] = mshift[8*tid+1]
          = mshift[8*tid+2] = mshift[8*tid+3]
          = mshift[8*tid+4] = 0;

        mlstack_max = mrstack_max = 8192;
        lstack = (int *) malloc(mlstack_max * sizeof(int));
        rstack = (int *) malloc(mrstack_max * sizeof(int));
        assert(lstack != NULL); assert(rstack != NULL);

#pragma omp barrier

#pragma omp for 
        for (i=m0; i<mi_end; i++){
          
          z2val = z2[i];
          if (z2val == HOLEVAL) {
            continue;
          }

          psitmp = z0[i];
          
#if SQRT_PRECOMPUTED
          r = psitmp;
#else
          r = sqrt(2.0*psitmp);
#endif          
          if ((r<a_nover_in && myrank_radiald > 0) ||
              (r>a_nover_out && myrank_radiald < (nradial_dom-1))) {
           
            mshift[8*tid]++;

            /* particle to move right */
            if (r > a_nover_out ) {
              rpos = mshift[8*tid+1]++;
              if (rpos == mrstack_max){
                /* double the size */
		//printf("radial mype=%d mrstack_max=%d\n",parallel_decomp->mype,mrstack_max);
                stack_tmp = (int *) 
                  malloc(2*mrstack_max*sizeof(int));
		assert(stack_tmp != NULL);
                memcpy(stack_tmp, rstack, mrstack_max*sizeof(int));
                free(rstack);
                rstack = stack_tmp;
                mrstack_max = 2 * mrstack_max;
              }
              rstack[rpos] = i;
            }
            else { /* particle to move left */
              lpos = mshift[8*tid+2]++;
              if (lpos == mlstack_max){
		//printf("radial mype=%d mlstack_max=%d\n",parallel_decomp->mype,mlstack_max);
                stack_tmp = (int *) 
                  malloc(2*mlstack_max*sizeof(int));
		assert(stack_tmp != NULL);
                memcpy(stack_tmp, lstack, mlstack_max*sizeof(int));
                free(lstack);
                lstack = stack_tmp;
                mlstack_max = 2 * mlstack_max;
              }
              lstack[lpos] = i;
            }

          }

        } // end i=m0, mi_end loop

        /* Merge partial arrays */
        if (tid == 0) {
          mshift[8*0+3] = 0;
          mshift[8*0+4] = 0;
          for (i=1; i<nthreads; i++) {
            mshift[8*i+3] = mshift[8*(i-1)+3] + mshift[8*(i-1)+1];
            mshift[8*i+4] = mshift[8*(i-1)+4] + mshift[8*(i-1)+2];
          }
          msendright  = mshift[8*(nthreads-1)+3] + mshift[8*(nthreads-1)+1];
          msendleft   = mshift[8*(nthreads-1)+4] + mshift[8*(nthreads-1)+2];
          msend = msendright + msendleft;
          if (msend >= mi) {                                             
            fprintf(stderr, "Error! mype %d, msend %d, left %d, right"
                    " %d, mi %d\n",                                    
                    parallel_decomp->mype, msend, msendleft, msendright, mi);
            MPI_Abort(MPI_COMM_WORLD, 1);
          }
        }

#pragma omp barrier

        memcpy(kzi+mshift[8*tid+3], rstack, mshift[8*tid+1]*sizeof(int));
        memcpy(kzi+msendright+mshift[8*tid+4], lstack,
               mshift[8*tid+2]*sizeof(int));
        free(lstack);
        free(rstack);
      } // end parallel 
      
      iteration++;
      if (iteration > 1) {
        mrecv = 0;
        MPI_Allreduce(&msend, &mrecv, 1, MPI_INT, MPI_SUM,
                      parallel_decomp->partd_comm );                                                          
        if (mrecv == 0) {                                                                        
          break;                                                                                   
        }
      }
      
      /* allocate mem for incoming particle data */     
      if (msendleft + msendright >= parallel_decomp->sendbuf_size) {
        fprintf(stderr, "Radial Error! PE %d, msendleft %d, msendright %d, "
                "sendbuf_size %d\n", parallel_decomp->mype,
		msendleft, msendright,
		parallel_decomp->sendbuf_size);
        MPI_Abort(MPI_COMM_WORLD, 1);
      } else {
        sendleft = parallel_decomp->sendbuf;
        sendright = parallel_decomp->sendbuf + nparam*msendleft;
      } 

      /* pack particle data */                                                                
#pragma omp parallel for private(pos)                                                           
      for (int i=0; i<msendright; i++) {                                                      
        pos = kzi[i];                                                                       
        sendright[nparam*i+0] = z0[pos];                                                    
        sendright[nparam*i+1] = z1[pos];                                                    
        sendright[nparam*i+2] = z2[pos];                                                    
        z2[pos] = HOLEVAL;                                                                  
        sendright[nparam*i+3] = z3[pos];                                                    
        sendright[nparam*i+4] = z4[pos];                                                    
        sendright[nparam*i+5] = z5[pos];                                                    
        sendright[nparam*i+6] = z00[pos];                                                   
        sendright[nparam*i+7] = z01[pos];                                                   
        sendright[nparam*i+8] = z02[pos];                                                   
        z02[pos] = HOLEVAL;                                                                 
        sendright[nparam*i+9] = z03[pos];                                                   
        sendright[nparam*i+10] = z04[pos];                                                  
        sendright[nparam*i+11] = z05[pos];                                                  
      }  

#pragma omp parallel for private(pos)                                                           
      for (int i=0; i<msendleft; i++) {                                                       
        pos = kzi[i+msendright];                                                            
        sendleft[nparam*i+0] = z0[pos];                                                     
        sendleft[nparam*i+1] = z1[pos];                                                     
        sendleft[nparam*i+2] = z2[pos];                                                     
        z2[pos] = HOLEVAL;                                                                  
        sendleft[nparam*i+3] = z3[pos];                                                     
        sendleft[nparam*i+4] = z4[pos];                                                     
        sendleft[nparam*i+5] = z5[pos];                                                     
        sendleft[nparam*i+6] = z00[pos];                                                    
        sendleft[nparam*i+7] = z01[pos];                                                    
        sendleft[nparam*i+8] = z02[pos];                                                    
        z02[pos] = HOLEVAL;                                                                 
        sendleft[nparam*i+9] = z03[pos];                                                    
        sendleft[nparam*i+10] = z04[pos];                                                   
        sendleft[nparam*i+11] = z05[pos];                                                   
      }     

      /* send # of particles to move right, get # to recv from left */                        
      isendtag = parallel_decomp->myrank_partd;                                        
      irecvtag = radial_decomp->left_radial_pe;                                                    
      mrecvleft = 0;                                               
      MPI_Sendrecv(&msendright, 1, MPI_INT, radial_decomp->right_radial_pe,                
                   isendtag, &mrecvleft, 1, MPI_INT, radial_decomp->left_radial_pe,
                   irecvtag, parallel_decomp->partd_comm, &istatus1);

      /* send # of particles to move left, get # to recv from right */                        
      irecvtag = radial_decomp->right_radial_pe;
      mrecvright = 0;                                                                         
      MPI_Sendrecv(&msendleft, 1, MPI_INT, radial_decomp->left_radial_pe,              
                   isendtag, &mrecvright, 1, MPI_INT, radial_decomp->right_radial_pe,      
                   irecvtag, parallel_decomp->partd_comm, &istatus2);

      /* allocate mem for incoming particle data */                                           
      if (mrecvleft + mrecvright >= parallel_decomp->recvbuf_size) {                          
        fprintf(stderr, "Radial Error! PE %d, mrecvleft %d, mrecvright %d, "                       
                "recvbuf_size %d\n", parallel_decomp->mype,                                 
                mrecvleft, mrecvright, 
                parallel_decomp->recvbuf_size);  
        MPI_Abort(MPI_COMM_WORLD, 1);
      }                                                                                       
      recvleft = parallel_decomp->recvbuf;                                                    
      recvright = parallel_decomp->recvbuf + nparam*mrecvleft;  
                                                                             
      /* send particles to right neighbor, recv from left */                                  
      irecvtag = radial_decomp->left_radial_pe;                                                
      MPI_Sendrecv(sendright, msendright*nparam, MPI_MYREAL,
                   radial_decomp->right_radial_pe,
                   isendtag, recvleft, mrecvleft*nparam, MPI_MYREAL,
                   radial_decomp->left_radial_pe,
                   irecvtag, parallel_decomp->partd_comm, &istatus1);
      /* send particles to left neighbor, recv from right */                                  
      irecvtag = radial_decomp->right_radial_pe;
      MPI_Sendrecv(sendleft, msendleft*nparam, MPI_MYREAL,
                   radial_decomp->left_radial_pe,
                   isendtag, recvright, mrecvright*nparam, MPI_MYREAL,
                   radial_decomp->right_radial_pe,
                   irecvtag, parallel_decomp->partd_comm, &istatus2);

      /* copy received data to particle arrays */                                             
      assert(mi_end + mrecvleft + mrecvright < mimax);

#pragma omp parallel for private(pos)
      for (int i=0; i<mrecvleft; i++) {                                                       
        pos = mi_end + i;                                                                   
        z0[pos]  = recvleft[nparam*i+0];                                                    
        z1[pos]  = recvleft[nparam*i+1];                                                    
        z2[pos]  = recvleft[nparam*i+2];                                                    
        z3[pos]  = recvleft[nparam*i+3];                                                    
        z4[pos]  = recvleft[nparam*i+4];                                                    
        z5[pos]  = recvleft[nparam*i+5];                                                    
        z00[pos] = recvleft[nparam*i+6];                                                    
        z01[pos] = recvleft[nparam*i+7];                                                    
        z02[pos] = recvleft[nparam*i+8];                                                    
        z03[pos] = recvleft[nparam*i+9];                                                    
        z04[pos] = recvleft[nparam*i+10];                                                   
        z05[pos] = recvleft[nparam*i+11];                                                   
      }                                                                                       

#pragma omp parallel for private(pos)
      for (int i=0; i<mrecvright; i++) {                                                      
        pos = mi_end + mrecvleft + i;                                                       
        z0[pos]  = recvright[nparam*i+0];                                                   
        z1[pos]  = recvright[nparam*i+1];                                                   
        z2[pos]  = recvright[nparam*i+2];                                                   
        z3[pos]  = recvright[nparam*i+3];                                                   
        z4[pos]  = recvright[nparam*i+4];                                                   
        z5[pos]  = recvright[nparam*i+5];                                                   
        z00[pos] = recvright[nparam*i+6];                                                   
        z01[pos] = recvright[nparam*i+7];                                                   
        z02[pos] = recvright[nparam*i+8];                                                   
        z03[pos] = recvright[nparam*i+9];                                                   
        z04[pos] = recvright[nparam*i+10];                                                  
        z05[pos] = recvright[nparam*i+11];                                                  
      }    

      /* update m0, m_end, and other global counts */                                         
      m0 = mi_end;                                                                            
      params->mi += mrecvleft + mrecvright;                                                   
      mi_end = params->mi;                                                                    
      params->holecount += msend;     
        
    } // end while section

  } // end while
                                                 
  if (iteration > (nradial_dom+1)) {                                                            
    fprintf(stderr, "Error! Endless particle sorting loop at PE %d\n",                      
            parallel_decomp->mype);                                                         
    MPI_Abort(MPI_COMM_WORLD, 1);
  }                                                                                           

  params->mi = mi_end;                                                 
                                                                    
  _mm_free(mshift);
    
  return 0;

}

int shifti_radial_a2a( gtc_bench_data_t *gtc_input ){

  gtc_global_params_t *params = &(gtc_input->global_params);
  gtc_radial_decomp_t *radial_decomp = &(gtc_input->radial_decomp);
  gtc_particle_data_t *particle_data = &(gtc_input->particle_data);
  gtc_particle_decomp_t *parallel_decomp  = &(gtc_input->parallel_decomp );
  gtc_aux_particle_data_t *aux_particle_data = &(gtc_input->aux_particle_data);  

  real *z0, *z1, *z2, *z3, *z4, *z5;
  z0 = particle_data->z0; z1 = particle_data->z1;
  z2 = particle_data->z2; z3 = particle_data->z3;
  z4 = particle_data->z4; z5 = particle_data->z5;

  real *z00, *z01, *z02, *z03, *z04, *z05, *ztmp;
  z00 = particle_data->z00; z01 = particle_data->z01;  
  z02 = particle_data->z02; z03 = particle_data->z03;
  z04 = particle_data->z04; z05 = particle_data->z05;
  ztmp = particle_data->ztmp;

  int m, p, ps, pr, q;
  int mi = params->mi;
  int mimax = params->mimax;
  MPI_Comm radiald_comm = radial_decomp->radiald_comm;
  int nradial_dom = radial_decomp->nradial_dom;
  int myrank_radiald = radial_decomp->myrank_radiald;
  if( nradial_dom==1 ) return 0;

  int *dest = aux_particle_data->kzi;
  int *sendcounts = (int*) malloc(  nradial_dom * sizeof(int) );
  int *recvcounts = (int*) malloc(  nradial_dom * sizeof(int) );
  int *sdispls    = (int*) malloc( (nradial_dom+1) * sizeof(int) );
  int *rdispls    = (int*) malloc( (nradial_dom+1) * sizeof(int) );
  assert( sendcounts != NULL );
  assert( recvcounts != NULL );
  assert( sdispls != NULL );
  assert( rdispls != NULL );
  for( p=0; p<nradial_dom; p++ ) sendcounts[p]=0;

  real psitmp, r, a0, a1;
  a0 = params->a0;
  a1 = params->a1;

  #ifdef _OPENMP
  #pragma omp parallel for private( psitmp, r, p )
  #endif
  for( m=0; m<mi; m++ ){
    
    psitmp = z0[m];
    #if SQRT_PRECOMPUTED
      r = psitmp;
    #else
      r = sqrt( 2.0 * psitmp );
    #endif
    p = (r*r-a0*a0)/(a1*a1-a0*a0)*nradial_dom;
    dest[m]=p;

    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    sendcounts[p]++;
  }
  
  // convert counts to offsets
  sdispls[0] = 0;
  for( p=1; p<=nradial_dom; p++ )
    sdispls[p] = sdispls[p-1] + sendcounts[p-1];

  //convert destination list to reorder list
  for( p=0; p<nradial_dom; p++ ) sendcounts[p]=0;
  #ifdef _OPENMP
  #pragma omp parallel for private(p,q)
  #endif
  for( m=0; m<mi; m++ ){

    p = dest[m];

    #ifdef _OPENMP
    #pragma omp atomic capture
    #endif
    q = sendcounts[p]++;

    dest[m] = sdispls[p] + q;    
  }
          
  MPI_Alltoall( sendcounts, 1, MPI_INT,
		recvcounts, 1, MPI_INT,
		radiald_comm );

  rdispls[0] = 0;
  for( p=1; p<=nradial_dom; p++ ) 
    rdispls[p] = rdispls[p-1] + recvcounts[p-1];
  assert( 0 <= rdispls[nradial_dom] && rdispls[nradial_dom] < mimax );

  int i, nparam=12;
  real *zsrc=NULL;
  for( i=0; i<nparam; i++ ){

    switch(i){
    case( 0): zsrc=z0;  break;
    case( 1): zsrc=z1;  break;
    case( 2): zsrc=z2;  break;
    case( 3): zsrc=z3;  break;
    case( 4): zsrc=z4;  break;
    case( 5): zsrc=z5;  break;
    case( 6): zsrc=z00; break;
    case( 7): zsrc=z01; break;
    case( 8): zsrc=z02; break;
    case( 9): zsrc=z03; break;
    case(10): zsrc=z04; break;
    case(11): zsrc=z05; break;
    default: assert(0);
    }

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for( m=0; m<mi; m++ ){
      ztmp[ dest[m] ] = zsrc[m];
    }

    MPI_Alltoallv( ztmp, sendcounts, sdispls, MPI_MYREAL,
		   zsrc, recvcounts, rdispls, MPI_MYREAL,
		   radiald_comm );
    //if( parallel_decomp->mype==0 ) 
    //printf("shifti_radial_buffered: completed a2a[%d] \n", i );
  }
  params->mi = mi = rdispls[ nradial_dom ];
  
  /*
  for( m=0; m<mi; m++ ){
    
    psitmp = z0[m];
    #if SQRT_PRECOMPUTED
      r = psitmp;
    #else
      r = sqrt( 2.0 * psitmp );
    #endif

   if( r < a_nover_in  && myrank_radiald != 0 )
     printf("particle in wrong radial domain (lower).\n");
   if( r > a_nover_out && myrank_radiald != nradial_dom-1 )
     printf("particle in wrong radial domain (upper) \n");
  }
  */
  
  //free( dest );
  free( sendcounts );
  free( recvcounts );
  free( sdispls );
  free( rdispls );
  
  return 0;
}//shifti_radial_a2a

#endif
