#include "bench_gtc.h"
#include "err_check.h"

int diagnosis(gtc_bench_data_t *gtc_input) {
  
  gtc_global_params_t *params;
  gtc_particle_decomp_t *parallel_decomp;
  gtc_diagnosis_data_t *diagnosis_data;
 
  params = &(gtc_input->global_params);
  parallel_decomp = &(gtc_input->parallel_decomp);
  diagnosis_data = &(gtc_input->diagnosis_data);

  int numberpe = parallel_decomp->numberpe;
  int mype = parallel_decomp->mype;

  int mflux = params->mflux;
  real a0 = params->a0;
  real a1 = params->a1;
  real rw = params->rw;
  real rc = params->rc;
  int nbound = params->nbound;
  real kappati = params->kappati;
  real gyroradius = params->gyroradius;
  real qion = params->qion;
  real aion = params->aion;
  real aelectron = params->aelectron;
  real tite = params->tite;

  real eradial = (diagnosis_data->scalar_data)[11]; 

  int mstepall = params->mstepall;
  int ndiag = params->ndiag;
  int mstep = params->mstep;


  real *scalar_data, *eflux, *rmarker;
  real r, rfac, kappa, vthi, vthe, tem_inv;
  real *xnormal = (real *) _mm_malloc(mflux*sizeof(real), IDEAL_ALIGNMENT);
  real *fdum = (real *) _mm_malloc((11+2*mflux)*sizeof(real), IDEAL_ALIGNMENT);
  real *adum = (real *) _mm_malloc((11+2*mflux)*sizeof(real), IDEAL_ALIGNMENT);
  assert(xnormal!=NULL);
  assert(fdum!=NULL);
  assert(adum!=NULL);
  for (int i=0; i<11+2*mflux; i++){
    fdum[i] = 0.0;
    adum[i] = 0.0;
  }
  
  scalar_data = diagnosis_data->scalar_data;
  eflux = diagnosis_data->eflux;
  rmarker = diagnosis_data->rmarker;

  if (mype==0){
      for (int i=0; i<mflux; i++){
        r = a0 + (a1-a0)*( (real)(i+1)-0.5 )/((real)mflux ); 
        rfac = rw*(r - rc);
        rfac = rfac*rfac;
        rfac = rfac*rfac*rfac;
        rfac = max(0.1, exp(-rfac));

        kappa = 1.0;
        if (nbound==0) kappa = 0.0;
        kappa = 1.0 - kappa + kappa*rfac;
       
        xnormal[i] = 1.0/(kappa*kappati*gyroradius);
	/* BMA: reduce verbosity for NERSC runs
        if (istep==ndiag)
           printf("kappa_T at radial_bin %d is %e gyroradius=%e\n", 
		  i, kappa*kappati, gyroradius);
	*/
      }
  }
   
#if USE_MPI 
     MPI_Bcast(&mstepall, 1, MPI_INT, 0, MPI_COMM_WORLD);
     params->mstepall = mstepall;
#endif
   
  // global sum of fluxes
  /* all these quantities come from summing up contributions from the 
     particles, sothe MPI_Reduce involves all the MPI processes */
  vthi = gyroradius*abs(qion)/aion;
  vthe = vthi*sqrt(aion/(aelectron*tite));
  tem_inv = 1.0/(aion*vthi*vthi);
  /*
  fdum[0] = diagnosis_data->efield;
  fdum[1] = diagnosis_data->entropyi;
  fdum[2] = diagnosis_data->entropye;
  fdum[3] = diagnosis_data->dflowi/vthi;
  fdum[4] = diagnosis_data->dflowe/vthe;
  fdum[5] = diagnosis_data->pfluxi/vthi; 
  fdum[6] = diagnosis_data->pfluxe/vthi;
  fdum[7] = diagnosis_data->efluxi*tem_inv/vthi;
  fdum[8] = diagnosis_data->efluxe*tem_inv/vthi;
  fdum[9] = diagnosis_data->particles_energy[0];     
  fdum[10] = diagnosis_data->particles_energy[1];
  for (int i=0; i<mflux; i++){
    fdum[11+i] = (diagnosis_data->eflux[i])*tem_inv/vthi;
    fdum[11+mflux+i] = diagnosis_data->rmarker[i];
    }
  */
  fdum[0] = scalar_data[10];
  fdum[1] = scalar_data[8];
  fdum[2] = scalar_data[9];
  fdum[3] = scalar_data[6]/vthi;
  fdum[4] = scalar_data[7]/vthe;
  fdum[5] = scalar_data[2]/vthi;
  fdum[6] = scalar_data[3]/vthi;
  fdum[7] = scalar_data[0]*tem_inv/vthi;
  fdum[8] = scalar_data[1]*tem_inv/vthi;

  fdum[9] = scalar_data[12];
  fdum[10] = scalar_data[13];

  for (int i=0; i<mflux; i++){
    fdum[11+i] = eflux[i]*tem_inv/vthi;
    fdum[11+mflux+i] = rmarker[i];
  }

#if USE_MPI
  MPI_Reduce(fdum, adum, 11+2*mflux, MPI_MYREAL, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

  if (mype==0){
#ifdef USE_MPI
    for (int i=0; i<11+2*mflux; i++){
      fdum[i] = adum[i];
    }
#endif

    // normalization
    real tmarker= 0.0;
    for (int i=0; i<mflux; i++){
      fdum[11+mflux+i] = max(1.0, fdum[11+mflux+i]);
      fdum[11+i] = fdum[11+i]*xnormal[i]/fdum[11+mflux+i];
      tmarker = tmarker + fdum[11+mflux+i];
    }

    fdum[0] = sqrt(fdum[0]/((real)numberpe));
    for (int i=1; i<9; i++){
        fdum[i] = fdum[i]/tmarker;
    }

    /* BMA: reduce verbosity and I/O of diagnostics for NERSC runs
    printf("istep+mstepall=%d efield=%e eradial=%e entropyi=%e, dflowi=%e, pfluxi=%e, efluxi=%e eflux[2]=%e rmarker[2]=%e particles_energy[0]=%e particles_energy[1]=%e\n",
           istep, fdum[0], eradial, fdum[1], fdum[3], fdum[5], fdum[7], fdum[11+mflux/2], fdum[11+mflux+mflux/2], fdum[9], fdum[10]);

    FILE *pFile;
    if (istep==1)
      pFile = fopen("diag_c.txt", "w");
    else
      pFile = fopen("diag_c.txt", "a");
    fprintf(pFile, "%d %e %e %e %e %e %e %e %e %e %e\n",  istep, fdum[0], eradial, fdum[1], fdum[3], fdum[5], fdum[7], fdum[11+mflux/2], fdum[11+mflux+mflux/2], fdum[9], fdum[10]);
    fclose(pFile);
    */

    if ( istep>=mstep )
      err_check_setFOM( fdum[0], fdum[20] );
  }  

  _mm_free(xnormal);
  _mm_free(fdum);
  _mm_free(adum);
  
    return 0;
}


