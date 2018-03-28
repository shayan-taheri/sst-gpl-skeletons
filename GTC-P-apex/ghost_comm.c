#include "bench_gtc.h"

#if USE_MPI
int sum_plane(gtc_bench_data_t *gtc_input){
  gtc_global_params_t     *params;
  gtc_field_data_t        *field_data;
  gtc_radial_decomp_t     *radial_decomp;

  params            = &(gtc_input->global_params);
  field_data        = &(gtc_input->field_data);
  radial_decomp     = &(gtc_input->radial_decomp);

  int mpsi = params->mpsi;
  int mzeta = params->mzeta;
  
  const int*  __restrict__ igrid;
  const int*  __restrict__ mtheta;

  igrid = field_data->igrid; 
  mtheta = field_data->mtheta;

  real *densityi = field_data->densityi;

  int igrid_in = radial_decomp->igrid_in;
  int *ghost_comm_list = radial_decomp->ghost_comm_list;
  int *ghost_start = radial_decomp->ghost_start;
  int *ghost_end = radial_decomp->ghost_end;
  int ghost_comm_num = radial_decomp->ghost_comm_num;
  real *ghost_sendrecvbuf = radial_decomp->ghost_sendrecvbuf;
  int ghost_bufsize = radial_decomp->ghost_bufsize;
  int *nghost_comm_list = radial_decomp->nghost_comm_list;
  int *nghost_start = radial_decomp->nghost_start;
  int *nghost_end = radial_decomp->nghost_end;
  real *nghost_sendrecvbuf = radial_decomp->nghost_sendrecvbuf;
  int nghost_comm_num = radial_decomp->nghost_comm_num;
  int nghost_bufsize = radial_decomp->nghost_bufsize;
  int nproc_radial_partd = radial_decomp->nproc_radial_partd;
  int ipsi_nover_in = radial_decomp->ipsi_nover_in;
  int ipsi_nover_out = radial_decomp->ipsi_nover_out;

  int proc_num, count, ipsi_s,ipsi_e, offset, offset2;
  MPI_Request *recv_request = (MPI_Request *) malloc(nghost_comm_num*sizeof(MPI_Request));
  MPI_Status  *recv_status = (MPI_Status *) malloc(nghost_comm_num*sizeof(MPI_Status));

  MPI_Request *send_request = (MPI_Request *) malloc(ghost_comm_num*sizeof(MPI_Request));
  MPI_Status  *send_status = (MPI_Status *) malloc(ghost_comm_num*sizeof(MPI_Status));

  MPI_Request *recv_request1 = (MPI_Request *) malloc(ghost_comm_num*sizeof(MPI_Request));
  MPI_Status  *recv_status1 = (MPI_Status *) malloc(ghost_comm_num*sizeof(MPI_Status));

  MPI_Request *send_request1 = (MPI_Request *) malloc(nghost_comm_num*sizeof(MPI_Request));
  MPI_Status  *send_status1 = (MPI_Status *) malloc(nghost_comm_num*sizeof(MPI_Status));

  //MPI_Barrier(radial_decomp->radiald_comm);
  offset = 0;
  for (int i=0; i<nghost_comm_num; i++){
    proc_num = nghost_comm_list[i];
    ipsi_s = nghost_start[i];
    ipsi_e = nghost_end[i];

    count = (igrid[ipsi_e] + mtheta[ipsi_e] - igrid[ipsi_s] + 1)*(1+mzeta);
    
    MPI_Irecv(nghost_sendrecvbuf+offset, count, MPI_MYREAL, proc_num, proc_num, 
	      radial_decomp->radiald_comm, &recv_request[i]);
    offset += count;
  }
  assert(offset == nghost_bufsize*(mzeta+1));

  for (int i=0; i<ghost_comm_num; i++){
    proc_num = ghost_comm_list[i];
    ipsi_s = ghost_start[i];
    ipsi_e = ghost_end[i];
    offset = (igrid[ipsi_s] - igrid_in)*(1+mzeta);

    count = (igrid[ipsi_e] + mtheta[ipsi_e] - igrid[ipsi_s] + 1)*(1+mzeta);
    MPI_Isend(densityi+offset, count, MPI_MYREAL, proc_num, radial_decomp->myrank_radiald, 
	      radial_decomp->radiald_comm, &send_request[i]);
  }

  MPI_Waitall(nghost_comm_num, recv_request, recv_status);
  MPI_Waitall(ghost_comm_num, send_request, send_status);
  
  offset = 0;
  for (int i=0; i<nghost_comm_num; i++){
    proc_num = nghost_comm_list[i];
    ipsi_s = nghost_start[i];
    ipsi_e = nghost_end[i];
    offset2 = (igrid[ipsi_s] - igrid_in)*(1+mzeta);

    count = (igrid[ipsi_e] + mtheta[ipsi_e] - igrid[ipsi_s] + 1)*(1+mzeta);

    //#pragma omp parallel for num_threads(NUM_THREADS)
#pragma omp parallel for
    for (int j=0; j<count; j++){
      densityi[j+offset2] += nghost_sendrecvbuf[j+offset];
    }
    offset += count;
  }
 
  offset = 0;
  for (int i=0; i<ghost_comm_num; i++){
    proc_num = ghost_comm_list[i];
    ipsi_s = ghost_start[i];
    ipsi_e = ghost_end[i];

    count = (igrid[ipsi_e] + mtheta[ipsi_e] - igrid[ipsi_s] + 1)*(1+mzeta);
    MPI_Irecv(ghost_sendrecvbuf+offset, count, MPI_MYREAL, proc_num, proc_num,
              radial_decomp->radiald_comm, &recv_request1[i]);
    offset += count;
  }
  assert(offset == ghost_bufsize*(mzeta+1));

  for (int i=0; i<nghost_comm_num; i++){
    proc_num = nghost_comm_list[i];
    ipsi_s = nghost_start[i];
    ipsi_e = nghost_end[i];
    offset = (igrid[ipsi_s] - igrid_in)*(1+mzeta);
    count = (igrid[ipsi_e] + mtheta[ipsi_e] - igrid[ipsi_s] + 1)*(1+mzeta);

    MPI_Isend(densityi+offset, count, MPI_MYREAL, proc_num, radial_decomp->myrank_radiald,
              radial_decomp->radiald_comm, &send_request1[i]);
  }
  
  MPI_Waitall(ghost_comm_num, recv_request1, recv_status1);
  MPI_Waitall(nghost_comm_num, send_request1, send_status1);

  offset = 0;
  for (int i=0; i<ghost_comm_num; i++){
    proc_num = ghost_comm_list[i];
    ipsi_s = ghost_start[i];
    ipsi_e = ghost_end[i];
    offset2 = (igrid[ipsi_s] - igrid_in)*(1+mzeta);
    count = (igrid[ipsi_e] + mtheta[ipsi_e] - igrid[ipsi_s] + 1)*(1+mzeta);
    //#pragma omp parallel for num_threads(NUM_THREADS)
#pragma omp parallel for
    for (int j=0; j<count; j++){
      densityi[j+offset2] = ghost_sendrecvbuf[j+offset];
    }
    offset += count;
  }
  assert(offset == ghost_bufsize*(mzeta+1));

  free(recv_request);
  free(recv_status);
  free(send_request);
  free(send_status);

  free(recv_request1);
  free(recv_status1);
  free(send_request1);
  free(send_status1);

  return 0;
}

int fix_radial_ghosts(gtc_bench_data_t *gtc_input, real *data, int mzeta, int dim){
  gtc_global_params_t     *params;
  gtc_field_data_t        *field_data;
  gtc_radial_decomp_t     *radial_decomp;

  params            = &(gtc_input->global_params);
  field_data        = &(gtc_input->field_data);
  radial_decomp     = &(gtc_input->radial_decomp);

  int mpsi = params->mpsi;

  const int*  __restrict__ igrid;
  const int*  __restrict__ mtheta;

  igrid = field_data->igrid;
  mtheta = field_data->mtheta;

  int igrid_in = radial_decomp->igrid_in;
  int *ghost_comm_list = radial_decomp->ghost_comm_list;
  int *ghost_start = radial_decomp->ghost_start;
  int *ghost_end = radial_decomp->ghost_end;
  int ghost_comm_num = radial_decomp->ghost_comm_num;
  real *ghost_sendrecvbuf = radial_decomp->ghost_sendrecvbuf;
  int ghost_bufsize = radial_decomp->ghost_bufsize;
  int *nghost_comm_list = radial_decomp->nghost_comm_list;
  int *nghost_start = radial_decomp->nghost_start;
  int *nghost_end = radial_decomp->nghost_end;
  real *nghost_sendrecvbuf = radial_decomp->nghost_sendrecvbuf;
  int nghost_comm_num = radial_decomp->nghost_comm_num;
  int nghost_bufsize = radial_decomp->nghost_bufsize;
  int nproc_radial_partd = radial_decomp->nproc_radial_partd;
  int ipsi_nover_in = radial_decomp->ipsi_nover_in;
  int ipsi_nover_out = radial_decomp->ipsi_nover_out;

  int proc_num, count, ipsi_s,ipsi_e, offset, offset2;
  MPI_Request *recv_request = (MPI_Request *) malloc(ghost_comm_num*sizeof(MPI_Request));
  MPI_Status  *recv_status = (MPI_Status *) malloc(ghost_comm_num*sizeof(MPI_Status));

  MPI_Request *send_request = (MPI_Request *) malloc(nghost_comm_num*sizeof(MPI_Request));
  MPI_Status  *send_status = (MPI_Status *) malloc(nghost_comm_num*sizeof(MPI_Status));

  offset = 0;
  for (int i=0; i<ghost_comm_num; i++){
    proc_num = ghost_comm_list[i];
    ipsi_s = ghost_start[i];
    ipsi_e = ghost_end[i];

    count = (igrid[ipsi_e] + mtheta[ipsi_e] - igrid[ipsi_s] + 1)*(1+mzeta)*dim;

    MPI_Irecv(ghost_sendrecvbuf+offset, count, MPI_MYREAL, proc_num, proc_num,
              radial_decomp->radiald_comm, &recv_request[i]);
    offset += count;
  }
  assert(offset == ghost_bufsize*(mzeta+1)*dim);

  for (int i=0; i<nghost_comm_num; i++){
    proc_num = nghost_comm_list[i];
    ipsi_s = nghost_start[i];
    ipsi_e = nghost_end[i];
    offset = (igrid[ipsi_s] - igrid_in)*(1+mzeta)*dim;
    count = (igrid[ipsi_e] + mtheta[ipsi_e] - igrid[ipsi_s] + 1)*(1+mzeta)*dim;

    MPI_Isend(data+offset, count, MPI_MYREAL, proc_num, radial_decomp->myrank_radiald,
              radial_decomp->radiald_comm, &send_request[i]);
  }

  MPI_Waitall(ghost_comm_num, recv_request, recv_status);
  MPI_Waitall(nghost_comm_num, send_request, send_status);

  offset = 0;
  for (int i=0; i<ghost_comm_num; i++){
    proc_num = ghost_comm_list[i];
    ipsi_s = ghost_start[i];
    ipsi_e = ghost_end[i];
    offset2 = (igrid[ipsi_s] - igrid_in)*(1+mzeta)*dim;
    count = (igrid[ipsi_e] + mtheta[ipsi_e] - igrid[ipsi_s] + 1)*(1+mzeta)*dim;
    //#pragma omp parallel for num_threads(NUM_THREADS)
#pragma omp parallel for
    for (int j=0; j<count; j++){
      data[j+offset2] = ghost_sendrecvbuf[j+offset];
    }
    offset += count;
  }
  assert(offset == ghost_bufsize*(mzeta+1)*dim);

  free(recv_request);
  free(recv_status);
  free(send_request);
  free(send_status);

  return 0;
}

#endif
