#include "pmi.h"
#include <sprockit/errors.h>
#include <sumi/message.h>
#include <sumi/transport.h>
#include <sstmac/software/api/api.h>
#include <sstmac/software/process/app.h>
#include <sstmac/software/process/operating_system.h>
#include <sstmac/software/process/thread.h>
#include <sstmac/libraries/sumi/sumi_transport.h>

extern sumi::transport* active_transport();

static sstmac::sw::app* current_app(){
  return sstmac::sw::operating_system::current_thread()->parent_app();
}

extern "C" int PMI_Init()
{
  return PMI_SUCCESS;
}

extern "C" int PMI_Finalize()
{
  return PMI_SUCCESS;
}

extern "C" int PMI_Get_rank(int *ret)
{
  *ret = active_transport()->rank();
  return PMI_SUCCESS;
}

extern "C" int PMI_Get_size(int* ret)
{
  *ret = active_transport()->nproc();
  return PMI_SUCCESS;
}

extern "C" int PMI_Allgather(void *in, void *out, int len)
{
  //for now, I know that any data from this will be completely ignored
  //so don't bother doing it
  return PMI_SUCCESS;
}


extern "C" int PMI_Barrier()
{
  auto api = active_transport();
  api->barrier(42);
  api->collective_block(sumi::collective::barrier, 42);
  return PMI_SUCCESS;
}

extern "C" int 
PMI_Get_numpes_on_smp(int* num)
{
  *num = 1; //for now
  return PMI_SUCCESS;
}

extern "C" int
PMI2_KVS_Put(const char key[], const char value[])
{
  return PMI_SUCCESS;
}

extern "C" int
PMI2_KVS_Get(const char *jobid, int src_pmi_id, const char key[], char value [], int maxvalue, int *vallen)
{
  return PMI_SUCCESS;
}

extern "C" int
PMI2_KVS_Fence(void)
{
  return PMI_SUCCESS;
}

extern "C" int
PMI2_Abort(void)
{
  sprockit::abort("unimplemented: PMI2_Abort");
  return PMI_SUCCESS;
}

extern "C" int
PMI2_Job_GetId(char jobid[], int jobid_size)
{
  auto thr = sstmac::sw::operating_system::current_thread();
  ::sprintf(jobid, "%d", thr->aid());
  return PMI_SUCCESS;
}

extern "C" int 
PMI2_Init(int *spawned, int *size, int *rank, int *appnum)
{
  auto api = active_transport();
  *size = api->nproc();
  *rank = api->rank();
  *appnum = 0; //for now, we can't do multiple mpiexec launch
  *spawned = 0; //I have no idea what this is - zero!
  return PMI_SUCCESS;
}


