#include "gni_pub.h"
#include "pmi.h"
#include <sprockit/errors.h>
#include <sumi/message.h>
#include <sumi/transport.h>
#include <sstmac/software/api/api.h>
#include <sstmac/software/process/operating_system.h>
#include <sstmac/software/process/thread.h>
#include <sstmac/libraries/sumi/sumi_transport.h>

MakeDebugSlot(ugni);

#define debug(tport, ...) \
  debug_printf(sprockit::dbg::ugni, "GNI Rank %d: %s", tport->rank(), \
    sprockit::printf(__VA_ARGS__).c_str())

class gni_message : public sumi::message {
 public:
  void* data() const {
    switch(payload_type()){
      case header:
      case eager_payload:
      case rdma_put:
        return local_buffer().ptr;
        break;
      case rdma_get:
        return remote_buffer().ptr;
        break;
      case software_ack:
      case nvram_get:
      case rdma_put_ack:
      case rdma_get_ack:
      case rdma_get_nack:
      case eager_payload_ack:
      case failure:
      case none:
        spkt_abort_printf("bad payload type %s in fetching gni message data", tostr(payload_type()));
        return nullptr;
    }
  }

  uint64_t inst_id() const {
    return inst_id_;
  }

  void set_inst_id(uint64_t id) {
    inst_id_ = id;
  }

  gni_post_type_t type() const {
    return type_;
  }

  void set_type(gni_post_type_t ty) {
    type_ = ty;
  }

  uint64_t msg_id() const {
    return msg_id_;
  }

  void set_msg_id(uint64_t msg_id) {
    msg_id_ = msg_id;
  }

  message* clone(payload_type_t ty) const override {
    gni_message* msg = new gni_message;
    clone_into(msg);
    return msg;
  }

 protected:
  void clone_into(gni_message* msg) const {
    msg->type_ = type_;
    msg->inst_id_ = inst_id_;
    msg->msg_id_ = msg_id_;
    sumi::message::clone_into(msg);
  }

 private:
  uint64_t inst_id_;
  gni_post_type type_;
  uint64_t msg_id_;
};

class gni_rdma_message : public gni_message {
 public:
  gni_rdma_message() : sync_flag_addr_(0), pd_(nullptr) {}

  void get_completed(gni_post_descriptor_t** descr){
    switch(payload_type()){
    case eager_payload_ack:
    case rdma_put_ack:
    case rdma_get_ack:
      //acks must return the original pd
      *descr = pd_;
      break;
    default: {
      gni_post_descriptor* pd = new gni_post_descriptor;
      pd->first_operand = first_op_;
      pd->second_operand = second_op_;
      pd->cqwrite_value = cqwrite_val_;
      pd->local_addr = (uint64_t) local_buffer().ptr;
      pd->remote_addr = (uint64_t) remote_buffer().ptr;
      pd->length = byte_length();
      pd->sync_flag_value = sync_flag_value_;
      pd->sync_flag_addr = sync_flag_addr_;
      *descr = pd;
      break;
    }
    }
  }

  ~gni_rdma_message(){}

  void write_sync_value() override {
    if (sync_flag_addr_){
      uint64_t* ptr = (uint64_t*) sync_flag_addr_;
      *ptr = sync_flag_value_;
    }
  }

  void set_post_descr(gni_post_descriptor_t* pd){
    pd_ = pd;
  }

  void set_first_operand(uint64_t op){
    first_op_ = op;
  }

  void set_second_operand(uint64_t op){
    second_op_ = op;
  }

  void set_sync_flag_value(uint64_t f){
    sync_flag_value_ = f;
  }

  void set_sync_flag_addr(uint64_t a){
    sync_flag_addr_ = a;
  }

  void set_cqwrite(uint64_t val){
    cqwrite_val_ = val;
  }

  message* clone(payload_type_t ty) const override {
    gni_rdma_message* msg = new gni_rdma_message;
    switch (ty){
    case eager_payload_ack:
    case rdma_put_ack:
    case rdma_get_ack:
      msg->pd_ = pd_;
      //acks only need the original post descriptor
      gni_message::clone_into(msg);
      break;
    default:
      clone_into(msg);
      break;
    }
    return msg;
  }

 protected:
  void clone_into(gni_rdma_message* msg) const {
    msg->first_op_ = first_op_;
    msg->second_op_ = second_op_;
    msg->cqwrite_val_ = cqwrite_val_;
    gni_message::clone_into(msg);
  }

 private:
  gni_post_descriptor_t* pd_;
  uint64_t first_op_;
  uint64_t second_op_;
  uint64_t sync_flag_value_;
  uint64_t sync_flag_addr_;
  uint64_t cqwrite_val_;
};

class gni_smsg_message : public gni_message {
 public:
  void set_tag(int tag){
    tag_ = tag;
  }

  int tag() const {
    return tag_;
  }

  message* clone(payload_type_t ty) const override {
    gni_smsg_message* msg = new gni_smsg_message;
    clone_into(msg);
    return msg;
  }

 protected:
  void clone_into(gni_smsg_message* msg) const {
    msg->tag_ = tag_;
    gni_message::clone_into(msg);
  }

 private:
  int tag_;
};

class ugni_transport : public sstmac::sumi_transport
{
  RegisterAPI("ugni", ugni_transport)
 public:
  ugni_transport(sprockit::sim_parameters* params,
                 sstmac::sw::software_id sid,
                 sstmac::sw::operating_system* os) :
    sstmac::sumi_transport(params, sid, os),
    pending_smsg_(4) //leave room for 4 cq for now, this can grow
  {
  }

  void init() override {
    sstmac::sumi_transport::init();
    //these are not ever used, but are needed for faking out real implementations
    setenv("PMI_GNI_LOC_ADDR", "0", 0);
    setenv("PMI_GNI_PTAG", "0", 0);
    setenv("PMI_GNI_COOKIE", "0", 0);
    setenv("PMI_GNI_DEV_ID", "0", 0);
  }

  int allocate_cq(){
    int new_cq = sstmac::sumi_transport::allocate_cq();
    if (new_cq >= pending_smsg_.size()){
      pending_smsg_.resize(new_cq+1);
    }
    return new_cq;
  }

  void set_pending_smsg(int cq_id, gni_smsg_message* smsg){
    pending_smsg_[cq_id] = smsg;
  }

  gni_smsg_message* take_pending_smsg(int cq_id){
    auto ret = pending_smsg_[cq_id];
    pending_smsg_[cq_id] = nullptr;
    return ret;
  }

  gni_return_t cg_get_event(gni_cq_handle_t cq_hndl, gni_cq_entry_t* event_data){
    sumi::message* msg = blocking_poll(cq_hndl);
    intptr_t msgPtr = (intptr_t) msg;
    *event_data = msgPtr;
    return GNI_RC_SUCCESS;
  }

 private:
  std::vector<gni_smsg_message*> pending_smsg_;
};

ugni_transport* sstmac_ugni()
{
  sstmac::sw::thread* t = sstmac::sw::operating_system::current_thread();
  return t->get_api<ugni_transport>();
}

sumi::transport* active_transport()
{
  return sstmac_ugni(); 
}

extern "C" gni_return_t
GNI_GetCompleted(
  gni_cq_handle_t cq_hndl,
  gni_cq_entry_t event_data,
  gni_post_descriptor_t **post_descr)
{
  debug(sstmac_ugni(), "GNI_CqGetCompleted(...)");
  gni_message* e = (gni_message*) event_data;
  auto rdma_msg = static_cast<gni_rdma_message*>(e);
  rdma_msg->get_completed(post_descr);
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t
GNI_CqGetEvent(
  gni_cq_handle_t cq_hndl,
  gni_cq_entry_t* event_data) {
  auto tport = sstmac_ugni();
  double timeout = 10e-3; //try for at least 10 ms
  sumi::message* msg = tport->poll(true, cq_hndl, timeout);
  if (msg){
    debug(tport, "GNI_CqGetEvent(...)");
    intptr_t msg_ptr = (intptr_t) static_cast<gni_message*>(msg);
    *event_data = msg_ptr;
    return GNI_RC_SUCCESS;
  } else {
    return GNI_RC_NOT_DONE;
  }
}

extern "C" uint64_t gni_cq_get_data(gni_cq_entry_t entry){
  gni_message* e = (gni_message*) entry;
  return (uint64_t) e->data();
}

extern "C" uint64_t gni_cq_get_source(gni_cq_entry_t entry){
  sprockit::abort("gni_cq_get_source: not implemented");
  return 0;
}

extern "C" uint64_t gni_cq_get_status(gni_cq_entry_t){
  sprockit::abort("gni_cq_get_status: not implemented");
  return 0;
}

extern "C" uint64_t gni_cq_get_info(gni_cq_entry_t){
  sprockit::abort("gni_cq_get_info: not implemented");
  return 0;
}

extern "C" uint64_t gni_cq_overrun(gni_cq_entry_t){
  //always return zero - for now we can't have queue overruns
  return 0;
}

extern "C" uint64_t gni_cq_rem_overrun(gni_cq_entry_t){
  sprockit::abort("gni_cq_rem_overrun: not implemented");
  return 0;
}

extern "C" uint64_t gni_cq_get_inst_id(gni_cq_entry_t entry){
  gni_message* e = (gni_message*) entry;
  return e->inst_id();
}

extern "C" uint64_t gni_cq_get_rem_inst_id(gni_cq_entry_t entry){
  gni_message* e = (gni_message*) entry;
  return e->inst_id();
}

extern "C" uint64_t gni_cq_get_tid(gni_cq_entry_t entry){
  sprockit::abort("gni_cq_get_tid: not implemented");
  return 0;
}

extern "C" uint64_t gni_cq_get_msg_id(gni_cq_entry_t entry){
  gni_message* e = (gni_message*) entry;
  return e->msg_id();
}

extern "C" uint64_t gni_cq_get_type(gni_cq_entry_t entry){
  gni_message* e = (gni_message*) entry;
  return e->type();
}

extern "C" uint64_t gni_cq_get_block_id(gni_cq_entry_t){
  sprockit::abort("gni_cq_get_block_id: not implemented");
  return 0;
}

extern "C" uint64_t gni_cq_get_unsuccessful_cnt(gni_cq_entry_t){
  sprockit::abort("gni_cq_get_unsuccessful_cnt: not implemented");
  return 0;
}

extern "C" uint64_t gni_cq_get_marker_id(gni_cq_entry_t){
  sprockit::abort("gni_cq_get_marker_id: not implemented");
  return 0;
}

extern "C" uint64_t gni_cq_get_failed_enqueue_cnt(gni_cq_entry_t){
  sprockit::abort("gni_cq_get_failed_enqueue_cnt: not implemented");
  return 0;
}

extern "C" uint64_t gni_cq_get_ce_id(gni_cq_entry_t){
  sprockit::abort("gni_cq_get_ce_id: not implemented");
  return 0;
}

extern "C" uint64_t gni_cq_get_reductn_id(gni_cq_entry_t){
  sprockit::abort("gni_cq_get_reductn_id: not implemented");
  return 0;
}

extern "C" uint64_t gni_cq_get_trans_type(gni_cq_entry_t){
  sprockit::abort("gni_cq_get_trans_type: not implemented");
  return 0;
}

extern "C" void     gni_cq_set_inst_id(gni_cq_entry_t* entry, uint64_t id){
  gni_message* e = (gni_message*) entry;
  e->set_inst_id(id);
}

extern "C" void     gni_cq_set_rem_inst_id(gni_cq_entry_t* entry, uint64_t id){
  gni_message* e = (gni_message*) entry;
  e->set_inst_id(id);
}

extern "C" void     gni_cq_set_tid(gni_cq_entry_t *,uint64_t){
  sprockit::abort("gni_cq_set_tid: not implemented");
}

extern "C" void     gni_cq_set_msg_id(gni_cq_entry_t* entry, uint64_t id){
  gni_message* e = (gni_message*) entry;
  e->set_msg_id(id);
}

extern "C" void     gni_cq_set_type(gni_cq_entry_t* entry, uint64_t ty){
  gni_message* e = (gni_message*) entry;
  e->set_type((gni_post_type_t)ty);
}

extern "C" void     gni_cq_clr_status(gni_cq_entry_t *){
  sprockit::abort("gni_cq_clr_status: not implemented");
}

extern "C" unsigned gni_cq_status_dla_overflow(gni_cq_entry_t){
  sprockit::abort("gni_cq_status_dla_overflow: not implemented");
  return 0;
}

extern "C" unsigned gni_cq_bte_enq_status(gni_cq_entry_t){
  sprockit::abort("gni_cq_bte_enq_status: not implemented");
  return 0;
}

extern "C" gni_return_t
GNI_CdmCreate(uint32_t inst_id, uint8_t ptag, uint32_t cookie,
              uint32_t modes, gni_cdm_handle_t *cdm_hndl){
  cdm_hndl->inst_id = inst_id;
  if (GNI_CDM_MODE_DUAL_EVENTS | modes){
    cdm_hndl->dual_events = true;
  } else {
    cdm_hndl->dual_events = false;
  }
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_EpBind(
  gni_ep_handle_t      ep_hndl,
  uint32_t             remote_addr,
  uint32_t             remote_id)
{
  ep_hndl->ep_id = remote_id;
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_EpUnbind(
  gni_ep_handle_t      ep_hndl)
{
  debug(sstmac_ugni(), "GNI_EpUnbind(ep=%d,...)", ep_hndl->ep_id);
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_EpIdle(gni_ep_handle_t ep_hndl)
{
  return GNI_RC_SUCCESS; //no op
}

extern "C" gni_return_t GNI_CdmHold(
  uint8_t ptag,
  uint32_t cookie,
  int *nic_fd_ret)
{
  *nic_fd_ret = 0; //basically a no-op
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t
GNI_CdmRelease(int nic_fd){
  return GNI_RC_SUCCESS;
} //no-op


extern "C" gni_return_t GNI_CdmAttach(
  gni_cdm_handle_t    cdm_hndl,
  uint32_t            device_id,
  uint32_t            *local_addr,
  gni_nic_handle_t    *nic_hndl)
{
 //no-op, nic handles do nothing in SST
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t
GNI_CdmDestroy(gni_cdm_handle_t cdm_hndl)
{
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t
GNI_CdmGetNicAddress(uint32_t device_id, uint32_t* address, uint32_t* cpu_id){
  *cpu_id = 0;
  //don't do physical NICs here, just do ranks
  *address = sstmac_ugni()->rank();
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_EpCreate(
  IN  gni_nic_handle_t    nic_hndl,
  IN  gni_cq_handle_t     src_cq_hndl,
  OUT gni_ep_handle_t     *ep_hndl){
  debug(sstmac_ugni(), "GNI_EpCreate(nic=%d,cq=%d)", nic_hndl.ep_id, src_cq_hndl);
  gni_ep_t* ep = new gni_ep_t;
  ep->cq_id = src_cq_hndl;
  ep->ep_id = nic_hndl.ep_id;
  *ep_hndl = ep;
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_EpDestroy(
  IN gni_ep_handle_t      ep_hndl){
  gni_ep_t* ep = ep_hndl;
  delete ep;
  return GNI_RC_SUCCESS;
}

extern "C" uint64_t gni_ce_res_get_status(gni_ce_result_t *){
  sprockit::abort("gni_ce_res_get_status: not implemented");
  return 0;
}

extern "C" uint64_t gni_ce_res_status_ok(gni_ce_result_t *){
  sprockit::abort("gni_ce_res_status_ok: not implemented");
  return 0;
}

extern "C" uint64_t gni_ce_res_get_fpe(gni_ce_result_t *){
  sprockit::abort("gni_ce_res_get_fpe: not implemented");
  return 0;
}

extern "C" uint64_t gni_ce_res_get_red_id(gni_ce_result_t *){
  sprockit::abort("gni_ce_res_get_red_id: not implemented");
  return 0;
}

extern "C" gni_return_t GNI_DqueueInit(
  gni_nic_handle_t nic_hndl,
  gni_dqueue_in_attr_t *in_attrs,
  gni_dqueue_out_attr_t *out_attrs,
  gni_dqueue_handle_t *dqueue) {
  spkt_abort_printf("unimplemented: GNI_DqueueInit()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_DqueueFini(
  gni_dqueue_handle_t dqueue) {
  spkt_abort_printf("unimplemented: GNI_DqueueFini()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_DqueueConnect(
  gni_dqueue_handle_t dqueue,
  gni_dqueue_out_attr_t *attrs,
  uint32_t nattrs) {
  spkt_abort_printf("unimplemented: GNI_DqueueConnect()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_DqueuePut(
  gni_dqueue_handle_t dqueue,
  void *data,
  uint32_t length,
  uint32_t pe_addr,
  uint32_t inst_id) {
  spkt_abort_printf("unimplemented: GNI_DqueuePut()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_DqueueProgress(
  gni_dqueue_handle_t dqueue) {
  spkt_abort_printf("unimplemented: GNI_DqueueProgress()");
  return GNI_RC_SUCCESS;
}

extern "C"  gni_return_t GNI_SuspendJob(
  uint32_t     device_id,
  uint64_t     job_id,
  uint8_t      ptag,
  uint32_t     cookie,
  uint32_t     timeout) {
  spkt_abort_printf("unimplemented: GNI_SuspendJob()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_ResumeJob(
  uint32_t     device_id,
  uint64_t     job_id,
  uint8_t      ptag,
  uint32_t     cookie) {
  spkt_abort_printf("unimplemented: GNI_ResumeJob()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_ConfigureNTT(
  int                         device_id,
  gni_ntt_descriptor_t        *ntt_desc,
  uint32_t                    *ntt_base) {
  spkt_abort_printf("unimplemented: GNI_ConfigureNTT()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_ConfigureNTTandJob(
  int                         device_id,
  uint64_t                    job_id,
  uint8_t                     ptag,
  uint32_t                    cookie,
  gni_job_limits_t            *limits,
  gni_ntt_descriptor_t        *ntt_desc,
  uint32_t                    *ntt_base) {
  spkt_abort_printf("unimplemented: GNI_ConfigureNTTandJob()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_GetCapabilities(
  gni_revision_info_t         *revision_info) {
  spkt_abort_printf("unimplemented: GNI_GetCapabilities()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_ValidateCapabilities(
  gni_revision_info_t         first_revision_info,
  gni_revision_info_t         second_revision_info) {
  spkt_abort_printf("unimplemented: GNI_ValidateCapabilities()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_EpPostData(
  gni_ep_handle_t      ep_hndl,
  void                 *in_data,
  uint16_t             data_len,
  void                 *out_buf,
  uint16_t             buf_size) {
  spkt_abort_printf("unimplemented: GNI_EpPostData()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_EpPostDataWId(
  gni_ep_handle_t      ep_hndl,
  void                 *in_data,
  uint16_t             data_len,
  void                 *out_buf,
  uint16_t             buf_size,
  uint64_t             datagram_id) {
  spkt_abort_printf("unimplemented: GNI_EpPostDataWId()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_EpPostDataTest(
  gni_ep_handle_t     ep_hndl,
  gni_post_state_t    *post_state,
  uint32_t            *remote_addr,
  uint32_t            *remote_id) {
  spkt_abort_printf("unimplemented: GNI_EpPostDataTest()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_EpPostDataTestById(
  gni_ep_handle_t     ep_hndl,
  uint64_t            datagram_id,
  gni_post_state_t    *post_state,
  uint32_t            *remote_addr,
  uint32_t            *remote_id) {
  spkt_abort_printf("unimplemented: GNI_EpPostDataTestById()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_EpPostDataWait(
  gni_ep_handle_t     ep_hndl,
  uint32_t            timeout,
  gni_post_state_t    *post_state,
  uint32_t            *remote_addr,
  uint32_t            *remote_id) {
  spkt_abort_printf("unimplemented: GNI_EpPostDataWait()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_EpPostDataWaitById(
  gni_ep_handle_t     ep_hndl,
  uint64_t            datagram_id,
  uint32_t            timeout,
  gni_post_state_t    *post_state,
  uint32_t            *remote_addr,
  uint32_t            *remote_id) {
  spkt_abort_printf("unimplemented: GNI_EpPostDataWaitById()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_PostDataProbe(
  gni_nic_handle_t    nic_hndl,
  uint32_t            *remote_addr,
  uint32_t            *remote_id) {
  spkt_abort_printf("unimplemented: GNI_PostDataProbe()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_PostDataProbeById(
  gni_nic_handle_t    nic_hndl,
  uint64_t            *datagram_id) {
  spkt_abort_printf("unimplemented: GNI_PostDataProbeById()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_PostdataProbeWaitById(
  gni_nic_handle_t    nic_hndl,
  uint32_t            timeout,
  uint64_t            *datagram_id) {
  spkt_abort_printf("unimplemented: GNI_PostdataProbeWaitById()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_EpPostDataCancel(
  gni_ep_handle_t      ep_hndl) {
  spkt_abort_printf("unimplemented: GNI_EpPostDataCancel()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_EpPostDataCancelById(
  gni_ep_handle_t      ep_hndl,
  uint64_t             datagram_id) {
  spkt_abort_printf("unimplemented: GNI_EpPostDataCancelById()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_MemRegister(
  gni_nic_handle_t  nic_hndl,
  uint64_t          address,
  uint64_t          length,
  gni_cq_handle_t   dst_cq_hndl,
  uint32_t          flags,
  uint32_t          vmdh_index,
  gni_mem_handle_t  *mem_hndl) {
  debug(sstmac_ugni(), "GNI_MemRegister(cq=%d,...)", dst_cq_hndl);
  mem_hndl->cq_handle = dst_cq_hndl;
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_MemRegisterSegments(
  gni_nic_handle_t  nic_hndl,
  gni_mem_segment_t *mem_segments,
  uint32_t          segments_cnt,
  gni_cq_handle_t   dst_cq_hndl,
  uint32_t          flags,
  uint32_t          vmdh_index,
  gni_mem_handle_t  *mem_hndl) {
  mem_hndl->cq_handle = dst_cq_hndl;
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_SetMddResources(
  gni_nic_handle_t     nic_hndl,
  uint32_t             num_entries) {
  spkt_abort_printf("unimplemented: GNI_SetMddResources()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_MemDeregister(
  gni_nic_handle_t     nic_hndl,
  gni_mem_handle_t     *mem_hndl) {
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_MemHndlQueryAttr(
  gni_mem_handle_t            *mem_hndl,
  gni_mem_handle_attr_t       attr,
  int                         *yesno) {
  spkt_abort_printf("unimplemented: GNI_MemHndlQueryAttr()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_RebuildMemHndl (
  gni_mem_handle_t    *src_mem_hndl,
  uint32_t            vmdh_index,
  gni_mem_handle_t    *dst_mem_hndl) {
  spkt_abort_printf("unimplemented: GNI_RebuildMemHndl ()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_MemQueryHndls(
  gni_nic_handle_t  nic_hndl,
  int               fd,
  gni_mem_handle_t *mem_hndl,
  uint64_t         *address,
  uint64_t         *length) {
  spkt_abort_printf("unimplemented: GNI_MemQueryHndls()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CqCreate(
  gni_nic_handle_t    nic_hndl,
  uint32_t            entry_count,
  uint32_t            delay_count,
  gni_cq_mode_t       mode,
  void                (*handler)(gni_cq_entry_t *,void *),
  void                *context,
  gni_cq_handle_t     *cq_hndl) {
  if (context || handler){
    spkt_abort_printf("cannot create completion queue with callback function");
  }
  auto tport = sstmac_ugni();
  *cq_hndl = tport->allocate_cq();
  debug(tport, "GNI_CqCreate(...)");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CqDestroy(
  gni_cq_handle_t      cq_hndl) {
  return GNI_RC_SUCCESS;
}

extern bool pmi_finalized;

static gni_return_t post_rdma_op(
  gni_ep_handle_t              ep_hndl,
  gni_post_descriptor_t        *post_descr) 
{
  ugni_transport* api = sstmac_ugni();
  gni_rdma_message* msg = new gni_rdma_message;
  msg->set_inst_id(api->rank());
  msg->set_post_descr(post_descr);
  msg->set_type(post_descr->type);
  msg->set_byte_length(post_descr->length);
  msg->set_first_operand(post_descr->first_operand);
  msg->set_second_operand(post_descr->second_operand);
  msg->set_cqwrite(post_descr->cqwrite_value);
  msg->remote_buffer().ptr = (void*) post_descr->remote_addr;
  msg->local_buffer().ptr = (void*) post_descr->local_addr;
  gni_cq_handle_t src_cq = post_descr->src_cq_hndl
      ? post_descr->src_cq_hndl
      : ep_hndl->cq_id;
  switch(post_descr->type)
  {
  case GNI_POST_FMA_GET:
  case GNI_POST_RDMA_GET:
  {
    bool recv_ack = (post_descr->cq_mode | GNI_CQMODE_GLOBAL_EVENT)
                    || (post_descr->cq_mode | GNI_CQMODE_LOCAL_EVENT);
    bool send_ack = post_descr->cq_mode | GNI_CQMODE_REMOTE_EVENT;
    api->rdma_get(ep_hndl->ep_id, msg,
                  send_ack ? post_descr->remote_mem_hndl.cq_handle : sumi::message::no_ack,
                  recv_ack ? src_cq : sumi::message::no_ack);
    return GNI_RC_SUCCESS;
  }
  case GNI_POST_FMA_PUT_W_SYNCFLAG:
    msg->set_sync_flag_value(post_descr->sync_flag_value);
    msg->set_sync_flag_addr(post_descr->sync_flag_addr);
  case GNI_POST_RDMA_PUT:
  case GNI_POST_FMA_PUT:
  {
    bool send_ack = (post_descr->cq_mode | GNI_CQMODE_GLOBAL_EVENT)
                    || (post_descr->cq_mode | GNI_CQMODE_LOCAL_EVENT);
    bool recv_ack = post_descr->cq_mode | GNI_CQMODE_REMOTE_EVENT;
    api->rdma_put(ep_hndl->ep_id, msg,
                  send_ack ? src_cq : sumi::message::no_ack,
                  recv_ack ? post_descr->remote_mem_hndl.cq_handle : sumi::message::no_ack);
    return GNI_RC_SUCCESS;
  } 
  default:
    return GNI_RC_INVALID_PARAM;
  }
}

extern "C" gni_return_t GNI_PostRdma(
  gni_ep_handle_t              ep_hndl,
  gni_post_descriptor_t        *post_descr) {
  gni_rdma_message* msg = new gni_rdma_message;
  debug(sstmac_ugni(), "GNI_PostRdma(...)");
  return post_rdma_op(ep_hndl, post_descr);
}



extern "C" gni_return_t GNI_PostFma(
  gni_ep_handle_t              ep_hndl,
  gni_post_descriptor_t        *post_descr) {
  gni_rdma_message* msg = new gni_rdma_message;
  debug(sstmac_ugni(), "GNI_PostFma(...)");
  return post_rdma_op(ep_hndl, post_descr);
}

extern "C" gni_return_t GNI_CtPostFma(
  gni_ep_handle_t              ep_hndl,
  gni_post_descriptor_t        *post_descr) {
  spkt_abort_printf("unimplemented: GNI_CtPostFma()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_PostCqWrite(
  gni_ep_handle_t              ep_hndl,
  gni_post_descriptor_t        *post_descr) {
  spkt_abort_printf("unimplemented: GNI_PostCqWrite()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CtPostCqWrite(
  gni_ep_handle_t              ep_hndl,
  gni_post_descriptor_t        *post_descr) {
  spkt_abort_printf("unimplemented: GNI_CtPostCqWrite()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CqWaitEvent(
  gni_cq_handle_t     cq_hndl,
  uint64_t            timeout,
  gni_cq_entry_t      *event_data) {
  spkt_abort_printf("unimplemented: GNI_CqWaitEvent()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CqVectorWaitEvent(
  gni_cq_handle_t     *cq_hndls,
  uint32_t            num_cqs,
  uint64_t            timeout,
  gni_cq_entry_t      *event_data,
  uint32_t            *which) {
  spkt_abort_printf("unimplemented: GNI_CqVectorWaitEvent()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CqVectorMonitor(
  gni_cq_handle_t     *cq_hndls,
  uint32_t            num_cqs,
  uint64_t            timeout,
  uint32_t            *which) {
  spkt_abort_printf("unimplemented: GNI_CqVectorMonitor()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CqInterruptMask(
  gni_cq_handle_t cq_hndl) {
  spkt_abort_printf("unimplemented: GNI_CqInterruptMask()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CqInterruptUnmask(
  gni_cq_handle_t cq_hndl) {
  spkt_abort_printf("unimplemented: GNI_CqInterruptUnmask()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CqTestEvent(
  gni_cq_handle_t      cq_hndl) {
  spkt_abort_printf("unimplemented: GNI_CqTestEvent()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CqErrorStr(
  gni_cq_entry_t      entry,
  void                *buffer,
  uint32_t            len) {
  spkt_abort_printf("unimplemented: GNI_CqErrorStr()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CqErrorRecoverable(
  gni_cq_entry_t      entry,
  uint32_t            *recoverable) {
  spkt_abort_printf("unimplemented: GNI_CqErrorRecoverable()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_SmsgBufferSizeNeeded(
  gni_smsg_attr_t     *smsg_attr,
  unsigned int        *size) {
  //don't do anything interesting with this
  //don't actually allocate any buffers
  *size = 8;
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_SmsgInit(
  gni_ep_handle_t      ep_hndl,
  gni_smsg_attr_t      *local_smsg_attr,
  gni_smsg_attr_t      *remote_smsg_attr) {
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_SmsgSetDeliveryMode(
  gni_nic_handle_t	nic_handle,
  uint16_t 		dlvr_mode) {
  //no-op for now
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_SmsgSend(
  gni_ep_handle_t      ep_hndl,
  void                 *header,
  uint32_t             header_length,
  void                 *data,
  uint32_t             data_length,
  uint32_t             msg_id) {
  return GNI_SmsgSendWTag(ep_hndl, header, header_length, data, data_length, msg_id, -1);
}

extern "C" gni_return_t GNI_SmsgSendWTag(
  gni_ep_handle_t      ep_hndl,
  void                 *header,
  uint32_t             header_length,
  void                 *data,
  uint32_t             data_length,
  uint32_t             msg_id,
  uint8_t              tag) {
  gni_smsg_message* msg = new gni_smsg_message;
  msg->set_byte_length(data_length + header_length);
  msg->set_tag(tag);
  if (header || data){
    char* data_buf = new char[header_length + data_length];
    if (header) ::memcpy(data_buf, header, header_length);
    if (data) ::memcpy(data_buf + header_length, data, data_length);
    msg->local_buffer().ptr = data_buf;
  }
  ugni_transport* api = sstmac_ugni();
  msg->set_msg_id(msg_id);
  msg->set_inst_id(api->rank());
  debug(api, "GNI_SmsgSendWTag(tag=%d,...)", tag);
  api->smsg_send(ep_hndl->ep_id, sumi::message::eager_payload, msg,
                 sumi::message::no_ack, ep_hndl->cq_id);
  return GNI_RC_SUCCESS;
}

static gni_smsg_message* poll_for_next_smsg(ugni_transport* api, int cq_id){
  double timeout = 10e-3; //block for at least 10 ms
  sumi::message* msg = api->poll(sumi::message::eager_payload, true, cq_id, timeout);
  if (msg->payload_type() == sumi::message::eager_payload){
    return dynamic_cast<gni_smsg_message*>(msg);
  } else {
    return nullptr;
  }
}

extern "C" gni_return_t GNI_SmsgGetNext(
  gni_ep_handle_t     ep_hndl,
  void                **header) {
  ugni_transport* api = sstmac_ugni();
  debug(api, "GNI_GetNext(...)")
  auto smsg = poll_for_next_smsg(api, ep_hndl->cq_id);
  if (smsg){
    *header = smsg->local_buffer().ptr;
    return GNI_RC_SUCCESS;
  } else {
    return GNI_RC_NOT_DONE;
  }
}

extern "C" gni_return_t GNI_SmsgGetNextWTag(
  gni_ep_handle_t     ep_hndl,
  void                **header,
  uint8_t             *tag) {
  ugni_transport* api = sstmac_ugni();
  debug(api, "GNI_GetNextWTag(tag=%d,...)", *tag)
  auto smsg = poll_for_next_smsg(api, ep_hndl->cq_id);
  if (smsg){
    if (*tag == GNI_SMSG_ANY_TAG){
      api->set_pending_smsg(ep_hndl->cq_id, smsg);
      *tag = smsg->tag();
      return GNI_RC_SUCCESS;
    } else if (*tag == smsg->tag()){
      api->set_pending_smsg(ep_hndl->cq_id, smsg);
      *tag = smsg->tag();
      return GNI_RC_SUCCESS;
    } else {
      return GNI_RC_NO_MATCH;
    }
  } else {
    return GNI_RC_NOT_DONE;
  }
}

extern "C" gni_return_t GNI_SmsgRelease(
  gni_ep_handle_t      ep_hndl) {
  auto api = sstmac_ugni();
  debug(api, "GNI_SmsgRelease(...)");
  gni_smsg_message* smsg = api->take_pending_smsg(ep_hndl->cq_id);
  if (!smsg){
    return GNI_RC_INVALID_STATE;
  } else {
    char* buffer = (char*) smsg->local_buffer().ptr;
    if (buffer) delete[] buffer;
    delete smsg;
  }
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_MsgqInit(
  gni_nic_handle_t            nic_hndl,
  gni_msgq_rcv_cb_func        *rcv_cb,
  void                        *cb_data,
  gni_cq_handle_t             snd_cq,
  gni_msgq_attr_t             *attrs,
  gni_msgq_handle_t           *msgq_hndl) {
  spkt_abort_printf("unimplemented: GNI_MsgqInit()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_MsgqRelease(
  gni_msgq_handle_t    msgq_hndl) {
  spkt_abort_printf("unimplemented: GNI_MsgqRelease()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_MsgqIdle(
  gni_msgq_handle_t    msgq_hndl) {
  spkt_abort_printf("unimplemented: GNI_MsgqIdle()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_MsgqGetConnAttrs(
  gni_msgq_handle_t   msgq_hndl,
  uint32_t            pe_addr,
  gni_msgq_ep_attr_t  *attrs,
  uint32_t            *attrs_size) {
  spkt_abort_printf("unimplemented: GNI_MsgqGetConnAttrs()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_MsgqConnect(
  gni_msgq_handle_t    msgq_hndl,
  uint32_t             pe_addr,
  gni_msgq_ep_attr_t   *attrs) {
  spkt_abort_printf("unimplemented: GNI_MsgqConnect()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_MsgqConnRelease(
  gni_msgq_handle_t    msgq_hndl,
  uint32_t             pe_addr) {
  spkt_abort_printf("unimplemented: GNI_MsgqConnRelease()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_MsgqSend(
  gni_msgq_handle_t    msgq_hndl,
  gni_ep_handle_t      ep_hndl,
  void                 *hdr,
  uint32_t             hdr_len,
  void                 *msg,
  uint32_t             msg_len,
  uint32_t             msg_id,
  uint8_t              msg_tag) {
  spkt_abort_printf("unimplemented: GNI_MsgqSend()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_MsgqProgress(
  gni_msgq_handle_t    msgq_hndl,
  uint32_t             timeout) {
  spkt_abort_printf("unimplemented: GNI_MsgqProgress()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_MsgqSize(
  gni_msgq_attr_t     *attrs,
  uint32_t            *size) {
  spkt_abort_printf("unimplemented: GNI_MsgqSize()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_SmsgSetMaxRetrans(
  gni_nic_handle_t     nic_handle,
  uint16_t             max_retrans) {
  spkt_abort_printf("unimplemented: GNI_SmsgSetMaxRetrans()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_SubscribeErrors(
  gni_nic_handle_t    nic_handle,
  uint32_t            device_id,
  gni_error_mask_t    mask,
  uint32_t            EEQ_size,
  gni_err_handle_t    *err_handle) {
  spkt_abort_printf("unimplemented: GNI_SubscribeErrors()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_ReleaseErrors(
  gni_err_handle_t     err_handle) {
  spkt_abort_printf("unimplemented: GNI_ReleaseErrors()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_GetErrorMask(
  gni_err_handle_t    err_handle,
  gni_error_mask_t    *mask) {
  spkt_abort_printf("unimplemented: GNI_GetErrorMask()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_SetErrorMask(
  gni_err_handle_t     err_handle,
  gni_error_mask_t     mask_in,
  gni_error_mask_t     *mask_out) {
  spkt_abort_printf("unimplemented: GNI_SetErrorMask()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_GetErrorEvent(
  gni_err_handle_t     err_handle,
  gni_error_event_t    *event) {
  spkt_abort_printf("unimplemented: GNI_GetErrorEvent()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_WaitErrorEvents(
  gni_err_handle_t    err_handle,
  gni_error_event_t   *events,
  uint32_t            events_size,
  uint32_t            timeout,
  uint32_t            *num_events) {
  spkt_abort_printf("unimplemented: GNI_WaitErrorEvents()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_SetErrorPtag(
  gni_err_handle_t     err_handle,
  uint8_t              ptag) {
  spkt_abort_printf("unimplemented: GNI_SetErrorPtag()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_GetNumLocalDevices(
  int *num_devices) {
  spkt_abort_printf("unimplemented: GNI_GetNumLocalDevices()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_GetLocalDeviceIds(
  int len,
  int *device_ids) {
  spkt_abort_printf("unimplemented: GNI_GetLocalDeviceIds()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_GetVersion(
  uint32_t    *version) {
  spkt_abort_printf("unimplemented: GNI_GetVersion()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_GetVersionInformation(
  gni_version_info_t  *version_info) {
  spkt_abort_printf("unimplemented: GNI_GetVersionInformation()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_GetDeviceType(
  gni_nic_device_t    *dev_type) {
  *dev_type = GNI_DEVICE_GEMINI;
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_GetDevResInfo(
  uint32_t            device_id,
  gni_dev_res_t       res_type,
  gni_dev_res_desc_t  *res_desc) {
  spkt_abort_printf("unimplemented: GNI_GetDevResInfo()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_GetJobResInfo(
  uint32_t            device_id,
  uint8_t             ptag,
  gni_job_res_t       res_type,
  gni_job_res_desc_t  *res_desc) {
  spkt_abort_printf("unimplemented: GNI_GetJobResInfo()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_SetJobResInfo(
  uint32_t      device_id,
  uint8_t       ptag,
  gni_job_res_t res_type,
  uint64_t      res_value) {
  spkt_abort_printf("unimplemented: GNI_SetJobResInfo()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_GetNttGran(
  uint32_t    device_id,
  uint32_t    *ntt_gran) {
  spkt_abort_printf("unimplemented: GNI_GetNttGran()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_GetPtag(
  uint32_t    device_id,
  uint32_t    cookie,
  uint8_t     *ptag) {
  spkt_abort_printf("unimplemented: GNI_GetPtag()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CeCreate(
  gni_nic_handle_t    nic_hndl,
  gni_ce_handle_t     *ce_hndl) {
  spkt_abort_printf("unimplemented: GNI_CeCreate()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CeGetId(
  gni_ce_handle_t     ce_hndl,
  uint32_t            *ce_id) {
  spkt_abort_printf("unimplemented: GNI_CeGetId()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_EpSetEventData(
  gni_ep_handle_t      ep_hndl,
  uint32_t             local_event,
  uint32_t             remote_event)
{
  gni_ep_t* ep = ep_hndl;
  ep->local_event = local_event;
  ep->remote_event = remote_event;
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_EpSetCeAttr(
  gni_ep_handle_t      ep_hndl,
  uint32_t             ce_id,
  uint32_t             child_id,
  gni_ce_child_t       child_type) {
  spkt_abort_printf("unimplemented: GNI_EpSetCeAttr()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CeConfigure(
  gni_ce_handle_t      ce_hndl,
  gni_ep_handle_t      *child_eps,
  uint32_t             num_child_eps,
  gni_ep_handle_t      parent_ep,
  gni_cq_handle_t      cq_hndl,
  uint32_t             modes) {
  spkt_abort_printf("unimplemented: GNI_CeConfigure()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CeCheckResult(
  gni_ce_result_t      *result,
  uint32_t             length) {
  spkt_abort_printf("unimplemented: GNI_CeCheckResult()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CeDestroy(
  gni_ce_handle_t      ce_hndl) {
  spkt_abort_printf("unimplemented: GNI_CeDestroy()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_SetBIConfig(
  uint32_t     device_id,
  uint16_t     bw,
  uint16_t     aot_bw,
  uint16_t     modes) {
  spkt_abort_printf("unimplemented: GNI_SetBIConfig()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_GetBIConfig(
  uint32_t             device_id,
  gni_bi_desc_t       *desc) {
  spkt_abort_printf("unimplemented: GNI_GetBIConfig()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_BISyncWait(
  uint32_t     device_id,
  uint32_t    timeout) {
  spkt_abort_printf("unimplemented: GNI_BISyncWait()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_GetNicStat(
  gni_nic_handle_t     nic_hndl,
  gni_statistic_t      stat,
  uint32_t            *value) {
  spkt_abort_printf("unimplemented: GNI_GetNicStat()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_ResetNicStat(
  gni_nic_handle_t     nic_hndl,
  gni_statistic_t      stat) {
  spkt_abort_printf("unimplemented: GNI_ResetNicStat()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CompChanCreate(
  gni_nic_handle_t             nic_hndl,
  gni_comp_chan_handle_t      *chan_hndl) {
  spkt_abort_printf("unimplemented: GNI_CompChanCreate()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CompChanDestroy(
  gni_comp_chan_handle_t chan_hndl) {
  spkt_abort_printf("unimplemented: GNI_CompChanDestroy()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CompChanFd(
  gni_comp_chan_handle_t       chan_hndl,
  int                         *comp_chan_fd) {
  spkt_abort_printf("unimplemented: GNI_CompChanFd()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CompChanGetEvent(
  gni_comp_chan_handle_t       chan_hndl,
  gni_cq_handle_t             *cq_hndl) {
  spkt_abort_printf("unimplemented: GNI_CompChanGetEvent()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CqAttachCompChan(
  gni_cq_handle_t              cq_hndl,
  gni_comp_chan_handle_t       chan_hndl) {
  spkt_abort_printf("unimplemented: GNI_CqAttachCompChan()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_CqArmCompChan(
  gni_cq_handle_t      *cq_hndls,
  uint32_t             num_cqs) {
  spkt_abort_printf("unimplemented: GNI_CqArmCompChan()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_SetDeviceOrbMMR(
  gni_nic_handle_t      nic_hndl,
  gni_dev_orb_mmr_t     orb_mmr_type,
  uint64_t              new_orb_mmr_value,
  uint64_t             *old_orb_mmr_value) {
  spkt_abort_printf("unimplemented: GNI_SetDeviceOrbMMR()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_ConfigureJob(
  uint32_t             device_id,
  uint64_t             job_id,
  uint8_t              ptag,
  uint32_t             cookie,
  gni_job_limits_t     *limits) {
  spkt_abort_printf("unimplemented: GNI_ConfigureJob()");
  return GNI_RC_SUCCESS;
}

extern "C" gni_return_t GNI_ConfigureJobFd(
  uint32_t             device_id,
  uint64_t             job_id,
  uint8_t              ptag,
  uint32_t             cookie,
  gni_job_limits_t     *limits,
  int                 *nic_fd) {
  spkt_abort_printf("unimplemented: GNI_ConfigureJobFd()");
  return GNI_RC_SUCCESS;
}


