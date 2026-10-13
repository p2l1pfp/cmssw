#ifndef L1TRIGGER_PHASE2L1PARTICLEFLOWS_L1TEMBEDDINGMODELEMU_H
#define L1TRIGGER_PHASE2L1PARTICLEFLOWS_L1TEMBEDDINGMODELEMU_H

#include "DataFormats/L1TParticleFlow/interface/layer1_emulator.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TEmbeddingModelServer.h"

namespace l1ct {

  struct L1TEmbeddingModelOutput {
    static constexpr size_t N_OUT = L1TEmbeddingModelServer::N_out;
    L1TEmbeddingModelServer::output_t embedding[N_OUT];
  };

  // Standalone "top function" emulation
  void L1TEmbeddingModelEmu_top(const std::vector<PuppiObjEmu> &particles, L1TEmbeddingModelServer::input_t *inputOut);

  std::vector<unsigned int> l1tEmbeddingModelSortOrder(const std::vector<PuppiObjEmu> &particles);

  void L1TEmbeddingModel_emu(const std::vector<PuppiObjEmu> &particles,
                             L1TEmbeddingModelOutput &output,
                             const L1TEmbeddingModelServer &model,
                             L1TEmbeddingModelServer::input_t *inputOut = nullptr);

}  // namespace l1ct

#endif
