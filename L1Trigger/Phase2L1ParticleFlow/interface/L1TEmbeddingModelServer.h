#ifndef L1TRIGGER_PHASE2L1PARTICLEFLOWS_L1TEMBEDDINGMODELSERVER_H
#define L1TRIGGER_PHASE2L1PARTICLEFLOWS_L1TEMBEDDINGMODELSERVER_H

#include <cstddef>
#include <memory>
#include <string>

#include "DataFormats/L1TParticleFlow/interface/datatypes.h"
#include "hls4ml/emulator.h"

class L1TEmbeddingModelServer {
public:
  using input_t = ap_fixed<11, 3>;   // model input feature (N_FEATURES=12 per particle)
  using output_t = ap_fixed<12, 2>;  // model output component (N_OUT=20 per event)

  explicit L1TEmbeddingModelServer(const std::string &modelName);

  L1TEmbeddingModelServer(const L1TEmbeddingModelServer &) = delete;
  L1TEmbeddingModelServer &operator=(const L1TEmbeddingModelServer &) = delete;

  void run(const input_t *input, output_t *output, size_t batchSize = 1, size_t nThreads = 1) const;

  static constexpr size_t N_inp = 1200;
  static constexpr size_t N_out = 20;

private:
  hls4mlEmulator::ModelLoader loader_;
  std::shared_ptr<hls4mlEmulator::Model> model_;
};

#endif
