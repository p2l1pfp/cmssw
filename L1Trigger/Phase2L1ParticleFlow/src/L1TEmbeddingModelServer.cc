#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TEmbeddingModelServer.h"

#include "ap_fixed.h"

#include <stdexcept>
#include <string>

using input_t = L1TEmbeddingModelServer::input_t;
using output_t = L1TEmbeddingModelServer::output_t;

struct ModelInputs {
  input_t *data;
};

struct ModelOutputs {
  output_t *data;
};

L1TEmbeddingModelServer::L1TEmbeddingModelServer(const std::string &modelName) : loader_(modelName) {
  try {
    model_ = loader_.load_model();
  } catch (std::runtime_error &e) {
    throw std::runtime_error("L1TEmbeddingModelServer: failed to load hls4ml model \"" + loader_.model_name() +
                             "\": " + e.what());
  }
}

void L1TEmbeddingModelServer::run(const input_t *input, output_t *output, size_t batchSize, size_t nThreads) const {
  (void)nThreads;
  for (size_t i = 0; i < batchSize; ++i) {
    ModelInputs modelInput = {const_cast<input_t *>(input + i * N_inp)};
    model_->prepare_input(modelInput);
    model_->predict();
    ModelOutputs modelOutput = {output + i * N_out};
    model_->read_result(&modelOutput);
  }
}
