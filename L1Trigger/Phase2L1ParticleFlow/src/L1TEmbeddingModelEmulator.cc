#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TEmbeddingModelEmulator.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>

namespace l1ct {

  namespace {
    constexpr unsigned int MAX_PUPPI = 100;
    constexpr unsigned int N_FEATURES = 12;
    constexpr unsigned int N_CONT = 4;

    constexpr float PT_LOG_MAX = 8.3179f;
    constexpr float ETA_MIN = -5.0309f;
    constexpr float ETA_MAX = 5.0353f;
    constexpr float NORM_EPS = 1e-8f;

    inline float normalizePt(float pt) {
      float ptLog = std::log1p(pt > 0.0f ? pt : 0.0f);
      return (ptLog - 0.0f) / (PT_LOG_MAX - 0.0f + NORM_EPS);
    }
    inline float normalizeEta(float eta) { return (eta - ETA_MIN) / (ETA_MAX - ETA_MIN + NORM_EPS); }
    inline float normalizePhi(float phi) { return (phi + M_PI) / (2.0f * M_PI); }

    inline int pidToCategory(int pdgId) {
      switch (pdgId) {
        case -211:
          return 0;
        case -13:
          return 1;
        case -11:
          return 2;
        case 11:
          return 3;
        case 13:
          return 4;
        case 22:
          return 5;
        case 130:
          return 6;
        case 211:
          return 7;
        default:
          return -1;
      }
    }
  }  // namespace

  std::vector<unsigned int> l1tEmbeddingModelSortOrder(const std::vector<PuppiObjEmu> &particles) {
    std::vector<unsigned int> order;
    order.reserve(particles.size());
    for (unsigned int j = 0; j < particles.size(); ++j)
      order.push_back(j);
    std::sort(order.begin(), order.end(), [&particles](unsigned int a, unsigned int b) {
      const PuppiObjEmu &pa = particles[a];
      const PuppiObjEmu &pb = particles[b];
      if (pa.hwPt != pb.hwPt)
        return pa.hwPt > pb.hwPt;
      if (pa.hwEta != pb.hwEta)
        return pa.hwEta < pb.hwEta;
      return a < b;
    });
    return order;
  }

  void L1TEmbeddingModelEmu_top(const std::vector<PuppiObjEmu> &particles, L1TEmbeddingModelServer::input_t *inputOut) {
    std::fill(inputOut, inputOut + L1TEmbeddingModelServer::N_inp, L1TEmbeddingModelServer::input_t(0));
    std::vector<unsigned int> order = l1tEmbeddingModelSortOrder(particles);
    const unsigned int n = std::min<unsigned int>(order.size(), MAX_PUPPI);
    for (unsigned int i = 0; i < n; ++i) {
      const PuppiObjEmu &p = particles[order[i]];
      unsigned int off = i * N_FEATURES;
      inputOut[off + 0] = L1TEmbeddingModelServer::input_t(normalizePt(p.floatPt()));
      inputOut[off + 1] = L1TEmbeddingModelServer::input_t(normalizeEta(p.floatEta()));
      inputOut[off + 2] = L1TEmbeddingModelServer::input_t(normalizePhi(p.floatPhi()));
      inputOut[off + 3] = L1TEmbeddingModelServer::input_t(p.hwId.charged() ? p.floatDxy() : 0.0f);
      int cat = pidToCategory(p.pdgId());
      if (cat >= 0)
        inputOut[off + 4 + cat] = L1TEmbeddingModelServer::input_t(1);
    }
  }

  void L1TEmbeddingModel_emu(const std::vector<PuppiObjEmu> &particles,
                             L1TEmbeddingModelOutput &output,
                             const L1TEmbeddingModelServer &model,
                             L1TEmbeddingModelServer::input_t *inputOut) {
    std::array<L1TEmbeddingModelServer::input_t, L1TEmbeddingModelServer::N_inp> input;
    L1TEmbeddingModelEmu_top(particles, input.data());

    if (inputOut)
      std::copy(input.begin(), input.end(), inputOut);

    std::array<L1TEmbeddingModelServer::output_t, L1TEmbeddingModelServer::N_out> result;
    model.run(input.data(), result.data(), 1, 1);

    std::memcpy(
        output.embedding, result.data(), L1TEmbeddingModelServer::N_out * sizeof(L1TEmbeddingModelServer::output_t));
  }

}  // namespace l1ct
