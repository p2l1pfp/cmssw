#ifndef L1Trigger_Phase2L1ParticleFlow_BitwisePFAlgo_h
#define L1Trigger_Phase2L1ParticleFlow_BitwisePFAlgo_h

#include "L1Trigger/Phase2L1ParticleFlow/interface/PFAlgoBase.h"

struct pfalgo_config;

namespace l1tpf_impl { 
class BitwisePFAlgo : public PFAlgoBase {
    public:
        BitwisePFAlgo( const edm::ParameterSet& ) ;
        ~BitwisePFAlgo() ;
        virtual void runPF(Region &r) const override;

    protected:
        enum AlgoChoice { algo3, algo2hgc } algo_;
        pfalgo_config * config_;
  };

} // end namespace

#endif
