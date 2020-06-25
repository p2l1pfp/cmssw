#ifndef L1Trigger_Phase2L1ParticleFlow_BitwisePuppiAlgo_h
#define L1Trigger_Phase2L1ParticleFlow_BitwisePuppiAlgo_h

#include "L1Trigger/Phase2L1ParticleFlow/interface/PUAlgoBase.h"

struct linpuppi_config;

namespace l1tpf_impl { 

  class BitwisePuppiAlgo : public PUAlgoBase {
    public:
        BitwisePuppiAlgo( const edm::ParameterSet& ) ;
        virtual ~BitwisePuppiAlgo() ;

        const std::vector<std::string> & puGlobalNames() const override ;
        void doPUGlobals(const std::vector<Region> &rs, float z0, float npu, std::vector<float> & globals) const override ;
        void runChargedPV(Region &r, float z0) const override ;
        void runNeutralsPU(Region &r, float z0, float npu, const std::vector<float> & globals) const override ;

    private:
        enum AlgoChoice { linpuppi_algo, linpuppi_flt_algo, fwdlinpuppi_algo, fwdlinpuppi_flt_algo } algo_;
        std::vector<linpuppi_config *> configs_; 
        bool debug_;      
  };

} // end namespace

#endif
