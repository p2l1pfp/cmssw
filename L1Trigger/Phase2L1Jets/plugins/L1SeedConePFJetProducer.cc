#include "L1Trigger/Phase2L1Jets/interface/L1SeedConePFJetProducer.hh"

l1t::PFJet makeJet(std::vector<l1t::PFCandidate> parts){

    l1t::PFCandidate seed = parts.at(0);

    auto sumpt = [](float a, const l1t::PFCandidate& b){
        return a + b.pt();
    };

    // Sum the pt
    float pt = std::accumulate(parts.begin(), parts.end(), 0., sumpt);

    // pt weighted d eta
    std::vector<float> pt_deta;
    pt_deta.resize(parts.size());
    std::transform(parts.begin(), parts.end(), pt_deta.begin(), [&seed, &pt](const l1t::PFCandidate &part){
            return (part.pt() / pt) * (part.eta() - seed.eta());});
    // Accumulate the pt weighted etas. Init to the seed eta, start accumulating at begin()+1 to skip seed
    float eta = std::accumulate(pt_deta.begin()+1, pt_deta.end(), seed.eta());

    // pt weighted d phi
    std::vector<float> pt_dphi;
    pt_dphi.resize(parts.size());
    std::transform(parts.begin(), parts.end(), pt_dphi.begin(), [&seed, &pt](const l1t::PFCandidate &part){
            return (part.pt() / pt) * (part.phi() - seed.phi());});
    // Accumulate the pt weighted phis. Init to the seed phi, start accumulating at begin()+1 to skip seed
    float phi = std::accumulate(pt_dphi.begin()+1, pt_dphi.end(), seed.phi());

    l1t::PFJet jet(pt, eta, phi);

    return jet;
}

L1SeedConePFJetProducer::L1SeedConePFJetProducer(const edm::ParameterSet& cfg) : 
    _coneSize( cfg.getParameter<double>("coneSize")),
    _nJets( cfg.getParameter<unsigned>("nJets")),
    _l1PFToken( consumes<vector<l1t::PFCandidate>>(cfg.getParameter<edm::InputTag>("L1PFObjects")))
{
    produces< l1t::PFJetCollection >();
}

void L1SeedConePFJetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    std::unique_ptr<l1t::PFJetCollection> newPFJetCollection(new l1t::PFJetCollection);

    edm::Handle<  l1t::PFCandidateCollection > l1PFCandidates;
    iEvent.getByToken( _l1PFToken, l1PFCandidates);
    l1t::PFCandidateCollection work;

    std::copy((*l1PFCandidates).begin(), (*l1PFCandidates).end(), std::back_inserter(work));

    std::sort(work.begin(), work.end(), [](l1t::PFCandidate i,l1t::PFCandidate j){return(i.pt() > j.pt());});   

    std::vector<l1t::PFJet> jets;
    jets.reserve(_nJets);

    while(!work.empty() && jets.size() < _nJets){
        // Take the first (highest pt) candidate as a seed
        l1t::PFCandidate seed = work.at(0);
        // Get the particles within a _coneSize of the seed
        std::vector<l1t::PFCandidate> particlesInCone;
        std::copy_if(work.begin(), work.end(), std::back_inserter(particlesInCone), [&](const l1t::PFCandidate &part){
                return reco::deltaR<l1t::PFCandidate,l1t::PFCandidate>(seed, part) <= _coneSize;});
        jets.push_back(makeJet(particlesInCone));
        // remove the clustered particles
        work.erase(std::remove_if(work.begin(), work.end(), [&](const l1t::PFCandidate &part){
                return reco::deltaR<l1t::PFCandidate,l1t::PFCandidate>(seed, part) <= _coneSize;}), work.end());
    } 
    std::sort(jets.begin(), jets.end(), [](l1t::PFJet i,l1t::PFJet j){return(i.pt() > j.pt());}); 
    newPFJetCollection->swap(jets);
    iEvent.put( std::move(newPFJetCollection) );
}

/////////////
// DESTRUCTOR
L1SeedConePFJetProducer::~L1SeedConePFJetProducer()
{
}  

//////////
// END JOB
void L1SeedConePFJetProducer::endRun(const edm::Run& run, const edm::EventSetup& iSetup)
{
}

////////////
// BEGIN JOB
void L1SeedConePFJetProducer::beginRun(const edm::Run& run, const edm::EventSetup& iSetup )
{
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1SeedConePFJetProducer);
