// input: L1TkTracks and Muon Stubs 
//
// 
//
// Author: Vladimir Rekovic
// Version: 1   Date: 10.Feb.2019
// Descritiption: 
//      Match L1Tracks to EndCap muon stubs.
//      Stubs are taken to be inputs to EMTF, EMTFHitCollection
//      Use DynamicWindow matching algorithm which assumes muon 
//      detector object (here Stubs) have coordinates at 2nd station.  
//      Match the two and produce a collection of TkMuon

// user include files
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/L1TCorrelator/interface/TkMuon.h"
#include "DataFormats/L1TCorrelator/interface/TkMuonFwd.h"
#include "L1Trigger/L1TMuon/interface/MicroGMTConfiguration.h"
#include "L1Trigger/L1TTrackMatch/interface/L1TkMuCorrDynamicWindows.h"
#include "L1Trigger/L1TMuonEndCap/interface/Common.h"

// system include files
#include <memory>
#include <string>

using namespace l1t;

class L1TkMuonStubProducer : public edm::EDProducer {
public:

  typedef TTTrack< Ref_Phase2TrackerDigi_ >  L1TTTrackType;
  typedef std::vector< L1TTTrackType > L1TTTrackCollectionType;

  struct PropState { //something simple, imagine it's hardware emulation
    PropState() :
      pt(-99),  eta(-99), phi(-99),
      sigmaPt(-99),  sigmaEta(-99), sigmaPhi(-99),
      valid(false) {}
    float pt;
    float eta;
    float phi;
    float sigmaPt;
    float sigmaEta;
    float sigmaPhi;
    bool valid;

  };

  enum AlgoType {
    kTP = 1,
    kDynamicWindows = 2
  };

  explicit L1TkMuonStubProducer(const edm::ParameterSet&);
  ~L1TkMuonStubProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  
  void makeMuonsME0Extended(const edm::Handle<EMTFHitCollection>& , TkMuonCollection& ) const;

  // algo for endcap regions using dynamic windows for making the match trk + muStub
  void runOnMuonHitCollection(const edm::Handle<EMTFHitCollection>&,
                          const edm::Handle<L1TTTrackCollectionType>&,
                          TkMuonCollection& tkMuons) const;

  void cleanStubs(const EMTFHitCollection &, EMTFHitCollection &) const;

  // int emtfMatchAlgoVersion_ ;         
  AlgoType emtfMatchAlgoVersion_ ;         

  std::unique_ptr<L1TkMuCorrDynamicWindows> dwcorr_;
  bool requireBX0_;
  int  mu_stub_station_;

  const edm::EDGetTokenT< EMTFTrackCollection >          emtfTCToken; // the track collection, directly from the EMTF and not formatted by GT
  const edm::EDGetTokenT< EMTFHitCollection >            emtfHCToken; // the hit collection, directly from the EMTF which stored the input Hits
  const edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > trackToken;
} ;


L1TkMuonStubProducer::L1TkMuonStubProducer(const edm::ParameterSet& iConfig) :
  emtfTCToken(consumes< EMTFTrackCollection >         (iConfig.getParameter<edm::InputTag>("L1EMTFTrackCollectionInputTag"))),
  emtfHCToken(consumes< EMTFHitCollection >           (iConfig.getParameter<edm::InputTag>("L1EMTFHitCollectionInputTag"))),
  trackToken (consumes< std::vector<TTTrack< Ref_Phase2TrackerDigi_> > > (iConfig.getParameter<edm::InputTag>("L1TrackInputTag")))
{
   
   // configuration of the EMTF algorithm type
   std::string emtfMatchAlgoVersionString = iConfig.getParameter<std::string>("emtfMatchAlgoVersion");
   std::transform(emtfMatchAlgoVersionString.begin(), emtfMatchAlgoVersionString.end(), emtfMatchAlgoVersionString.begin(), ::tolower); // make lowercase

   mu_stub_station_ = iConfig.getParameter<int>("mu_stub_station");
   requireBX0_ = iConfig.getParameter<bool>("require_BX0");

   if (emtfMatchAlgoVersionString == "dynamicwindows")
      emtfMatchAlgoVersion_ = kDynamicWindows;
   else
    throw cms::Exception("TkMuAlgoConfig") << "the ID of the EMTF algo matcher passed is invalid\n";
   

   produces<TkMuonCollection>();
   produces<TkMuonCollection>("ME0Ext");

   // initializations
   if (emtfMatchAlgoVersion_ == kDynamicWindows)
   {
      // FIXME: to merge eventually into an unique file with bith phi and theta boundaries
      std::string fIn_bounds_name = iConfig.getParameter<edm::FileInPath>("emtfcorr_boundaries").fullPath();
      std::string fIn_theta_name  = iConfig.getParameter<edm::FileInPath>("emtfcorr_theta_windows").fullPath();
      std::string fIn_phi_name    = iConfig.getParameter<edm::FileInPath>("emtfcorr_phi_windows").fullPath();
      std::string fIn_S1_theta_name  = iConfig.getParameter<edm::FileInPath>("emtfcorr_S1_theta_windows").fullPath();
      std::string fIn_S1_phi_name    = iConfig.getParameter<edm::FileInPath>("emtfcorr_S1_phi_windows").fullPath();
      auto bounds = L1TkMuCorrDynamicWindows::prepare_corr_bounds(fIn_bounds_name.c_str(), "h_dphi_l");
      TFile* fIn_theta = TFile::Open (fIn_theta_name.c_str());
      TFile* fIn_phi   = TFile::Open (fIn_phi_name.c_str());
      TFile* fIn_S1_theta = TFile::Open (fIn_theta_name.c_str());
      TFile* fIn_S1_phi   = TFile::Open (fIn_phi_name.c_str());
      dwcorr_ = std::unique_ptr<L1TkMuCorrDynamicWindows> (new L1TkMuCorrDynamicWindows(bounds, fIn_theta, fIn_phi, fIn_S1_theta, fIn_S1_phi));

      // files can be closed since the correlator code clones the TF1s
      fIn_theta->Close();
      fIn_phi->Close();

      // FIXME: more initialisation using the parameters passed from the cfg
      dwcorr_->set_safety_factor  (iConfig.getParameter<double>("final_window_factor"));
      dwcorr_->set_sf_initialrelax(iConfig.getParameter<double>("initial_window_factor"));
      
      dwcorr_->set_relaxation_pattern(
        iConfig.getParameter<double>("pt_start_relax"),
        iConfig.getParameter<double>("pt_end_relax")
        );
      
      dwcorr_->set_do_relax_factor(iConfig.getParameter<bool>("do_relax_factors"));

      //
      dwcorr_ -> set_n_trk_par       (iConfig.getParameter<int>("n_trk_par"));
      dwcorr_ -> set_min_trk_p       (iConfig.getParameter<double>("min_trk_p"));
      dwcorr_ -> set_max_trk_aeta    (iConfig.getParameter<double>("max_trk_aeta"));
      dwcorr_ -> set_max_trk_chi2    (iConfig.getParameter<double>("max_trk_chi2"));
      dwcorr_ -> set_min_trk_nstubs  (iConfig.getParameter<int>("min_trk_nstubs"));
      dwcorr_ -> set_do_trk_qual_presel(true);
   }
}

L1TkMuonStubProducer::~L1TkMuonStubProducer() {
}

// ------------ method called to produce the data  ------------
void
L1TkMuonStubProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // the L1Mu objects
  edm::Handle<EMTFTrackCollection> l1emtfTCH;
  edm::Handle<EMTFHitCollection> l1emtfHCH;

  iEvent.getByToken(emtfTCToken, l1emtfTCH);
  iEvent.getByToken(emtfHCToken, l1emtfHCH);

  // the L1Tracks
  edm::Handle<L1TTTrackCollectionType> l1tksH;
  iEvent.getByToken(trackToken, l1tksH);

  TkMuonCollection oc_endcap_tkmuonStub;
  TkMuonCollection oc_me0Extended_tkmuonStub;

  // process each of the MTF collections separately! -- we don't want to filter the muons
  //if (emtfMatchAlgoVersion_ == kDynamicWindows) 
    runOnMuonHitCollection(l1emtfHCH, l1tksH, oc_endcap_tkmuonStub);
  //else
    //throw cms::Exception("TkMuAlgoConfig") << "trying to run an invalid algorithm (this should never happen)\n";

  // now combine all trk muons into a single output collection!
  std::unique_ptr<TkMuonCollection> oc_tkmuon(new TkMuonCollection());
  for (const auto& p : {oc_endcap_tkmuonStub}){
    oc_tkmuon->insert(oc_tkmuon->end(), p.begin(), p.end());
  }

   makeMuonsME0Extended(l1emtfHCH,oc_me0Extended_tkmuonStub);

  // now combine all trk muons into a single output collection!
  std::unique_ptr<TkMuonCollection> oc_me0Ext(new TkMuonCollection());
  for (const auto& p : {oc_me0Extended_tkmuonStub}){
    oc_me0Ext->insert(oc_me0Ext->end(), p.begin(), p.end());
  }
   
  // put the new track+muon objects in the event!
  iEvent.put( std::move(oc_tkmuon));
  // put the new me0  objects in the event!
  iEvent.put( std::move(oc_me0Ext),"ME0Ext");
};


void
L1TkMuonStubProducer::makeMuonsME0Extended(const edm::Handle<EMTFHitCollection>& muonStubH, TkMuonCollection& tkMuons) const

{
  const EMTFHitCollection& l1muStubs = (*muonStubH.product());

  // collection for cleaned stubs
  EMTFHitCollection cleanedStubs;

  // reserve to size of original coolection
  cleanedStubs.reserve(l1muStubs.size());

  // fill collection of cleaned stubs
  cleanStubs(l1muStubs, cleanedStubs);


  for (auto muStub : cleanedStubs) {

    //if(muStub.Subsystem() != EMTFHit::kME0) continue;
    continue;

    
    float stubBend = muStub.Bend();
    stubBend *= 3.63/1000.; // in rad
    float pt = 0.25 * 1.0 / abs(stubBend); // in GeV (scale by ~1/4)

    float eta = muStub.Eta_sim();
    float phi = muStub.Phi_sim() * TMath::Pi()/180.;

    if(abs(eta) < 2.4) continue;

    if(pt < 5.0) continue;
    if (pt > 80) pt=80;  // saturate to remove inf pt due to limited bend angle resolution
    
    // build a L1Candidate
    reco::Candidate::PolarLorentzVector muonLV(pt,eta,phi,0);
    TkMuon muon(muonLV);

    int charge = 1;
    if(stubBend > 0) charge = -1;

    int quality = muStub.Quality();

    
    muon.setCharge(charge);
    muon.setQuality(quality);

    // store 
    //cout << "Storing ME0 muon pt = " << pt << endl;
    tkMuons.push_back( muon);

  }

  return;
}

void
L1TkMuonStubProducer::runOnMuonHitCollection(const edm::Handle<EMTFHitCollection>& muonStubH,
                                     const edm::Handle<L1TTTrackCollectionType>& l1tksH,
                                     TkMuonCollection& tkMuons) const

{
  const EMTFHitCollection& l1muStubs = (*muonStubH.product());

  // collection for cleaned stubs
  EMTFHitCollection cleanedStubs;

  // reserve to size of original coolection
  cleanedStubs.reserve(l1muStubs.size());
  
  // fill collection of cleaned stubs
  cleanStubs(l1muStubs, cleanedStubs);

  const L1TTTrackCollectionType& l1trks = (*l1tksH.product());
  auto corr_muStub_idxs = dwcorr_->find_match_stub(l1muStubs, l1trks, mu_stub_station_, requireBX0_);
  // it's a vector with as many entries as the L1TT vector.
  // >= 0 : the idx in the muStub vector of matched muStubs
  // < 0: no match
  
  // sanity check
  if (corr_muStub_idxs.size() != l1trks.size())
    throw cms::Exception("TkMuAlgoOutput") << "the size of tkmu indices does not match the size of input trk collection\n";
  
  for (uint il1ttrack = 0; il1ttrack < corr_muStub_idxs.size(); ++il1ttrack)
  {
    int muStub_idx = corr_muStub_idxs.at(il1ttrack);
    if (muStub_idx < 0)
      continue;

    const L1TTTrackType& matchTk = l1trks[il1ttrack];
    const auto& p3 = matchTk.momentum();
    const auto& tkv3 = matchTk.POCA();
    const auto& curve = matchTk.rInv();
    float p4e = sqrt(0.105658369*0.105658369 + p3.mag2() );
    math::XYZTLorentzVector l1tkp4(p3.x(), p3.y(), p3.z(), p4e);

    edm::Ref< RegionalMuonCandBxCollection > l1muRef; // FIXME! The reference to the muon is null 
    edm::Ptr< L1TTTrackType > l1tkPtr(l1tksH, il1ttrack);
    float trkisol = -999; // FIXME: now doing as in the TP algo
    TkMuon l1tkmu(l1tkp4, l1muRef, l1tkPtr, trkisol);
    l1tkmu.setTrkzVtx( (float)tkv3.z() );
    
    // Set curvature
    l1tkmu.setTrackCurvature(curve);
    //Set charge
    int charge =  (curve > 0) ? 1 : -1;
    l1tkmu.setCharge(charge);

    tkMuons.push_back(l1tkmu);
  }


  return;
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TkMuonStubProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void 
L1TkMuonStubProducer::cleanStubs(const EMTFHitCollection &  muStubs, EMTFHitCollection & cleanedStubs) const {

    // if empty collection don't do anything
    if(muStubs.size() == 0) return;
    
    // copy the first stub in the new collection
    const EMTFHit & muStub = muStubs[0];
    cleanedStubs.push_back(muStub);

    for(uint ms = 0; ms < muStubs.size(); ms++) {

      const EMTFHit & muStub = muStubs[ms];

      int n_duplicate = 0;

      for(uint i = 0; i <cleanedStubs.size(); i++) {

        const EMTFHit & cStub = cleanedStubs[i];
        int dSubsystem = cStub.Subsystem() - muStub.Subsystem();
        int dStation = cStub.Station() - muStub.Station();
        int dChamber = cStub.Chamber() - muStub.Chamber();
        int dBend = cStub.Bend() - muStub.Bend();
        float aDeltaPhi = abs(cStub.Phi_sim() * TMath::Pi()/180. - muStub.Phi_sim() * TMath::Pi()/180.);
        //cout << "   aDeltaPhi = " << aDeltaPhi << endl;

        // duplicate stubs defined as having same phi
        if(aDeltaPhi < 0.0001 && dSubsystem == 0 && dStation == 0 && dChamber <= 1 && dBend == 0) {

          n_duplicate++;

        } // end if

      } // end for cleaned

      // did not find any duplicate stubs and the stub is not Neighbor
      if(n_duplicate == 0 && muStub.Neighbor() == 0) {
        cleanedStubs.push_back(muStub);
      }

    } // end for muStubs

    /*
    printf("----------------- Cleaned stubs ---------------------------------------------------- \n");
    for(uint i = 0; i <cleanedStubs.size(); i++) {
        const EMTFHit & cStub = cleanedStubs[i];
        printStub(cStub);
    }
    printf("------------------------------------------------------------------------------- \n");
    */

}


//define this as a plug-in
DEFINE_FWK_MODULE(L1TkMuonStubProducer);



