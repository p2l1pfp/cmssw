// -*- C++ -*-
//
// Package:    PrintL1TkObjects
// Class:      PrintL1TkObjects
// 
/**\class PrintL1TkObjects PrintL1TkObjects.cc SLHCUpgradeSimulations/PrintL1TkObjects/src/PrintL1TkObjects.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Emmanuelle Perez,40 1-A28,+41227671915,
//         Created:  Thu Nov 14 11:22:13 CET 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"


// Gen-level stuff:
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/L1TCorrelator/interface/TkPrimaryVertex.h"
#include "DataFormats/L1TCorrelator/interface/TkEtMiss.h"
#include "DataFormats/L1TCorrelator/interface/TkEtMissFwd.h"
#include "DataFormats/L1TCorrelator/interface/TkEm.h"
#include "DataFormats/L1TCorrelator/interface/TkEmFwd.h"
#include "DataFormats/L1TCorrelator/interface/TkElectron.h"
#include "DataFormats/L1TCorrelator/interface/TkElectronFwd.h"
#include "DataFormats/L1TCorrelator/interface/TkJet.h"
#include "DataFormats/L1TCorrelator/interface/TkJetFwd.h"
#include "DataFormats/L1TCorrelator/interface/TkHTMiss.h"
#include "DataFormats/L1TCorrelator/interface/TkHTMissFwd.h"
#include "DataFormats/L1TCorrelator/interface/TkMuon.h"
#include "DataFormats/L1TCorrelator/interface/TkMuonFwd.h"
#include "DataFormats/L1TCorrelator/interface/TkTau.h"
#include "DataFormats/L1TCorrelator/interface/TkTauFwd.h"

#include "TFile.h"
#include "TH1F.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace l1t;



//
// class declaration
//

class PrintL1TkObjects : public edm::EDAnalyzer {
public:

  //  typedef TTTrack< Ref_PixelDigi_ >  L1TkTrackType;
  //  typedef std::vector< L1TkTrackType >  L1TkTrackCollectionType;
  typedef TTTrack< Ref_Phase2TrackerDigi_ >  L1TTTrackType;
  typedef std::vector< L1TTTrackType > L1TTTrackCollectionType;
  
  explicit PrintL1TkObjects(const edm::ParameterSet&);
  ~PrintL1TkObjects();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  //virtual void endRun(edm::Run const&, edm::EventSetup const&);
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  // ----------member data ---------------------------
  
  
  // HEPMC token
  const edm::EDGetTokenT<edm::HepMCProduct> hepMCToken;

  // to test the TkPrimaryVertex :
  const edm::EDGetTokenT< TkPrimaryVertexCollection > tkPvToken;

  // for L1TrackEtmiss:
  const edm::EDGetTokenT< TkEtMissCollection >  tkEtMissToken;

  // for TkEms
  const edm::EDGetTokenT< TkJetCollection >  tkJetToken;

  // for L1TkHTMiss
  const edm::EDGetTokenT< TkHTMissCollection > tkHTMissToken;

  // for TkJets
  const edm::EDGetTokenT< TkEmCollection > tkPhotonToken;

  // for L1TkElectron
  const edm::EDGetTokenT< TkElectronCollection > tkElecToken;

  // for L1TkMuon
  const edm::EDGetTokenT< TkMuonCollection > tkMuonToken;

  // for L1TkTau
  const edm::EDGetTokenT< TkTauCollection > tkTauToken;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PrintL1TkObjects::PrintL1TkObjects(const edm::ParameterSet& iConfig) :
   hepMCToken(consumes<edm::HepMCProduct>(edm::InputTag("generator",""))),
   tkPvToken(consumes< TkPrimaryVertexCollection >(iConfig.getParameter<edm::InputTag>("L1TkVtxInputTag"))),
   tkEtMissToken(consumes< TkEtMissCollection> (iConfig.getParameter<edm::InputTag>("TkEtMissInputTag"))),
   tkJetToken(consumes<TkJetCollection>(iConfig.getParameter<edm::InputTag>("TkJetsInputTag"))),
   tkHTMissToken(consumes<TkHTMissCollection>(iConfig.getParameter<edm::InputTag>("TkHTMInputTag"))),
   tkPhotonToken(consumes<TkEmCollection>(iConfig.getParameter<edm::InputTag>("TkPhotonsInputTag"))),
   tkElecToken(consumes<TkElectronCollection>(iConfig.getParameter<edm::InputTag>("TkElectronsInputTag"))),
   tkMuonToken(consumes<TkMuonCollection>(iConfig.getParameter<edm::InputTag>("TkMuonsInputTag"))),
   tkTauToken(consumes<TkTauCollection>(iConfig.getParameter<edm::InputTag>("TkTausInputTag")))
{
   //now do what ever initialization is needed

}


PrintL1TkObjects::~PrintL1TkObjects()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PrintL1TkObjects::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   std::cout << " ----  a new event ----- " << std::endl;

	// First, retrieve the generated primary vertex


  edm::Handle<edm::HepMCProduct> HepMCEvt;
  iEvent.getByToken(hepMCToken,HepMCEvt);

     float zvtx_gen = -999;

 if (HepMCEvt.isValid() ) {

     const HepMC::GenEvent* MCEvt = HepMCEvt->GetEvent();
     const double mm=0.1;

     for ( HepMC::GenEvent::vertex_const_iterator ivertex = MCEvt->vertices_begin(); ivertex != MCEvt->vertices_end(); ++ivertex )
         {
             bool hasParentVertex = false;
 
             // Loop over the parents looking to see if they are coming from a production vertex
             for (
                 HepMC::GenVertex::particle_iterator iparent = (*ivertex)->particles_begin(HepMC::parents);
                 iparent != (*ivertex)->particles_end(HepMC::parents);
                 ++iparent
             )
                 if ( (*iparent)->production_vertex() )
                 {
                     hasParentVertex = true;
                     break;
                }
 
             // Reject those vertices with parent vertices
             if (hasParentVertex) continue;
 
             // Get the position of the vertex
             HepMC::FourVector pos = (*ivertex)->position();
	     zvtx_gen = pos.z()*mm; 

	     break;  // there should be one single primary vertex

          }  // end loop over gen vertices

     std::cout << " Generated zvertex : " << zvtx_gen << std::endl;
  }


	// ----------------------------------------------------------------------
	//
        // retrieve the L1 vertices
        
 edm::Handle<TkPrimaryVertexCollection> L1VertexHandle;
 iEvent.getByToken(tkPvToken,L1VertexHandle);
 std::vector<TkPrimaryVertex>::const_iterator vtxIter;
 
 if ( L1VertexHandle.isValid() ) {
     std::cout << " -----  TkPrimaryVertex objects   ----- " << std::endl;
     vtxIter = L1VertexHandle->begin();     // only one algorithm is run in the TkPrimaryVertexProducer
					    // (in contrast to earlier, under-dev, versions of the code)
        float z = vtxIter -> zvertex();
        float sum = vtxIter -> sum();
        std::cout << " a vertex with  zvtx " << z << " (cm) and SumPT " << sum << " (GeV) " << std::endl;
 }

	//
        // ----------------------------------------------------------------------
	// retrieve the EtMiss objects
	//

 edm::Handle<TkEtMissCollection> L1TkEtMissHandle;
 iEvent.getByToken(tkEtMissToken, L1TkEtMissHandle);
 std::vector<TkEtMiss>::const_iterator etmIter;

 if (L1TkEtMissHandle.isValid() ) {
    std::cout << " -----  L1TkEtMiss objects  -----  " << std::endl; 
    etmIter = L1TkEtMissHandle -> begin();	// idem: only one TrkMET now.
	float etmis = etmIter -> et();
	const edm::Ref< TkPrimaryVertexCollection > vtxRef = etmIter -> vtxRef();
	float zvtx = vtxRef -> zvertex();
        float etMissPU = etmIter -> etMissPU();
	std::cout << " ETmiss = " << etmis << " for zvtx = " << zvtx << " and ETmiss from PU = " << etMissPU << std::endl;
 }


        //
        // ----------------------------------------------------------------------
        // retrieve the TkJet objects
        //
        
 edm::Handle<TkJetCollection> L1TkJetsHandle;
 iEvent.getByToken(tkJetToken, L1TkJetsHandle);
 std::vector<TkJet>::const_iterator jetIter ;

 if ( L1TkJetsHandle.isValid() ) {
    std::cout << " -----   TkJet  objects -----  " << std::endl;
    for (jetIter = L1TkJetsHandle -> begin(); jetIter != L1TkJetsHandle->end(); ++jetIter) {
        float et = jetIter -> pt();
        float phi = jetIter -> phi();
        float eta = jetIter -> eta();
	float jetvtx = jetIter -> jetVtx();
        const edm::Ref< JetBxCollection > Jetref = jetIter -> jetRef();
        float et_L1Jet = Jetref -> et();
	const std::vector< edm::Ptr< L1TTTrackType > > trkPtrs = jetIter -> trkPtrs() ;

        std::cout << " a Jet candidate ET eta phi zvertex " << et << " " << eta << " " << phi << " " << jetvtx  << std::endl;
        std::cout << "                Calo  ET, typ " << et_L1Jet << std::endl;
	std::cout << "      Tracks associated with the jet : " << std::endl;
        for (int it=0; it < (int)trkPtrs.size(); it++) {
            edm::Ptr< L1TTTrackType > aTrack = trkPtrs.at( it );
            std::cout << "             a track PT " << aTrack -> momentum().perp() << " eta " << 
		   aTrack -> momentum().eta() << " phi " << aTrack -> momentum().phi() << std::endl;
        }

    }
 }

        //
        // ----------------------------------------------------------------------
        // retrieve HT and HTM
	//

 edm::Handle<TkHTMissCollection> L1TkHTMHandle;
 iEvent.getByToken(tkHTMissToken, L1TkHTMHandle);

 if ( L1TkHTMHandle.isValid() ) {
	std::cout << " -----  TkHTMiss: size (should be 1) = " << L1TkHTMHandle -> size() << std::endl;
	std::vector<TkHTMiss>::const_iterator HTMIter = L1TkHTMHandle -> begin();
	float HTT = HTMIter -> etTotal();
	float HTM = HTMIter -> EtMiss();
	//float HTM_the_same = HTMIter -> et();

	// phi of the HTM vector :
	float phi = HTMIter -> phi();
	std::cout << " HTT = " << HTT << " HTM = " << HTM << " " << "phi(HTM) = " << phi << std::endl;
	 
	// access the L1TkJets used to build HT and HTM :
	const edm::RefProd< TkJetCollection > jetCollRef = HTMIter -> jetCollectionRef();
 	std::vector<TkJet>::const_iterator jet = jetCollRef -> begin();
	std::cout << " ET of the first L1TkJet = " << jet -> et() << std::endl;
 }


        //
        // ----------------------------------------------------------------------
        // retrieve the TkEm objects
	//

 edm::Handle<TkEmCollection> L1TkPhotonsHandle;
 iEvent.getByToken(tkPhotonToken, L1TkPhotonsHandle);
 std::vector<TkEm>::const_iterator phoIter ;

 if ( L1TkPhotonsHandle.isValid() ) {
    std::cout << " -----   TkEm  objects -----  " << std::endl;
    for (phoIter = L1TkPhotonsHandle -> begin(); phoIter != L1TkPhotonsHandle->end(); ++phoIter) {
	float et = phoIter -> pt();
	float phi = phoIter -> phi();
	float eta = phoIter -> eta();
        float trkisol = phoIter -> trkIsol() ;
	const edm::Ref< EGammaBxCollection > EGref = phoIter -> EGRef();
	float et_L1Calo = EGref -> et();
	float eta_calo = EGref -> eta();
	float phi_calo = EGref -> phi();

	std::cout << " a photon candidate ET eta phi trkisol " << et << " " << eta << " " << phi << " " << trkisol << std::endl;
	std::cout << "                Calo  ET eta phi " << et_L1Calo << " " << eta_calo << " " << phi_calo << std::endl; 
    }
 }


        //
        // ----------------------------------------------------------------------
        // retrieve the TkElectron objects
        //

 edm::Handle<TkElectronCollection> L1TkElectronsHandle;
 iEvent.getByToken(tkElecToken, L1TkElectronsHandle);
 std::vector<TkElectron>::const_iterator eleIter ;

 if ( L1TkElectronsHandle.isValid() ) {
    std::cout << " -----   TkElectron  objects -----  " << std::endl;
    for (eleIter = L1TkElectronsHandle -> begin(); eleIter != L1TkElectronsHandle->end(); ++eleIter) {
        float et = eleIter -> pt();
        float phi = eleIter -> phi();
        float eta = eleIter -> eta();
    	float trkisol = eleIter -> trkIsol() ;
	float ztr = eleIter -> trkzVtx() ;
        std::cout << "an electron candidate ET eta phi trkisol ztr " << et << " " << eta << " " << phi <<  " " << trkisol << " " << ztr << std::endl;

        const edm::Ref< EGammaBxCollection > EGref = eleIter -> EGRef();
        if ( EGref.isNonnull() ) {
           float et_L1Calo = EGref -> et();
           float eta_calo = EGref -> eta();
           float phi_calo = EGref -> phi();
           std::cout << "                Calo  ET eta phi " << et_L1Calo << " " << eta_calo << " " << phi_calo << std::endl;
	}
	else {
	    std::cout << " .... edm::Ref to EGamma is unvalid !! ?  " << std::endl;
	}

        const edm::Ptr< L1TTTrackType > TrkRef = eleIter -> trkPtr();
	if ( TrkRef.isNonnull() ) {
            float pt_track = TrkRef -> momentum().perp();
            float phi_track = TrkRef -> momentum().phi();
            float eta_track = TrkRef -> momentum().eta();
            float ztrack = TrkRef -> POCA().z();
            std::cout << "                Track PT eta phi ztr " << pt_track << " " << eta_track << " " << phi_track << " " << ztrack << std::endl;
	}
	else {
	    std::cout << " ... edm::Ptr to L1Tracks is unvalid (e.g. electron was matched to stubs) " << std::endl;
	}
    }
 }

        //  
        // ----------------------------------------------------------------------
        // retrieve the L1TkMuons
        //

 edm::Handle<TkMuonCollection> L1TkMuonsHandle;
 iEvent.getByToken(tkMuonToken, L1TkMuonsHandle);
 std::vector<TkMuon>::const_iterator muIter;

 if ( L1TkMuonsHandle.isValid() ) {
    std::cout << " -----   L1TkMuonPaticle objects ---- " << std::endl;
    for (muIter = L1TkMuonsHandle -> begin(); muIter != L1TkMuonsHandle->end(); ++muIter) {
        float pt = muIter -> pt();
        float eta = muIter -> eta();
        float phi = muIter -> phi();
        float zvtx = muIter -> trkzVtx();
        std::cout << " a muon candidate pt eta phi " << pt << " " << eta << " " << phi << " zvertex = " << zvtx << std::endl;
    }
 }

        //  
        // ----------------------------------------------------------------------
        // retrieve the L1TkMuons
        //

 edm::Handle<TkTauCollection> L1TkTausHandle;
 iEvent.getByToken(tkTauToken,L1TkTausHandle);
 std::vector<TkTau>::const_iterator tauIter;

 if ( L1TkTausHandle.isValid() ) {
    std::cout << " -----   L1TkTauPaticle objects ---- " << std::endl;
    for (tauIter = L1TkTausHandle -> begin(); tauIter != L1TkTausHandle->end(); ++tauIter) {
        float pt = tauIter -> pt();
        float eta = tauIter -> eta();
        float phi = tauIter -> phi();
        float zvtx = tauIter -> trkzVtx();
        std::cout << " a tai candidate pt eta phi " << pt << " " << eta << " " << phi << " zvertex = " << zvtx << std::endl;

    }
 }


}


// ------------ method called once each job just before starting event loop  ------------
void 
PrintL1TkObjects::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PrintL1TkObjects::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
PrintL1TkObjects::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
PrintL1TkObjects::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
PrintL1TkObjects::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
PrintL1TkObjects::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PrintL1TkObjects::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PrintL1TkObjects);
