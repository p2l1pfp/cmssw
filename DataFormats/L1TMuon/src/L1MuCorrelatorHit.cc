#include "DataFormats/L1TMuon/interface/L1MuCorrelatorHit.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
L1MuCorrelatorHit::L1MuCorrelatorHit() :
  etaRegion_(0),phiRegion_(0),depthRegion_(0),phi_(0), phiB_(0),id_(0), quality_(-1), bxNum_(17),
  eta_(0),etaQuality_(-1),type_(0) {}

L1MuCorrelatorHit::L1MuCorrelatorHit(int etaRegion,int phiRegion,int depthRegion,uint tfLayer,int phi,int phiB,int id,int bx,int quality,int eta,int etaQuality,int type):
  etaRegion_(etaRegion),
  phiRegion_(phiRegion),
  depthRegion_(depthRegion),
  tfLayer_(tfLayer),
  phi_(phi),
  phiB_(phiB),
  id_(id),
  quality_(quality),
  bxNum_(bx),
  eta_(eta),
  etaQuality_(etaQuality),
  type_(type)
{

}

L1MuCorrelatorHit::~L1MuCorrelatorHit() {}


bool L1MuCorrelatorHit::operator==(const L1MuCorrelatorHit& id) const {

  if ( etaRegion_             != id.etaRegion_ )               return false;
  if ( phiRegion_             != id.phiRegion_ )               return false;
  if ( depthRegion_           != id.depthRegion_)               return false;
  if ( id_                    != id.id_ )                      return false;
  if ( phi_                   != id.phi_ )                     return false;
  if ( phiB_                  != id.phiB_)                     return false;
  if ( quality_               != id.quality_ )                 return false;
  if ( bxNum_                 != id.bxNum_ )                   return false;
  if ( eta_                   != id.eta_ )                     return false;
  if ( etaQuality_            != id.etaQuality_ )              return false;
  if ( type_                  != id.type_ )                    return false;
  return true;
}

//
// output stream operator for phi track segments
//
ostream& operator<<(ostream& s, const L1MuCorrelatorHit& id) {

  s.setf(ios::right,ios::adjustfield);
  s << "BX: "              << setw(5) << id.bxNum_  << " "
  << "etaregion:"          << setw(5) << id.etaRegion_  << " "
  << "phiRegion: "         << setw(5) << id.phiRegion_  << " "
  << "depthRegion: "       << setw(5) << id.depthRegion_  << " "
  << "stub ID: "           << setw(5) << id.id_  << " "
  << "phi: "               << setw(5) << id.phi_  << " "
  << "phiB: "              << setw(4) << id.phiB_ << " "
  << "quality: "           << setw(4) << id.quality_ << " "
  << "eta:"                << setw(4) << id.eta_ << " "
  << "qeta1:"              << setw(4) << id.etaQuality_ << " "
  << "type:"               << setw(4) <<id.type_;
  return s;

}
