//-------------------------------------------------
//
/**  \class L1MuCorrelatorHit
 *
 *  Class that creates a super-primitive for all chambers
 *
 *
 *   M.Bachtis (UCLA)
 */
//
//--------------------------------------------------
#ifndef L1MU_CORRELATORHIT_H
#define L1MU_CORRELATORHIT_H
//---------------
// C++ Headers --
//---------------

#include <iosfwd>
#include <vector>

//----------------------
// Base Class Headers --
//----------------------

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

#include "DataFormats/L1TMuon/interface/BMTF/L1MuBMTrackSegLoc.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/Common/interface/Ref.h"

//              ---------------------
//              -- Class Interface --
//              ---------------------

class L1MuCorrelatorHit;

typedef std::vector<L1MuCorrelatorHit> L1MuCorrelatorHitCollection;
typedef edm::Ref<L1MuCorrelatorHitCollection > L1MuCorrelatorHitRef;
typedef std::vector<edm::Ref<L1MuCorrelatorHitCollection> > L1MuCorrelatorHitRefVector;

class L1MuCorrelatorHit {

  public:

    /// default constructor
    L1MuCorrelatorHit();

    /// constructor
    L1MuCorrelatorHit(int etaRegion,int phiRegion,int depthRegion,uint tfLayer,int phi,int phiB,int id,int bx,int quality,int eta=0,int etaQuality=-1,int type=0);
    ~L1MuCorrelatorHit();
    /// return wheel
    inline int etaRegion() const { return etaRegion_; }
    /// return sector
    inline int phiRegion() const { return phiRegion_; }
    /// return station
    inline int depthRegion() const { return depthRegion_; }
    /// return track finder layer
    inline uint tfLayer() const { return tfLayer_; }
    /// return phi
    inline int phi() const { return phi_; }
    /// return phib
    inline int phiB() const { return phiB_; }
    /// return quality code
    inline int quality() const { return quality_; }
    /// return tag (second TS tag)
    inline int id() const { return id_; }
    /// return bunch crossing
    inline int bxNum() const { return bxNum_; }

    /// return first eta
    inline int eta() const { return eta_; }
    /// return first eta quality
    inline int etaQuality() const { return etaQuality_; }
    //return type
    inline int type() const { return type_; }


    inline bool isDT()  const { return (type_==0); }
    inline bool isRPCBarrel() const { return (type_==1); }
    inline bool isCSC() const { return (type_==2); }
    inline bool isRPCEndcap() const { return (type_==3); }


    void setEta(int eta, int etaQ) {
      eta_ = eta;
      etaQuality_ = etaQ;
    }

    /// assignment operator
    L1MuCorrelatorHit& operator=(const L1MuCorrelatorHit&);
    /// equal operator
    bool operator==(const L1MuCorrelatorHit&) const;
    /// unequal operator
    bool operator!=(const L1MuCorrelatorHit&) const;

    /// overload output stream operator for phi track segments
    friend std::ostream& operator<<(std::ostream&, const L1MuCorrelatorHit&);

  private:

    int               etaRegion_;           //In the barrel this is wheel. In the endcap it is 6-ring       
    int               phiRegion_;           //In the barrel it is sector. In the endcap it is chamber
    int               depthRegion_;         //Station
    uint               tfLayer_;             //TF Layer
    int               phi_;                 // global position angle in units of 30 degrees/2048
    int               phiB_;                // bending angle  only in barrel for now
    int               id_;                  // stub id in case of more stubs per chamber
    int               quality_;             // 
    int               bxNum_;               // bunch crossing identifier
    int               eta_;                 // eta coordinate - in units of 3.0/512.
    int               etaQuality_;          // quality of the eta information
    int               type_;                //Type: 0 TwinMux or DT, 1 RPC Barrel, 2 CSC, 3 RPC endcap
};

#endif
