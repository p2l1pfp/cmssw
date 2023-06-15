#ifndef DataFormats_L1THGCal_HGCalCluster_HW_h
#define DataFormats_L1THGCal_HGCalCluster_HW_h

#include "DataFormats/L1TParticleFlow/interface/bit_encoding.h"

#include <array>
// #include <cstdint>

namespace l1thgcfirmware {

  typedef ap_uint<14> e_t; // E, E_EM
  typedef ap_uint<4> gctbits_t;
  typedef ap_uint<8> eFraction_t;
  typedef ap_uint<6> nLayer_t;

  typedef ap_uint<10> eta_t;
  typedef ap_int<9> phi_t;
  typedef ap_uint<12> z_t;
  typedef ap_uint<1> spare_w2_t;
  typedef ap_uint<10> nTC_t;
  typedef ap_uint<7> qualFlags_t;

  typedef ap_uint<7> sigma_E_t;
  typedef ap_uint<13> spare_w3_t;
  typedef ap_uint<7> sigma_z_t;
  typedef ap_uint<7> sigma_phi_t;
  typedef ap_uint<5> sigma_eta_t;
  typedef ap_uint<7> sigma_roz_t;

  static constexpr int wordLength = 64;
  static constexpr int nWordsPerCluster = 4;
  typedef uint64_t ClusterWord;
  typedef std::array<ClusterWord, nWordsPerCluster> ClusterWords;
//   namespace Scales {
    
//   }

  struct HGCalCluster_HW {

    // First word
    e_t e;
    e_t e_em;
    gctbits_t gctBits;
    eFraction_t fractionInCE_E;
    eFraction_t fractionInCoreCE_E;
    eFraction_t fractionInEarlyCE_E;
    nLayer_t firstLayer;

    // Second word
    eta_t w_eta;
    phi_t w_phi;
    z_t w_z;
    spare_w2_t spare_w2 = 0;
    nTC_t nTC;
    qualFlags_t qualFlags;

    // Third word
    sigma_E_t sigma_E;
    nLayer_t lastLayer;
    nLayer_t showerLength;
    spare_w3_t spare_w3 = 0;
    sigma_z_t sigma_z;
    sigma_phi_t sigma_phi;
    nLayer_t coreShowerLength;
    sigma_eta_t sigma_eta;
    sigma_roz_t sigma_roz;

    inline bool operator==(const HGCalCluster_HW &other) const {
      return e == other.e && e_em == other.e_em && w_eta == other.w_eta && w_phi == other.w_phi && w_z == other.w_z  && sigma_roz == other.sigma_roz;
    }

    inline bool operator!=(const HGCalCluster_HW &other) const {
      return !(*this == other);
    }

    inline void clear() {
      e = 0;
      e_em = 0;
      gctBits = 0;
      fractionInCE_E = 0;
      fractionInCoreCE_E = 0;
      fractionInEarlyCE_E = 0;
      firstLayer = 0;
      w_eta = 0;
      w_phi = 0;
      w_z = 0;
      nTC = 0;
      qualFlags = 0;
      sigma_E = 0;
      lastLayer = 0;
      showerLength = 0;
      spare_w3 = 0;
      sigma_z = 0;
      sigma_phi = 0;
      coreShowerLength = 0;
      sigma_eta = 0;
      sigma_roz = 0;
    }

    // float floatE() const { return Scales::floatE(e); }
    // float floatE_EM() const { return Scales::floatE(e_em); }
    // float floatFraction() const { return Scales::floatEta(hwEta); }

    static const int BITWIDTH_FIRSTWORD = e_t::width + e_t::width + gctbits_t::width + eFraction_t::width + eFraction_t::width + eFraction_t::width + nLayer_t::width;
    inline ap_uint<BITWIDTH_FIRSTWORD> pack_firstWord() const {
      ap_uint<BITWIDTH_FIRSTWORD> ret;
      unsigned int start = 0;
      pack_into_bits(ret, start, e);
      pack_into_bits(ret, start, e_em);
      pack_into_bits(ret, start, gctBits);
      pack_into_bits(ret, start, fractionInCE_E);
      pack_into_bits(ret, start, fractionInCoreCE_E);
      pack_into_bits(ret, start, fractionInEarlyCE_E);
      pack_into_bits(ret, start, firstLayer);
      return ret;
    }

    static const int BITWIDTH_SECONDWORD = eta_t::width + phi_t::width + z_t::width + spare_w2_t::width + nTC_t::width + qualFlags_t::width;
    inline ap_uint<BITWIDTH_SECONDWORD> pack_secondWord() const {
      ap_uint<BITWIDTH_SECONDWORD> ret;
      unsigned int start = 0;
      pack_into_bits(ret, start, w_eta);
      pack_into_bits(ret, start, w_phi);
      pack_into_bits(ret, start, w_z);
      pack_into_bits(ret, start, spare_w2);
      pack_into_bits(ret, start, nTC);
      pack_into_bits(ret, start, qualFlags);
      return ret;
    }

    static const int BITWIDTH_THIRDWORD = sigma_E_t::width + nLayer_t::width + nLayer_t::width + spare_w3_t::width + sigma_z_t::width + sigma_phi_t::width + nLayer_t::width + sigma_eta_t::width + sigma_roz_t::width;
    inline ap_uint<BITWIDTH_THIRDWORD> pack_thirdWord() const {
      ap_uint<BITWIDTH_THIRDWORD> ret;
      unsigned int start = 0;
      pack_into_bits(ret, start, sigma_E);
      pack_into_bits(ret, start, lastLayer);
      pack_into_bits(ret, start, showerLength);
      pack_into_bits(ret, start, spare_w3);
      pack_into_bits(ret, start, sigma_z);
      pack_into_bits(ret, start, sigma_phi);
      pack_into_bits(ret, start, coreShowerLength);
      pack_into_bits(ret, start, sigma_eta);
      pack_into_bits(ret, start, sigma_roz);
      return ret;
    }

    inline ClusterWords pack() const {
      ClusterWords packed;
      ap_uint<BITWIDTH_FIRSTWORD> firstWord = this->pack_firstWord();
      ap_uint<BITWIDTH_SECONDWORD> secondWord = this->pack_secondWord();
      ap_uint<BITWIDTH_THIRDWORD> thirdWord = this->pack_thirdWord();
      packed[0] = firstWord;
      packed[1] = secondWord;
      packed[2] = thirdWord;
      packed[3] = 0;
      return packed;
    }

    inline static void unpack_firstWord(const ap_uint<BITWIDTH_FIRSTWORD> &src, HGCalCluster_HW& cluster) {
      cluster.initFromBits_firstWord(src);
      return;
    }

    inline void initFromBits_firstWord(const ap_uint<BITWIDTH_FIRSTWORD> &src ) {
      unsigned int start = 0;
      unpack_from_bits(src, start, e);
      unpack_from_bits(src, start, e_em);
      unpack_from_bits(src, start, gctBits);
      unpack_from_bits(src, start, fractionInCE_E);
      unpack_from_bits(src, start, fractionInCoreCE_E);
      unpack_from_bits(src, start, fractionInEarlyCE_E);
      unpack_from_bits(src, start, firstLayer);
    }

    inline static void unpack_secondWord(const ap_uint<BITWIDTH_SECONDWORD> &src, HGCalCluster_HW& cluster) {
      cluster.initFromBits_secondWord(src);
      return;
    }

    inline void initFromBits_secondWord(const ap_uint<BITWIDTH_SECONDWORD> &src ) {
      unsigned int start = 0;
      unpack_from_bits(src, start, w_eta);
      unpack_from_bits(src, start, w_phi);
      unpack_from_bits(src, start, w_z);
      unpack_from_bits(src, start, spare_w2);
      unpack_from_bits(src, start, nTC);
      unpack_from_bits(src, start, qualFlags);

    }

    inline static void unpack_thirdWord(const ap_uint<BITWIDTH_THIRDWORD> &src, HGCalCluster_HW& cluster) {
      cluster.initFromBits_thirdWord(src);
      return;
    }

    inline void initFromBits_thirdWord(const ap_uint<BITWIDTH_THIRDWORD> &src ) {
      unsigned int start = 0;
      unpack_from_bits(src, start, sigma_E);
      unpack_from_bits(src, start, lastLayer);
      unpack_from_bits(src, start, showerLength);
      unpack_from_bits(src, start, spare_w3);
      unpack_from_bits(src, start, sigma_z);
      unpack_from_bits(src, start, sigma_phi);
      unpack_from_bits(src, start, coreShowerLength);
      unpack_from_bits(src, start, sigma_eta);
      unpack_from_bits(src, start, sigma_roz);
    }


    inline static HGCalCluster_HW unpack(const ClusterWords &src) {
        HGCalCluster_HW cluster;
        ap_uint<BITWIDTH_FIRSTWORD> firstWord = src[0];
        unpack_firstWord( firstWord, cluster );
        ap_uint<BITWIDTH_SECONDWORD> secondWord = src[1];
        unpack_secondWord( secondWord, cluster );
        ap_uint<BITWIDTH_THIRDWORD> thirdWord = src[2];
        unpack_thirdWord( thirdWord, cluster );
        return cluster;
    }
  };

  inline void clear(HGCalCluster_HW &c) { c.clear(); }

}  // namespace l1thgcfirmware

#endif