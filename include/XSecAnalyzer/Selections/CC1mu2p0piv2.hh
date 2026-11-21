#pragma once

// ROOT includes
#include "TVector3.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/SelectionBase.hh"
#include "XSecAnalyzer/KICalculator.hh"

class CC1mu2p0piv2 : public SelectionBase {

public:

  CC1mu2p0piv2( const std::string& sel_name = "CC1mu2p0piv2" );

  virtual std::string categorize_event( AnalysisEvent& ev ) override;
  virtual bool is_selected( AnalysisEvent& ev ) override;

protected:

  // Pandora results of the most recent call to is_selected(), for use by
  // derived selections that combine Pandora with other reconstructions
  struct PandoraResult {
    bool selected = false;
    bool nu_mu_cc = false;
    bool reco_vertex_in_FV = false;
    bool has_vertex = false; // a reconstructed (raw and SCE) vertex exists
    bool has_muon_candidate = false;
    bool muon_contained = false;
    bool muon_passed_mom_cuts = false;
    bool muon_quality_if_contained = false;
    bool muon_ng2_mip = false;
    bool no_reco_showers = false;
    bool two_p_candidates = false;
    bool no_extra_pfps = false;
    bool protons_contained = false;
    bool protons_near_vertex = false;
    bool protons_no_descendents = false;
    double lead_p_ng2_hipfrac = -9999.;
    double sub_p_ng2_hipfrac = -9999.;
    double p_mu_range = -9999.;
    double p_mu_mcs = -9999.;
    double cos_gamma_pp = -9999.;
    // reco 3-momenta of the muon and proton candidates (valid if has_p3)
    bool has_p3 = false;
    bool has_lead_p = false, has_sub_p = false;
    TVector3 p3mu, p3_lead_p, p3_sub_p;
    // Sideband topologies (see is_selected): the selection's muon plus
    // exactly three good protons (sb_3p), one or two good protons and one
    // pion-like track at the vertex (sb_npi), or exactly one good proton
    // (sb_1p), with no showers or other direct daughters
    bool sb_3p = false, sb_npi = false, sb_1p = false;
  };
  PandoraResult pd_;

  // Computes the reco observables from the muon, leading and subleading
  // proton 3-momenta and stores them in the output tree (branches named
  // prefix + observable). If ok is false, all observables are set to BOGUS,
  // except, if lead_only_ok is true, those of the muon and the leading
  // proton alone (single-particle and leading-proton imbalance variables).
  void fill_observables( TreeMap::Wrapper& out, bool ok, const TVector3& p3mu,
    const TVector3& p3_lead_p, const TVector3& p3_sub_p,
    const std::string& prefix = "reco_", bool lead_only_ok = false ) const;

  KICalculator::CalcType ki_calc_type_;
};
