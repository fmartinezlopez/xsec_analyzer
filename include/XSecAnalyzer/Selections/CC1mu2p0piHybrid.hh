#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/CC1mu2p0piv2.hh"

// nu_mu CC0pi2p selection combining three reconstructions: an event is
// selected when at least two of Pandora (a looser version of CC1mu2p0piv2),
// Wire-Cell and LANTERN find the 1mu2p topology, plus the CC1mu2p0piv2
// events with nearly collinear protons, which the other two cannot confirm.
// The kinematics come from Pandora whenever it voted, otherwise from
// LANTERN (protons) with the aligned Pandora or Wire-Cell muon.
// Signal definition and event categories are those of CC1mu2p0piv2.
// See studies/stage7_combined/README.md for the tuning and performance.
class CC1mu2p0piHybrid : public CC1mu2p0piv2 {

public:

  CC1mu2p0piHybrid();

  virtual bool is_selected( AnalysisEvent& ev ) override final;

  // Events in one of the sidebands are also kept by ProcessNTuples
  virtual bool last_keep_extra() const override { return in_sideband_; }

private:

  bool in_sideband_ = false;
};
