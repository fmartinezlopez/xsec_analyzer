#pragma once

// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/BinSchemeBase.hh"

class ProtonMultiplicityBinScheme : public BinSchemeBase {

  public:

    ProtonMultiplicityBinScheme();
    virtual void DefineBlocks() override;
};
