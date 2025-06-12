#pragma once

#include "Selection.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
/// \brief    Tournament Selection approach
////////////////////////////////////////////////////////////////////////
class CTournamentSelection : public CSelection
{
public:

    virtual void ResetSelection(const CPopulation& population) override;

    virtual  unsigned int Select(const CPopulation& population, default_random_engine& randomEngineu, unsigned int skip=UINT_MAX) override;

private:
    uniform_int_distribution<unsigned int>  mRandomGenerator;

};

