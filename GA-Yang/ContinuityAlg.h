#pragma once
#include "PostProcessingAlg.h"

////////////////////////////////////////////////////////////////////////
/// \brief    Continuity algorithm
////////////////////////////////////////////////////////////////////////
class CContinuityAlg :
    public CPostProcessingAlg
{
public:
    CContinuityAlg() : CPostProcessingAlg() {}
    virtual ~CContinuityAlg() {}

    virtual void Update(const std::vector<CChromosome>& iElitism) override;
    virtual bool DerivedExecute(const CGenotypeData& ioGenotypeData, int iOrder) override;
    virtual void Print(const CGenotypeData& ioGenotypeData, const std::string& iBaseFileName, int iOrder) const override;
    virtual void UpdateGenotypeData(CGenotypeData& ioGenotypeData) override {}

private:

    vector<vector<CChromosome>>  mBestChromosomesPerRun;

    vector<vector<pair<CChromosome,CChromosome>>>  mAlgorithmResults;

    void PrintAllBests(const CGenotypeData& ioGenotypeData, std::ofstream& oStream) const;

};

