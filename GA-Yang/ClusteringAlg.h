#pragma once
#include "PostProcessingAlg.h"

////////////////////////////////////////////////////////////////////////
/// \brief    Clustering algorithm
////////////////////////////////////////////////////////////////////////
class CClusteringAlg : public CPostProcessingAlg
{
public:
    CClusteringAlg(CConfigParams::EPostProcessingAlgorithm iAlgorithm) : CPostProcessingAlg(), mCurrentOrder(-1), mAlgorithm(iAlgorithm) {}
	virtual ~CClusteringAlg();

    virtual void Update(const std::vector<CChromosome>& iElitism) override;
    virtual bool DerivedExecute(const CGenotypeData& ioGenotypeData, int iOrder) override;
    virtual void Print(const CGenotypeData& ioGenotypeData, const std::string& iBaseFileName, int iOrder) const override;
    virtual void UpdateGenotypeData(CGenotypeData& ioGenotypeData) override;

private:

    // Histogram for each order from 1 to [Barcode length - 1]. Each element in the histogram is a Chromosome and the number of times it appears
    vector<vector<CChromosome>>  mBestSnpsPerRun;

    // The best Chromosome of current run
    CChromosome  mBestChromosome;

    // Best Chromosomes of each run  (as string + fitness) 
    vector<pair <string, double> >  mBestChromosomePerRun;

    // Barcode of the best combination of each run
    vector<pair <string, double> >  mBestBarcodePerRun;

    vector<vector<pair<CChromosome, unsigned int>>>  mAlgorithmResults;

    int mCurrentOrder;

    CConfigParams::EPostProcessingAlgorithm mAlgorithm;

	void CalcBestBarcode(unsigned int iOrder, const CChromosome& iChromsome, unsigned int iCurrIndex, CChromosome*& oBestChromosome, const CGenotypeData& ioGenotypeData) const;

};

