
#include "Population.h"
#include "ChromosomeYang.h"

#include <boost/math/special_functions/factorials.hpp>

CPopulation::CPopulation(default_random_engine& randomEngine, CGenotypeData& genotypeData, unsigned int numberOfSNPs, unsigned int numberOfGeneotype) : 
    mRandomEngine(randomEngine), 
    mNumberOfSNPs(numberOfSNPs), mSize(0),
    mSnpHistogram(NULL), mGenotypeHistogram(NULL), mGenotypeData(genotypeData), mSelectionAlg(NULL),
    mSnpDistribution(0, numberOfSNPs-1), 
    mGenotypDistribution(0, CGenotypeData::NUMBER_OF_GENOTYPES - 1),
    mChormosomeDistribution(0, CChromosome::RequiredLength()-1),
    mTrappedIterations(0)
{
    mSelectionAlg = CSelection::CreateConcrete();

	mSize = numberOfSNPs / CChromosome::RequiredLength(); 

#ifndef ONLY_2_GENOTYPE
	mSize *= numberOfGeneotype;
#endif

    // Size must be odd
    mSize = mSize / 2 * 2;

	// Maximum is hard read from command line
	mSize = min(mSize, CConfigParams::GetTheInstance().PopulationSize);

	// Minimum is 100 (or half of Number of SNPs)
	mSize = max(int(mSize), min (100, int(numberOfSNPs / 2) ));

	// AMIAD- 2025 - unmark in order to get the size from configuraiton file 
	// mSize = CConfigParams::GetTheInstance().PopulationSize;

    // crossover is done for all chromosome in population
    mCrossoverDistribution = uniform_int_distribution<unsigned int>(0, mSize-1);
    mCrossoverRate = CConfigParams::GetTheInstance().CrossoverRate == 0 ? 0 :
    		max((unsigned int)(CConfigParams::GetTheInstance().CrossoverRate * mSize + 0.5),(unsigned int)(1));

    // Mutation is done for all element in chromosome
    mMutationDistribution = uniform_int_distribution<unsigned int>(0, (mSize * CChromosome::RequiredLength() * 2 - 1));
    mMutationRate = CConfigParams::GetTheInstance().MutationRate == 0 ? 0 :
    		max((unsigned int)(CConfigParams::GetTheInstance().MutationRate * mSize * CChromosome::RequiredLength() * 2 + 0.5), (unsigned int)(1));

    mMaxTrappedIterations = CConfigParams::GetTheInstance().TrapRatio == 0 ? 0 :
    		max((unsigned int)(CConfigParams::GetTheInstance().TrapRatio * CConfigParams::GetTheInstance().NumOfIterations + 0.5) , (unsigned int)(1));
    mVibratedElements = CConfigParams::GetTheInstance().VibrationRate == 0 ? 0 :
    		max((unsigned int)(CConfigParams::GetTheInstance().VibrationRate * mSize + 0.5), (unsigned int)(1));
}

CPopulation::~CPopulation()
{
    delete [] mSnpHistogram;
    delete [] mGenotypeHistogram;

    delete mSelectionAlg;
}
