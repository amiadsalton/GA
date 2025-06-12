
#include "ChromosomeMooney.h"
#include "ConfigParams.h"
#include "GenotypeData.h"

#ifdef _MSC_VER
#include <windows.h>
#endif

unsigned int CChromosomeMooney::sChi2TableSize=0;
boost::math::chi_squared CChromosomeMooney::sBoostDist(2);
unsigned int* CChromosomeMooney::sTableEntryWeight=nullptr;

//====================
// CChromosomeMooney::CChromosomeMooney(const CChromosomeMooney& input) :
// CChromosome (input)
// {
//     mObserved = new int[sChi2TableSize];
//     mExpected = new double[sChi2TableSize];
// 
//     memcpy(mObserved, input.mObserved, sChi2TableSize);
//     memcpy(mExpected, input.mExpected, sChi2TableSize);
// }


//====================
// Method: Build
// Responsibility: builds random (and consistent) solution 
//====================
void CChromosomeMooney::Build(default_random_engine& randomEngine, 
                        uniform_int_distribution<unsigned int>& snpDistribution,
                        uniform_int_distribution<unsigned int>& genotypeDistribution)
{
    for (unsigned int i=0; i<mLength; i++)
    {
        mSNPs[i] = snpDistribution(randomEngine);
    }
    Fix(randomEngine, snpDistribution);
    mFitenssValid = false;
}

//====================
// Method: Fix
// Responsibility: Confirm that there are not 2 identical SNPs. If there are then replace one of them
//====================
void CChromosomeMooney::Fix(default_random_engine& randomEngine, 
                          uniform_int_distribution<unsigned int>& snpDistribution)
{
    for (unsigned int i=0; i < mLength; i++)
    {
        bool exist;
        do 
        {
            exist = false;
            for (unsigned int j=0; j < i; j++)
            {
                if (mSNPs[i] == mSNPs[j])
                {
                    // Need to replace i 
                    mSNPs[i] = snpDistribution(randomEngine);
                    exist = true;
                    mFitenssValid = false;
                    break;
                }
            }
        } while (exist != false);
    }

    // Sort to make oprerator== more efficient
//     if (mFitenssValid == false)
//     {
        Sort();
//    }
}

//====================
// Method: Crossover
// Responsibility: Perform crossover with input solution
// Method: Define a random index and replace the chromosomes from this index with input's chromosome
//====================
void CChromosomeMooney::CrossOver(CChromosomeMooney& mate,
                            default_random_engine& randomEngine,
                            uniform_int_distribution<unsigned int>& chormosomeDistribution)
{
    // Non crossable regions (between p1 and p2 for SNP and between p3 and p4 for genotypes)
    unsigned int p1=chormosomeDistribution(randomEngine);
    unsigned int p2=chormosomeDistribution(randomEngine);

    unsigned int end=min(p1+p2, mLength);
    for (unsigned int i=p1; i<end ; i++)
    {
        ReplaceSNP(mSNPs[i], mate.mSNPs[i]);
    }

    mFitenssValid = false;
    mate.mFitenssValid = false;
}

//====================
// Method: Mutation
// Responsibility: Performs mutation
// Method: Get a random index, where to mutate the chromosome and replace it with its
//         subsequent individual of the same type.
//====================
void CChromosomeMooney::Mutation(unsigned int mutationRate,
                           default_random_engine& randomEngine,
                           uniform_int_distribution<unsigned int>& mutationDistribution,
                           uniform_int_distribution<unsigned int>& snpDistribution,
                           uniform_int_distribution<unsigned int>& genotypeDistribution)
{
    for (unsigned int i=0 ;i<mLength ;i++)
    {
        if (mutationDistribution(randomEngine) < mutationRate)
        {
            mSNPs[i] = snpDistribution(randomEngine);
            mFitenssValid = false;
        }
    }
}

//====================
// Method: Fitness
// Responsibility: Calculate the fitness of current solution
// Comment: Since we are deal with find minimum problem then the evaluation value is reduce from
//          constant which is the maximum possible fitness.
//====================
double CChromosomeMooney::CalculateFitness(const vector<vector<BYTE>>& caseData,
                                 const vector<vector<BYTE>>& controlData)
{
    for (unsigned int i=0; i< sChi2TableSize; i++)
    {
        mExpected[i] = 1;
        mObserved[i] = 1;
    }

    // Find Expected and Observed
    for (unsigned int i=0; i<caseData.size() ; i++)
    {
//         unsigned int entry = sTableEntryWeight[0] * caseData[i][mSNPs[0]];
//         entry += sTableEntryWeight[1] * caseData[i][mSNPs[1]];
//         entry += sTableEntryWeight[2] * caseData[i][mSNPs[2]];
//         entry += sTableEntryWeight[3] * caseData[i][mSNPs[3]];
        
        unsigned int entry = 0;
        for (unsigned int j=0; j<mLength ;j++)
        {
            entry += (sTableEntryWeight[j] * caseData[i][mSNPs[j]]);
        }    
         mObserved[entry]++;

//         entry = sTableEntryWeight[0] * controlData[i][mSNPs[0]];;
//         entry += sTableEntryWeight[1] * controlData[i][mSNPs[1]];
//         entry += sTableEntryWeight[2] * controlData[i][mSNPs[2]];
//         entry += sTableEntryWeight[3] * caseData[i][mSNPs[3]];

        entry = 0;
        for (unsigned int j=0; j<mLength ;j++)
        {
            entry += (sTableEntryWeight[j] * controlData[i][mSNPs[j]]);
        }    
        mExpected[entry]++;
    }

    double chi2=0;
    for (unsigned int i=0; i<sChi2TableSize ;i++)
    {
        chi2 += pow(mExpected[i]-mObserved[i],2) / mExpected[i];
    }

//     if (chi2 > 2)
//     {
//         chi2 = log(chi2);
//     }

    mFitness = chi2;//boost::math::cdf(sBoostDist,chi2);

    mFitenssValid = true;

    return mFitness;
}

void CChromosomeMooney::Init()
{
    sChi2TableSize = (unsigned int)(pow(3,sLength));
    sBoostDist = boost::math::chi_squared(sChi2TableSize-1);

    delete [] sTableEntryWeight;
    sTableEntryWeight = new uint32_t[sLength];
    for (unsigned int j=0; j<sLength ;j++)
    {
        sTableEntryWeight[j] = (unsigned int)(pow(3,(sLength - 1 - j)));
    }  
}

double CChromosomeMooney::Fitness(const vector<vector<BYTE>>& caseData, const vector<vector<BYTE>>& controlData)
{
	if (mFitenssValid == true)
	{
		return mFitness;
	}
#pragma omp critical
	{
		sFitenssCalculationCounter++;
	}
	return CalculateFitness(caseData, controlData);
}

