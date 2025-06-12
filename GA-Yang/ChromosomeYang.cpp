
#include "ChromosomeYang.h"
#include "ConfigParams.h"
#include "GenotypeData.h"

#ifdef _MSC_VER
#include <windows.h>
#endif
#include <algorithm>
 
#ifdef ONLY_2_GENOTYPE
const unsigned int CChromosomeYang::NUMBER_OF_GENOTYPES=1;
#else
const unsigned int CChromosomeYang::NUMBER_OF_GENOTYPES=3;
#endif

//====================
CChromosomeYang& CChromosomeYang::operator=(const CChromosomeYang& input)
{
    if (this == &input)
    {
        return *this;
    }

    CChromosome::operator=(input);

    return *this;
}

//====================
// Method: Build
// Responsibility: builds random (and consistent) solution 
//====================
void CChromosomeYang::Build(default_random_engine& randomEngine, 
                        uniform_int_distribution<unsigned int>& snpDistribution,
                        uniform_int_distribution<unsigned int>& genotypeDistribution)
{
    for (unsigned int i=0; i<mLength; i++)
    {
        mSNPs[i] = snpDistribution(randomEngine);
    }
    Fix(randomEngine, snpDistribution);

    for (unsigned int i=0; i<mLength; i++)
    {
#ifdef ONLY_2_GENOTYPE
        mGenotypes[i] = 2;
#else
        mGenotypes[i] = genotypeDistribution(randomEngine);
#endif
    }
    mFitenssValid = false;
}

//====================
// Method: Fix
// Responsibility: Confirm that there are not 2 identical SNPs. If there are then replace one of them
//====================
void CChromosomeYang::Fix(default_random_engine& randomEngine, 
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
//    if (mFitenssValid == false)
//    {
        Sort();
//    }
}

//====================
// Method: Crossover
// Responsibility: Perform crossover with input solution
// Method: Define a random index and replace the chromosomes from this index with input's chromosome
//====================
void CChromosomeYang::CrossOver(CChromosomeYang& mate,
                            default_random_engine& randomEngine,
                            uniform_int_distribution<unsigned int>& chormosomeDistribution)
{
    // Non crossable regions (between p1 and p2 for SNP and between p3 and p4 for genotypes)
    unsigned int p1=chormosomeDistribution(randomEngine);
    unsigned int p2=chormosomeDistribution(randomEngine);

#ifndef ONLY_2_GENOTYPE  
     unsigned int p3=chormosomeDistribution(randomEngine);
     unsigned int p4=chormosomeDistribution(randomEngine);
#endif

    unsigned int end=min(p1+p2, mLength);
    for (unsigned int i=p1; i<end ; i++)
    {
        ReplaceSNP(mSNPs[i], mate.mSNPs[i]);
    }

#ifndef ONLY_2_GENOTYPE  
     end=min(p3+p4, mLength);
     for (unsigned int i=p3; i<end ; i++)
     {
         ReplaceGenotype(mGenotypes[i], mate.mGenotypes[i]);
     }
#endif

    mFitenssValid = false;
    mate.mFitenssValid = false;
}

//====================
// Method: Mutation
// Responsibility: Performs mutation
// Method: Get a random index, where to mutate the chromosome and replace it with its
//         subsequent individual of the same type.
//====================
void CChromosomeYang::Mutation(unsigned int mutationRate,
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

#ifndef ONLY_2_GENOTYPE  
     for (unsigned int i=0 ;i<mLength ;i++)
     {
         if (mutationDistribution(randomEngine) < mutationRate)
         {
             mGenotypes[i] = genotypeDistribution(randomEngine);
             mFitenssValid = false;
         }
     }
#endif
}

double CChromosomeYang::Fitness(const vector<vector<BYTE>>& caseData, const vector<vector<BYTE>>& controlData)
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

//====================
// Method: Fitness
// Responsibility: Calculate the fitness of current solution
// Comment: Since we are deal with find minimum problem then the evaluation value is reduce from
//          constant which is the maximum possible fitness.
//====================
double CChromosomeYang::CalculateFitness(const vector<vector<BYTE>>& caseData,
                            const vector<vector<BYTE>>& controlData)
{

//     mSNPs[0] = 0;
//     mSNPs[1] = 497;
// 
//     for (int kk = 0 ; kk<3 ;kk++)
//     {
//         for (int jj = 0 ; jj<3 ;jj++)
//         {
// 
//         mGenotypes[0] = kk;
//         mGenotypes[1] = jj;

    unsigned int totalCase(0);

    // Iterates over all items (genotype of one people) in case data set, and for each item check 
    // if it has match to current chromosome's genotype
    for (size_t i=0; i<caseData.size(); i++)
    {
        bool match(true);
        for (unsigned int j=0; j<mLength; j++)
        {
            if (caseData[i][mSNPs[j]] != mGenotypes[j])
            {
                match = false;
                break;
            }
        }
        if (match)
        {
            totalCase++;
        }
    }


    // Get the number of occurrences of all genotypes in case data set
    unsigned int totalControl(0);
    for (size_t i=0; i<controlData.size(); i++)
    {
        bool match(true);
        for (unsigned int j=0; j<mLength; j++)
        {
            if (controlData[i][mSNPs[j]] != mGenotypes[j])
            {
                match = false;
                break;
            }
        }
        if (match)
        {
            totalControl++;
        }
    }

#ifdef ORIGINAL_FITNESS
     mFitness = int(totalControl) - int(totalCase);
#else
     mFitness = (totalCase+1.0f) / (totalControl+1.0f);
     if (mFitness < 1 && mFitness > FLT_EPSILON)
     {
         mFitness = - (1.0f / mFitness);
     }
#endif

     mAdditionalFitnessData.clear();
     mAdditionalFitnessData.push_back(to_string(totalCase));
     mAdditionalFitnessData.push_back(to_string(totalControl));
//         }
//     }

    mFitenssValid = true;

    return mFitness;
}



