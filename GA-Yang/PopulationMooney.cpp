
#include "PopulationMooney.h"
#include "ConfigParams.h"
#include "GenotypeData.h"
#include "Selection.h"
#include "ChromosomeMooney.h"

#include <math.h>
#include <algorithm>
#include <omp.h>

CPopulationMooney::CPopulationMooney(default_random_engine& randomEngine, CGenotypeData& genotypeData, unsigned int numberOfSNPs) : 
    CPopulation(randomEngine, genotypeData, numberOfSNPs, 1), mCurrentGeneration(), mNextGeneration()
{
}

CPopulationMooney::~CPopulationMooney()
{
    delete [] mCurrentGeneration;
    delete [] mNextGeneration;

//     for (size_t i=0; i<mSize ;i++)
//     {
//          delete mCurrentGeneration[i];
//          delete mNextGeneration[i];
//     }
}

void CPopulationMooney::Init()
{    
    mCurrentGeneration = new CChromosomeMooney[mSize];
    mNextGeneration = new CChromosomeMooney[mSize];

    // Build chromosomes
    for (unsigned int i=0; i<mSize; i++) 
    {
        //mNextGeneration[i] = new CChromosomeMooney();
        //mCurrentGeneration[i] = new CChromosomeMooney();
        mCurrentGeneration[i].Build(mRandomEngine, mSnpDistribution, mGenotypDistribution);

        // If it has already been created then drop it (by decreasing the counter)
        for (unsigned j=0; j<i ; j++)
        {
            if ((mCurrentGeneration[j] == mCurrentGeneration[i]))
            {
                i--;
                break;
            }
        }
    }

    Fiteness();
}

// Return whether first element is greater than the second
// bool ComparePred( CChromosomeMooney elem1, CChromosomeMooney elem2 )
// {
//     return elem1 < elem2;
// }

void CPopulationMooney::Fiteness()
{
    // Calculate fitness 
    mFitnessMax = FLT_MIN;
    mFitnessMin = FLT_MAX;
    mFitnessTotal = 0;
    #pragma omp parallel for
    for (int i=0; i<int(mSize); i++) 
    {
        double fitness=mCurrentGeneration[i].Fitness(mGenotypeData.CaseData(), mGenotypeData.ControlData());
        #pragma omp critical
        {
            mFitnessMin = min(fitness, mFitnessMin);
            mFitnessMax = max(fitness, mFitnessMax);
        }
        mFitnessTotal += fitness;
    }
    mFitnessMean = float(mFitnessTotal) / mSize;
    mFitnessTotalNormalized = (mFitnessMean - mFitnessMin) * mSize;

    mFitenssStd = 0;
    #pragma omp parallel for
    for (int i=0; i<int(mSize) ; i++)
    {
        mFitenssStd += pow((mCurrentGeneration[i].Fitness() - mFitnessMean), 2);
    }
    mFitenssStd = sqrt(mFitenssStd) / mSize;

    sort(mCurrentGeneration, mCurrentGeneration + mSize);//, ComparePred);

    // Initiate selection algorithm before starting selection
    mSelectionAlg->ResetSelection(*this);

#ifdef _STATISTICAL_DATA
    if (mSnpHistogram == NULL)
    {
        mSnpHistogram = new unsigned int[mNumberOfSNPs];
        mGenotypeHistogram = new unsigned int[CGenotypeData::NUMBER_OF_GENOTYPES];
    }

    memset(mSnpHistogram, 0 , mNumberOfSNPs*sizeof(unsigned int));
    memset(mGenotypeHistogram, 0 , CGenotypeData::NUMBER_OF_GENOTYPES*sizeof(unsigned int));
    for (unsigned int i=0; i<mSize; i++) 
    {
        for (unsigned int j=0; j<mCurrentGeneration[i].Length(); j++)
        {
            mSnpHistogram[mCurrentGeneration[i].SNP(j)]++; 
            mGenotypeHistogram[mCurrentGeneration[i].Genotype(j)]++;  
        }
    }

    // Print statistics
    cout << "\n\nSNPs distribution: \n";
    for (unsigned int i=0; i<mNumberOfSNPs; i++)
    {
        cout << i+1 << " = " << mSnpHistogram[i] << endl;
    }

    cout << "\n\nGenotypes distribution: \n";
    for (int i=0; i<CGenotypeData::NUMBER_OF_GENOTYPES; i++)
    {
        cout << i+1 << " = " << mGenotypeHistogram[i] << endl;
    }
#endif
}


// Method: NextGeneration
// Algorithm:
//          1. Elitism - copy best solutions (ratio can be controlled)
//          2. Until there are enough solutions do:
//              a. Select 2 solutions 
//              b. In user defined probability perform crossover between them and select the first to be a new solution
//              c. Perform mutation for the selected solution.
//              d. Calculates solution's evaluation value
//          3. Copy the new generation into current generation
//          4. Calculates the fitness of current generation and sort it
//          5. Vibrate population in case of trap in local maximum
void CPopulationMooney::NextGeneration()
{
    double bestFitness = mFitnessMax;

    // Phase 1.
    unsigned int elitismIndex = Elitism();

    //omp_set_num_threads(4);
    #pragma omp parallel for
    for (int i = elitismIndex; i < int(mSize); i+=2)
    {
        // Stage 2. a.
        unsigned int firtIndex=mSelectionAlg->Select(*this, mRandomEngine);
        unsigned int secondIndex=mSelectionAlg->Select(*this, mRandomEngine,firtIndex);
        mNextGeneration[i] = mCurrentGeneration[firtIndex];
        mNextGeneration[i+1] = mCurrentGeneration[secondIndex]; 

        // Stage 2. b.
        if (mCrossoverDistribution(mRandomEngine) < mCrossoverRate)
        {
            // Perform crossover only for chromosomes according to required rate
            mNextGeneration[i].CrossOver(mNextGeneration[i+1], mRandomEngine, mChormosomeDistribution);
        }

        mNextGeneration[i].Mutation(mMutationRate, mRandomEngine, mMutationDistribution, mSnpDistribution, mGenotypDistribution);
        mNextGeneration[i+1].Mutation(mMutationRate, mRandomEngine, mMutationDistribution, mSnpDistribution, mGenotypDistribution);

        // Fix chromosomes if they became invalid by crossover or mutation
        mNextGeneration[i].Fix(mRandomEngine, mSnpDistribution);
        mNextGeneration[i+1].Fix(mRandomEngine, mSnpDistribution);

//         if (CheckForDuplication(mNextGeneration[i], i))
//         {
//             mNextGeneration[i].Build(mRandomEngine, mSnpDistribution, mGenotypDistribution);
//         }
//         if (CheckForDuplication(mNextGeneration[i+1], i+1))
//         {
//             mNextGeneration[i+1].Build(mRandomEngine, mSnpDistribution, mGenotypDistribution);
//         }
    }

    // 3. Replace current generation
    for (unsigned int i = 0; i < mSize ; i++)
    {
        mCurrentGeneration[i] = mNextGeneration[i];
    }

    // 4. Calculate the fitness of population
    Fiteness();

    // 5. Vibrate
    Vibrate(bestFitness);
}

bool CPopulationMooney::CheckForDuplication(const CChromosomeMooney& input, unsigned int currentIndex)
{
    // If it has already been created then drop it (by decreasing the counter)
    for (unsigned int j=0; j < currentIndex ; j++)
    {
        // Search for identical element in next generation group
        if (mNextGeneration[j] == input)
        {
            return true;
        }
    }

    return false;
}

void CPopulationMooney::DifferentChromosomes(vector<const CChromosome*>& oChromosomes, int iNumberOfSolutions) const 
{
    if (iNumberOfSolutions == -1)
    {
        iNumberOfSolutions = mSize;
    }
    int currentIndex = mSize-1;

    oChromosomes.clear();
    while (oChromosomes.size() < iNumberOfSolutions && currentIndex > 0)
    {
        const CChromosome& tmp=mCurrentGeneration[currentIndex--];
        int i=0;
        for (; i<oChromosomes.size() ;i++)
        {
            if (*oChromosomes[i] == tmp)
            {
                break;
            }
        }
        if (i == oChromosomes.size())
        {
            oChromosomes.push_back(&tmp);
        }
    }
}

unsigned int CPopulationMooney::Elitism()
{
    unsigned int elitismNum = ElitismGroupSize();

    // Take the best different "elitismNum" solutions from current generation
    vector<const CChromosome*>  differentChromosomes;
    DifferentChromosomes(differentChromosomes, elitismNum);
    size_t i=0;
    unsigned int nextGenIndex(0);
    for (; i < differentChromosomes.size() ;i++)
    {
        mNextGeneration[nextGenIndex++] = *dynamic_cast<const CChromosomeMooney*>(differentChromosomes[i]);
    }

    for (; i < elitismNum ;i++)
    {
        // No more different solutions for elitism, so find new solutions   
        mNextGeneration[nextGenIndex++].Build(mRandomEngine, mSnpDistribution, mGenotypDistribution);
    }

    return elitismNum;
}

void CPopulationMooney::Vibrate(double lastBestFitness)
{
    if (mVibratedElements == 0)
    {
        // No point to continue, because no element should be replaced by vibrate
        return;
    }

    if (abs(lastBestFitness - mFitnessMax) < FLT_EPSILON)
    {
        mTrappedIterations++;
    }
    else
    {
        mTrappedIterations = 0;
    }

    if (mTrappedIterations < mMaxTrappedIterations)
    {
        return;  // No need to vibrate
    }

    // Copy the elitism elements first 
    unsigned int index = Elitism();
    for (unsigned int i=0; i < index; i++)
    {
        mCurrentGeneration[mSize - 1 - i] = mNextGeneration[i];
    }

    mTrappedIterations = 0;
    for (int i=mVibratedElements-1; i>=0 ; i--)
    {
        mCurrentGeneration[i].Build(mRandomEngine, mSnpDistribution, mGenotypDistribution);

        // If it has already been created then drop it (by decreasing the counter)
//         for (unsigned j=i+1; j<mSize ; j++)
//         {
//             if (mCurrentGeneration[j] == mCurrentGeneration[i])
//             {
//                 i++;
//                 break;
//             }
//         } 
    }

    Fiteness();
}