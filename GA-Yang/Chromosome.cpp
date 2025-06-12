
#include "Chromosome.h"
#include "GenotypeData.h"
#include "ChromosomeMooney.h"

#ifdef _MSC_VER
#include <windows.h>
#endif
#include <algorithm>

unsigned int CChromosome::sLength=0;
uint64_t CChromosome::sFitenssCalculationCounter=0;

void CChromosome::Init(unsigned int iLength)
{
    sLength = iLength;
    sFitenssCalculationCounter = 0;
    CChromosomeMooney::Init();
}

CChromosome& CChromosome::operator=(const CChromosome& input)
{
    if (this == &input)
    {
        return *this;
    }

    if (mLength != input.mLength)
    {
        mLength = input.mLength;

        delete [] mSNPs;
        delete [] mGenotypes;

        mSNPs = new unsigned int[mLength];
        mGenotypes = new BYTE[mLength];
    }

    // No need to reallocate
#ifdef _MSC_VER
    memcpy_s(mSNPs, mLength*sizeof(unsigned int), input.mSNPs, mLength*sizeof(unsigned int));
    memcpy_s(mGenotypes, mLength, input.mGenotypes, mLength);
#else
    memcpy(mSNPs, input.mSNPs, mLength*sizeof(unsigned int));
    memcpy(mGenotypes, input.mGenotypes, mLength);
#endif
    
    mFitness = input.mFitness;
    mFitenssValid = input.mFitenssValid;
    mAdditionalFitnessData = input.mAdditionalFitnessData;

    return *this;
}


