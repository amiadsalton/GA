
#pragma once 

#include "GenotypeData.h"
#include "ConfigParams.h"

#include <random>
#include <float.h>
#include <assert.h>
#include <string.h>

using namespace std;

////////////////////////////////////////////////////////////////////////
/// \brief    CChromosome is one element in GA population
///
/// Role: Represents one solution (chromosome of bits). 
///
/// Responsibilities: 
///      1. Access to SNPs that compose this chromosome and their genotypes.
///      2. Knows how to build random and consistent Chromosome
///      3. Interface for genetic operators: mutation and crossover
///      4. Interface for calculating the Fitness of this chromosome
////////////////////////////////////////////////////////////////////////
class CChromosome 
{
 
public:

    /// \brief    Empty constructor.
    CChromosome():
        mLength(sLength),
        mFitness(DBL_MAX), 
        mFitenssValid(false),
        mAdditionalFitnessData()
    {
        mSNPs = new unsigned int[mLength];
        mGenotypes = new BYTE[mLength];
        memset(mSNPs, 0, mLength*sizeof(unsigned int));
        memset(mGenotypes, 2, mLength);
    }

    /// \brief    predefined chromosome
    CChromosome(unsigned int iLength, const unsigned int* iSNPs, const BYTE* iGenotypes, double iFitness = DBL_MAX, const vector<string>& iAdditionaFitnesslInfo = vector<string>()):
        mLength(iLength),
        mFitness(iFitness),
        mAdditionalFitnessData(iAdditionaFitnesslInfo)
    {
        mSNPs = new unsigned int[mLength];
        mGenotypes = new BYTE[mLength];

        if (iSNPs != nullptr)
        {
#ifdef _MSC_VER
            memcpy_s(mSNPs, iLength*sizeof(unsigned int), iSNPs, iLength*sizeof(unsigned int));
#else
            memcpy(mSNPs, iSNPs, iLength*sizeof(unsigned int));
#endif
        }
        if (iGenotypes != nullptr)
        {
#ifdef _MSC_VER
            memcpy_s(mGenotypes, iLength, iGenotypes, iLength);
#else
            memcpy(mGenotypes, iGenotypes, iLength);
#endif
        }
        mFitenssValid = (mFitness != DBL_MAX);
    }

    /// \brief    predefined chromosome of length 1 
    CChromosome(unsigned int iSNP, BYTE iGenotypes = 2, double iFitness = DBL_MAX, const vector<string>& iAdditionaFitnesslInfo = vector<string>()):
        mLength(1),
        mFitness(iFitness),
        mAdditionalFitnessData(iAdditionaFitnesslInfo)
    {
        mSNPs = new unsigned int[mLength];
        mSNPs[0] = iSNP;
        mGenotypes = new BYTE[mLength];
        mGenotypes[0] = iGenotypes;
        mFitenssValid = (mFitness != DBL_MAX);
    }

    /// \brief    Copy constructor.
    CChromosome(const CChromosome& input) :
        mLength(input.mLength), mFitness(input.mFitness), mFitenssValid(input.mFitenssValid),
        mSNPs(new unsigned int[input.mLength]),  mGenotypes(new BYTE[input.mLength]),
        mAdditionalFitnessData(input.mAdditionalFitnessData)
    {
#ifdef _MSC_VER
        memcpy_s(mSNPs, mLength*sizeof(unsigned int), input.mSNPs, mLength*sizeof(unsigned int));
        memcpy_s(mGenotypes, mLength, input.mGenotypes, mLength);
#else
        memcpy(mSNPs, input.mSNPs, mLength*sizeof(unsigned int));
        memcpy(mGenotypes, input.mGenotypes, mLength);
#endif
    }

    /// \brief    Destructor.
    virtual ~CChromosome()
    {
        delete [] mSNPs;
        delete [] mGenotypes;
    }

    /// \brief    Assignment operator.
    CChromosome& operator=(const CChromosome& input);

    //virtual void Clone(const CChromosome& input) { CChromosome(input); }

    /// \brief    Assignment operator.
    virtual CChromosome& Assignment(const CChromosome& input) { return *this = input; }

    static void Init(unsigned int iLength);

    /// \brief    Comparison operator.
    bool operator==(const CChromosome& input) const;
    bool operator==(const vector<unsigned int>& input) const;

    bool operator<(const CChromosome& input) const;
    bool operator>(const CChromosome& input) const;

    double Fitness()  const { assert(mFitenssValid == true); return mFitness; }

    // Used in post processing
    void Fitness(double input)  { mFitness = input; mFitenssValid = true; }

	/// \brief    Fitness function. 
	/// \param[in] caseData   Genotype distributions of all analyzed SNPs in case dataset
	/// \param[in] controlData   Genotype distributions of all analyzed SNPs in control dataset
	/// \return int  Object's fitness
	virtual double Fitness(const vector<vector<BYTE>>& caseData,
		const vector<vector<BYTE>>& controlData) { return 0; }
	
	//-------------------
    // Getters
    //-------------------
    static unsigned int RequiredLength()  { return sLength; }
    
    static uint64_t FitnessCalculationCounter()  { return sFitenssCalculationCounter; }

    unsigned int Length() const { return mLength; }

    unsigned int SNP(unsigned int index) const { assert(index < mLength);  return mSNPs[index]; }
    unsigned int Genotype(unsigned int index) const { assert(index < mLength);  return mGenotypes[index]; }

	void SNP(unsigned int index, unsigned int iValue) { assert(index < mLength);   mSNPs[index] = iValue; mFitenssValid = false;  }
    void Genotype(unsigned int index, BYTE iValue) { assert(index < mLength);   mGenotypes[index] = iValue; mFitenssValid = false;	}

    const unsigned int* SNP() const {  return mSNPs; }
    BYTE* Genotype() const { return mGenotypes; }

    const vector<string>& AdditionalFitnessData() const { return mAdditionalFitnessData; }

protected:

		
    /// \brief    Chromosome required length (number of SNPs)
    static unsigned int sLength;
	
    /// \brief    Chromosome current length (number of SNPs)
    unsigned int mLength;

    /// \brief    Chromosome's SNPs. 
    unsigned int* mSNPs;

    /// \brief    Genotype of each SNP {1,2,3}. 
    BYTE* mGenotypes;

    /// \brief    Fitness of current object.
    double mFitness;

    /// \brief    true iff the object's contents has been changed from the last fitness calculation.
    bool mFitenssValid;

    /// \brief Additional  Fitness information
    vector<string> mAdditionalFitnessData;

    static uint64_t sFitenssCalculationCounter;

    /// \brief    Replace the value of 2 SNPs.
    /// \param[out] c1    first SNP to replace.
    /// \param[out] c2    second SNP to replace.
    static void ReplaceSNP(unsigned int& c1, unsigned int& c2);

    /// \brief    Replace the value of 2 Genotypes.
    /// \param[out] c1    first Genotype to replace.
    /// \param[out] c2    second Genotype to replace.
    static void ReplaceGenotype(BYTE& c1, BYTE& c2);
};

////////////////////////////////////////////////////////////////////////
/// \brief    CChromosome serializer 
////////////////////////////////////////////////////////////////////////
struct CChromosomeSerializer
{
    CChromosomeSerializer(const  CChromosome& iChromosome, const CGenotypeData& iData) : mChromosome(iChromosome), mData(iData) {}
    const  CChromosome& mChromosome;
    const CGenotypeData& mData;
};

inline void CChromosome::ReplaceSNP(unsigned int& c1, unsigned int& c2)
{
    unsigned int snp=c1;
    c1 = c2;
    c2 = snp;
}

inline void CChromosome::ReplaceGenotype(BYTE& c1, BYTE& c2)
{
    BYTE genotype=c1;
    c1 = c2;
    c2 = genotype;
}

inline bool CChromosome::operator<(const CChromosome& input) const
{
    return mFitness < input.mFitness;
}

inline bool CChromosome::operator>(const CChromosome& input) const
{
    return mFitness > input.mFitness;
}

inline bool CChromosome::operator==(const CChromosome& input) const 
{
    return ((memcmp(mSNPs, input.mSNPs, mLength * sizeof(int)) == 0) && (memcmp(mGenotypes, input.mGenotypes, mLength) == 0));
}

inline bool CChromosome::operator==(const vector<unsigned int>& input) const 
{
    for (size_t j=0; j<mLength ;j++)
    {
        if (mSNPs[j] != input[j] )
        {
            return false;
        }
        if (mGenotypes[j] != 2)
        {
            return false;
        }
    }
    return true;
}

inline std::ostream& operator<<(std::ostream& out, const CChromosomeSerializer& iSerializer)
{
    for (unsigned int i=0; i < iSerializer.mChromosome.Length() ;i++)
    { 
        out <<  "[" << iSerializer.mData.TotalCase()[iSerializer.mChromosome.SNP(i)].mIndex << "-" << iSerializer.mChromosome.Genotype(i) << "]";    
    }
    out << ";" << iSerializer.mChromosome.Fitness() << ";" ;

    if (!iSerializer.mChromosome.AdditionalFitnessData().empty())
    {
        for (auto elem : iSerializer.mChromosome.AdditionalFitnessData())
        {
            out << elem << ";";
        }
    }

    return out;
}

