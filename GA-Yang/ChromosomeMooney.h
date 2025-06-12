
#pragma once 

#include "Chromosome.h"
#include <boost/math/distributions/chi_squared.hpp>

using namespace std;

////////////////////////////////////////////////////////////////////////
/// \brief    CChromosomeMooney is one element in GA population
///
/// Role: Represents chromosome that executes Mooney's fitness function 
///
/// Responsibilities: 
///         1.	Implements the interface that is defined by the base class CChromosome.
///         2.	Knows how to calculate Mooney’s fitness function.
///         
////////////////////////////////////////////////////////////////////////
class CChromosomeMooney : public CChromosome 
{
 
public:

    /// \brief    Empty constructor.
    CChromosomeMooney() : CChromosome(), mObserved(new int[sChi2TableSize]), mExpected(new int[sChi2TableSize])
    {
    }

    /// \brief    predefined chromosome
    CChromosomeMooney(unsigned int iLength, const unsigned int* iSNPs, const BYTE* iGenotypes) : 
        CChromosome(iLength, iSNPs, iGenotypes),
        mObserved(new int[sChi2TableSize]),
        mExpected(new int[sChi2TableSize]) {}

    /// \brief    Copy constructor.
    CChromosomeMooney(const CChromosomeMooney& input) : CChromosome(input), mObserved(new int[sChi2TableSize]), mExpected(new int[sChi2TableSize]) {}

    /// \brief    Copy constructor.
    CChromosomeMooney(const CChromosome& input) : CChromosome(input), mObserved(new int[sChi2TableSize]), mExpected(new int[sChi2TableSize]) {}

    /// \brief    Destructor.
    virtual  ~CChromosomeMooney()
    {
        delete [] mObserved;
        delete [] mExpected;
    }

    /// \brief    Assignment operator.
    CChromosomeMooney& operator=(const CChromosomeMooney& input);

    //virtual void Clone(const CChromosomeMooney& input) { ::CChromosomeMooney(input); }

    /// \brief    Assignment operator.
    CChromosomeMooney& Assignment(const CChromosome& input) override { return *this = dynamic_cast<const CChromosomeMooney&>(input); }

    static void Init();

    /// \brief    Builds random (and consistent) chromosome 
    /// \param[in] randomEngine    Pseudo random generator which is used to for initialization.
    /// \param[in] snpDistribution    Distribution of SNPs (used for random generations)
    /// \param[in] genotypeDistribution    Distribution of Genotypes (used for random generations)
    void Build(default_random_engine& randomEngine,  
               uniform_int_distribution<unsigned int>& snpDistribution,
               uniform_int_distribution<unsigned int>& genotypeDistribution);

    /// \brief    Confirm that there are no identical SNPs in barcode  
    /// \param[in] randomEngine    Pseudo random generator which is used to for initialization.
    /// \param[in] snpDistribution    Distribution of SNPs (used for random generations)
    void Fix(default_random_engine& randomEngine, 
             uniform_int_distribution<unsigned int>& snpDistribution);

    /// \brief    Fitness function.
    /// \param[in] caseData   Genotype distributions of all analyzed SNPs in case dataset
    /// \param[in] controlData   Genotype distributions of all analyzed SNPs in control dataset
    /// \return int  Object's fitness
    virtual double Fitness(const vector<vector<BYTE>>& caseData,
						   const vector<vector<BYTE>>& controlData) override;

	using CChromosome::Fitness;

    //-------------------
    // Genetic operators
    //-------------------

    /// \brief    Crossover operator. 
    /// \param[in] mate    mate for crossover (is changed as well)
    /// \param[in] randomEngine    Pseudo random generator which is used to for snpDistribution.
    /// \param[in] snpDistribution    Distribution of chromosome (used for random selection of crossover points)
    void CrossOver(CChromosomeMooney& mate, 
                           default_random_engine& randomEngine,
                           uniform_int_distribution<unsigned int>& chormosomeDistribution);

    /// \brief    Mutation operator. 
    /// \param[in] mutationRate    Mutation rate which is used to decide whether to perform mutation or not    
    /// \param[in] randomEngine    Pseudo random generator which is used to for mutationDistribution.
    /// \param[in] mutationDistribution    Distribution of mutation (used for random selection of mutation points)
    /// \param[in] snpDistribution    Distribution of SNPs (used to select new SNP in case of mutation)
    /// \param[in] genotypeDistribution    Distribution of Genotypes ((used to select new genotype in case of mutation)
    void Mutation(unsigned int mutationRate,
                          default_random_engine& randomEngine,
                          uniform_int_distribution<unsigned int>& mutationDistribution,
                          uniform_int_distribution<unsigned int>& snpDistribution,
                          uniform_int_distribution<unsigned int>& genotypeDistribution);


protected:

    // Chi square data
    static unsigned int sChi2TableSize;
    static boost::math::chi_squared sBoostDist;

    int* mObserved;
    int* mExpected;

    // Used for Chi square table entry 
    static unsigned int* sTableEntryWeight;

private:

    void Sort() { std::sort(mSNPs, mSNPs + mLength); }

    /// \brief    Calculate the fitness value. 
    /// \param[in] caseData   Genotype distributions of all analyzed SNPs in case dataset
    /// \param[in] controlData   Genotype distributions of all analyzed SNPs in control dataset
    /// \return int  Object's fitness
    double CalculateFitness(const vector<vector<BYTE>>& caseData,
                            const vector<vector<BYTE>>& controlData);
};

inline CChromosomeMooney& CChromosomeMooney::operator=(const CChromosomeMooney& input)
{
    if (this == &input)
    {
        return *this;
    }

    CChromosome::operator=(input);

    return *this;
}
