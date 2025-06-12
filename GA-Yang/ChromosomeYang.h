
#pragma once 

#include "GenotypeData.h"
#include "Chromosome.h"

#include <random>
#include <float.h>
#include <omp.h>

using namespace std;

//#define  ONLY_2_GENOTYPE
//#define ORIGINAL_FITNESS

////////////////////////////////////////////////////////////////////////
/// \brief    CChromosomeYang is one element in GA population
///
/// Role: Represents chromosome that executes Yang’s fitness function 
///
/// Responsibilities: 
///         1.	Implements the interface that is defined by the base class CChromosome.
///         2.	Knows how to calculate Yang’s fitness function.
///         
////////////////////////////////////////////////////////////////////////
class CChromosomeYang : public CChromosome
{
 
public:

    /// \brief    Empty constructor.
    CChromosomeYang() : CChromosome() {}

    /// \brief    predefined chromosome
    CChromosomeYang(unsigned int iLength, const unsigned int* iSNPs, const BYTE* iGenotypes) : CChromosome(iLength, iSNPs, iGenotypes) {}

    /// \brief    Copy constructor.
    CChromosomeYang(const CChromosomeYang& input) : CChromosome(input) {}

    /// \brief    Copy constructor.
    CChromosomeYang(const CChromosome& input) : CChromosome(input) {}

    /// \brief    Destructor.
    virtual ~CChromosomeYang() {}

    /// \brief    Assignment operator.
    CChromosomeYang& operator=(const CChromosomeYang& input);

    /// \brief    Assignment operator.
    CChromosomeYang& Assignment(const CChromosome& input) override { return *this = dynamic_cast<const CChromosomeYang&>(input); }

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
    void CrossOver(CChromosomeYang& mate, 
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

    static const unsigned int NUMBER_OF_GENOTYPES;

private:

    void Sort();

    /// \brief    Calculate the fitness value. 
    /// \param[in] caseData   Genotype distributions of all analyzed SNPs in case dataset
    /// \param[in] controlData   Genotype distributions of all analyzed SNPs in control dataset
    /// \return int  Object's fitness
    double CalculateFitness(const vector<vector<BYTE>>& caseData,
                            const vector<vector<BYTE>>& controlData);
};

inline void CChromosomeYang::Sort()
{
    for (unsigned i=0; i<mLength ;i++)
    {
        for (unsigned int j=0; j<mLength ;j++)
        {
            if (mSNPs[i] < mSNPs[j])
            {
                ReplaceSNP(mSNPs[i], mSNPs[j]);
                ReplaceGenotype(mGenotypes[i], mGenotypes[j]);
            }
        }
    }
}