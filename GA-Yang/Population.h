#pragma once


#include "GenotypeData.h"
#include "Chromosome.h"
#include "Selection.h"

#include <random>


using namespace std;

////////////////////////////////////////////////////////////////////////
/// \brief    Population of one generation
///
/// Role: Represents the population (chromosomes collection) of one generation 
///
/// Responsibilities: 
///      1. Holds all Chromosomes that belong to this generation
///      2. Knows how to build the first generation
///      3. Knows how to build the next generation by:
///         a. Selects n/2 pairs of chromosomes to be the parents of the next generation
///         b. Reproduces each pair 
///         c. Knows how to execute the mutation and crossover in reproduction process
///         d. Executes the genetic operators Elitism and Vibration next generation
///      4. Calculates the fitness and statistical data of current population
/// 
////////////////////////////////////////////////////////////////////////
class CPopulation
{

public:

    /// \brief    Empty constructor.
    /// \param[in] randomEngine    Pseudo random generator which is used to for initialization and evolution.
    /// \param[in] genotypeData    input dataset.
    /// \param[in] numberOfSNPs    Number of analyzed SNPs.
    CPopulation(default_random_engine& randomEngine, CGenotypeData& genotypeData, unsigned int numberOfSNPs, unsigned int numberOfGeneotype);

    virtual ~CPopulation();

    /// \brief    Created the initial Population
    virtual void Init() = 0;

    /// \brief    Prepare the next generation from current generation
    virtual void NextGeneration() = 0;

    //-------------------
    // Getters
    //-------------------
    const unsigned int* SnpHistogram() const { return mSnpHistogram; }
    const unsigned int* GenotypeHistogram() const { return mGenotypeHistogram; }
    double FitnessMin() const { return mFitnessMin; }
    double FitnessMax() const { return mFitnessMax; }
    double FitnessTotal() const { return mFitnessTotal; }
    double FitnessTotalNormalized() const { return mFitnessTotalNormalized; }
    double FitnessMean() const { return mFitnessMean; }
    double FitenssStd() const { return mFitenssStd; }
    unsigned int Size() const { return mSize; }

    virtual const CChromosome& operator[](size_t index) const = 0;

    /// \brief    brief Find all different chromosome in current generation
    /// \param[out] oChromosomes    list of all different chromosomes.
    /// \param[in] iNumberOfSolutions    Number of different items to find, if  =-1 then all different items are retrieved.
    virtual void DifferentChromosomes(vector<const CChromosome*>& oChromosomes, int ioNumberOfSolutions = -1) const = 0;

    // The population size depends on chromosome length and number of SNPs
    unsigned int PopulationSize() { return mSize; }

    /// \brief    Size of elitism group
    unsigned int ElitismGroupSize() const;

protected:

    /// \brief     Population  size  
    unsigned int mSize;

    /// \brief      Pseudo random generator which is used to for initialization and evolution
    default_random_engine&   mRandomEngine;

    /// \brief      Number of analyzed SNPs
    unsigned int mNumberOfSNPs;

    /// \brief      Reference to input genotypes data 
    CGenotypeData& mGenotypeData;

    CSelection* mSelectionAlg;

    /// \brief    Used for random selection of SNPs.
    uniform_int_distribution<unsigned int> mSnpDistribution;

    /// \brief    Used for random selection of Genotype.
    uniform_int_distribution<unsigned int> mGenotypDistribution;

    /// \brief    Used for random selection point for crossover and mutation
    uniform_int_distribution<unsigned int> mChormosomeDistribution;

    /// \brief    Used for random selection point for crossover
    uniform_int_distribution<unsigned int> mCrossoverDistribution;

    /// \brief    Crossover rate (translated to population size)
    unsigned int mCrossoverRate;

    /// \brief    Used for random selection point for mutation
    uniform_int_distribution<unsigned int> mMutationDistribution;

    /// \brief    Mutation rate (translated to number of elements in population)
    unsigned int mMutationRate;

    /// \brief    Counts the number of iterations that the best fit doesn't change
    unsigned int mTrappedIterations;

    /// \brief    Maximum number of iterations that the algorithm can be trapped before vibration
    unsigned int mMaxTrappedIterations;

    /// \brief    Number of elements to replace when vibration is required
    unsigned int mVibratedElements;

    // Statistical data
    //-----------------
    unsigned int* mSnpHistogram;
    unsigned int* mGenotypeHistogram;

    double mFitnessMin;
    double mFitnessMax;
    double mFitnessMean;
    double mFitenssStd; 
    double mFitnessTotal;

    /// \brief    Normalized relates to FitnessMin
    double mFitnessTotalNormalized;
};

inline unsigned int CPopulation::ElitismGroupSize() const
{
    unsigned int elitismNum=(unsigned int)(CConfigParams::GetTheInstance().ElitismRate * mSize + 0.5);
    elitismNum = (elitismNum / 2) * 2;   // Make it a even
    return max(elitismNum,(unsigned int)2);
}
