#pragma once

#include "ChromosomeMooney.h"
#include "Population.h"
#include "GenotypeData.h"


//#define YANG 1
#ifdef YANG
//typedef CChromosomeYang CChromosomeMooney;
#else
//typedef CChromosomeMooney CChromosomeMooney;
#endif


using namespace std;

class CSelection;

////////////////////////////////////////////////////////////////////////
/// \brief    CChromosome is one element in GA population
///
/// Role: Represents one generation of Yang chromosomes 
///
/// Responsibilities: 
///      1. Efficient runt time implementation of CPopulation methods to the specific fitness function (Mooney’s)
///
////////////////////////////////////////////////////////////////////////
class CPopulationMooney : public CPopulation
{

public:

    /// \brief    Empty constructor.
    /// \param[in] randomEngine    Pseudo random generator which is used to for initialization and evolution.
    /// \param[in] genotypeData    input dataset.
    /// \param[in] numberOfSNPs    Number of analyzed SNPs.
    CPopulationMooney(default_random_engine& randomEngine, CGenotypeData& genotypeData, unsigned int numberOfSNPs);

    virtual ~CPopulationMooney();

    /// \brief    Created the initial Population
    virtual void Init() override;

    /// \brief    Prepare the next generation from current generation
    virtual void NextGeneration() override;


    virtual  const CChromosome& operator[](size_t index) const override { return mCurrentGeneration[index]; }

    virtual void DifferentChromosomes(vector<const CChromosome*>& oChromosomes, int ioNumberOfSolutions = -1) const override;

private:

    /// \brief     All elements of current generation (mSize = Population size)      
    CChromosomeMooney* mCurrentGeneration;

    /// \brief     Temporary storage of next generation      
    CChromosomeMooney* mNextGeneration;

    /// \brief    Calculate the fitness of current generation and sort 
    void Fiteness();

    /// \brief    Elitism - select the best elements from current generation and transfer them to next generation 
    /// \return unsigned int    number of selected elements by Elitism 
    unsigned int  Elitism();

    /// \brief    Vibrate - if the algorithm is trap in local extreme point then it replaces some of the low fitness chromosomes 
    /// \param[in] lastBestFitness    Best fitness of previous iteration.
    void Vibrate(double lastBestFitness);

    /// \brief    brief check if input items already exist in current generation
    /// \param[in] first    first chromosome to check.
    /// \param[out] currentIndex    index in current generation.
    /// \return bool    true if there is duplication.
    bool CheckForDuplication(const CChromosomeMooney& input, unsigned int currentIndex);


};