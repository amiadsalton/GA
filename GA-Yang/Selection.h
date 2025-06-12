#pragma once

#include "ConfigParams.h"

#include <random>
#include <limits.h>

using namespace std;

class CPopulation;
class CChromosome;

////////////////////////////////////////////////////////////////////////
/// \brief    Interface for selection algorithm
///
/// Role: Represents GA selection algorithm 
///
/// Responsibilities: 
///      1. Interface for selecting one chromosome out of current generation
///      2. Factory of concrete algorithm
////////////////////////////////////////////////////////////////////////
class CSelection
{
public:

    enum ESelectionAlgorithm
    {
        eRoulette,
        eTournament
    };

    /// \brief    Access to Singleton object
    static CSelection* CreateConcrete();

    /// \brief    Initiate algorithm before sequence of selections
    /// \param[in] population   population to select from
    virtual void ResetSelection(const CPopulation& population) = 0;

    /// \brief    Load from file
    /// \param[in] population   population to select from
    /// \param[in] skip   skip this item 
    /// \return const CChromosome&     Selected chromosome 
    virtual unsigned int Select(const CPopulation& population,  default_random_engine& randomEngine, unsigned int skip=UINT_MAX) = 0;

protected:
    
};

