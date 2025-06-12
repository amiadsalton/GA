
#include "TournamentSelection.h"
#include "Population.h"

void CTournamentSelection::ResetSelection(const CPopulation& population) 
{
    mRandomGenerator = uniform_int_distribution<unsigned int>(0 , population.Size()-1);
}

 unsigned int CTournamentSelection::Select(const CPopulation& population, default_random_engine& randomEngine, unsigned int skip)
{
    unsigned int selected= mRandomGenerator(randomEngine);
    while (selected == skip)
    {
        selected = mRandomGenerator(randomEngine);
    }


    unsigned int mate = mRandomGenerator(randomEngine);
    while (mate == skip)
    {
        mate = mRandomGenerator(randomEngine);
    }

    if (population[mate] > population[selected])
    {
        selected = mate;
    }


    return selected;
}

