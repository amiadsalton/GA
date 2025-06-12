
#include "RouletteWheelSelection.h"
#include "Population.h"
#include <omp.h>

void CRouletteWheelSelection::ResetSelection(const CPopulation& population) 
{
    mRandomGenerator = uniform_real_distribution<double>(0,population.FitnessTotalNormalized());
}

 unsigned int CRouletteWheelSelection::Select(const CPopulation& population, default_random_engine& randomEngine, unsigned int skip)
{

    double wheelLocation;
    #pragma omp critical  
    {
        wheelLocation = mRandomGenerator(randomEngine);
    }

    int index = population.Size() - 1;
    if (index == skip)
    {
        index--;
    }
    double currSum = population[index].Fitness() - population.FitnessMin();
    while ( (currSum < wheelLocation) && (index > 0) )
    {
        index--;
        if (index == skip)
        {
            index--;
            if (index < 0)
            {
                break;
            }
        }
        currSum += population[index].Fitness() - population.FitnessMin();
    }

    return index;
}

