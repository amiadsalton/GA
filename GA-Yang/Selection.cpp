
#include "Selection.h"
#include "RouletteWheelSelection.h"
#include "TournamentSelection.h"

 CSelection* CSelection::CreateConcrete()
 {
     CConfigParams::ESelectionAlgorithm type=CConfigParams::GetTheInstance().SelectionAlgorithm;
     if (type == eRoulette)
         return new CRouletteWheelSelection();
     return new CTournamentSelection();
 }
