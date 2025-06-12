#include "PostProcessingAlg.h"

#include "ClusteringAlg.h"
#include "ContinuityAlg.h"

CPostProcessingAlg* CPostProcessingAlg::CreateConcrete()
{
    CConfigParams::EPostProcessingAlgorithm type=CConfigParams::GetTheInstance().PostProcessingAlg;
    if (type == CConfigParams::eContinuity)
    {
        return new CContinuityAlg();
    }
    return new CClusteringAlg(type);
}

