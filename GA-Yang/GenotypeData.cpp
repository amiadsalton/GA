
#include "ConfigParams.h"

#include "GenotypeData.h"
#include "HapSampleData.h"
#include "SyntheticData.h"
#include "YangData.h"

 CGenotypeData* CGenotypeData::CreateConcrete()
 {
     CConfigParams::EGenotypeDataProvider  type=CConfigParams::GetTheInstance().Provider;
     if (type == CConfigParams::eHapMap)
     {
         return new CHapSampleData();
     }
     else if(type == CConfigParams::eSynthetic)
     {
         return new CSyntheticData();
     }
     else if (type == CConfigParams::eYangData)
     {
         return new CYangData(true);
     }
     else //if (type == CConfigParams::eCrohnData)
     {
         return new CYangData(false);
     }
 }

 void CGenotypeData::RemoveSNP(const vector<int>& iSnpsToRemove)
 {
     vector<int>::const_iterator iter = iSnpsToRemove.begin();
     while (iter != iSnpsToRemove.end() )
     {
         mTotalCase.erase(mTotalCase.begin() + *iter);
         mTotalControl.erase(mTotalControl.begin() + *iter);
         iter++;
     }

     for (size_t i=0; i<mCaseData.size() ; i++)
     {
         vector<int>::const_iterator iter = iSnpsToRemove.begin();
         while (iter != iSnpsToRemove.end() )
         {
             mCaseData[i].erase(mCaseData[i].begin() + *iter);
             mControlData[i].erase(mControlData[i].begin() + *iter);
             iter++;
         }
     }
 }
