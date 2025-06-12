
#include "SyntheticData.h"

#include <time.h>
#include <random>
#include <assert.h>

const uint32_t DATA_SET_SIZE=5000;
const uint32_t NUMBER_OF_SNP=500;

const uint32_t SNP_FIRST=0;
const uint32_t SNP_SECOND=150;

bool CSyntheticData::Load()
{
    bool result=true;
    
    mCaseData.resize(DATA_SET_SIZE);
    mControlData.resize(DATA_SET_SIZE);
    for (size_t i=0; i<DATA_SET_SIZE ; i++)
    {
        mCaseData[i].resize(NUMBER_OF_SNP);
        mControlData[i].resize(NUMBER_OF_SNP);
    }

    mTotalCase.resize(NUMBER_OF_SNP);
    mTotalControl.resize(NUMBER_OF_SNP);
    for (size_t i=0; i<NUMBER_OF_SNP ;i++)
    {
        memset(mTotalCase[i].mEvents, 0, NUMBER_OF_GENOTYPES);
        memset(mTotalControl[i].mEvents, 0, NUMBER_OF_GENOTYPES);
    }

    default_random_engine randomEngine((uint32_t)time(0));
    uniform_int_distribution<uint32_t>  snpIndex(0,NUMBER_OF_SNP-1);
    uniform_int_distribution<uint32_t>  dataIndex(0,DATA_SET_SIZE-1);

    unique_ptr<bool>  caseAllocatted(new bool[DATA_SET_SIZE]);
    unique_ptr<bool>  controlAllocatted(new bool[DATA_SET_SIZE]);

    CSnpData snpCase,snpControl;

    vector<pair<uint32_t, uint32_t>> barcode;

    for (uint32_t i=0; i<NUMBER_OF_SNP ;i++)
    {
        mTotalCase[i].mName =  mTotalControl[i].mName = snpCase.mName = snpControl.mName = ("rs" + std::to_string(i));

        if ((i == SNP_FIRST) || (i == SNP_SECOND))
        {
            continue;
//             snpCase.mEvents[0] = 580;
//             snpCase.mEvents[1] = 320;
//             snpCase.mEvents[2] = 100;
// 
//             snpControl.mEvents[0] = 700;
//             snpControl.mEvents[1] = 280;
//             snpControl.mEvents[2] = 20;
        }

        else
        {
            // Set Case distribution
            uniform_int_distribution<uint32_t>  genRatio(uint32_t(DATA_SET_SIZE*0.4f + 0.5f),uint32_t(DATA_SET_SIZE*0.6f + 0.5f));
            snpCase.mEvents[0] = genRatio(randomEngine);
            genRatio = uniform_int_distribution<uint32_t>(uint32_t((DATA_SET_SIZE-snpCase.mEvents[0])*0.4f + 0.5f) , uint32_t((DATA_SET_SIZE-snpCase.mEvents[0])*0.6f+ 0.5f));
            snpCase.mEvents[1] = genRatio(randomEngine);
            snpCase.mEvents[2] = DATA_SET_SIZE - snpCase.mEvents[0] - snpCase.mEvents[1];
            assert(snpCase.mEvents[2] < DATA_SET_SIZE);

            // Set the Control distribution
            genRatio = uniform_int_distribution<uint32_t>(snpCase.mEvents[0]-uint32_t(0.02*DATA_SET_SIZE),
                                                          snpCase.mEvents[0]+uint32_t(0.02*DATA_SET_SIZE));
            snpControl.mEvents[0] = genRatio(randomEngine);
            genRatio = uniform_int_distribution<uint32_t>(snpCase.mEvents[1]-uint32_t(0.02*DATA_SET_SIZE),
                                                          snpCase.mEvents[1]+uint32_t(0.02*DATA_SET_SIZE));
            snpControl.mEvents[1] = genRatio(randomEngine);
            snpControl.mEvents[2] = DATA_SET_SIZE - snpControl.mEvents[0] - snpControl.mEvents[1];
            assert(snpControl.mEvents[2] < DATA_SET_SIZE);
        }

        memset(caseAllocatted.get(), 0, DATA_SET_SIZE);
        memset(controlAllocatted.get(), 0, DATA_SET_SIZE);

        for (uint32_t j=0; j<NUMBER_OF_GENOTYPES ;j++ )
        {
            barcode.clear();
            barcode.push_back(make_pair(i,j));
            GenerateGenotype(snpCase.mEvents[j], barcode, caseAllocatted, mCaseData, mTotalCase, randomEngine, dataIndex );
            GenerateGenotype(snpControl.mEvents[j], barcode, controlAllocatted, mControlData, mTotalControl, randomEngine, dataIndex );
        }

        if (!(snpControl == mTotalControl[i]) || !(snpCase == mTotalCase[i]))
        {
            break;
        }
    }

    // Set selected SNPs
    memset(caseAllocatted.get(), 0, DATA_SET_SIZE);
    memset(controlAllocatted.get(), 0, DATA_SET_SIZE);

    //   100 75 80      410  200 20
    //   75 150 120      200  100 20
    //   80 120 200     20   20 10
    barcode.clear();
    barcode.push_back(make_pair(SNP_FIRST,0));
    barcode.push_back(make_pair(SNP_SECOND,0));
    GenerateGenotype(100, barcode, caseAllocatted, mCaseData, mTotalCase, randomEngine, dataIndex );
    GenerateGenotype(410, barcode, controlAllocatted, mControlData, mTotalControl, randomEngine, dataIndex );

    barcode.clear();
    barcode.push_back(make_pair(SNP_FIRST,1));
    barcode.push_back(make_pair(SNP_SECOND,0));
    GenerateGenotype(75, barcode, caseAllocatted, mCaseData, mTotalCase, randomEngine, dataIndex );
    GenerateGenotype(200, barcode, controlAllocatted, mControlData, mTotalControl, randomEngine, dataIndex );

    barcode.clear();
    barcode.push_back(make_pair(SNP_FIRST,2));
    barcode.push_back(make_pair(SNP_SECOND,0));
    GenerateGenotype(80, barcode, caseAllocatted, mCaseData, mTotalCase, randomEngine, dataIndex );
    GenerateGenotype(20, barcode, controlAllocatted, mControlData, mTotalControl, randomEngine, dataIndex );

    barcode.clear();
    barcode.push_back(make_pair(SNP_FIRST,0));
    barcode.push_back(make_pair(SNP_SECOND,1));
    GenerateGenotype(75, barcode, caseAllocatted, mCaseData, mTotalCase, randomEngine, dataIndex );
    GenerateGenotype(200, barcode, controlAllocatted, mControlData, mTotalControl, randomEngine, dataIndex );

    barcode.clear();
    barcode.push_back(make_pair(SNP_FIRST,1));
    barcode.push_back(make_pair(SNP_SECOND,1));
    GenerateGenotype(150, barcode, caseAllocatted, mCaseData, mTotalCase, randomEngine, dataIndex );
    GenerateGenotype(100, barcode, controlAllocatted, mControlData, mTotalControl, randomEngine, dataIndex );

    barcode.clear();
    barcode.push_back(make_pair(SNP_FIRST,2));
    barcode.push_back(make_pair(SNP_SECOND,1));
    GenerateGenotype(120, barcode, caseAllocatted, mCaseData, mTotalCase, randomEngine, dataIndex );
    GenerateGenotype(20, barcode, controlAllocatted, mControlData, mTotalControl, randomEngine, dataIndex );

    barcode.clear();
    barcode.push_back(make_pair(SNP_FIRST,0));
    barcode.push_back(make_pair(SNP_SECOND,2));
    GenerateGenotype(80, barcode, caseAllocatted, mCaseData, mTotalCase, randomEngine, dataIndex );
    GenerateGenotype(20, barcode, controlAllocatted, mControlData, mTotalControl, randomEngine, dataIndex );

    barcode.clear();
    barcode.push_back(make_pair(SNP_FIRST,1));
    barcode.push_back(make_pair(SNP_SECOND,2));
    GenerateGenotype(120, barcode, caseAllocatted, mCaseData, mTotalCase, randomEngine, dataIndex );
    GenerateGenotype(20, barcode, controlAllocatted, mControlData, mTotalControl, randomEngine, dataIndex );

    barcode.clear();
    barcode.push_back(make_pair(SNP_FIRST,2));
    barcode.push_back(make_pair(SNP_SECOND,2));
    GenerateGenotype(200, barcode, caseAllocatted, mCaseData, mTotalCase, randomEngine, dataIndex );
    GenerateGenotype(10, barcode, controlAllocatted, mControlData, mTotalControl, randomEngine, dataIndex );

    if (mTotalCase[SNP_FIRST].mEvents[0] != 255 || mTotalCase[SNP_FIRST].mEvents[1] != 345 || mTotalCase[SNP_FIRST].mEvents[2] != 400)
    {
        result = false;
    }
    if (mTotalCase[SNP_SECOND].mEvents[0] != 255 || mTotalCase[SNP_SECOND].mEvents[1] != 345 || mTotalCase[SNP_SECOND].mEvents[2] != 400)
    {
        result = false;
    }

    if (mTotalControl[SNP_FIRST].mEvents[0] != 630 || mTotalControl[SNP_FIRST].mEvents[1] != 320 || mTotalControl[SNP_FIRST].mEvents[2] != 50)
    {
        result = false;
    }
    if (mTotalControl[SNP_SECOND].mEvents[0] != 630 || mTotalControl[SNP_SECOND].mEvents[1] != 320 || mTotalControl[SNP_SECOND].mEvents[2] != 50)
    {
        result = false;
    }

    return result;
}

void CSyntheticData::GenerateGenotype(uint32_t iEventsNum, const vector<pair<uint32_t,uint32_t>>& iBarcode,  
                                      unique_ptr<bool>& ioAllocated, vector<vector<BYTE>>& ioData, vector<CSnpData>& iTotal,
                                      default_random_engine& iRandomEngine, uniform_int_distribution<uint32_t>& iDataDistribution )
{
	for (uint32_t j=0; j<iEventsNum ;j++)
	{
	    uint32_t index=iDataDistribution(iRandomEngine);
	    while  (ioAllocated.get()[index] == true)
	    {
	        if (++index == DATA_SET_SIZE)
	        {
	            index = 0;
	        }
	    }
	    ioAllocated.get()[index] = true;
        for(size_t i=0; i<iBarcode.size() ;i++)
        {
	        ioData[index][iBarcode[i].first] = iBarcode[i].second;

            iTotal[iBarcode[i].first].mEvents[iBarcode[i].second]++;
        }
	}
}
