// SyntheticGenerator.cpp : Defines the entry point for the console application.
//


#include <sstream>
#include <time.h>
#include <random>
#include <vector>
#include <boost/filesystem.hpp>

#include "../GA-Yang/GenotypeData.h"
#include "ConfigParams.h"

using namespace std;

typedef unsigned char BYTE;


/// \brief    Data of each item in Case population. Each item has genotype definition for each SNP. 
vector<vector<BYTE>> mCaseData;

/// \brief    Data of each item in Control population. Each item has genotype definition for each SNP. 
vector<vector<BYTE>> mControlData;

/// \brief    Genotype distributions of all analyzed SNPs in case dataset
vector<CGenotypeData::CSnpData> mTotalCase;

/// \brief    Genotype distributions of all analyzed SNPs in control dataset
vector<CGenotypeData::CSnpData> mTotalControl;

void GenerateGenotype(unsigned int iEventsNum, const vector<pair<unsigned int,unsigned int>>& iBarcode,  
    unique_ptr<bool>& ioAllocated, vector<vector<BYTE>>& ioData, vector<CGenotypeData::CSnpData>& iTotal,
    default_random_engine& iRandomEngine, uniform_int_distribution<unsigned int>& iDataDistribution )
{
    const CConfigParams& config = CConfigParams::GetTheInstance();   

    for (unsigned int j=0; j<iEventsNum ;j++)
    {
        unsigned int index=iDataDistribution(iRandomEngine);
        while  (ioAllocated.get()[index] == true)
        {
            if (++index == config.PopulationSize())
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
static int counter=0;
void Recursion(size_t size, int* distribution,unsigned int& ioTotal)
{
    const CConfigParams& config = CConfigParams::GetTheInstance();   
    if (size == config.SelectedSNPs().size()-1)
    {
        float toatlDist=1.0f;
        for (size_t i=0 ; i<config.SelectedSNPs().size()-1 ; i++ )
        {
            toatlDist *= config.CaseDistribution()[distribution[i]];
        }

        for (int i=CGenotypeData::NUMBER_OF_GENOTYPES-1; i>=0 ;i--)
        {
            size_t numOfEvents = size_t(config.PopulationSize() * toatlDist * config.CaseDistribution()[i] + 0.5f);
            for (size_t j=0; j < numOfEvents && ioTotal < config.PopulationSize() ; j++, ioTotal++)
            {
                for (size_t k=0 ; k<config.SelectedSNPs().size()-1 ; k++)
                {
                    mCaseData[ioTotal][config.SelectedSNPs()[k]] = distribution[k]; 
                }
                mCaseData[ioTotal][config.SelectedSNPs()[config.SelectedSNPs().size()-1]] = i;
                counter++;
            }
        }     
    }
    else
    {
        for (int i=CGenotypeData::NUMBER_OF_GENOTYPES-1; i>=0 ;i--)
        {
            distribution[size] = i;
            Recursion(size+1, distribution, ioTotal);
        }
    }
}

void GenerateSelected()
{
    const CConfigParams& config = CConfigParams::GetTheInstance();   

    int* distribution = new int[config.SelectedSNPs().size()];
    memset(distribution, 0, sizeof(int) * config.SelectedSNPs().size());

    unsigned int total=0;
    Recursion(0, distribution, total);
}


bool Load(string dirName)
{
    bool result=true;

    const CConfigParams& config = CConfigParams::GetTheInstance();   

    mCaseData.resize(config.PopulationSize());
    mControlData.resize(config.PopulationSize());
    for (size_t i=0; i<config.PopulationSize() ; i++)
    {
        mCaseData[i].resize(config.NumberOfSNPs());
        mControlData[i].resize(config.NumberOfSNPs());
    }

    mTotalCase.resize(config.NumberOfSNPs());
    mTotalControl.resize(config.NumberOfSNPs());
    for (size_t i=0; i<config.NumberOfSNPs() ;i++)
    {
        memset(mTotalCase[i].mEvents, 0, CGenotypeData::NUMBER_OF_GENOTYPES);
        memset(mTotalControl[i].mEvents, 0, CGenotypeData::NUMBER_OF_GENOTYPES);
    }

    default_random_engine randomEngine((unsigned int)time(0));
    uniform_int_distribution<unsigned int>  snpIndex(0,config.NumberOfSNPs()-1);
    uniform_int_distribution<unsigned int>  dataIndex(0,config.PopulationSize()-1);

    unique_ptr<bool>  caseAllocatted(new bool[config.PopulationSize()]);
    unique_ptr<bool>  controlAllocatted(new bool[config.PopulationSize()]);

    CGenotypeData::CSnpData snpCase,snpControl;

    vector<pair<unsigned int, unsigned int>> barcode;

    for (unsigned int i=0; i<config.NumberOfSNPs() ;i++)
    {
        ostringstream  tmp;
        tmp << "rs" << i;
        mTotalCase[i].mName =  mTotalControl[i].mName = snpCase.mName = snpControl.mName = tmp.str();

        size_t j=0;
        uniform_int_distribution<unsigned int>  genRatio(unsigned int(config.PopulationSize()*0.4f + 0.5f),unsigned int(config.PopulationSize()*0.6f + 0.5f));
        for (; j < config.SelectedSNPs().size() ;j++) 
        {
            if (config.SelectedSNPs()[j] == i)
            {
                // Set selected SNP
                snpCase.mEvents[0] = unsigned int(config.CaseDistribution()[0] * config.PopulationSize() + 0.5f); //genRatio(randomEngine);
                //genRatio = uniform_int_distribution<unsigned int>(unsigned int((config.PopulationSize()-snpCase.mEvents[0])*0.4f + 0.5f) , unsigned int((config.PopulationSize()-snpCase.mEvents[0])*0.6f+ 0.5f));
                snpCase.mEvents[1] = unsigned int(config.CaseDistribution()[1] * config.PopulationSize() + 0.5f); // genRatio(randomEngine);
                snpCase.mEvents[2] = config.PopulationSize() - snpCase.mEvents[0] - snpCase.mEvents[1];
                _ASSERT(snpCase.mEvents[2] < config.PopulationSize());
                break;
            }
        }

        if (j == config.SelectedSNPs().size())
        {
            // Set Case distribution (default)
            snpCase.mEvents[0] = unsigned int(config.DefaultDistribution()[0] * config.PopulationSize() + 0.5f); //genRatio(randomEngine);
            //genRatio = uniform_int_distribution<unsigned int>(unsigned int((config.PopulationSize()-snpCase.mEvents[0])*0.4f + 0.5f) , unsigned int((config.PopulationSize()-snpCase.mEvents[0])*0.6f+ 0.5f));
            snpCase.mEvents[1] = unsigned int(config.DefaultDistribution()[1] * config.PopulationSize() + 0.5f); // genRatio(randomEngine);
            snpCase.mEvents[2] = config.PopulationSize() - snpCase.mEvents[0] - snpCase.mEvents[1];
            _ASSERT(snpCase.mEvents[2] < config.PopulationSize());
        }
        else
        {
            continue;
        }

        // Set the Control distribution
        //genRatio = uniform_int_distribution<unsigned int>(snpCase.mEvents[0]-unsigned int(0.02*config.PopulationSize()),
        //                                                  snpCase.mEvents[0]+unsigned int(0.02*config.PopulationSize()));
        snpControl.mEvents[0] = unsigned int(config.DefaultDistribution()[0] * config.PopulationSize() + 0.5f); // genRatio(randomEngine);
        //genRatio = uniform_int_distribution<unsigned int>(snpCase.mEvents[1]-unsigned int(0.02*config.PopulationSize()),
        //                                                  snpCase.mEvents[1]+unsigned int(0.02*config.PopulationSize()));
        snpControl.mEvents[1] = unsigned int(config.DefaultDistribution()[1] * config.PopulationSize() + 0.5f); //genRatio(randomEngine);
        snpControl.mEvents[2] =  config.PopulationSize() - snpControl.mEvents[0] - snpControl.mEvents[1];
        _ASSERT(snpControl.mEvents[2] < config.PopulationSize());

        memset(caseAllocatted.get(), 0, config.PopulationSize());
        memset(controlAllocatted.get(), 0, config.PopulationSize());

        for (unsigned int j=0; j<CGenotypeData::NUMBER_OF_GENOTYPES ;j++ )
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

    GenerateSelected();

    return result;
}


int main(int argc, char* argv[])
{
    const CConfigParams& config = CConfigParams::GetTheInstance();   

    string dirName("..\\");
    dirName +=config.DirectoryName();

    if (Load(dirName) == false)
    {
        cout << "Failed to generate data" << endl;
        char ch;
        cin >> ch;
        exit(-1);
    }

    if (!boost::filesystem::exists(dirName))
    {
	        boost::filesystem::create_directory(dirName);
    }
    if (boost::filesystem::exists(dirName + "\\Summary.txt"))
    {
        boost::filesystem::remove(dirName + "\\Summary.txt");
    }
    boost::filesystem::copy_file(CConfigParams::GetTheInstance().FILE_NAME, dirName + "\\Summary.txt");
  

    ofstream  caseFile(dirName + "\\case_genotypes.dat");
    ofstream  controlFile(dirName + "\\anticase_genotypes.dat");
    for (unsigned int i=0; i<config.NumberOfSNPs() ;i++)
    {
        caseFile << "1\t" << mTotalCase[i].mName << "\t2222\t0.00000\t";
        controlFile << "1\t" << mTotalCase[i].mName << "\t2222\t0.00000\t";
        for (unsigned j=0; j<mCaseData.size()-1 ; j++)
        {
            caseFile << int(mCaseData[j][i]) << "\t";
            controlFile << int(mControlData[j][i]) << "\t";
        }
        caseFile << int(mCaseData[mCaseData.size()-1][i]);
        controlFile << int(mControlData[mCaseData.size()-1][i]);
        if (i < config.NumberOfSNPs() - 1)
        {
            caseFile << endl;
            controlFile << endl;
        }
    }
	return 0;
}

