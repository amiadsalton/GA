#include "GeneticAlgorithm.h"
#include "ConfigParams.h"
#include "GenotypeData.h"
#include "PopulationYang.h"
#include "PopulationMooney.h"
#include "PostProcessingAlg.h"

#include <time.h>
#include <sstream>
#include <fstream>
#include <numeric>

#ifndef _MSC_VER
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#else
#include <boost/filesystem.hpp>
#endif



inline bool EqualFunctor(const pair<CChromosome, uint32_t>& first, const pair<CChromosome, uint32_t>& second)
{
    return (first == second);
}

bool CGeneticAlgorithm::Load()
{
    mData = CGenotypeData::CreateConcrete(); 
    if (mData == NULL)
    {
        return false;
    }

    bool ret = mData->Load();
    if (ret == false)
    {
        cout << "Failed to load data files" << endl;
        return false;
    }

    if (mData->CaseData().size() != mData->ControlData().size())
    {
        cout << "Mismatch between number of cases and number of controls" << endl;
        return false;
    }
    mNumberOfSNPs = mData->NumberOfSNPs();

    // Define mSelectedBarcode
    vector<unsigned int> selectedBarcodeSNPs;
    if (!FindSelectedBarcoseInDataSet(selectedBarcodeSNPs))
    {
        return false;
    }

    if (!selectedBarcodeSNPs.empty())
    {
        vector<BYTE> genotype(selectedBarcodeSNPs.size(), 2);

        CChromosome::Init((unsigned int)(selectedBarcodeSNPs.size()));
        if (CConfigParams::GetTheInstance().Algorithm == CConfigParams::eMooney)
        {
            mSelectedBarcode = new CChromosomeMooney((unsigned int)(selectedBarcodeSNPs.size()), &selectedBarcodeSNPs[0], &genotype[0]);
            dynamic_cast<CChromosomeMooney*>(mSelectedBarcode)->Fitness(mData->CaseData(), mData->ControlData());
        }
        else
        {
            mSelectedBarcode = new CChromosomeYang((unsigned int)(selectedBarcodeSNPs.size()), &selectedBarcodeSNPs[0],  &genotype[0]);
            dynamic_cast<CChromosomeYang*>(mSelectedBarcode)->Fitness(mData->CaseData(), mData->ControlData());
        }

        cout << "Selected barcode fitness: " << mSelectedBarcode->Fitness() << endl;
    }
    else
    {
        cout << "No Selected barcode" << endl;
    }

    return true;
}


bool CGeneticAlgorithm::Init(const string& iOutputFileName)
{
    bool ret=true;

    boost::system::error_code ec;
    boost::filesystem::path dir(iOutputFileName);
    if (boost::filesystem::exists(dir))
    {
        cout << "Directory " << dir.string() << " already exits" << endl;
        return false;
    }
    boost::filesystem::create_directory(iOutputFileName,ec);
    if (ec != boost::system::errc::success)
    {
        cout << "Failed to create directory " << dir.string() << endl;
        return false;
    }

    boost::filesystem::copy_file(CConfigParams::FILE_NAME,
    						 	 dir / CConfigParams::FILE_NAME,
								 boost::filesystem::copy_option::overwrite_if_exists,
								 ec);
    if (ec != boost::system::errc::success)
    {
    	cout << "Failed to copy config file: " << ec.message() << endl;
        return false;
    }

    // Open all files
    dir /= string("Run");
    mOutputFileName = dir.string();

    // Print SNPs information
    boost::filesystem::path tmp = boost::filesystem::path(mOutputFileName);
    tmp = tmp.parent_path();
    tmp = tmp.parent_path();
    tmp /= "SnpNames.csv";
    ofstream snpName(tmp.string());
    for (size_t i=0 ; i<mData->TotalCase().size() ; i++)
    {
        snpName << mData->TotalCase()[i].mName << ",";
        for (int j=0 ; j < CGenotypeData::NUMBER_OF_GENOTYPES ; j++ )
        {
            snpName <<  mData->TotalCase()[i].mEvents[j] << "," << mData->TotalControl()[i].mEvents[j] << ",";
        }
        snpName << endl;
    }
    snpName.close();

    if (mSelectedBarcode != nullptr)
    {
        DataSetAnalysis();
    }

    // Create post processing algorithm
    mPostProcessing = CPostProcessingAlg::CreateConcrete();

	mGenerationsCounter = 0;

#ifdef _STATISTICAL_DATA
    ostringstream mOutputFileName;
    mOutputFileName << CConfigParams::GetTheInstance().DirectoryName() << "/SNPSummary" << ".csv";
    ofstream mOutputFile(mOutputFileName.str());
    mOutputFile << "Id,name,0,1,2,name,0,1,2" << endl;
    for (unsigned int i=0; i < mNumberOfSNPs ;i++)
    {
        mOutputFile << i+1 << "," << mData->TotalCase()[i].mName;
        for (int j=0; j<3 ; j++)
        {
            mOutputFile  << "," << mData->TotalCase()[i].mEvents[j];
        }
        mOutputFile << "," << mData->TotalControl()[i].mName;
        for (int j=0; j<3 ; j++)
        {
            mOutputFile  << "," << mData->TotalControl()[i].mEvents[j];
        }  
        mOutputFile << endl;
    }
    mOutputFile.close();
#endif

    return true;
}

bool CGeneticAlgorithm::RunAllOrders()
{
    bool result(true);

    mNumberOfIterations = CConfigParams::GetTheInstance().NumOfIterations;

    unsigned int length = CConfigParams::GetTheInstance().ChromosomeLength;
    if  (length > 0)
    {
        // The order is pre-defined, need to run the algorithm only on for this order. 
        CChromosome::Init(length);


        bool result = Run(length, mOutputFileName +  "_" + to_string(length));

        mPostProcessing->Execute(*mData, length);
        mPostProcessing->Print(*mData, mOutputFileName, length);

        mUserInput = eNone;
        
		cout << "Number of executed generations: " << mGenerationsCounter << endl 
			 << "Number of fitness calculations: " << CChromosome::FitnessCalculationCounter() << endl;

        return result;
    }


    // Run the algorithm for varies chromosomes length
    for (unsigned int i=2; i<=CConfigParams::GetTheInstance().NumberOfOrders ;i++)
    {
        cout << "\n=============\nRun for length: " << i  << "\n=============\n";

        CChromosome::Init(i);
 
        //mNumberOfIterations = (unsigned int) (CConfigParams::GetTheInstance().NumOfIterations * (double(i) / CConfigParams::GetTheInstance().NumberOfRuns)) ;

        if (!Run(i, mOutputFileName +  "_" + to_string(i)))
        {
            result = false;
            break;
        }
        else
        {
            cout <<  "\nNumber of Fitness calculations: " << CChromosome::FitnessCalculationCounter() << endl;
            if (!mPostProcessing->Execute(*mData, i))
            {
                mPostProcessing->Print(*mData, mOutputFileName, i);
                cout <<  "No more solutions after Run: " << i;
                break;
            }
            mPostProcessing->Print(*mData, mOutputFileName, i);

            // Update data set: Remove all SNPs that are not relevant any more
            mPostProcessing->UpdateGenotypeData(*mData);

            // Update the number  of SNPs that should be checked
            mNumberOfSNPs = mData->NumberOfSNPs();

            // Update the SNPs indexes of mSelectedBarcode
            if (mSelectedBarcode != nullptr)
            {
	            vector<unsigned int> selectedBarcodeSNPs;
	            if (!FindSelectedBarcoseInDataSet(selectedBarcodeSNPs))
	            {
	                result = false;
	                break;
	            }
	            CChromosome* tmp = mSelectedBarcode;
	            mSelectedBarcode = new CChromosome(tmp->Length(), &(selectedBarcodeSNPs[0]), tmp->Genotype(), tmp->Fitness(), tmp->AdditionalFitnessData());
	            delete tmp;
            }
        }

        if (mUserInput == eExit)
        {
            break;
        }

        mUserInput = eNone;
    }

    mUserInput = eNone;

	cout << "Number of executed generations: " << mGenerationsCounter << endl;

    return result;
}

bool CGeneticAlgorithm::Run(unsigned int iLength, const string& iOutputFileName)
{
    mOutputFile.open(iOutputFileName + ".csv");
    mOutputFile << "TimeElapsed, Max, Min, Mean, Std, Number of items.,";
    if (mSelectedBarcode != nullptr)
    {
        mOutputFile << " Selected: " << mSelectedBarcode->Fitness();
    }
    mOutputFile <<  ",Top"  << endl;

    mBestFile.open(iOutputFileName + "_best.txt");

    for (unsigned int i=0; i<CConfigParams::GetTheInstance().NumberOfExecutions ;i++)
    {
        cout << "\nExecution no. " << i  << endl;
        mCurrentExecution = i;
        if (Execution() == false)
        {
            return false;
        }

        // Analyze user input
        if (mUserInput == eExit || mUserInput == eSkipExecution)
        {
            break;
        }
    }

    mOutputFile.close();
    mBestFile.close();

    if (mUserInput == eSkipExecution)
    {
        mUserInput = eNone;
    }

    return true;
}


bool CGeneticAlgorithm::Execution()
{

    // Build population
    default_random_engine randomEngine((unsigned int)time(0));

    CPopulation*  population=NULL;
    if (CConfigParams::GetTheInstance().Algorithm == CConfigParams::eMooney)
    {
        population = new CPopulationMooney(randomEngine, *mData, mNumberOfSNPs);
    }
    else
    {
        population = new CPopulationYang(randomEngine, *mData, mNumberOfSNPs);

    }
     
    if ( /*(population->Size() < mData->TotalCase().size() && population->Size() < 100) || */ population->Size() < 2)
    { 
        cout << "Population size is too small " << " Population size: " << population->Size() << " Number of SNPs: " << mData->TotalCase().size() << endl;
        return false;
    }

    population->Init();

    time_t startTime;
    time(&startTime);
 
    vector<CChromosome> lastElitism;

    unsigned int stuckCriteria = 0;

    unsigned int elitismNum = population->ElitismGroupSize();

    // Loop until number of iteration (or stuck in local extreme) 
    for (unsigned int i=0; i<mNumberOfIterations; i++, mGenerationsCounter++)
    {
        population->NextGeneration(); 

        vector<const CChromosome*>  differentChromosomes;
        population->DifferentChromosomes(differentChromosomes, -1);

        PrintResults(*population, differentChromosomes, startTime, elitismNum);

        if (EqualElitism(lastElitism,differentChromosomes) == elitismNum)
        {
            if (++stuckCriteria == CConfigParams::GetTheInstance().ExecutionIsStcuk)
            {
                cout << "\nExecution is stuck after " << i << " Generations " << endl;
                break;
            }
        }
        else
        { 
            stuckCriteria = 0;
            lastElitism.clear();
            for (unsigned int i = 0; i<elitismNum ; i++)
            {
                lastElitism.push_back(*differentChromosomes[i]);
            }
        }

        // Prepare result file
        mOutputFile.flush();
        mBestFile.flush();

        cout << i << " ";

        // Analyze user input
        while (mUserInput == ePause)
        {
#ifdef _MSC_VER
        	Sleep(1000);
#else
        	sleep(1);
#endif
        }
        if (mUserInput == eExit || mUserInput == eSkipGeneration || mUserInput == eSkipExecution)
        {
            break;
        }

        mUserInput = eNone;
    }

    mPostProcessing->Update(lastElitism);

    delete population;

    if (mUserInput == eSkipGeneration)
    {
        mUserInput = eNone;
    }

    return true;
}


void CGeneticAlgorithm::PrintResults(const CPopulation& population, const vector<const CChromosome*>& iDifferentChromosmes,  time_t startTime, uint32_t iElitismNum) 
{
    // General information
    // ------------------
    // Time elapsed
    time_t currTime;
    time(&currTime);
    mOutputFile << currTime - startTime << ",";

    mOutputFile.precision(12);
    mOutputFile << population.FitnessMax() << "," << population.FitnessMin() << ","
               /*<< population.FitnessTotal() << "," */<< population.FitnessMean() << "," 
               << population.FitenssStd() << "," ;

    mOutputFile << iDifferentChromosmes.size() << ",";

    // Selected Barcode:
    // ------------------
    // Find the index of selected SNP
    if (mSelectedBarcode != nullptr)
    {
	    int barcodeIndex=-1;
	    double fitness=0;
	    for (int i=0; i<iDifferentChromosmes.size() ;i++)
	    {
	        if (*iDifferentChromosmes[i] == *mSelectedBarcode)
	        {
	            barcodeIndex = i;
	            fitness = iDifferentChromosmes[i]->Fitness();
	            break;
	        }
	    }
	
	    mOutputFile << barcodeIndex;
	    mOutputFile << CChromosomeSerializer(*mSelectedBarcode, *mData);
	    mOutputFile << ",";
    }

    // 3 Best Chromosomes [1, length-1]:
    // ------------------
    for (unsigned int j=0; j < min(iDifferentChromosmes.size(),size_t(3)) ;j++)
    {
        if (j>0)
        {
            mOutputFile << ",";
        }
        mOutputFile <<  CChromosomeSerializer(*iDifferentChromosmes[j], *mData);
    } 
    mOutputFile << endl;

    // Best Chromosomes:
    // ------------------
    for (unsigned int j=0; j < min(size_t(iElitismNum),iDifferentChromosmes.size()) ;j++)
    {
        if (j>0)
        {
            mBestFile << ",";
        }
        mBestFile << CChromosomeSerializer(*iDifferentChromosmes[j], *mData);
    } 
    mBestFile << endl;
}


int CGeneticAlgorithm::NumberOfCreatedItems(vector<int>& iCurrentGenotypes, const vector<vector<BYTE>>& iDataSet) const
{
    int numberOfItems(0);
    for (size_t i=0; i<iDataSet.size(); i++)
    {
        bool match(true);
        for (unsigned int j=0; j<iCurrentGenotypes.size(); j++)
        {
            if (iDataSet[i][mSelectedBarcode->SNP()[j]] != iCurrentGenotypes[j])
            {
                match = false;
                break;
            }
        }
        if (match)
        {
            numberOfItems++;
        }
    }
    return  numberOfItems;
}

 void CGeneticAlgorithm::OneCombination(ofstream& oOutputFile, ifstream& iInputProbFile,
                           vector<int>& oDiff, vector<int>& oDiffRatio, 
                           const vector<vector<BYTE>>& iDataSet,
                           vector<int>& ioCurrentGenotypes, unsigned int iSNP) const
{
    for (int i=0; i<3; i++)
    {
        ioCurrentGenotypes[iSNP] = i;
        if (iSNP > 0)
        {
            OneCombination(oOutputFile,iInputProbFile, oDiff, oDiffRatio, iDataSet, ioCurrentGenotypes, iSNP - 1);
        }
        else
        {
            double tmp;
            iInputProbFile >> tmp;
            int tmpRequested=int(tmp*iDataSet.size()+0.5);
            int tmpCreated(NumberOfCreatedItems(ioCurrentGenotypes, iDataSet));
            oDiff.push_back(tmpRequested - tmpCreated);
            oDiffRatio.push_back(int((double(tmpRequested - tmpCreated) / tmpRequested)*100.0 + 0.5));
            for (unsigned int j=0 ; j<mSelectedBarcode->Length() ; j++)
            {
                oOutputFile << ioCurrentGenotypes[j] << ",";
            }
            oOutputFile << tmpRequested << "," << tmpCreated << "," << (tmpRequested - tmpCreated) << "," << oDiffRatio[oDiffRatio.size()-1] << endl; 
        }
    }

}

void CGeneticAlgorithm::OneDataSetAnalysis(const vector<vector<BYTE>>& iDataSet, const string& iDataType) const
{
    boost::filesystem::path fileName(mOutputFileName);
    fileName = fileName.parent_path();
    fileName = fileName.parent_path();
    fileName /= string("disease_model_" + iDataType + ".txt");
	
	ifstream inputProbFile(fileName.string());
	if (inputProbFile.bad())
	{
	    cout << "\nCan't find input probabilities files " << endl;
	    return;
	}
	
	ofstream analysisFile(mOutputFileName + "_compare_" + iDataType + ".csv");
	if (analysisFile.bad())
	{
	    cout << "\nCan't generate comparison file " << endl;
	    return;
	}
    for (unsigned int  i=0; i<mSelectedBarcode->Length() ;i++)
    {
        analysisFile << mSelectedBarcode->SNP()[i] <<",";
    }
	analysisFile << "Requested,Created,Diff,Diff%" << endl;

    int numberOfSNPs;
    inputProbFile >> numberOfSNPs;
    if (numberOfSNPs != mSelectedBarcode->Length())
    {
        cout << "\nMismatch between input data and defined chromosome in config file " << endl;
        return;
    }
    string line;
	for (unsigned int i=0; i<mSelectedBarcode->Length(); i++)
	{
	    inputProbFile >> line;
	}
	
	vector<int> diff, diffRatio;
	
	vector<int>  currentGenotypes(mSelectedBarcode->Length(),0);
	OneCombination(analysisFile, inputProbFile, diff, diffRatio, iDataSet, currentGenotypes,mSelectedBarcode->Length()-1);

    int sum = std::accumulate(diff.begin(), diff.end(), 0);

    analysisFile << "\nNUM,SUM,MAX,MIN" << endl;
    analysisFile << diff.size() - std::count(diff.begin(), diff.end(), 0) << "," << sum << "," << *std::max_element(diff.begin(), diff.end()) << "," << *std::min_element(diff.begin(), diff.end()) << endl;
    analysisFile << diffRatio.size() - std::count(diffRatio.begin(), diffRatio.end(), 0) << "," << sum << "," << *std::max_element(diffRatio.begin(), diffRatio.end()) << "," << *std::min_element(diffRatio.begin(), diffRatio.end()) << endl;

}

void CGeneticAlgorithm::DataSetAnalysis() const
{
    OneDataSetAnalysis(mData->CaseData(), "case");
    OneDataSetAnalysis(mData->ControlData(), "anticase");
}


bool CGeneticAlgorithm::FindSelectedBarcoseInDataSet(vector<unsigned int>& selectedBarcodeSNPs)
{
    const vector<string>& selectedBarcode=CConfigParams::GetTheInstance().SelectedBarcode;

    if (selectedBarcode.empty())
    {
        return true;
    }

    for (size_t i=0; i<selectedBarcode.size() ;i++)
    {
        for (size_t j=0; j < mData->TotalCase().size() ;j++)
        {
            if (mData->TotalCase()[j].mName == selectedBarcode[i]) 
            {
                selectedBarcodeSNPs.push_back((unsigned int)(j));
                break;
            }
        }
    }

    if (selectedBarcodeSNPs.size() != selectedBarcode.size())
    {
        cout << "\nMismatch between input data chromosome size " << selectedBarcodeSNPs.size() << " and defined chromosome in config file " << selectedBarcode.size() << endl;
        for (auto elem : selectedBarcodeSNPs)
        {
            cout << mData->TotalCase()[elem].mName << ",";
        }
        cout << endl;
        return false;
    }

    sort(selectedBarcodeSNPs.begin(), selectedBarcodeSNPs.end());

    return true;
}