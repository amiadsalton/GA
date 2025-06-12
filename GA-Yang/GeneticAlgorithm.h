#pragma once

#include <string>
#include <vector>
#include <map>
#include <stdint.h>

#include "Chromosome.h"

using namespace std;

class CGenotypeData;
class CPopulation;
class CPostProcessingAlg;


////////////////////////////////////////////////////////////////////////
/// \brief    Execution of Genetic Algorithm for GWAS
///
/// Role: Controls the execution of the Genetic Algorithm and reports the results  
///
/// Responsibilities: 
///      1. Loads the input dataset
///      2. Run the algorithm for all requested order 
///      3. Run the algorithm for one order
///      4. Knows how many generations should be executed
///      5. Execute all Generations
///      6. Generate the next generation population
///      7. Knows how to print the algorithm results. 
/// 
////////////////////////////////////////////////////////////////////////
class CGeneticAlgorithm
{
public:

    enum EUserInput 
    {
        eNone = 0,
        eExit,
        ePause,
        eResume,
        eSkipExecution,
        eSkipGeneration,
        eLast,
    };
    
    CGeneticAlgorithm() : mSelectedBarcode(nullptr) , mGenerationsCounter(0) {}
    ~CGeneticAlgorithm()  { delete mSelectedBarcode; }

    bool Load();

    /// \brief    Initialize algorithm.
    /// \return bool    false in case of failure.
    bool Init(const string& iOutputFileName);

    /// \brief    Execute genetic algorithm.
    /// 
    /// It runs the algorithm for variant SNPs combination length (order), and then it analyzes the results in order to find the max combination length that
    /// has correlation with the disease. 
    /// It calls for the function Run with combination length that varies from 2 to user defined maximum order.
    /// \return bool    false in case of failure.
    bool RunAllOrders();

    EUserInput GetUserInput() const { return mUserInput; }
    void SetUserInput(EUserInput iRequest) { mUserInput = iRequest; }

private:
    
    /// \brief    Number of iteration in one GA run
    unsigned int mNumberOfIterations;

    /// \brief    Number of analyzed SNPs
    unsigned int mNumberOfSNPs;

    /// \brief    Reference to input genotype data of cases and controls
    CGenotypeData*  mData;

    /// \brief    The barcode that we are looking for
    CChromosome* mSelectedBarcode;

    CPostProcessingAlg* mPostProcessing;

    unsigned int mCurrentExecution;

	/// \brief    Counts the total number of Generations
	unsigned int mGenerationsCounter;

    string mOutputFileName;

    ofstream mOutputFile;
    ofstream mHistFile;
    ofstream mBestFile;
    ofstream mSummaryFile;

    EUserInput mUserInput; 

    /// \brief    Execute genetic algorithm.
    ///
    /// For input SNPs combination length, this function looks for the SNPs that have best correlation with disease. 
    /// It call for Execution configurable number of times
    /// \param[in] iLength    GA chromosome length which is the length of SNPs combination which the algorithm is looking for 
    /// \param[in] iOutputFileName    Filename which is used to record algorithm results
    /// \return bool    false in case of failure.
    bool Run(unsigned int iLength, const string& iOutputFileName);

    /// \brief    One execution of Genetic Algorithm.
    /// 
    /// This is one round of GA algorithm: evolution of chromosomes population for configurable number of generations 
    /// \return void    description.
    bool Execution();

    /// \brief    Execute genetic algorithm.
    /// \param[in] population    reference to current generation.
    /// \param[in] outputFile    reference to a file which algorithm results should be written into.
    /// \param[in] startTime    start algorithm time
    void PrintResults(const CPopulation& population, const vector<const CChromosome*>& iDifferentChromosmes,  time_t startTime, uint32_t iElitismNum);

    unsigned int EqualElitism(const vector<CChromosome>& ioLastGeneration, const vector<const CChromosome*>& iCurrGeneration) const;

    // Compares the generated data set with the requested probabilities and print results into csv file. 
    void DataSetAnalysis() const;
    void OneDataSetAnalysis(const vector<vector<BYTE>>& iDataSet, const string& iDataType) const;
    void OneCombination(ofstream& oOutputFile, 
                        ifstream& iInputProbFile,
                        vector<int>& oDiff, 
                        vector<int>& oDiffRatio,
                        const vector<vector<BYTE>>& iDataSet,
                        vector<int>& ioCurrentGenotypes,
                        unsigned int iSNP) const;
    int NumberOfCreatedItems(vector<int>& iCurrentGenotypes, const vector<vector<BYTE>>& iDataSet) const;

    bool FindSelectedBarcoseInDataSet(vector<unsigned int>& selectedBarcodeSNPs);
};

inline unsigned int CGeneticAlgorithm::EqualElitism(const vector<CChromosome>& ioLastGeneration, const vector<const CChromosome*>& iCurrGeneration) const
{
    unsigned int result = 0;

    unsigned int size =(unsigned int)(min(ioLastGeneration.size(), iCurrGeneration.size()));

    if (size > 0)
    {
        for ( ; result<size ; result++)
        {
            if (ioLastGeneration[result].Fitness() < iCurrGeneration[result]->Fitness())
            {
                break;
            }
        }
    }

    return result;
}
