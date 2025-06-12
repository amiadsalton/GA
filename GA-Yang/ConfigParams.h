#pragma once

#include <memory>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////
/// \brief    Configurable parameters
///
/// Role: Represents a collection of all configurable parameters
///
/// Responsibilities: 
///      1. Singleton
///      2. Access to all configurable parameters
///      3. Knows how to load this object from configuration file
///
////////////////////////////////////////////////////////////////////////
class CConfigParams
{
public:

    /// \brief    Genotype data provider
    enum EGenotypeDataProvider { eHapMap = 0, eSynthetic, eYangData, eCrohnData };

    /// \brief    Selection algorithm
    enum ESelectionAlgorithm { eRouletteWheel = 0 };

    /// \brief    Selection algorithm
    enum EPostProcessingAlgorithm { eClusteringFitenss = 0, eClusteringPosition, eClusteringSum , eContinuity, eNumberOfPostProcessingAlgorithms };

    /// \brief    GA algorithm
    enum EAlgorithm { eYang = 0, eMooney };

    /// \brief    Access to Singleton object
    static CConfigParams& GetTheInstance();

    /// \brief    Load from file
    virtual void Load();


    /// \brief  FILE_NAME  Name of config file
    static const string  FILE_NAME;

    // Configurable parameters
    string DirectoryName;
    string CaseFileName;
    string ControlFileName;
    unsigned int ChromosomeLength;
    unsigned int NumOfIterations;
    unsigned int PopulationSize;
    float ElitismRate;
    float CrossoverRate;
    float MutationRate;
    float TrapRatio;
    float VibrationRate;
    EGenotypeDataProvider Provider;
    float DisplayRatio;
    vector<string> SelectedBarcode;
    float HomogeneousRatio;
    bool IgnoreGenotype2;

    ESelectionAlgorithm SelectionAlgorithm;

    EAlgorithm Algorithm;

    unsigned int NumberOfExecutions;
    unsigned int ExecutionIsStcuk;

    unsigned int NumberOfOrders;

    EPostProcessingAlgorithm  PostProcessingAlg;

    unsigned int HaltCriteria;

private: 

    /// \brief    Singleton object
    static unique_ptr<CConfigParams> sTheInstance;

    /// \brief    Object can be created only by Singleton function
    CConfigParams();

    ///\brief    Create the singleton instance
    /// 
    /// Create the singleton object (this function comes to make sure the get instance will be inlined) 
    static void  CreateTheInstance();

    /// \brief    reads the value of input keyword from config file
    /// \param[in] file    reference to configuration file.
    /// \param[in] keyword    keyword to retrieve.
    /// \return std::string    retrieved value of input keyword.
    string ReadNextToken(ifstream& file, const string& keyword);
};

inline CConfigParams&   CConfigParams::GetTheInstance() 
{
    // This statement is required in case of more than one thread tries to create the object at the same time
    if (sTheInstance == NULL)      
    {
        // Do that in separate function in order to make it inline and save run time
        CreateTheInstance();
    }

    return *sTheInstance;
}

