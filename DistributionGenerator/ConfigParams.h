#pragma once

#include <memory>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

////////////////////////////////////////////////////////////////////////
/// \brief    Configurable parameters
///
/// Role: Access to all configurable parameters
///
/// Responsibilities: 
///      1. Singleton
///      2. Knows how to load object from file
////////////////////////////////////////////////////////////////////////
class CConfigParams
{
public:

    /// \brief    Access to Singleton object
    static CConfigParams& GetTheInstance();

    /// \brief    Load from file
    virtual void Load();

    // Getters
    const string& DirectoryName() const { return mDirectoryName; }
    unsigned int ChromosomeLength() const { return mLength; }

    unsigned int NominalMin() const { return mNominalMin; }
    unsigned int NominalMax() const { return mNominalMax; }
    unsigned int PopulationSize() const { return mPopulationSize; }

    double RationalFactor() const { return mRationalFactor; }

    const vector<string>& SelectedSNPs() const { return mSelectedSNPs; }

    const vector<double>& MAF() const { return mMAF; }

    //const vector<double>& Factors() const { return mFactors; }

    /// \brief  FILE_NAME  Name of config file
    static const string  FILE_NAME;

private:

    // Configurable parameters
    string mDirectoryName;
    unsigned int mLength;
    
    unsigned int mNominalMin;
    unsigned int mNominalMax;
    unsigned int mPopulationSize;

    double mRationalFactor;

    vector<string> mSelectedSNPs;

    vector<double> mMAF;

    vector<double> mFactors;

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

