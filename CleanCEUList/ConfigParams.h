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
    unsigned int NumberOfSNPs() const { return mNumberOfSNPs; }
    unsigned int StartChr() const { return mStartChr; }
    unsigned int EndChr() const { return mEndChr; }
    const vector<string>& SelectedBarcode() const { return mSelectedBarcode; }
    const vector<int>& SelectedChr() const { return mSelectedChr; }

    unsigned int BacrcodeLength() const { return mBarcodeLength; }
    float MinMaf() const {  return mMinMaf; }
    float MaxMaf() const {  return mMaxMaf; }

    /// \brief  FILE_NAME  Name of config file
    static const string  FILE_NAME;

private:

    // Configurable parameters

    unsigned int mNumberOfSNPs;
    unsigned int mStartChr;
    unsigned int mEndChr;
    vector<string> mSelectedBarcode;
    vector<int> mSelectedChr;
    unsigned int mBarcodeLength;
    float mMinMaf;
    float mMaxMaf;

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

