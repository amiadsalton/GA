#pragma once

#include <string>
#include <vector>

#ifdef _MSC_VER
#include <windows.h>
#else
#include <string.h>
typedef unsigned char       BYTE;
#endif

using namespace std;

////////////////////////////////////////////////////////////////////////
/// \brief    Interface for case and control dataset
///
/// Role: Represents input data set of case and controls genotypes.
///
/// Responsibilities: 
///      1. Access to Case and Control datasets 
///      2. Common interface for all concrete data set
///      3. Factory of concrete data set
////////////////////////////////////////////////////////////////////////
class CGenotypeData
{
public:

    /// \brief    Number of possible genotypes: AA, Aa or aa
    static const int NUMBER_OF_GENOTYPES = 3;

    /// \brief    Input data for one SNP
    struct CSnpData
    {
        CSnpData() : mName() { memset(mEvents, 0, NUMBER_OF_GENOTYPES*sizeof(unsigned int)); }

        /// \brief    SNP name as appears in input data set
        string mName;

        /// \brief    SNP index as appears in input data set
        unsigned int mIndex;

        /// \brief    Occurrences of each genotype AA, Aa or aa
        unsigned int mEvents[NUMBER_OF_GENOTYPES];

        bool operator==(const CSnpData& input) { return mName == input.mName && memcmp(mEvents, input.mEvents, NUMBER_OF_GENOTYPES) == 0;}
    };

    /// \brief    Access to Singleton object.
    static CGenotypeData* CreateConcrete();

    /// \brief    Load from file
    virtual bool Load() = 0;

    /// \brief    Retrieve the number of analyzed SNPs
    unsigned int NumberOfSNPs() const { return (unsigned int)(mCaseData[0].size()); }

    /// \brief    Remove the i-th SNP and all it associated data. 
    void RemoveSNP(const vector<int>& iSnpsToRemove);

    //-------------------
    // Getters
    //-------------------
    const vector<vector<BYTE>>& CaseData() const { return mCaseData; }
    const vector<vector<BYTE>>& ControlData() const { return mControlData; }
    const vector<CSnpData>& TotalCase() const { return mTotalCase; }
    const vector<CSnpData>& TotalControl() const { return mTotalControl; }

protected:


    /// \brief    Data of each item in Case population. Each item has genotype definition for each SNP. 
    vector<vector<BYTE>> mCaseData;

    /// \brief    Data of each item in Control population. Each item has genotype definition for each SNP. 
    vector<vector<BYTE>> mControlData;

    /// \brief    Genotype distributions of all analyzed SNPs in case dataset
    vector<CSnpData> mTotalCase;

    /// \brief    Genotype distributions of all analyzed SNPs in control dataset
    vector<CSnpData> mTotalControl;
};

