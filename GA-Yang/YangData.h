#pragma once

#include "GenotypeData.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
/// \brief    Concrete data set of case and controls data
///
/// Role: Concrete data set generated by http://www.hapsample.org/
///
/// Responsibilities: 
///      1. Knows how to load case and controls data from file
////////////////////////////////////////////////////////////////////////
class CYangData : public CGenotypeData
{
public:

    CYangData(bool iYangRange) : mYangRange(iYangRange) {}

    virtual bool Load() override;

private:

    /// \brief    Load from file SNPs data into input collection 
    /// \param[in] fileName   File that contains cases data
    /// \param[out] dataSet   Data set of all individuals in population
    /// \param[out] totalSnp   All SNPs data
    bool LoadLine(const string& line, vector<CSnpData>& snpSet, vector<vector<BYTE>>& dataSet);

    /// \brief    For  Yang's data range is 1-3 instead of 0 -2 
    bool mYangRange;

};

