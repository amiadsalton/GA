
#include "YangData.h"
#include "ConfigParams.h"

#include <sstream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

bool CYangData::Load()
{
    bool result(true);

    const std::string target_path(CConfigParams::GetTheInstance().DirectoryName);
    const boost::regex my_filter(".(dat|sln)$" );

    string caseFileName, controlFileName;
    boost::filesystem::directory_iterator end_itr; // Default ctor yields past-the-end
    for( boost::filesystem::directory_iterator i( target_path ); i != end_itr; ++i )
    {
        // Skip if not a file
        if( !boost::filesystem::is_regular_file( i->status() ) ) 
            continue;

        // Skip if no match
        if( !boost::regex_search(i->path().leaf().string(),my_filter) ) 
            continue;

        // File matches, open it
        string fileName = i->path().leaf().string();
        static const string CASE_FILE_NAME("case_genotypes");
        static const string ANTICASE_FILE_NAME("anticase_genotypes");
        if (fileName.substr(0,CASE_FILE_NAME.size()) == CASE_FILE_NAME)
        {
            caseFileName = i->path().string();
        }
        else if (fileName.substr(0,ANTICASE_FILE_NAME.size()) == ANTICASE_FILE_NAME)
        {
            controlFileName = i->path().string();
        } 
    }

    float homogeneousRatio  = CConfigParams::GetTheInstance().HomogeneousRatio;
    if (!caseFileName.empty() && !controlFileName.empty())
    {
        ifstream  caseFile(caseFileName), controlFile(controlFileName);
        if (!caseFile.good() || caseFile.eof() || !controlFile.good() || controlFile.eof() )
        {
            return false;
        }

        while (caseFile.good() && !caseFile.eof() && controlFile.good() && !controlFile.eof())
        {
            CSnpData snpCase, snpControl;

            // Read each line until EOF
            static const int BUFFER_SIZE=2<<20;
            static char tmp[BUFFER_SIZE];
            {
                memset(tmp, 0, BUFFER_SIZE);
                caseFile.getline(tmp, BUFFER_SIZE);               
                string tmpStr(tmp);
                if (LoadLine(tmpStr, mTotalCase, mCaseData) == false)
                {
                    // EOF
                    break;
                }
            }

            {
                memset(tmp, 0, BUFFER_SIZE);
                controlFile.getline(tmp, BUFFER_SIZE);
                string  tmpStr(tmp);
                if (LoadLine(tmpStr, mTotalControl, mControlData) == false)
                {
                    result = false;
                    break;
                }
            }
        }
    }

    return result;
}

bool CYangData::LoadLine(const string& line, vector<CSnpData>& snpSet, vector<vector<BYTE>>& dataSet)
{
    // Each line contains data for one item in case/control group
    if (!line.empty())
    {
        // Retrieve data from one line
        istringstream  stream(line);

        dataSet.push_back(vector<BYTE>());

        // For  Yang's data range is 1-3 instead of 0 -2 
        int sub = mYangRange ? 1 : 0;

        int index(0);
        do 
        {
            // Count the number of occurrences of each genotype
            int elem;
            stream >> elem;
            elem -= sub;

            if (snpSet.size() <= index)
            {
                CSnpData snp;
                snp.mIndex = index + 1;
                snp.mName = "SNP" + to_string(snp.mIndex);
                snpSet.push_back(snp);
            }
            snpSet[index].mEvents[elem]++;

            dataSet[dataSet.size()-1].push_back(BYTE(elem));

            index++;
        } while (stream.good() && !stream.eof());
    }
    return !line.empty();
}
