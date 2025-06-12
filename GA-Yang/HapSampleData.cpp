
#include "HapSampleData.h"
#include "ConfigParams.h"

#include <sstream>
#include <fstream>
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

#undef min
#undef max 

//#define  BLACK_LIST
#define  BLACK_LIST_SIZE 10

bool CHapSampleData::Load()
{
    bool result(true);

    const std::string target_path(CConfigParams::GetTheInstance().DirectoryName);
    boost::system::error_code ec;
    if (!boost::filesystem::exists(target_path,ec) || ec != boost::system::errc::success)
    {
    	cout << "Can't open directory: " << target_path << endl;
    	return false;
    }

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

    vector<string> blackList = BuildBlackList(caseFileName);

    float homogeneousRatio  = CConfigParams::GetTheInstance().HomogeneousRatio;
    if (!caseFileName.empty() && !controlFileName.empty())
    {
        ifstream  caseFile(caseFileName), controlFile(controlFileName);
        if (!caseFile.good() || caseFile.eof() || !controlFile.good() || controlFile.eof() )
        {
            return false;
        }

        unsigned int snpIndex = 1;
        while (caseFile.good() && !caseFile.eof() && controlFile.good() && !controlFile.eof())
        {
            CSnpData snpCase, snpControl;
            snpCase.mIndex = snpControl.mIndex = snpIndex++;

            // Read each line until EOF
            static const int BUFFER_SIZE=2<<20;
            static char tmp[BUFFER_SIZE];
            {
	            memset(tmp, 0, BUFFER_SIZE);
	            caseFile.getline(tmp, BUFFER_SIZE);               
	            string tmpStr(tmp);
	            if (LoadSNP(tmpStr, snpCase, mCaseData, blackList) == false)
	            {
	                // EOF
	                break;
	            }
            }

            {
	            memset(tmp, 0, BUFFER_SIZE);
	            controlFile.getline(tmp, BUFFER_SIZE);
	            string  tmpStr(tmp);
	            if (LoadSNP(tmpStr, snpControl, mControlData, blackList) == false)
	            {
	                result = false;
	                break;
	            }
            }

            // Don't push empty elements or elements that have Genotype 2 in control than in case.  
            bool remove=false;  
            if ( (CConfigParams::GetTheInstance().IgnoreGenotype2 == true) && (snpControl.mEvents[2] > snpCase.mEvents[2]) )
            {
                remove = true;
            }
            else
            {
                for (size_t i=0; i<NUMBER_OF_GENOTYPES ;i++)
                {
                    if ( ((float(snpCase.mEvents[i]) / mCaseData.size()) > homogeneousRatio) &&
                        ((float(snpControl.mEvents[i]) / mControlData.size()) > homogeneousRatio) )
                    {
                        remove = true;
                        break;
                    }
                }
            }


            //if (snpCase.mEvents[2] == 0)
            //{
            //    remove = true;
            //}

            if (remove)
            {
                for (size_t j=0; j<mCaseData.size() ;j++)
                {
                    mCaseData[j].pop_back();
                    mControlData[j].pop_back();
                }
            }
            else  if (find(blackList.begin(), blackList.end(), snpCase.mName) == blackList.end())
            {
                mTotalCase.push_back(snpCase);
                mTotalControl.push_back(snpControl);
            }
        }
    }

    return result;
}

bool CHapSampleData::LoadSNP(const string& line, CSnpData& snp, vector<vector<BYTE>>& dataSet, const vector<string>& iBlackList)
{
    if (!line.empty())
    {
        // Retrieve data from one line
        istringstream  stream(line);
        string dummy1, dummy2, dummy3;
        stream >> dummy1 >> snp.mName >> dummy2 >> dummy3;

        //            file << dummy1 << " " << snp.mName << " " << dummy2 << " " << dummy3<< " " ;

        if (find(iBlackList.begin(), iBlackList.end(), snp.mName) != iBlackList.end())
        {
            return true;
        }

        int index(0);
        do 
        {
            // Count the number of occurrences of each genotype
            int elem;
            stream >> elem;
            snp.mEvents[elem]++;
            //                file << elem << " ";

            // There is no struct for that item yet
            if (dataSet.size() <= index)
            {
                dataSet.push_back(vector<BYTE>());
            }
            dataSet[index].push_back(BYTE(elem));
            index++;
        } while (stream.good() && !stream.eof());
    }
    return !line.empty();
}

bool CHapSampleData::LoadFile(const string& fileName, vector<vector<BYTE>>& dataSet, vector<CSnpData>& totalSnp)
{

    ifstream  casesFile(fileName);
    if (!casesFile.good() || casesFile.eof())
    {
        return false;
    }

//    ofstream file("kkk.txt");

    float homogeneousRatio  = CConfigParams::GetTheInstance().HomogeneousRatio;

    while (casesFile.good() && !casesFile.eof())
    {
        // Read each line until EOF
        static char tmp[8192];
        memset(tmp, 0, 8192);
        casesFile.getline(tmp, 8192);
        string line(tmp);
        
        if (!line.empty())
        {
            // Retrieve data from one line
            istringstream  stream(line);
            string dummy1, dummy2, dummy3;
            CSnpData snp;
            stream >> dummy1 >> snp.mName >> dummy2 >> dummy3;

//            file << dummy1 << " " << snp.mName << " " << dummy2 << " " << dummy3<< " " ;

            int index(0);
            do 
            {
                // Count the number of occurrences of each genotype
                int elem;
                stream >> elem;
                snp.mEvents[elem]++;
//                file << elem << " ";

                // There is no struct for that item yet
                if (dataSet.size() <= index)
                {
                    dataSet.push_back(vector<BYTE>());
                }
                dataSet[index].push_back(BYTE(elem));
                index++;
            } while (stream.good() && !stream.eof());
//            file << endl;

            // Don't push empty elements 
            bool remove=false;
          
            for (size_t i=0; i<NUMBER_OF_GENOTYPES ;i++)
            {
                if ((float(snp.mEvents[i]) / dataSet.size()) > homogeneousRatio)
                {
                    remove = true;
                    break;
                }
			}
           
            if (remove)
            {
                for (size_t j=0; j<dataSet.size() ;j++)
                {
                    dataSet[j].pop_back();
                }
            }
            else
            {
                totalSnp.push_back(snp);
            }
        }
    }

    return true;
}

void CHapSampleData::Reorganize(const vector<pair<string,BYTE>>& barCode)
{
    // At first find the genotype with lowest number events 
    size_t lowestNumberOfEvents=mCaseData.size();
    size_t lowestIndex=0;

    for (size_t i=0 ; i<barCode.size() ; i++)
    {
        size_t j=0;
        for (; j<mTotalCase.size() ; j++)
        {
            if (mTotalCase[j].mName == barCode[i].first)
            {
                break;
            }
        }
        if (mTotalCase[j].mEvents[barCode[i].second] < lowestNumberOfEvents)
        {
            lowestNumberOfEvents = mTotalCase[j].mEvents[barCode[i].second];
            lowestIndex = j;
        }
    }

    // Go over CaseData and reshuffle it, so there will be perfect match 
}

vector<string> CHapSampleData::BuildBlackList( const string &caseFileName ) const
{
    vector<string> blackList;

#ifdef BLACK_LIST
    // First pass - find selected Barcode and add close SNPs to a black list:
    vector<string> snpList;
    if (!caseFileName.empty())
    {
        ifstream  caseFile(caseFileName);
        if (!caseFile.good() || caseFile.eof() )
        {
            return blackList;
        }
        unsigned int snpIndex = 1;
        while (caseFile.good() && !caseFile.eof())
        {
            static const int BUFFER_SIZE=2<<20;
            static char tmp[BUFFER_SIZE];
            {
                memset(tmp, 0, BUFFER_SIZE);
                caseFile.getline(tmp, BUFFER_SIZE); 
                istringstream  stream(tmp);
                string dummy1, snpNamp;
                stream >> dummy1 >> snpNamp;
                snpList.push_back(snpNamp);
            }
        }

        vector<string> selected = CConfigParams::GetTheInstance().SelectedBarcode;
        for (size_t i=0; i<snpList.size() ; i++)
        {
            for (size_t j=0; j<selected.size() ; j++)
            {
                if (snpList[i] == selected[j])
                {
                    for (int k = std::max(0,int(i-BLACK_LIST_SIZE)) ; k < std::min(snpList.size(), i + BLACK_LIST_SIZE) ; k++)
                    {
                        if (k!=i)
                        {
                            blackList.push_back(snpList[k]);
                        }
                    }
                }
            }
        }
    }
#endif

    return blackList;
}
