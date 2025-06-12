
#include "ConfigParams.h"
#include <sstream>


unique_ptr<CConfigParams>  CConfigParams::sTheInstance;
const string CConfigParams::FILE_NAME="ConfigParams.txt";

CConfigParams::CConfigParams() 
{
}

void   CConfigParams::CreateTheInstance() 
{
    // This statement is required in case of more than one thread tries to create the object at the same time
    if (sTheInstance == NULL)      
    {
        try {
            sTheInstance = unique_ptr<CConfigParams>(new CConfigParams()); 
            sTheInstance->Load();                      
        } catch (...) {
            // This message will be written to log file anyway
            cout << "Failed to load config file: " + CConfigParams::FILE_NAME + ". Default parameters are used";
            sTheInstance = unique_ptr<CConfigParams>(new CConfigParams());   
        }
    }
}

void   CConfigParams::Load()
{
    ifstream file(FILE_NAME);
    if (file.good())
    {
        mDirectoryName = ReadNextToken(file, "Dir");

        string token=ReadNextToken(file, "Size");
        istringstream stream(token);
        stream >> mPopulationSize;

        token=ReadNextToken(file, "NumberOfSNPs");
        stream = istringstream(token);
        stream >> mNumberOfSNPs;

        token=ReadNextToken(file, "Length");
        stream = istringstream(token);
        stream >> mLength;

        float f;
        char ch;

        token=ReadNextToken(file, "Distribution");
        stream = istringstream(token);
        for (size_t i=0;  i<3 ;i++)
        {
            stream >> f >> ch;
            mDefaultDistribution.push_back(f);
        }

        token=ReadNextToken(file, "Selected");
        stream = istringstream(token);
        for (size_t i=0;  i<mLength ;i++)
        {
            int tmp;
            stream >> tmp >> ch;
            mSelectedSNPs.push_back(tmp);
        };

        token=ReadNextToken(file, "CaseDistribution");
        stream = istringstream(token);
        for (size_t i=0;  i<3 ;i++)
        {
            float f;
            stream >> f >> ch;
            mCaseDistribution.push_back(f);
        }
    }
}

inline string CConfigParams::ReadNextToken(ifstream& file, const string& keyword)
{
    static char tmp[2048];
    do 
    {
        file.getline(tmp,  2048);
    } while (file.good() && (tmp[0] == '#'));

    if (keyword.empty())
    {
        return string(tmp);
    }

    string line(tmp);
    if (line.substr(0, keyword.length()) == keyword)
    {
        return line.substr(keyword.length() + 1, line.length() - keyword.length() -1);
    }    

    cout << "Can't find " << keyword << " in config file\n";
    return "";
}
