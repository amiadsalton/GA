
#include "ConfigParams.h"
#include <sstream>
#include <boost/token_functions.hpp>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>


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

        string token=ReadNextToken(file, "Length");
        istringstream stream = istringstream(token);
        stream >> mLength;
       
        token=ReadNextToken(file, "Selected");
        boost::char_separator<char> sep(",");
        boost::tokenizer<boost::char_separator<char>> tokensBoost(token, sep);
        BOOST_FOREACH (const string& t, tokensBoost)
        {
            mSelectedSNPs.push_back(t);
        }

        token=ReadNextToken(file, "MAF");
        stream = istringstream(token);
        boost::tokenizer<boost::char_separator<char>> tokensBoost1(token, sep);
        BOOST_FOREACH (const string& t, tokensBoost1)
        {
            mMAF.push_back(std::stod(t)); 
        }

        token=ReadNextToken(file, "RationalFactor");
        stream = istringstream(token);
        stream >> mRationalFactor;

        token=ReadNextToken(file, "NominalMin");
        stream = istringstream(token);
        stream >> mNominalMin;

        token=ReadNextToken(file, "NominalMax");
        stream = istringstream(token);
        stream >> mNominalMax;

        token=ReadNextToken(file, "PopulationSize");
        stream = istringstream(token);
        stream >> mPopulationSize;

//         token=ReadNextToken(file, "NominalExponent");
//         stream = istringstream(token);
//         stream >> mNominalExponent;

// 
//         token=ReadNextToken(file, "Factors");
//         stream = istringstream(token);
//         boost::tokenizer<boost::char_separator<char>> tokensBoost2(token, sep);
//         BOOST_FOREACH (const string& t, tokensBoost2)
//         {
//             mFactors.push_back(std::stod(t)); 
//         }
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
