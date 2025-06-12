
#include "ConfigParams.h"
#include <sstream>
#include <boost/token_functions.hpp>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>


unique_ptr<CConfigParams>  CConfigParams::sTheInstance;
const string CConfigParams::FILE_NAME="ConfigParams.txt";

CConfigParams::CConfigParams() : mBarcodeLength(0)
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
        string token=ReadNextToken(file, "Number");
        istringstream stream(token);
        stream >> mNumberOfSNPs;

        token=ReadNextToken(file, "Start");
        stream = istringstream(token);
        stream >> mStartChr;

        token=ReadNextToken(file, "End");
        stream = istringstream(token);
        stream >> mEndChr;

        token=ReadNextToken(file, "Barcode");
        stream = istringstream(token);
        boost::char_separator<char> sep(",");
        boost::tokenizer<boost::char_separator<char>> tokensBoost(token, sep);
        BOOST_FOREACH (const string& t, tokensBoost)
        {
            mSelectedBarcode.push_back(t);
//             stream = istringstream(t);
//             stream >> tmpInt;
//             mSelectedBarcode.push_back(tmpInt-1);
        }

        token=ReadNextToken(file, "BarcodeChr");
        stream = istringstream(token);
        boost::tokenizer<boost::char_separator<char>> tokensBoost1(token, sep);
        BOOST_FOREACH (const string& t, tokensBoost1)
        {
            mSelectedChr.push_back(std::stoul(t)); 
        }

        token=ReadNextToken(file, "BarcodeSize");
        if (token.empty())
        {
            return;
        }

        stream = istringstream(token);
        stream >> mBarcodeLength;

        token=ReadNextToken(file, "MinMAF");
        stream = istringstream(token);
        stream >> mMinMaf;

        token=ReadNextToken(file, "MaxMAF");
        stream = istringstream(token);
        stream >> mMaxMaf;
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
