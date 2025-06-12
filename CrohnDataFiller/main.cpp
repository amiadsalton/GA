// SyntheticGenerator.cpp : Defines the entry point for the console application.
//

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <random>

#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>

#include <windows.h>
#include <time.h>

using namespace std;

struct CSnpData
{
    unsigned int mEvents[3];
    int mNumberOfNA;
    vector<int> mFiller;
    string mName;
    CSnpData(const string& iName) : mName(iName), mNumberOfNA(0) { memset(mEvents, 0 , 3* sizeof(int)); }
};

vector<CSnpData>  sCaseData;
vector<CSnpData>  sControlData;

bool Load(const string& iFileName)
{
    bool result  = true;

    ifstream  inputFile(iFileName);
    if (!inputFile.good() )
    {
        return false;
    }

    static const int BUFFER_SIZE=2<<20;
    static char tmp[BUFFER_SIZE];

    bool isFirstLine(true);
    int numberOfCases(0), numberOfControls(0);

    while (inputFile.good())
    {
        size_t currentSnp(0);

        memset(tmp, 0, BUFFER_SIZE);

        inputFile.getline(tmp, BUFFER_SIZE);   
        string line(tmp);
        
        vector<std::string> tokens;
        boost::char_separator<char> sep(", ");

        boost::tokenizer<boost::char_separator<char>> tokensBoost(line, sep);
        BOOST_FOREACH (const string& t, tokensBoost)
        {
            tokens.push_back(t);
        }

        if (tokens.empty())
        {
            break;
        }

        if (isFirstLine)
        {
            for (size_t i=1; i<tokens.size() ; i++)
            {
                sCaseData.push_back(CSnpData(tokens[i]));
                sControlData.push_back(CSnpData(tokens[i]));
            }
            isFirstLine = false;
            continue;
        }

        vector<CSnpData>* currSnpVector;
        if (tokens[0].substr(0,2) == "CD")
        {
            currSnpVector = &sCaseData;
            numberOfCases++;
        }
        else
        {
            currSnpVector = &sControlData;
            numberOfControls++;
        }
        currentSnp++;

        if (numberOfControls > numberOfCases)
        {
            break;
        }

        // Ignore the instance name
        for (size_t i=1; i<tokens.size() ; i++)
        {
            if (tokens[i] == "NA")
            {
                currSnpVector->operator[](i-1).mNumberOfNA++; 
            }
            else
            {
                currSnpVector->operator[](i-1).mEvents[stoi(tokens[i])]++; 
            }
        }
    }

    return true;
}


bool FillAndSave(const string& iFileName, int iNumberOfElements)
{
    // Read the input for each line and 
    bool result  = true;

    ifstream  inputFile(iFileName);
    ofstream  outputFile("Filled_" + iFileName);

    if (!inputFile.good() || !outputFile.good())
    {
        return false;
    }

    static const int BUFFER_SIZE=2<<20;
    static char tmp[BUFFER_SIZE];

    // Copy the first line
    memset(tmp, 0, BUFFER_SIZE);
    inputFile.getline(tmp, BUFFER_SIZE); 
    string line(tmp);
    line += '\n';
    outputFile.write(line.c_str(), line.size());

    default_random_engine randomEngine((unsigned int)time(0));

    for (int i=0 ; i< 2*iNumberOfElements ; i++)
    {
        size_t currentSnp(0);

        memset(tmp, 0, BUFFER_SIZE);

        inputFile.getline(tmp, BUFFER_SIZE);   
        if (!inputFile.good())
        {
            return false;
        }

        line = tmp;

        vector<std::string> tokens;
        boost::char_separator<char> sep(", ");

        boost::tokenizer<boost::char_separator<char>> tokensBoost(line, sep);
        BOOST_FOREACH (const string& t, tokensBoost)
        {
            tokens.push_back(t);
        }

        vector<CSnpData>* currSnpVector;
        if (tokens[0].substr(0,2) == "CD")
        {
            currSnpVector = &sCaseData;
        }
        else
        {
            currSnpVector = &sControlData;
        }

        for (size_t i=1; i<tokens.size() ; i++)
        {
            if (tokens[i] == "NA")
            {
                CSnpData& elem = currSnpVector->operator[](i-1);
                uniform_int_distribution<unsigned int> mFillterDistribution(0,elem.mFiller.size()-1);
                size_t random = mFillterDistribution(randomEngine);

                tokens[i] = to_string(elem.mFiller[random]);

                elem.mFiller.erase(elem.mFiller.begin() + random);
            }
        }

        // Save
        for (size_t i=0; i<tokens.size() - 1 ; i++)
        {
            outputFile << tokens[i] << ",";
        }
        outputFile << tokens[tokens.size() - 1] << endl;

        if (!outputFile.good())
        {
            return false;
        }
    }

    return true;
}


void  CreateFillers(vector<CSnpData>& ioData, int iNumberOfElements)
{
    for (CSnpData& snp : ioData)
    {
        if (snp.mNumberOfNA > 0)
        {
            int numberOfFilled0 = int(snp.mNumberOfNA * double(snp.mEvents[0]) / (iNumberOfElements - snp.mNumberOfNA) + 0.5);
            for (int i=0; i < numberOfFilled0 ; i++)
            {
                snp.mFiller.push_back(0);
            }
            int numberOfFilled1 = int(snp.mNumberOfNA * double(snp.mEvents[1]) / (iNumberOfElements - snp.mNumberOfNA) + 0.5);
            for (int i=0; i < numberOfFilled1 ; i++)
            {
                snp.mFiller.push_back(1);
            } 
            for (int i=0; i < snp.mNumberOfNA-(numberOfFilled0 + numberOfFilled1) ; i++)
            {
                snp.mFiller.push_back(2);
            } 
        }
    }
}

bool  ValidateData(const vector<CSnpData>& iData, int numberOfElements)
{
    for (size_t i = 0; i<iData.size() ;i++)
    {
        if ((iData[i].mNumberOfNA + iData[i].mEvents[0] + iData[i].mEvents[1] + iData[i].mEvents[2]) !=  numberOfElements)
        {
            cout << "Illegal number of events for SNP: " << i << endl;
            return false;
        }
    }

    return true;
}


int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout << "Missing input file name\n";
    }

    Load(argv[1]);

    int numberOfCases(sCaseData[0].mNumberOfNA + sCaseData[0].mEvents[0] + sCaseData[0].mEvents[1] + sCaseData[0].mEvents[2]);

    // Check validity
    if (ValidateData(sCaseData, numberOfCases) && ValidateData(sControlData, numberOfCases))
    {
        CreateFillers(sCaseData , numberOfCases);
        CreateFillers(sControlData , numberOfCases);
    }

    if (FillAndSave(argv[1], numberOfCases))
    {
        cout << "Date was successfully filled\n";
    }
    else
    {
        cout << "Error while trying to save data\n";
    }

    cout << "Enter any key to exit ==>";
    char ch;
    cin >> ch;

	return 0;
}

