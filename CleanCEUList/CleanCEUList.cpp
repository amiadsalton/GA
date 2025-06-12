// CleanCEUList.cpp : Defines the entry point for the console application.
//

#include "ConfigParams.h"

#include <fstream>
#include <string>
#include <algorithm>
#include <map>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/token_functions.hpp>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>

using namespace std;

static const int MIN_LENGTH_FROM_SELECTED_SNP = 500;

void SelectedBarcodeAndChr(const vector<boost::filesystem::directory_entry>& iAllChrFiles, vector<string>& oBarcode, vector<int>& oChromosomes)
{
    if (CConfigParams::GetTheInstance().BacrcodeLength() == 0)
    {
        oBarcode = CConfigParams::GetTheInstance().SelectedBarcode();
        oChromosomes = CConfigParams::GetTheInstance().SelectedChr();
        return;
    }

    cout << "Finding " << CConfigParams::GetTheInstance().BacrcodeLength() << " SNPs by chance. Results in file SelectedSNP.txt" << endl;

    char buffer[4096];
    boost::char_separator<char> sep(" ");
    map <int , pair<string,float> >  SnpsSelected;

    for (unsigned int i=0; i < CConfigParams::GetTheInstance().BacrcodeLength() ;i++)
    {
        int randomChr;
        while (true)
        {
            randomChr = int(float(rand ()) / RAND_MAX * 22);
            bool found(false);
            for (auto elem : SnpsSelected)
            {
                if (elem.first == randomChr)
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                SnpsSelected[randomChr] = make_pair("",0.0f);
                break;
            }
        }

        cout << "Finding SNP in chromosome " <<  randomChr + 1 << endl;

        ifstream  inputFile(iAllChrFiles[randomChr].path().string());
        // Ignore the header 
        inputFile.getline(buffer, 4096);

        vector<pair<string,float>>  SnpsVec;

        // Read all SNPs into input buffer
        while (inputFile.good())
        {
            inputFile.getline(buffer, 4096);

            string  snp(buffer);
            boost::tokenizer<boost::char_separator<char>> tokensBoost1(snp, sep);
            if (tokensBoost1.end() == tokensBoost1.begin())
            {
                break;
            }
            boost::tokenizer<boost::char_separator<char>>::iterator it1 = tokensBoost1.begin();
            boost::tokenizer<boost::char_separator<char>>::iterator it2 = it1;
            std::advance(it2,14);
            SnpsVec.push_back(make_pair(*it1, std::stof(*it2)));
        }

        while (true)
        {
            int randomSnp = int(float(rand ()) / RAND_MAX * SnpsVec.size());
            if (randomSnp > MIN_LENGTH_FROM_SELECTED_SNP*2 && randomSnp < (SnpsVec.size() - MIN_LENGTH_FROM_SELECTED_SNP*2))
            {
                if (SnpsVec[randomSnp].second > CConfigParams::GetTheInstance().MinMaf() && SnpsVec[randomSnp].second < CConfigParams::GetTheInstance().MaxMaf())
                {
                    SnpsSelected[randomChr] = SnpsVec[randomSnp];
                    break;
                }
            }
        }
    }

    ofstream selectedFile("SelectedSNP.txt");
    selectedFile << "Chr=";
    for (auto elem : SnpsSelected)
    {
        selectedFile << elem.first + 1 << ",";
        oChromosomes.push_back(elem.first + 1);
    }


    selectedFile << "\nBarcode=";
    for (auto elem : SnpsSelected)
    {
        selectedFile << elem.second.first << ",";
        oBarcode.push_back(elem.second.first);
    }

    selectedFile << "\nMAF=";
    for (auto elem : SnpsSelected)
    {
        selectedFile << elem.second.second << ",";
    }
}

int main(int argc, char* argv[])
{
    char ch;

    time_t currentTIme;
    time(&currentTIme);
    srand (unsigned int(currentTIme ));

    vector<boost::filesystem::directory_entry> allChr;
    copy(boost::filesystem::directory_iterator("input"), boost::filesystem::directory_iterator(), back_inserter(allChr));
    sort (allChr.begin(), allChr.end());

    vector<string> barcode;
    vector<int> sourceChr;
    SelectedBarcodeAndChr(allChr, barcode, sourceChr);

    ofstream  outputFile("SNP_" + to_string(CConfigParams::GetTheInstance().NumberOfSNPs()) + ".txt");

    char buffer[4096];

    unsigned int numerOfSnpPerChr = CConfigParams::GetTheInstance().NumberOfSNPs() / (CConfigParams::GetTheInstance().EndChr() - CConfigParams::GetTheInstance().StartChr() + 1);
        
    for (unsigned int i=CConfigParams::GetTheInstance().StartChr(); i<=CConfigParams::GetTheInstance().EndChr() ; i++ )
    {
        unsigned int numerOfSnpInThisChr = numerOfSnpPerChr;
        cout << "Chromosome " << i << endl;
        ifstream  inputFile(allChr[i-1].path().string());
        // Ignore the header 
        inputFile.getline(buffer, 4096);

        // Read all SNPs into input buffer
        vector<string> SnpsVec, SelectedVec;

        while (inputFile.good())
        {
            inputFile.getline(buffer, 4096);
            string  snp(buffer);
            SnpsVec.push_back(snp.substr(0, snp.find_first_of(' ')));
        }

        int indexOfSnpsInThisChr=-1;
        for (size_t j=0; j < sourceChr.size() ; j++)
        {
            if (sourceChr[j] == i)
            {
                size_t k=0;
                for ( ; k<SnpsVec.size() ;k++)
                {
                    if (SnpsVec[k] == barcode[j])
                    {
                        // chromosome was found 
                        indexOfSnpsInThisChr = int(j);
                        numerOfSnpInThisChr--;
                        outputFile << barcode[j]  << endl;
                        SelectedVec.push_back(barcode[j]);
                        break;
                    }
                }
                if (k == SnpsVec.size())
                {
                    cout << "Chromosome: " << j+1 << " was not found" << endl;
                    cout << "Please enter any key to exit ...\n";
                    cin >> ch;
                    return 0;
                }
            }
        }

        unsigned int counter = 0;
        while (counter < numerOfSnpInThisChr)
        {
            // Random SNPs
            double r, x;
            x = rand ();
            r = x / RAND_MAX; 

            int index = int(r * (SnpsVec.size() - 1));

            // If it is already selected, skip it 
            size_t j=0;
            for (; j < SelectedVec.size() ; j++)
            {
                if (SelectedVec[j] == SnpsVec[index])
                {
                    break;
                }
            }
            if (j < SelectedVec.size() )
            {
                continue;
            }

            // if it is close to selected SNPs (by 10) then skip it as well
            if (indexOfSnpsInThisChr != -1)
            {
                int k=max(0,index - MIN_LENGTH_FROM_SELECTED_SNP);
                for ( ;  k < min(int(SnpsVec.size()), index + MIN_LENGTH_FROM_SELECTED_SNP) ; k++)
                {
                    if (SnpsVec[k] == barcode[indexOfSnpsInThisChr])
                    {
                        break;
                    }
                }
                if (  k < min(int(SnpsVec.size()), index + MIN_LENGTH_FROM_SELECTED_SNP) )
                {
                    continue;
                }
            }

            // Legal SNPs , add it to collection
            SelectedVec.push_back(SnpsVec[index]);
            outputFile << SnpsVec[index]  << endl;
            counter++;
        }
    }

    cout << "Please enter any key to exit ...\n";
    cin >> ch;

    return 0;
}


