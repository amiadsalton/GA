// SyntheticGenerator.cpp : Defines the entry point for the console application.
//


#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <boost/filesystem.hpp>

#include "ConfigParams.h"

using namespace std;

typedef unsigned char BYTE;

// Generate Data for one SNP and continue in  recursion to other chromosomes
void GenerateData(vector<double>& oResults, vector<int>& ioGenotypes, int iChoromosomeIndex)
{
    if (iChoromosomeIndex < 0)
    {
        int countNumOf2 = 0;
        for (int i=0; i<ioGenotypes.size() ;i++)
        {
            if (ioGenotypes[i] == 2)
            {
                countNumOf2++;
            }
        }

        vector<double> MAF(CConfigParams::GetTheInstance().MAF());

        double result=1.0;
        for (int i=0; i<ioGenotypes.size() ;i++)
        {
            if (ioGenotypes[i] == 0)
            {
                result *= pow(1-MAF[i],2);
            }
            else if (ioGenotypes[i] == 1)
            {
                result *= 2 * (1-MAF[i]) * MAF[i];
            }
            else
            {
                result *= pow(MAF[i],2);
            }
        }

        // Enrich all 2 Genotypes
        if (countNumOf2 > 0)
        {
            result *= pow(CConfigParams::GetTheInstance().RationalFactor(), countNumOf2);

            double base = pow(CConfigParams::GetTheInstance().NominalMax() - CConfigParams::GetTheInstance().NominalMin() , 1.0 / CConfigParams::GetTheInstance().ChromosomeLength());
            double nominalGradient = (CConfigParams::GetTheInstance().NominalMin() + pow(base , countNumOf2 )) / CConfigParams::GetTheInstance().PopulationSize();

            result = max(result, nominalGradient);
        }

        oResults.push_back(result);
    }
    else
    {
        for (int j=0; j<=2 ; j++)
        {
            GenerateData(oResults, ioGenotypes, iChoromosomeIndex-1);
            ioGenotypes[iChoromosomeIndex]++;
        } 
        ioGenotypes[iChoromosomeIndex] = 0;
    }
}

// Generate Data for one SNP and continue in  recursion to other chromosomes
void Fix(vector<double>& ioResults, size_t&  ioIndex, vector<int>& ioGenotypes, double& oSumOfEnriched, int iChoromosomeIndex, double correction)
{
    if (iChoromosomeIndex < 0)
    {
        int countNumOf2 = 0;
        for each (auto item in ioGenotypes)
        {
            if (item == 2)
            {
                countNumOf2++;
            }
        }

        if (countNumOf2 == 0)
        {
            ioResults[ioIndex] *= correction;
        }
        else
        {
            oSumOfEnriched += ioResults[ioIndex];
        }
        ioIndex++;
    }
    else
    {
		// fix for both 0 and 1 indexes. 
        for (int j=0; j<=2 ; j++)
        {
            Fix(ioResults, ioIndex, ioGenotypes, oSumOfEnriched, iChoromosomeIndex-1, correction);
            ioGenotypes[iChoromosomeIndex]++;
        } 
        ioGenotypes[iChoromosomeIndex] = 0;
    }
}

int main(int argc, char* argv[])
{
    const CConfigParams& config = CConfigParams::GetTheInstance();   

    boost::filesystem::path dirName = config.DirectoryName();
    if (!dirName.string().empty() && !boost::filesystem::exists(dirName))
    {
	    boost::filesystem::create_directory(dirName);
    }

    boost::filesystem::path fileName = dirName;
    fileName /= (string("Model_") + to_string(config.ChromosomeLength()) + ".txt");
    if (boost::filesystem::exists(fileName))
    {
        boost::filesystem::remove(fileName);
    } 
   
    // Generate data
    vector<int> genotypes(config.ChromosomeLength(), 0);
    vector<double> results;
    GenerateData(results, genotypes, config.ChromosomeLength()-1);

    bool success(true);
    do 
    {
        double sum=0;
        for each (auto item in results)
        {
            sum += item;
        }

        if (abs(sum - 1) < FLT_EPSILON)
        {
            break;
        }
        
        double correction = 1.0 / sum;
        double sumOfEnriched=0;
        size_t index=0;
        Fix(results, index, genotypes, sumOfEnriched, config.ChromosomeLength()-1, correction);

        if (sumOfEnriched > 1)
        {
            success = false;
            cout << "Total probabilities of enriched combinations is bigger than 1 : " << sumOfEnriched << endl;
            break;
        }

    } while (true);

    if (success)
    {
	    ofstream  outputFile(fileName.string());
	    if (!outputFile.good())
	    {
	        cout << "Failed to open output file: " << fileName.string() << endl;
	    }
	    
	    outputFile << config.ChromosomeLength() << endl;
	    const vector<string>& selected  = config.SelectedSNPs();
	    for each (auto item in selected)
	    {
	        outputFile << item << endl << flush;
	    }
	
	    outputFile << "0.001\nAG\n";
	
	    outputFile.precision(15);
	    for each (auto item in results)
	    {
	        outputFile << item << endl << flush;
	    }
	
	    outputFile.close();
	    if (!outputFile.good())
	    {
	        cout << "Failed to save output file: " << fileName.string() << endl;
	    }
    }

    char ch;
    cout << "created file: " <<  fileName.string() << "\nPress any key to exit...\n";
    cin >> ch;

	return 0; 
}

