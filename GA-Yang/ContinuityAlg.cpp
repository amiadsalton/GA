#include "ContinuityAlg.h"
#include "Chromosome.h"

#include <algorithm>

using namespace std;

void CContinuityAlg::Update(const std::vector<CChromosome>& iElitism)
{
    if (iElitism.empty())
    {
        return;
    }

    // Allocate new entry for current order
    if (mBestChromosomesPerRun.empty() || 
        (!mBestChromosomesPerRun[mBestChromosomesPerRun.size()-1].empty() && iElitism[0].Length() > mBestChromosomesPerRun[mBestChromosomesPerRun.size()-1].at(0).Length()) )
    {
        mBestChromosomesPerRun.push_back(vector<CChromosome>());
    }


    // Update Best Chromosomes
    vector<CChromosome>& bests = mBestChromosomesPerRun[mBestChromosomesPerRun.size()-1];
    for (auto elem : iElitism)
    {
        if (find(bests.begin(), bests.end(), elem) == bests.end())
        {
            bests.push_back(elem);
        }
    }
    sort(bests.begin(), bests.end());
    reverse(bests.begin(), bests.end());
}


bool CContinuityAlg::DerivedExecute(const CGenotypeData& ioGenotypeData, int iOrder)
{
    // Analyze results in order to find the SNPs combination  with higher correlation with input disease:
    // Find only chromosomes that are subsequence in higher order:

    if (mBestChromosomesPerRun.size() == 1)
    {
        return true;
    }

    mAlgorithmResults.push_back(vector<pair<CChromosome,CChromosome>>());
    vector<pair<CChromosome,CChromosome>>& currResult=*mAlgorithmResults.rbegin();

    // Best chromosome in last iteration
    const CChromosome& best = *mBestChromosomesPerRun.rbegin()->begin();

    // Find if it has sub chromosome in previous bests
    for (auto chromosome : mBestChromosomesPerRun[mBestChromosomesPerRun.size() -2])
    {
        unsigned int k=0;
        for (; k < chromosome.Length() ; k++)
        {
            unsigned j=0;
            for (j; j<best.Length() ; j++)
            {
                if ( (chromosome.SNP(k) == best.SNP(j)) && (chromosome.Genotype(k) == best.Genotype(j)) )
                {
                    // SNP[k] matches between best and chromosome
                    break;
                }
            }
            if (j == best.Length())
            {
                // No match
                break;
            }
        }

        if (k == chromosome.Length())
        {
            // "best" is subsequence of "chromosome", is inserted into foundInHigherOrder
            currResult.push_back(make_pair(chromosome,best));
        }
    }

    return !currResult.empty();
}

void CContinuityAlg::Print(const CGenotypeData& ioGenotypeData, const std::string& iBaseFileName, int iOrder) const
{
    // Print all Bests
    ofstream summaryFile(iBaseFileName + "_" + to_string(iOrder) + "_Summary.txt");
    if (!summaryFile.good())
    {
        cout << "Failed to open file: " << endl;
        return;
    }

    summaryFile <<  "\nNumber of Fitness calculations: " << CChromosome::FitnessCalculationCounter() << endl;
    summaryFile <<  "\nTotal Number of Fitness calculations: " << mNumberOfFitnessCalcualtions << endl;

    PrintAllBests(ioGenotypeData, summaryFile);

    int order(iOrder);
    if (mAlgorithmResults.empty())
    {
        summaryFile << "\n=========\nNo Best Results for order: " << order << endl;
    }
    for (auto elem : mAlgorithmResults)
	{
        if (elem.empty())
        {
            summaryFile << "\n=========\nNo Best Results for order: " << ++order << endl;
        }
        else
        {
            order = elem[0].second.Length();
	        // Need to add 1 since the first element (0) in mBestChromosomesPerRun is of order 2. 
	        summaryFile << "\n=========\nOrder: " << order << endl;
	        for  (auto  chrPair : elem)
	        {
	            summaryFile <<  CChromosomeSerializer(chrPair.first, ioGenotypeData); ;
	            summaryFile << " in ";
                summaryFile <<  CChromosomeSerializer(chrPair.second, ioGenotypeData); ;
	            summaryFile << endl;
	        }
	        summaryFile << endl;
        }
	}

    summaryFile.close();
}

void CContinuityAlg::PrintAllBests(const CGenotypeData& ioGenotypeData, std::ofstream& oStream ) const
{
    for (auto best : mBestChromosomesPerRun)
    {
        if (best.size() == 0)
        {
            oStream << "No Bests for this order " <<  endl;
        }
        else
        {
            oStream << "Number of Bests of order " << best[0].Length()  << ": " << best.size() << endl;
            for (size_t i=0;  i < min(best.size(), size_t(30)) ; i++)
            {
                oStream << CChromosomeSerializer(best[i],ioGenotypeData);
            }
            oStream << endl;
        }
    }
}
