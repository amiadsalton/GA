#include "ClusteringAlg.h"
#include "ChromosomeMooney.h"
#include "ChromosomeYang.h"

#include <algorithm>
#include <map>
#include <boost/filesystem.hpp>

using namespace std;

CClusteringAlg::~CClusteringAlg() 
{
}

void CClusteringAlg::Update( const std::vector<CChromosome>& iElitism )
{
    if (iElitism.empty())
    {
        return;
    }

    ostringstream  bestChrStr;
    if (mCurrentOrder == -1 || (iElitism[0].Length() != (unsigned int)(mCurrentOrder)))
    {
        mBestSnpsPerRun.push_back(vector<CChromosome>());
        mBestChromosome = iElitism[0];
        mCurrentOrder = iElitism[0].Length();
    }

    if (iElitism[0].Fitness() > mBestChromosome.Fitness())
    {
        mBestChromosome = iElitism[0];
    }

    vector<CChromosome>& currHistogram = *mBestSnpsPerRun.rbegin();

    for (size_t i=0; i<iElitism.size() ;i++)
    {
        //  Find all different combinations in this element
        for (unsigned int j=0; j<iElitism[i].Length() ; j++)
        {
            size_t foundElem=0;
            for (; foundElem < currHistogram.size() ;foundElem++)
            {
                if (currHistogram[foundElem].SNP(0) == iElitism[i].SNP(j))
					// &&  currHistogram[foundElem].Genotype(0) == iElitism[i].Genotype(j))
                {
                    break;
                }
            }
            if (foundElem == currHistogram.size())
            {
                currHistogram.push_back(CChromosome(iElitism[i].SNP(j),iElitism[i].Genotype(j),0,vector<string>()));
                foundElem = currHistogram.size() - 1;
            }

            // We don't just count the number of time that each SNP appears  in elitism group, but also give bonus according to his position in the list. 
            if (mAlgorithm == CConfigParams::eClusteringFitenss)
            {
                currHistogram[foundElem].Fitness( currHistogram[foundElem].Fitness() + iElitism[i].Fitness());
            }
            else if (mAlgorithm == CConfigParams::eClusteringPosition)
            {
                currHistogram[foundElem].Fitness( currHistogram[foundElem].Fitness() + (iElitism.size() - j));
            }
            else  // Sum
            {
                currHistogram[foundElem].Fitness( currHistogram[foundElem].Fitness() + 1);
            }

        }
    }

    sort(currHistogram.begin(), currHistogram.end());
    reverse(currHistogram.begin(), currHistogram.end());
}

bool CClusteringAlg::DerivedExecute(const CGenotypeData& ioGenotypeData, int iOrder)
{
    mAlgorithmResults.push_back(vector<pair<CChromosome, unsigned int>>());
    vector<pair<CChromosome, unsigned int>>& currResult = *mAlgorithmResults.rbegin();

    vector<unsigned int> chr;

    if (mAlgorithmResults.size() == 1)
    {
        for  (auto elem : *(mBestSnpsPerRun.rbegin()))
        {
            currResult.push_back(make_pair(elem, ioGenotypeData.TotalCase()[elem.SNP(0)].mIndex));
            chr.push_back(elem.SNP(0));
        }
        if (currResult.empty())
        {
            return false;
        }
    }
    else
    {
	    const vector<pair<CChromosome, unsigned int>>& lastResult = mAlgorithmResults[mAlgorithmResults.size()-2];
	
	    // Analyze results in order to find the SNPs combination  with higher correlation with input disease:
	    if (mBestSnpsPerRun.empty())
	    {
	        return false;
	    }
	
	    const vector<CChromosome>& currBests = *mBestSnpsPerRun.rbegin();
	    if   (currBests.size() < iOrder)
	    {
	        return false;
	    }
	
	    for (size_t i=0; i<iOrder ;i++)
	    {
	        currResult.push_back((make_pair(currBests[i], ioGenotypeData.TotalCase()[currBests[i].SNP(0)].mIndex)));
	        chr.push_back(currBests[i].SNP(0));
	    }
    }

    ostringstream  bestChrStr;
    bestChrStr << CChromosomeSerializer(mBestChromosome,ioGenotypeData);
    mBestChromosomePerRun.push_back(make_pair(bestChrStr.str(), mBestChromosome.Fitness()));


    // Calculate the Fitness of the best combination
    vector<BYTE> genotype(iOrder, 2);
    ostringstream  bestBarcodeStr;
    double bestChrFitenss(0);
	CChromosome*  bestChr = nullptr;
    if (CConfigParams::GetTheInstance().Algorithm == CConfigParams::eMooney)
    {
		bestChr = new CChromosomeMooney(iOrder , &chr[0], &genotype[0]);
    }
    else
    {
		bestChr = new CChromosomeYang(iOrder , &chr[0], &genotype[0]);
    }
	bestChrFitenss = bestChr->Fitness(ioGenotypeData.CaseData(), ioGenotypeData.ControlData());
	bestBarcodeStr << "Best SNPs: " << CChromosomeSerializer(*bestChr, ioGenotypeData) << endl;

	// Calculate the best barcode of bests SNPs
	CChromosome* bestChromosome = nullptr;
	CalcBestBarcode(iOrder, *bestChr, 0, bestChromosome, ioGenotypeData);
	bestBarcodeStr << "Best Barcode of of Best SNPs: " << CChromosomeSerializer(*bestChromosome, ioGenotypeData) << endl;

	delete bestChromosome;
	mBestBarcodePerRun.push_back(make_pair(bestBarcodeStr.str(), bestChrFitenss));

    if (mBestBarcodePerRun.size() > CConfigParams::GetTheInstance().HaltCriteria)
    {
        if ( (mBestBarcodePerRun[mBestBarcodePerRun.size() - 1].second < mBestBarcodePerRun[mBestBarcodePerRun.size() - 2].second) &&
             (mBestBarcodePerRun[mBestBarcodePerRun.size() - 2].second < mBestBarcodePerRun[mBestBarcodePerRun.size() - 3].second) &&
             (mBestChromosomePerRun[mBestChromosomePerRun.size() - 1].second < mBestChromosomePerRun[mBestChromosomePerRun.size() - 2].second) &&
             (mBestChromosomePerRun[mBestChromosomePerRun.size() - 2].second < mBestChromosomePerRun[mBestChromosomePerRun.size() - 3].second)    )
        {
            // "Halt" criteria: Regression in Fitness  for the last 2 Runs. 
            cout << "Execution is stopped: " 
                 << "Bests Barcode[last]-" << mBestBarcodePerRun[mBestBarcodePerRun.size() - 1].second
                 << "Bests Barcode[last-1]-" << mBestBarcodePerRun[mBestBarcodePerRun.size() - 2].second
                 << "Bests Barcode[last-2]-" << mBestBarcodePerRun[mBestBarcodePerRun.size() - 3].second
                 << "Best Chromosome[last]-" << mBestChromosomePerRun[mBestChromosomePerRun.size() - 1].second
                 << "Best Chromosome[last-1]-" << mBestChromosomePerRun[mBestChromosomePerRun.size() - 2].second
                 << "Best Chromosome[last-2]-" << mBestChromosomePerRun[mBestChromosomePerRun.size() - 3].second
                 << endl;
            return false;
        }
    }

    return true;
}

void CClusteringAlg::Print(const CGenotypeData& ioGenotypeData, const std::string& iBaseFileName, int iOrder) const
{
    // Print all Bests
    ofstream summaryFile(iBaseFileName + "_" + to_string(iOrder) + "_Summary.txt");
    if (!summaryFile.good())
    {
        cout << "Failed to open file: " << endl;
        return;
    }

    summaryFile <<  "\nNumber of Fitness calculations: " << CChromosome::FitnessCalculationCounter() << endl;
    summaryFile <<  "Total Number of Fitness calculations: " << mNumberOfFitnessCalcualtions << endl << endl;

    summaryFile << "(Total SNPs - " << ioGenotypeData.TotalCase().size() << ")" << endl << endl;

    // Print all best results
    auto bests = mBestSnpsPerRun.rbegin();
    if (bests->size() == 0)
    {
        summaryFile << "No best results for order " << iOrder << endl;
    }
    else
    {
        summaryFile << "Number of Bests of order " << iOrder  << ": " << bests->size() << endl;
        for (auto best :  *bests)
        {
            summaryFile << CChromosomeSerializer(best,ioGenotypeData) << ",";
        }
        summaryFile << endl;
    }

    
    // Print algorithm results
    unsigned int order = iOrder;
    if (mAlgorithmResults.size() > 1)
    {
        order = 2;
    }
    for (size_t i=0; i < mAlgorithmResults.size() ; i++, order++)
    {
        summaryFile << "\n=========\nBest Length: " << order << endl; 

        for (auto elem : mAlgorithmResults[i])
        {
            summaryFile << elem.second << "-" << elem.first.Fitness() << endl;
        }

		summaryFile << "Best Chromosome of this Run:" << mBestChromosomePerRun[i].first << endl;
		cout << endl << "Best Chromosome of this Run:" << mBestChromosomePerRun[i].first << endl;
		summaryFile << mBestBarcodePerRun[i].first << endl;
		cout << mBestBarcodePerRun[i].first << endl;
	}


    // Print all best SNPs
    boost::filesystem::path  allFilePath(iBaseFileName);
    allFilePath += "_";
    allFilePath += "_All.txt";

    boost::system::error_code ec;
    if (!boost::filesystem::exists(allFilePath, ec) || ec != boost::system::errc::success)
    {
        // Create the file
        ofstream allFile(allFilePath.string());
    }

    ofstream allFile(allFilePath.string(),ios::app);
      
    // Collect all bests in a map (ordered by index)
    map<int, double>  orderedBests;
    const vector<CChromosome>& currBests = *mBestSnpsPerRun.rbegin();
    for (auto elem : currBests)
    {
        orderedBests[ioGenotypeData.TotalCase()[elem.SNP(0)].mIndex] = elem.Fitness();
    }

    allFile << "\n==========\nOrder:" << iOrder << "\t" << orderedBests.size() << endl;
    for (auto elem : orderedBests)
    {
        allFile << elem.first << "\t" << elem.second << endl;
    }
}

void CClusteringAlg::UpdateGenotypeData(CGenotypeData& ioGenotypeData) 
{
    vector<CChromosome>& bests = *mBestSnpsPerRun.rbegin();

    const vector<CGenotypeData::CSnpData>& allSnps = ioGenotypeData.TotalCase();
   
    vector<int>  snpsToRemove;
    // Must start from the end in order to keep the same index as was used in the algorithm
    for (int i=int(allSnps.size())-1 ; i>=0 ; i--)
    {
        if (find(bests.begin(), bests.end(), CChromosome(i, 0)) == bests.end() &&
			find(bests.begin(), bests.end(), CChromosome(i, 1)) == bests.end() &&
			find(bests.begin(), bests.end(), CChromosome(i, 2)) == bests.end()   )
        {
            snpsToRemove.push_back(i);
        }
    }
    
	cout << "Removed SNPs:" << endl;
	for (auto snp : snpsToRemove)
	{
		cout << ioGenotypeData.TotalCase()[snp].mIndex << ",";
	}
	cout << endl;
	
	ioGenotypeData.RemoveSNP(snpsToRemove);
}

void CClusteringAlg::CalcBestBarcode(unsigned int iOrder, const CChromosome& iChromsome, unsigned int iCurrIndex, CChromosome*& oBestChromosome, const CGenotypeData& ioGenotypeData) const
{
	CChromosome*  tmpChr = nullptr;
	if (CConfigParams::GetTheInstance().Algorithm == CConfigParams::eMooney)
	{
		tmpChr = new CChromosomeMooney(iChromsome);
	}
	else
	{
		tmpChr = new CChromosomeYang(iChromsome);
	}
	
	for (unsigned int j = 0; j < CGenotypeData::NUMBER_OF_GENOTYPES; j++)
	{
		tmpChr->Genotype(iCurrIndex, j);
		if (iCurrIndex == iOrder - 1)
		{
			double fitness = tmpChr->Fitness(ioGenotypeData.CaseData(), ioGenotypeData.ControlData());
			if (oBestChromosome == nullptr || fitness > oBestChromosome->Fitness())
			{
				delete oBestChromosome;
				if (CConfigParams::GetTheInstance().Algorithm == CConfigParams::eMooney)
				{
					oBestChromosome = new CChromosomeMooney(*tmpChr);
				}
				else
				{
					oBestChromosome = new CChromosomeYang(*tmpChr);
				}
			}
		}
		else
		{
			CalcBestBarcode(iOrder, *tmpChr, iCurrIndex + 1, oBestChromosome, ioGenotypeData);
		}
	}

	delete tmpChr;
}
