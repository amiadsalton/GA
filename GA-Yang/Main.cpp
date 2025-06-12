#pragma once 

#include "GeneticAlgorithm.h"
#include "ConfigParams.h"

#include <iostream>
#include <fstream>
#ifdef _MSC_VER
#include <windows.h>
#else
#include <unistd.h>
#endif
#include <time.h>
#include <omp.h>
#include <thread>
#include <boost/filesystem.hpp>

using namespace std;

CGeneticAlgorithm  sAlgo;
bool sMainThreadIsRunning(true);

void RunOneAlgo(int iAlgo, time_t iTimeStamp)
{
    cout << "\n\n*********************\nAlgorithm " << iAlgo << " Start Running\n*********************\n\n";

    CConfigParams::GetTheInstance().PostProcessingAlg = CConfigParams::EPostProcessingAlgorithm(iAlgo);

    boost::filesystem::path outputPath = CConfigParams::GetTheInstance().DirectoryName;
    outputPath /= "Run_";
    outputPath += to_string(iTimeStamp).c_str();
    outputPath += "_";
    outputPath += to_string(int(CConfigParams::GetTheInstance().Algorithm)).c_str();
    outputPath += "_";
    outputPath += to_string(int(CConfigParams::GetTheInstance().PostProcessingAlg)).c_str();

    if (sAlgo.Init(outputPath.string()))
    {
        // Copy the parameters file
        if (sAlgo.RunAllOrders())
        {           
            cout << "Algorithm execution successfully completed.";
        }
        else
        {
            cout << "Algorithm execution failed.";
        }
    }
    else
    {
        cout << "Algorithm initialization failed";
    }

}

void MainThread()
{
    time_t timeStamp;
    time(&timeStamp);
#ifdef _MSC_VER
        	Sleep(1000);
#else
        	sleep(1);
#endif

    if (sAlgo.Load())
    {
        if (CConfigParams::GetTheInstance().PostProcessingAlg == CConfigParams::eNumberOfPostProcessingAlgorithms)
        {
			RunOneAlgo(CConfigParams::eClusteringFitenss, timeStamp);
			RunOneAlgo(CConfigParams::eClusteringPosition, timeStamp);
			RunOneAlgo(CConfigParams::eClusteringSum, timeStamp);
			RunOneAlgo(CConfigParams::eContinuity, timeStamp);
        }
        else 
        {
            RunOneAlgo(CConfigParams::GetTheInstance().PostProcessingAlg, timeStamp);
        }
    }

    cout << "\nEnter any key to exit  ==>" << endl;

    sMainThreadIsRunning = false;
}

int main(int argV, char** argC)
{
    cout << "Genetic Algorithm execution\n\n";

    //omp_set_num_threads(4);


    thread mainThread(MainThread);

    unsigned int ch(CGeneticAlgorithm::eNone);
    do
    {
        // Wait until previous request was performed
        cout << "Please enter:\n1 - Exit\n2 - Pause\n3 - Resume\n4 - Stop Execution for current Run\n5 - Stop Generation for current Execution\n\n";
        cout << "==>\n";
        cin >> ch;
        if (ch >= (unsigned int)(CGeneticAlgorithm::eLast))
        {
            ch = CGeneticAlgorithm::eNone;
        }
        if (sMainThreadIsRunning == false)
        {
            break;
        }
        sAlgo.SetUserInput(CGeneticAlgorithm::EUserInput(ch));
        while (sAlgo.GetUserInput() != CGeneticAlgorithm::eNone && sAlgo.GetUserInput() != CGeneticAlgorithm::ePause)
        {
            // Wait until sAlgo analyzes user request
            cout << ".";
#ifdef _MSC_VER
        	Sleep(1000);
#else
        	sleep(1);
#endif
        }
    }
    while(ch != CGeneticAlgorithm::eExit);

    mainThread.join();

    return 1;
}