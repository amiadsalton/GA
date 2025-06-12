
#include "ConfigParams.h"
#include <sstream>
#include <boost/token_functions.hpp>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>


unique_ptr<CConfigParams>  CConfigParams::sTheInstance;
const string CConfigParams::FILE_NAME="ConfigParams.txt";

CConfigParams::CConfigParams() :
    ChromosomeLength(3),
    NumOfIterations(100),
    SelectedBarcode(),
    HaltCriteria(0)
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
        DirectoryName = ReadNextToken(file, "Dir");
        CaseFileName = ReadNextToken(file, "Case");
        ControlFileName = ReadNextToken(file, "Control");

        ChromosomeLength = std::stoul(ReadNextToken(file, "Length"));

        NumOfIterations = std::stoul(ReadNextToken(file, "Iterations"));

        PopulationSize  = std::stoul(ReadNextToken(file, "Size"));

        Provider = EGenotypeDataProvider(std::stoul(ReadNextToken(file, "Provider")));

        SelectionAlgorithm = ESelectionAlgorithm(std::stoul(ReadNextToken(file, "SelectionAlgorithm")));

        ElitismRate = std::stof(ReadNextToken(file, "Elitism"));

        CrossoverRate = std::stof(ReadNextToken(file, "Crossover"));

        MutationRate = std::stof(ReadNextToken(file, "Mutation"));

        TrapRatio = std::stof(ReadNextToken(file, "TrapRatio"));

        VibrationRate = std::stof(ReadNextToken(file, "Vibration"));

        Algorithm = EAlgorithm(std::stoul(ReadNextToken(file, "Algorithm")));

        DisplayRatio = std::stof(ReadNextToken(file, "DisplayRatio"));

        string token=ReadNextToken(file, "Barcode");
        boost::char_separator<char> sep(",");
        boost::tokenizer<boost::char_separator<char>> tokensBoost(token, sep);
        BOOST_FOREACH (const string& t, tokensBoost)
        {
            SelectedBarcode.push_back(t);
        }

        HomogeneousRatio = std::stof(ReadNextToken(file, "HomogeneousRatio"));
        
        IgnoreGenotype2 = std::stoul(ReadNextToken(file, "IgnoreGenotype2")) == 1 ? true : false;

        NumberOfExecutions = std::stoul(ReadNextToken(file, "NumberOfExecutions"));

        ExecutionIsStcuk = std::stoul(ReadNextToken(file, "ExecutionIsStcuk"));

        NumberOfOrders = std::stoul(ReadNextToken(file, "NumberOfOrders"));

        PostProcessingAlg = EPostProcessingAlgorithm(std::stoul(ReadNextToken(file, "PostProcessing")));

        HaltCriteria  = std::stoul(ReadNextToken(file, "HaltCriteria"));
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
        string  value = line.substr(keyword.length() + 1, line.length() - keyword.length() -1);
        if (!value.empty() && value[value.size()-1] == '\r')
        {
        	value.resize(value.size()-1);
        }
        return value;
    }    

    cout << "Can't find " << keyword << " in config file\n";
    return "";
}
