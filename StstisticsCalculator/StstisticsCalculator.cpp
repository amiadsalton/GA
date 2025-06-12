// StstisticsCalculator.cpp : Defines the entry point for the console application.
//

#include <windows.h>

#include <string>
#include <iostream>
#include <math.h>

inline double Factorial(int min, int max, int divider)
{
    double result = 1.0;
    for  (int i = min+1 ; i <= max ; i++)
    {
        result = result *  i / divider;
    }    
    return result;
}

int main(int argc, char* argv[])
{
    int n = std::stoul(argv[1], nullptr, 10);
    int m = std::stoul(argv[2], nullptr, 10);
    int f = std::stoul(argv[3], nullptr, 10);
    int p = std::stof(argv[4], nullptr) * n;

    double result = Factorial((p-f) , p ,n);
    double tmp = Factorial(n-p-(p-f), n-p,  n);
    result *= tmp;

//     double result = 1.0;
//     for  (int i=0 ;i<f ;i++)
//     {
//         result = result *  ((double(p)-i-1) / n);
//     }
// 
//     for  (int i=0 ;i<(p-f) ;i++)
//     {
//         result = result *  ((double(n)-p-(p-i-1)) / n);
//     }

    std::cout << " Result: " << result << std::endl;

//     char  ch;
//     std::cin >> ch;

	return 0;
}

