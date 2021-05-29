#include <iostream>

#include "boost/numeric/interval.hpp"
#include "boost/numeric/interval/io.hpp"
#include "Knapsack.h"
#include "ProblemFactory.h"
#include "BenchmarkManager.h"




void KnapsackPreprocess(std::vector<int>& weights, std::vector<int>& values)
{
    std::vector<BnB::Knapsack::Detail::item> rem;
    for(int i = 0; i < weights.size(); i++){
        BnB::Knapsack::Detail::item it;
        it.wt = weights[i];
        it.cost = values[i];
        it.ratio = static_cast<double>(it.cost) / static_cast<double>(it.wt);
        rem.push_back(it);
    }
    std::sort(rem.begin(),rem.end(),BnB::Knapsack::Detail::comp); // decreasing order of ratio of cost to weight;

    weights.clear();
    values.clear();

    for(int i = 0; i < rem.size(); i++)
    {
        weights.push_back(rem[i].wt);
        values.push_back(rem[i].cost);
    }
}




void RunSerialMeasurementsRandom(int argc, char* argv[], int mainIter, int subIter){
    BenchmarkManager manager;
    std::vector<double> runtimes;
    std::vector<double> standardDeviation;
    for(int i = 0;i < mainIter; i++)
    {
        auto [weights, values, W] = std::string(argv[argc - 2]) == "Corr" ? BnB::Knapsack::GenerateCorrelatedProblemConstants(std::pair(1,100),atoi(argv[argc - 1]), 55555,0.10,  0.5)
                                                                          : BnB::Knapsack::GenerateRandomProblemConstants(std::pair(1,100),atoi(argv[argc - 1]), 6413, 0.5);
        KnapsackPreprocess(weights, values);
        manager.SetKConsts(weights, values, W);
        auto result = manager.Benchmark(subIter, ProblemClass::KNAPSACK, argv, argc);
        std::cout << "it took " << result.first << " seconds" << std::endl;
        runtimes.push_back(result.first);
        standardDeviation.push_back(result.second);
    }
    if(mainIter  > 1)
        std::cout << "runtime = " << mean(runtimes) << " +- " << standardDev(runtimes, mean(runtimes)) << std::endl;
    else
        std::cout << "runtime = " << runtimes[0] << "+-" << standardDeviation[0] << std::endl;
}



void RunSerialMeasurements(int argc, char* argv[], std::string FileName, int mainIter, int subIter){
    BenchmarkManager manager;
    std::vector<double> runtimes;
    std::vector<double> standardDeviation;
    std::ifstream file(FileName);
    for(int i = 0;i < mainIter; i++)
    {
        if(file.peek() == EOF) break;
        auto [weights, values, W] = ReadKnapsackProblem(file);
        KnapsackPreprocess(weights, values);
        manager.SetKConsts(weights, values, W);
        auto result = manager.Benchmark(subIter, ProblemClass::KNAPSACK, argv, argc);
	    runtimes.push_back(result.first);
        std::cout << "it took " << result.first << " seconds iter: " << i << std::endl;
	    standardDeviation.push_back(result.second);
    }
    if(mainIter  > 1)
    	std::cout << "runtime = " << mean(runtimes) << " +- " << standardDev(runtimes, mean(runtimes)) << std::endl;
    else
	std::cout << "runtime = " << runtimes[0] << "+-" << standardDeviation[0] << std::endl;
}

void GatherTreeData(int argc, char* argv[], int seed){
    if(argv[argc - 2] == "Corr")
        std::cout << "we are solving a correlated problem" << std::endl;
    else if(argv[argc - 2] == "UnCorr")
        std::cout << "we are solving a uncorrelated problem" << std::endl;

    BenchmarkManager manager;
    auto [weights, values, W] = std::string(argv[argc - 2]) == "Corr" ? BnB::Knapsack::GenerateCorrelatedProblemConstants(std::pair(1,100),atoi(argv[argc - 1]), seed, 0.05)
                                                                      : BnB::Knapsack::GenerateRandomProblemConstants(std::pair(1,100),atoi(argv[argc - 1]), seed);
    KnapsackPreprocess(weights, values);
    manager.SetKConsts({30,50,40,10,40,30,10},{60,60,40,10,20,10,3} ,100);
    manager.GatherTreeData(ProblemClass::KNAPSACK, argv, argc);
}


void GenerateLongProblems(int argc, char* argv[], std::string filename){
    BenchmarkManager manager;
    for(int i = 0;i < 1000; i++){
        auto [weights, values, W] = std::string(argv[argc - 2]) == "Corr" ? BnB::Knapsack::GenerateCorrelatedProblemConstants(std::pair(1,100),atoi(argv[argc - 1]), 12341213,0.05,  0.5)
                                                                          : BnB::Knapsack::GenerateRandomProblemConstants(std::pair(1,100),atoi(argv[argc - 1]), 6413, 0.0);
        KnapsackPreprocess(weights, values);
        manager.SetKConsts(weights, values, W);
        double result = manager.Benchmark(1, ProblemClass::KNAPSACK, argv, argc).first;
        std::cout << "it took " << result << " seconds" << std::endl;
        if(result > 20 && result < 50 )
           SaveKnapsackProblem(weights,values,W, filename);
    }
}





//****************************** main ******************************//
int main(int argc, char* argv[]) {
    ///RunSerialMeasurements(argc, argv, "Knapsack95Corr.txt", 1, 1);
    //RunSerialMeasurementsRandom(argc, argv, 100, 1);
    //GeerateLongProblems(argc, argv, "Knapsack95Corr.txt");
    //GatherTreeData(argc, argv, 5122111);


    std::vector<interval> domain;
    SetProblemParams(domain);

    BenchmarkManager manager;

    ///auto [weights, values, W] =  BnB::Knapsack::GenerateCorrelatedProblemConstants(std::pair(1,100), 120, 1, 0.05);
    manager.SetKConsts({30,50,40,10,40,30,10},{60,60,40,10,20,10,3},100);
	///for(int i = 0; i < 1; i++){
    manager.SetIConsts(domain);
    ///manager.SetTConsts(BnB::TSP::GenerateRandomMatrix(12, {10,20}, i*10));
    //manager.SetLConsts(A,b,c,v);

    auto result = manager.Benchmark(1,ProblemClass::INTERVAL, argv, argc);
    std::cout << "avarage runtime " << result.first << " +- " << result.second<< std::endl;

    auto result2 = manager.Benchmark(1,ProblemClass::KNAPSACK, argv, argc);
    std::cout << "avarage runtime " << result2.first << " +- " << result2.second<< std::endl;


	//}
    ///manager.GatherTreeData(ProblemClass::INTERVAL, argv, argc);


  return 0;
}
