#pragma once

#include <chrono>
#include <Knapsack.h>
#include <numeric>
#include <cmath>
#include <BnB_MPI_Solver.h>
#include <TSP.h>
#include "BnB_OMP_Solver.h"
#include "BnB_Serial_Solver.h"
#include "Base.h"
#include "function.h"
#include "ProblemFactory.h"
#include "ArgumentParser.h"
#include "TreeGenerator.h"
#include "TreeGeneratorKnap.h"


class BenchmarkManager{
public:
    ~BenchmarkManager(){if(IsMPIInitialized) MPI_Finalize();}
    // benchmarks run until we reach a stable std
    std::pair<double,double> Benchmark(int iterations, ProblemClass problem, char** argc, int argv);
    void GatherTreeData(ProblemClass problem, char** argc, int argv);

    void SetPercentageOfSTD(float percent) { PercentageOfStandardDeviation = percent;}
    void SetKConsts(std::vector<int> weights, std::vector<int> costs, int Limit);
    void SetTConsts(const std::vector<std::vector<int>>& Matrix);
    void SetIConsts(std::vector<interval> domain);

private:
    void ParseArguments(char** args, int num);
    void PrintConfig(ProblemClass);

    template<typename P, typename S, typename Type> std::pair<double, double> BenchmarkLoop(int iterations, std::shared_ptr<BnB::Solver<P, S, Type>> solver, BnB::Problem_Definition<P,S,Type> Def, P Consts);
    template<typename P, typename S, typename Type> std::shared_ptr<BnB::Solver<P, S, Type>> GetSolver(Parallelization parallel);

    float PercentageOfStandardDeviation = 0.1f;
    Parameters params;
    BnB::Knapsack::Consts KConsts;
    BnB::TSP::Consts TConsts;
    Problems::IntervalConsts IConsts;

    bool IsMPIInitialized = false;
};

template<typename T>
T mean(const std::vector<T>& numbers){ return std::accumulate(numbers.begin(), numbers.end(), (T)0)/ numbers.size();}

template<typename T>
double standardDev(const std::vector<T>& numbers, T mean){
    if(numbers.size() == 1)
        return 1000000.0;
    T result = 0;
    for(const auto& num : numbers)
    {
        result += std::pow(num - mean, 2);
    }
    return std::sqrt(result/(numbers.size()-1));
}


void BenchmarkManager::SetKConsts(std::vector<int> weights, std::vector<int> costs, int Limit){
    std::get<0>(KConsts) = weights;
    std::get<1>(KConsts) = costs;
    std::get<2>(KConsts) = Limit;
}
void BenchmarkManager::SetIConsts(std::vector<interval> domain){
    std::get<0>(IConsts) = domain;
}


void BenchmarkManager::SetTConsts(const std::vector<std::vector<int>> &Matrix){
    std::get<0>(TConsts) = Matrix;
}

template<typename T, typename P>
int GetTotalValue(T& t, const P& p){};

template<>
int GetTotalValue(BnB::Knapsack::Params& t,const BnB::Knapsack::Consts& C)
{
    /*int i = std::get<0>(t).size();
    if(std::get<1>(t) + *(std::get<0>(C).end()-1) <= std::get<2>(C))
    {
        std::get<1>(t) += *(std::get<0>(C).end() - 1);
        std::get<2>(t) += *(std::get<1>(C).end()-1);
        std::get<0>(t).push_back(1);
    }else{
        std::get<0>(t).push_back(0);
    }*/

    // perform greedy to finish it
    int w = std::get<1>(t);
    int cost = std::get<2>(t);

    float bound = cost;
    for (int j = std::get<0>(t).size(); j < std::get<0>(C).size(); j++) {
        if (w == std::get<2>(C))
            break;
        else if (w + std::get<0>(C)[j] <= std::get<2>(C)) {
            bound += std::get<1>(C)[j];
            w += std::get<0>(C)[j];
        }
    }
    return std::get<3>(t);
}



template<typename P, typename S, typename Type>
std::pair<double, double> BenchmarkManager::BenchmarkLoop(int iterations, std::shared_ptr<BnB::Solver<P,S,Type>> solver, BnB::Problem_Definition<P,S,Type> Def, P Consts)
{
    using namespace std::chrono;
    std::vector<double> runtimes;
    double TotalRuntime = 0;
    double StdMeanRatio = 1.0;
    time_point<system_clock> start, end;
    int counter = 0;
    while(counter < iterations)
    {
        counter++;
        start = system_clock::now();

        if(params.goal == BnB::Goal::MAX)
        {
            auto result = solver->Maximize(Def, Consts);
            if(params.ParallelMode == Parallelization::MPI)
            {
                int rank;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                if(rank != 0)
                    continue;
            }
        }
        else{
            auto result = solver->Minimize(Def, Consts);
            std::cout << "corrected result " << Def.GetContainedUpperBound(Consts, result) << std::endl;
        }
        end = system_clock::now();

        runtimes.push_back(duration_cast<milliseconds>(end - start).count()/1000.0);
        TotalRuntime += duration_cast<milliseconds>(end - start).count()/1000.0;
        StdMeanRatio = static_cast<double>(standardDev(runtimes, mean(runtimes))) / static_cast<double>(mean(runtimes));
    }
    return std::pair<double, double>(mean(runtimes), standardDev(runtimes, mean(runtimes)));
}


template<typename P, typename S, typename Type>
std::shared_ptr<BnB::Solver<P,S,Type>> BenchmarkManager::GetSolver(Parallelization parallel)
{
    switch(parallel)
    {
        case Parallelization::MPI:
        {
            auto solver = std::make_shared<BnB::Solver_MPI<P, S, Type>>();
            solver->SetScheduler(params.scheduler_MPI);
            solver->SetSchedulerParameters()
                    ->CommFrequency(params.Communication_Frequency)
                    ->Eps(params.eps)
                    ->TraversMode(params.mode)
                    ->MaximalPackageSize(params.PackSize);
            static_cast<BnB::MPI_Scheduler_Hybrid<P,S,Type>*>(solver->SetSchedulerParameters())->Threads(params.OpenMPThreads);
            static_cast<BnB::MPI_Scheduler_WorkerOnly<P,S,Type>*>(solver->SetSchedulerParameters())->TermCheckFrequency(params.TermFreq);
            return solver;
        }
        case Parallelization::OMP:
        {
            auto solver = std::make_shared<BnB::Solver_OMP<P, S, Type>>();
            solver->SetNumThreads(params.OpenMPThreads);
            solver->SetScheduler(params.scheduler_OMP);
            solver->SetSchedulerParameters()
                    ->Eps(params.eps)
                    ->Traversal(params.mode);
            return solver;
        }
        case Parallelization::SERIAL:
            auto solver = std::make_shared<BnB::Solver_Serial<P, S, Type>>();
            solver->SetSchedulerParameters()
                    ->Eps(params.eps)
                    ->TraversMode(params.mode);
            return solver;
    }
}


std::pair<double, double> BenchmarkManager::Benchmark(int iterations, ProblemClass problem, char** argc, int argv)
{
    ParseArguments(argc, argv);
    if(params.ParallelMode == Parallelization::MPI)
    {
        if(!IsMPIInitialized){
            int provided;
            MPI_Init_thread(0,NULL, MPI_THREAD_SERIALIZED, &provided);
            IsMPIInitialized = true;
        }
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(rank == 0)
            PrintConfig(problem);
    } else
        PrintConfig(problem);


    if(problem == ProblemClass::KNAPSACK)
    {
        auto ProblemDef = BnB::Knapsack::GenerateFasterProblem(this->params.mode == BnB::TraversalMode::BFS);
        auto solver = GetSolver<BnB::Knapsack::Consts , BnB::Knapsack::Params ,float>(params.ParallelMode);
        return BenchmarkLoop<BnB::Knapsack::Consts, BnB::Knapsack::Params, float>(iterations, solver, ProblemDef, KConsts);
    }
    else if(problem == ProblemClass::INTERVAL)
    {
        auto ProblemDef =  Problems::GenerateIntervalProblem(params.goal);
        auto solver = GetSolver<Problems::IntervalConsts , Problems::IntervalParams , double>(params.ParallelMode);
        return BenchmarkLoop<Problems::IntervalConsts , Problems::IntervalParams , double>(iterations, solver, ProblemDef,IConsts);
    }else if(problem == ProblemClass::TSP)
    {
        auto ProblemDef = BnB::TSP::GenerateToyProblem();
        auto solver = GetSolver<BnB::TSP::Consts, BnB::TSP::Params, int>(params.ParallelMode);
        return BenchmarkLoop<BnB::TSP::Consts, BnB::TSP::Params, int>(iterations, solver, ProblemDef, TConsts);
    }
}


void BenchmarkManager::GatherTreeData(ProblemClass problem, char** argc, int argv)
{
    ParseArguments(argc, argv);
    if(params.ParallelMode == Parallelization::MPI)
    {
        int provided;
        MPI_Init_thread(0,NULL, MPI_THREAD_SERIALIZED, &provided);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(rank == 0)
            PrintConfig(problem);
    } else
        PrintConfig(problem);


    if(problem == ProblemClass::KNAPSACK)
    {
        auto ProblemDef = TreeGenerator::Knapsack::GenerateFasterProblem();
        auto solver = GetSolver<TreeGenerator::Knapsack::Consts, TreeGenerator::Knapsack::Params, float>(params.ParallelMode);
        if(params.goal == BnB::Goal::MAX)
            auto result = solver->Maximize(ProblemDef, KConsts);
        else
            auto result = solver->Minimize(ProblemDef, KConsts);
    }
    else if(problem == ProblemClass::INTERVAL)
    {
        auto ProblemDef =  TreeGenerator::Interval::GenerateIntervalProblem(params.goal);
        auto solver = GetSolver<TreeGenerator::Interval::Consts , TreeGenerator::Interval::Params , double>(params.ParallelMode);
        if(params.goal == BnB::Goal::MAX)
            auto result = solver->Maximize(ProblemDef, IConsts);
        else
            auto result = solver->Minimize(ProblemDef, IConsts);
    }

}


void BenchmarkManager::ParseArguments(char** args, int num)
{
    if(num < 3)
    {
        DisplayHelpInfo();
        assert(false);
    }

    for(int i = 1; i < num - 1; i+= 2)
    {
        std::string flag = args[i];
        std::string arg = args[i + 1];
        if(flag == "--Scheduler"){
           ParseScheduler(params, arg);
        }else if(flag == "--Eps"){
            params.eps = std::stod(arg);
        }else if(flag == "--Opt"){
            ParseOpt(params, arg);
        }else if(flag == "--Freq") {
            params.Communication_Frequency = std::stoi(arg);
        }else if(flag == "--Threads"){
            params.OpenMPThreads = std::stoi(arg);
        }else if(flag == "--Traversal"){
            ParseTraversal(params, arg);
        }else if(flag == "--Parallel"){
            ParseParallelMode(params, arg);
        }else if(flag == "--TermFreq") {
            params.TermFreq = std::stoi(arg);
        }else if(flag == "--PackSize") {
            params.PackSize = std::stoi(arg);
        }else if(flag == "--Other"){
                return;
        }else{
            DisplayHelpInfo();
            assert(false);
        }
    }
}


void BenchmarkManager::PrintConfig(ProblemClass p)
{
    std::string opt = (bool)params.goal ? "maximizing" : "minimizing";
    std::string prob;
    if(p == ProblemClass::KNAPSACK)
       prob = "Knapsack";
    else if(p == ProblemClass::INTERVAL)
       prob = "Interval";
    else if(p == ProblemClass::TSP)
        prob = "Traveling salesman problem";
    std::string para;
    if(params.ParallelMode== Parallelization::MPI)
        para = "parallelizing using MPI";
    else if(params.ParallelMode == Parallelization::OMP)
        para = "parallelizing using OMP";
    else if(params.ParallelMode == Parallelization::SERIAL)
        para = "not parallelizing";
    std::string search = params.mode == BnB::TraversalMode::DFS ? "depth first search" : "breadth first search";


    std::cout << "--------------------------------------------------------------------------" << std::endl;
    std::cout << "We are " << opt << " a " << prob << " problem with following parameters" << std::endl;
    std::cout << "We are " << para << " using " << search << std::endl;
    std::cout << "The precision of the optimization is " << params.eps << " and communcation frequency is " << params.Communication_Frequency << std::endl;
    std::cout << "--------------------------------------------------------------------------" << std::endl;
}
