#pragma once

#include "BenchmarkManager.h"


enum class ProblemClass {KNAPSACK, INTERVAL, TSP};
enum class Parallelization {MPI, OMP, SERIAL};

struct Parameters{
    BnB::MPI_Scheduler_Type scheduler_MPI;
    BnB::OMP_Scheduler_Type scheduler_OMP;
    double eps;
    BnB::TraversalMode mode = BnB::TraversalMode::DFS;
    BnB::Goal goal = BnB::Goal::MAX;
    int Communication_Frequency;
    int OpenMPThreads;
    int PackSize;
    int TermFreq;
    Parallelization ParallelMode;
};

void ParseParallelMode(Parameters& params, std::string arg)
{
    if(arg == "MPI")
        params.ParallelMode = Parallelization::MPI;
    else if(arg == "OMP")
        params.ParallelMode = Parallelization::OMP;
    else if(arg == "Serial")
        params.ParallelMode = Parallelization::SERIAL;

}

void ParseScheduler(Parameters& params, std::string arg)
{
    if (arg == "Prio")
        params.scheduler_MPI= BnB::MPI_Scheduler_Type::PRIORITY;
    else if( arg == "Hybrid")
        params.scheduler_MPI = BnB::MPI_Scheduler_Type::HYBRID;
    else if(arg == "One")
        params.scheduler_MPI = BnB::MPI_Scheduler_Type::ONESIDED;
    else if(arg == "Queue")
        params.scheduler_OMP = BnB::OMP_Scheduler_Type::QUEUE;
    else if(arg == "Task")
        params.scheduler_OMP = BnB::OMP_Scheduler_Type::TASKING;
    else
        assert(false);
}

void ParseOpt(Parameters& params, std::string arg)
{
    if(arg == "Max")
        params.goal = BnB::Goal::MAX;
    else if(arg == "Min")
        params.goal = BnB::Goal::MIN;
    else
        assert(false);
}

void ParseTraversal(Parameters& params, std::string arg)
{
    if(arg == "DFS")
        params.mode = BnB::TraversalMode::DFS;
    else if(arg == "BFS")
        params.mode = BnB::TraversalMode::BFS;
    else if(arg =="auto")
        params.mode = BnB::TraversalMode::AUTOMATIC;
    else
        assert(false);
}


void DisplayHelpInfo()
{
    std::cout << "The following parameters can be used" << std::endl;
    std::cout << "--Schedueler [Prio]/[Hybrid]/[One]/[Task]/[Queue]" << std::endl;
    std::cout << "--Traversal  [DFS]/[BFS]/[auto]" << std::endl;
    std::cout << "--Parallel   [MPI]/[OMP]/[Serial]" << std::endl;
    std::cout << "--Threads    [integer value]" << std::endl;
    std::cout << "--TermFreq   [integer value]" << std::endl;
    std::cout << "--PackSize   [integer value]" << std::endl;
    std::cout << "--Freq       [integer value]" << std::endl;
    std::cout << "--Eps        [numeric value]" << std::endl;
    std::cout << "--Opt        [Max, Min]" << std::endl;
}
