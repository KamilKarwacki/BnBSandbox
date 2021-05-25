#pragma once
#include "Base.h"
#include <bitset>
#include <fstream>
#include "boost/numeric/interval.hpp"
#include "boost/numeric/interval/io.hpp"
#include "function.h"

///#define GET_RANK(id) MPI_Comm_rank(MPI_COMM_WORLD, &id);
#define GET_RANK(id) id = omp_get_thread_num();


// functionalities to write tree data out while solving

// ------------------------------------------------- USER STUFF --------------------------------------------------------
namespace TreeGenerator::Interval {
    std::ofstream TreeDataOutput;
    int counter = 1;
    using Consts = BnB::Problem_Constants<
            std::vector<interval>>;

    using Params = BnB::Subproblem_Parameters<
            std::vector<interval>, // domain
            int, // parent
            int, // id
            int>;// depth



    Consts GenerateConsts() {
        Consts consts;
        SetProblemParams(std::get<0>(consts));
        return consts;
    }
// ----------------------------------------------

    std::function<std::vector<Params>(const Consts &, const Params &)>
    GeneratePrioritySplitFunction() {
        return [](const Consts &consts, const Params &params) {
            int id;
            GET_RANK(id);
            std::vector<Params> ret;
            auto SubDomain = std::get<0>(params);

            interval y;
            size_t n_split = 0;
            std::vector<interval> new_domain(SubDomain.size());
            std::vector<std::vector<interval> > split_space(SubDomain.size(), std::vector<interval>(2));
            for (size_t i = 0; i < SubDomain.size(); i++) { // fill split space
                n_split++;
                split_space[i][0].set(SubDomain[i].lower(), median(SubDomain[i]));
                split_space[i][1].set(median(SubDomain[i]), SubDomain[i].upper());
            }
            if (n_split > 0) {
                for (size_t c = 0; c < pow(2, n_split); c++) { // go through all possible splits
                    std::bitset<64> split_comb(c);
                    size_t i_split = 0;
                    for (size_t i = 0; i < SubDomain.size(); i++) {
                        new_domain[i] = split_space[i][split_comb[i_split]];
                        i_split++;
                    }
                    Params sol;
                    std::get<0>(sol) = new_domain;
                    std::get<1>(sol) = std::get<2>(params);  // parent
                    std::get<3>(sol) = counter++; // unique id
                    std::get<2>(sol) = 10000*id + std::get<3>(sol); // id scaled by proc
                    ret.push_back(sol);
                }
            }
            return ret;
        };
    }

    BnB::Problem_Definition<Consts, Params, double> GenerateIntervalProblem(BnB::Goal goal) {
        BnB::Problem_Definition<Consts, Params, double > Definition;
        Definition.SplitSolution = GeneratePrioritySplitFunction();

        if (goal == BnB::Goal::MIN) {
            Definition.GetEstimateForBounds = [](const Consts &consts, const Params &params) {
                int id;
                GET_RANK(id);
                std::ofstream TreeDataOutput;
                TreeDataOutput.open("output" + std::to_string(id) + ".txt", std::ios_base::app);
                TreeDataOutput << std::get<2>(params) <<";" << std::get<1>(params) << ";" << id << ";" << 0 << std::endl;
                interval inter = objective(std::get<0>(params));
                set_nan_to_inf(inter);
                TreeDataOutput.close();
                return std::tuple<double,double>(inter.lower(), inter.upper());
            };
        } else {
            Definition.GetEstimateForBounds = [](const Consts &consts, const Params &params) {
                int id;
                GET_RANK(id);
                std::ofstream TreeDataOutput;
                TreeDataOutput.open("output" + std::to_string(id) + ".txt", std::ios_base::app);
                TreeDataOutput << std::get<2>(params) <<";" << std::get<1>(params) << ";" << id << ";" << 0 << std::endl;
                interval inter = objective(std::get<0>(params));
                set_nan_to_inf(inter);
                TreeDataOutput.close();
                return std::tuple<double,double>(inter.upper(), inter.lower());
            };
        }

        Definition.GetContainedUpperBound = [](const Consts &consts, const Params &params) {
            auto SubDomain = std::get<0>(params);
            std::vector<double> x_mid(SubDomain.size());
            for (size_t i = 0; i < SubDomain.size(); i++)
                x_mid[i] = median(SubDomain[i]);
            return objective(x_mid);
        };

        Definition.IsFeasible = [](const Consts &consts, const Params &params) {
            return BnB::FEASIBILITY::Full;
        };

        Definition.PrintSolution = [](const Params &params) {
            std::cout << "The solution is in the interval:" << std::endl;
            auto Domain = std::get<0>(params);
            for (const auto &interval : Domain) {
                std::cout << "[" << std::setprecision(10) << interval.lower() << ", " << std::setprecision(10) << interval.upper() << "]" << std::endl;
            }
        };

        Definition.GetInitialSubproblem = [](const Consts &prob) {
            Params Initial;
            std::get<0>(Initial) = std::get<0>(prob);
            return Initial;
        };
        return Definition;
    }
}
