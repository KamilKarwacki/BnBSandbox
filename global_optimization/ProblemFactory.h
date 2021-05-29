#pragma once
#include "Base.h"
#include <bitset>
#include "boost/numeric/interval.hpp"
#include "boost/numeric/interval/io.hpp"
#include "function.h"
#include <iomanip>

// ------------------------------------------------- USER STUFF --------------------------------------------------------
namespace Problems {
    using IntervalConsts = BnB::Problem_Constants<
            std::vector<interval>>; // initial domain

    using IntervalParams = BnB::Subproblem_Parameters<
            std::vector<interval>>;


    IntervalConsts GenerateIntervalConsts() {
        IntervalConsts consts;
        SetProblemParams(std::get<0>(consts));
        return consts;
    }


// ----------------------------------------------
    BnB::Problem_Definition<IntervalConsts, IntervalParams, double> GenerateIntervalProblem(BnB::Goal goal) {
        BnB::Problem_Definition<IntervalConsts, IntervalParams, double> Definition;


        Definition.SplitSolution = [](const IntervalConsts &consts, const IntervalParams &params) {
            std::vector<IntervalParams> ret;
            auto SubDomain = std::get<0>(params);
            //const auto eps_x = std::get<1>(consts);
            std::vector<bool> x_mask(SubDomain.size(), true);
            interval y;
            //for (size_t i = 0; i < SubDomain.size(); i++)
             //   if (width(SubDomain[i]) <= eps_x)
              //      x_mask[i] = false;

            size_t n_split = 0;
            std::vector<interval> new_domain(SubDomain.size());
            std::vector<std::vector<interval> > split_space(SubDomain.size(), std::vector<interval>(2));
            for (size_t i = 0; i < SubDomain.size(); i++) { // fill split space
                if (x_mask[i]) {
                    n_split++;
                    split_space[i][0].set(SubDomain[i].lower(), median(SubDomain[i]));
                    split_space[i][1].set(median(SubDomain[i]), SubDomain[i].upper());
                } else {
                    split_space[i].resize(1);
                    new_domain[i] = SubDomain[i];
                }
            }
            if (n_split > 0) {
                for (size_t c = 0; c < pow(2, n_split); c++) { // go through all possible splits
                    std::bitset<64> split_comb(c);
                    size_t i_split = 0;
                    for (size_t i = 0; i < SubDomain.size(); i++) {
                        if (x_mask[i]) {
                            new_domain[i] = split_space[i][split_comb[i_split]];
                            i_split++;
                        }
                    }
                    IntervalParams sol;
                    std::get<0>(sol) = new_domain;
                    ret.push_back(sol);
                }
            }
            return ret;
        };

        if (goal == BnB::Goal::MIN) {
            Definition.GetEstimateForBounds = [](const IntervalConsts &consts, const IntervalParams &params)
            {
                interval inter = objective(std::get<0>(params));
                set_nan_to_inf(inter);
                return std::tuple<double, double>(inter.lower(), inter.upper());
            };
        } else {
            Definition.GetEstimateForBounds = [](const IntervalConsts &consts, const IntervalParams &params)
            {
                interval inter = objective(std::get<0>(params));
                set_nan_to_inf(inter);
                return std::tuple<double, double>(inter.upper(), inter.lower());
            };
        }

        Definition.GetContainedUpperBound = [](const IntervalConsts &consts, const IntervalParams &params) {
            auto SubDomain = std::get<0>(params);
            std::vector<double> x_mid(SubDomain.size());
            for (size_t i = 0; i < SubDomain.size(); i++)
                x_mid[i] = median(SubDomain[i]); // TODO this was a little different in the orginal solver

            return objective(x_mid);
        };

        Definition.IsFeasible = [](const IntervalConsts& consts, const IntervalParams &params){
            /*std::vector<interval> domain = std::get<0>(params);
            auto xlow = domain[0].lower();
            auto xup = domain[0].upper();

            auto ylow = domain[0].lower();
            auto yup = domain[0].upper();

            if(xup*xup + yup*yup <= 2.0)
                return BnB::FEASIBILITY::Full;
            else if(xlow*xlow + ylow*ylow <= 2.0)
                return BnB::FEASIBILITY ::PARTIAL;
            else
                return BnB::FEASIBILITY::NONE;*/
            return BnB::FEASIBILITY::Full;
        };

        Definition.PrintSolution = [](const IntervalParams &params) {
            std::cout << "The solution is in the interval:" << std::endl;
            auto Domain = std::get<0>(params);
            for (const auto &interval : Domain) {
                std::cout << "[" <<  std::setprecision(16) << interval.lower() << ", " << interval.upper() << "]" << std::endl;
            }
        };

        Definition.GetInitialSubproblem = [](const IntervalConsts &prob) {
            IntervalParams Initial;
            std::get<0>(Initial) = std::get<0>(prob);
            return Initial;
        };
        return Definition;
    }
}
