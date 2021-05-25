#pragma once
#include "Base.h"

namespace TreeGenerator::Knapsack {
    int counter = 1;
    ///std::ofstream TreeDataOutput;

    namespace Detail {
        struct item {
            int wt, cost;
            double ratio;
        };

        bool comp(item a, item b) {
            return a.ratio > b.ratio;
        }
    }

    using Consts = BnB::Problem_Constants<
            std::vector<int>, // Weights
            std::vector<int>, // Costs
            int>;             // Capacity of the Bag

    using Params = BnB::Subproblem_Parameters<
            std::vector<int>, // Taken objects
            int, // current weight
            int, // current value
            float,
            int, // parent
            int, // id
            int>;// depth

    std::tuple<float, float> NewBoundAsInPaper(const Consts &consts, const Params &params) {
        // calculate the cost as if all items were to fit in the sack
        int ObjectIndex = std::get<0>(params).size();
        int sumWeight = std::get<1>(params);
        int maxWeight = std::get<2>(consts);
        float sumCost = std::get<2>(params);
        int lastIndex = ObjectIndex;

        for (int i = ObjectIndex;
             i < std::get<0>(consts).size() and sumWeight + std::get<0>(consts)[i] <= maxWeight; i++) {
            sumCost += std::get<1>(consts)[i];
            sumWeight += std::get<0>(consts)[i];
            lastIndex = i;
        }

        if (lastIndex + 1 != std::get<0>(consts).size()) {
            float factor =
                    static_cast<float>(maxWeight - sumWeight) / static_cast<float>(std::get<0>(consts)[lastIndex + 1]);
            sumCost += std::get<1>(consts)[lastIndex + 1] * factor;
            sumCost += 1;
        }

        return std::tuple<float, float>(sumCost, -1);
    }


    BnB::Problem_Definition<Consts, Params, float> GenerateFasterProblem() {
        int id;
        GET_RANK(id);
        ///TreeDataOutput.open("output" + std::to_string(id) + ".txt");
        using namespace Detail;
        // object that will hold the functions called by the solver
        BnB::Problem_Definition<Consts, Params, float> KnapsackProblem;

        // you can assume that the subproblem is not feasible yet so it can be split
        KnapsackProblem.SplitSolution = [](const Consts &consts, const Params &params) {
            int id;
            GET_RANK(id);
            if (std::get<0>(params).size() >= std::get<0>(consts).size()) return std::vector<Params>();

            std::vector<Params> ret;
            int ObjectIndex = std::get<0>(params).size();
            Params s1 = params;
            Params s2 = params;
            std::get<0>(s1).push_back(0);
            std::get<3>(s1) = std::get<0>(NewBoundAsInPaper(consts, s1));
            std::get<4>(s1) = std::get<5>(params);  // parent
            std::get<6>(s1) = counter++; // unique id
            std::get<5>(s1) = 10000 * id + std::get<6>(s1); // id scaled by proc
            ret.push_back(s1); // without taking object at index pos
            if (std::get<0>(consts)[ObjectIndex] + std::get<1>(params) <= std::get<2>(consts)) {
                std::get<1>(s2) += std::get<0>(consts)[ObjectIndex];
                std::get<2>(s2) += std::get<1>(consts)[ObjectIndex];
                std::get<0>(s2).push_back(1);
                std::get<3>(s2) = std::get<0>(NewBoundAsInPaper(consts, s2));
                std::get<4>(s2) = std::get<5>(params);  // parent
                std::get<6>(s2) = counter++; // unique id
                std::get<5>(s2) = 10000 * id + std::get<6>(s2); // id scaled by proc

                ret.push_back(s2); // taking object at index pos
            }
            return ret;
        };

        KnapsackProblem.GetEstimateForBounds = [](const Consts &consts, const Params &params) {
            int id;
            GET_RANK(id);
            std::ofstream TreeDataOutput;
            TreeDataOutput.open("output" + std::to_string(id) + ".txt", std::ios_base::app);
            TreeDataOutput << std::get<5>(params) << ";" << std::get<4>(params) << ";" << id << ";" << 0 << std::endl;
            TreeDataOutput.close();
            return std::tuple(std::get<3>(params), -1);
        };


        KnapsackProblem.GetContainedUpperBound = [](const Consts &consts, const Params &params) {
            int w = std::get<1>(params);
            int cost = std::get<2>(params);

            return static_cast<float>(std::get<2>(params));
        };

        KnapsackProblem.IsFeasible = [](const Consts &consts, const Params &params) {
            return std::get<1>(params) > std::get<2>(consts) ? BnB::FEASIBILITY::NONE : BnB::FEASIBILITY::Full;
        };


        KnapsackProblem.PrintSolution = [](const Params &params) {

        };

        KnapsackProblem.GetInitialSubproblem = [](const Consts &prob) {
            Params Initial;
            std::get<3>(Initial) = std::get<0>(NewBoundAsInPaper(prob, Params()));
            return Initial;
        };

        return KnapsackProblem;
    }
}
