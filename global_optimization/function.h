#pragma once
#include "boost/numeric/interval.hpp"
#include "boost/numeric/interval/io.hpp"
#include "Objectives/objective.h"



template<typename T>
void SetProblemParams(std::vector<T>& domain) {
    size_t n = 3;
    //domain.resize(n,T(-512.0,512.0));
    domain.resize(n,T(-5.0, 5.1));
    //domain.push_back({-0.2, 0.2});
    //domain.push_back({0.5,1.5});
    //domain.push_back({-4.2, 4});
    //domain.push_back({0.1, 5.0});
}

template <typename T>
using checking = boost::numeric::interval_lib::checking_base<T>;

template <typename T>
using rounding = boost::numeric::interval_lib::rounded_transc_std<T>;

template <typename T>
using boost_interval_transc_t = boost::numeric::interval
        <T,
                boost::numeric::interval_lib::policies<
                        boost::numeric::interval_lib::save_state<rounding<T>>,
                        checking<T>
                >
        >;

using interval = boost_interval_transc_t<double>;

void set_nan_to_inf(interval &v) {
    using checking = interval::traits_type::checking;
    if (checking::is_nan(v.lower())) v.set(checking::neg_inf(), v.upper());
    if (checking::is_nan(v.upper())) v.set(v.lower(), checking::pos_inf());
}

template<>
void encodeParam(std::stringstream &ss, const std::vector<interval> &p) {
    size_t s = p.size();
    ss << s << " ";
    for (const auto &i : p)
        ss << i.lower() << " " << i.upper() << " ";
}

template<>
void decodeParam(std::stringstream &ss, std::vector<interval> &p) {
    size_t size;
    ss >> size;
    p.resize(size);
    for (size_t i = 0; i < size; i++) {
        double u, l;
        ss >> l >> u;
        p[i].set(l, u);
    }
}
