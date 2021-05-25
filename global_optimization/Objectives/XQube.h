#pragma once

template<typename T>
T objective(const std::vector<T> &x) {
    T y = T(0.0);
    for (size_t i=0; i<x.size(); i++) {
       y += pow(x[i],3) + pow(x[i], 2) - x[i]*10.0; // index i benutzen, x[i] nur einmal
    }
    return y;
}
