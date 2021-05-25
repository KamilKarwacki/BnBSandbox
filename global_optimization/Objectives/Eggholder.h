#pragma once

template<typename T>
T objective(const std::vector<T> &x) {
    T y = T(0);
    for(int i = 0; i < x.size() - 1; i+=2)
    {
    	y += -(x[i+1] + 47.0)*sin(sqrt(abs(x[i]*0.5 + (x[i+1] + 47.0)))) - x[i]*sin(sqrt(abs(x[i] - (x[i+1] + 47.0))));
    }
    return y;
}

