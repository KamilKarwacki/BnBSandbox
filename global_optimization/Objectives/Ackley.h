#pragma once

const double PI = 3.1415926535897;
const double e =  2.7182818284590 + 20;

template<typename T>
T objective(const std::vector<T> &x) {
    T y = T(0);
    for(int i = 0; i < x.size() - 1; i+=2)
    {
    	y += exp(-0.2*sqrt(0.5*(pow(x[i],2) + pow(x[i + 1],2))))*(-20.0)
		-exp(0.5*(cos(x[i]*(PI*2.0)) + cos(x[i+1])*(PI*2.0))) + e;
    }

    return y;
}

