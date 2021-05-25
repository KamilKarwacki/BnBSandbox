#pragma once


template<typename T>
T objective(const std::vector<T> &x) {
    T y = T(0);
    for(int i = 0; i < x.size() - 1; i+=2)
    {
    	y += pow(pow(x[i],2) + x[i+1] - 11.0,2) + pow(x[i] + pow(x[i+1],2) -7.0,2) - 0.0001*(x[i] +x[i+1]);
		
    }

    return y;
}

