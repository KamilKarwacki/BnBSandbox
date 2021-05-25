#pragma once

template<typename T>
T objective(const std::vector<T> &x) {
    T y = T(0.0);
    for (size_t i=0; i<x.size(); i++) {
	   if(i%2==0)
       	y += pow(x[i],2) + 0.0001*x[i];
	   else
		y += pow(x[i],2) - 0.0001*x[i];
    }
    return y;
}
