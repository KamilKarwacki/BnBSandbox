#pragma once

const double PI = 3.1415926535897;

template<typename T>
T objective(const std::vector<T> &x) {
    T y = T(0);
	int size = x.size();
    for(int i = 0; i < x.size(); i++)
    {
    	y +=10*size +  pow(x[i],2) - 10.0*cos(2*PI*x[i]);	
    }

    return y;
}

