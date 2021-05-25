#pragma once

template<typename T>
T objective(const std::vector<T> &x) {
    int i = 0;	
    T y = (1.0 + pow(x[i] + x[i+1] + 1,2)*(19.0-14.0*x[i]+3.0*pow(x[i],2) - 14.0*x[i+1] + 6.0*x[i]*x[i+1] + 3.0*pos(x[i+1],2))*(30.0 + pow(2.0*x[i]- 3.0*x[i+1],2)*(18.0-32.0*x[i] + 12.0*pow(x[i],2),48.0*x[i+2] - 36.0*x[i]*x[i+1] + 27.0*pow(x[i+1],2)));
	
    return y;
}

