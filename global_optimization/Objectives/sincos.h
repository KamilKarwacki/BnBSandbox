#pragma once

template<typename T>
T objective(const std::vector<T> &x) {
	T y = T(0);
	for(int i = 0; i < x.size() - 1; i+=2)
		y += sin(x[i]) * cos(x[i + 1]);
    return y;
}

