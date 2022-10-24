#ifndef UTILS_H
#define UTILS_H

#include <complex>

#ifndef cdouble
#define cdouble std::complex<double>
#endif
template <class T>
T* alloc1d(int& ni);

template <class T>
T** alloc2d(int& ni, int& nj);

template <class T>
void free1d(T** arr, int& ni);

template <class T>
void free2d(T*** arr, int& ni, int& nj);

template <class T>
void linspace(T* arr, T& xi, T& xf, int& n);

template <class T>
void write_array1d(T* arr, int& n, std::string &path);

template <class T>
void write_array2d(T** arr, int& ni, int& nj, std::string &path);

void write_array2d_complex(cdouble** arr, int& ni, int& nj, std::string &path);

int factorial(int n);
#endif
