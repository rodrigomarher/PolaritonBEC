#include <cmath>
#include <fstream>
#include <string>
#include <complex>
#include <iomanip>
#include <iostream>
#include "utils.h"

template <class T>
T* alloc1d(int& ni){
    T *arr;
    arr = new T[ni];
    return arr;
}

template <class T>
T** alloc2d(int& ni, int& nj){
    T** arr;
    arr = new T*[ni];
    for (int i=0; i<ni; i++){
        arr[i] = new T[nj];
    }
    return arr;
}

template <class T>
void free1d(T** arr, int& ni){
    delete[] (*arr);
}

template <class T>
void free2d(T*** arr, int& ni, int& nj){
    for(int i=0; i<ni;i++){
        delete[] (*arr)[i];
    }
    delete[] (*arr);
}


template <class T>
void linspace(T* arr, T& xi, T& xf, int& n){
    double dx = abs(xf-xi)/(double)(n-1);
    if (n==0){}
    else if (n==1){arr[0] = xi;}
    else{
        for(int i=0; i<n; i++){
            arr[i] = xi + i*dx;
        }
    }
}

template <class T>
void write_array1d(T* arr, int& n, std::string &path){
    std::ofstream outfile;
    outfile.open(path);
    if(outfile.is_open()){
        for(int i=0; i<n; i++){
            std::ostringstream doubleStr;
            doubleStr<<std::fixed<<std::setprecision(12);
            doubleStr<<arr[i];
            outfile<<doubleStr.str()<<" ";
        }
    }
    outfile.close();
}

template <class T>
void write_array2d(T** arr, int& ni, int& nj, std::string &path){
    std::ofstream outfile;
    outfile.open(path);
    if(outfile.is_open()){
        for(int i=0; i<ni; i++){
            for(int j=0;j<nj;j++){
                std::ostringstream doubleStr;
                doubleStr<<std::fixed<<std::setprecision(12);
                doubleStr<<arr[i][j];
                outfile<<doubleStr.str()<<" ";
            }
            outfile<<std::endl;
        }
    }
    outfile.close();
}

void write_array2d_complex(cdouble** arr, int& ni, int& nj, std::string &path){
    std::ofstream outfile;
    outfile.open(path);
    std::string sign;
    if(outfile.is_open()){
        for(int i=0; i<ni; i++){
            for(int j=0;j<nj;j++){
                sign = arr[i][j].imag()<0 ? "-": "+";
                std::ostringstream doubleStr;
                doubleStr<<std::fixed<<std::setprecision(12);
                doubleStr<<arr[i][j].real()<<sign<<std::abs(arr[i][j].imag())<<"j";
                outfile<<doubleStr.str()<<" ";
            }
            outfile<<std::endl;
        }
        outfile.close();
    }
}

int factorial(int n){
    if(n<2){
        return 1;
    }
    return n*factorial(n-1);
}

template int* alloc1d<int>(int& ni);
template double* alloc1d<double>(int& ni);
template cdouble* alloc1d<cdouble>(int& ni);
template int** alloc2d<int>(int& ni, int& nj);
template double** alloc2d<double>(int& ni, int& nj);
template cdouble** alloc2d<cdouble>(int& ni, int& nj);
template void free1d(int** arr,  int& ni);
template void free1d(double** arr,  int& ni);
template void free1d(cdouble** arr,  int& ni);
template void free2d(int*** arr,  int& ni, int& nj);
template void free2d(double*** arr,  int& ni, int& nj);
template void free2d(cdouble*** arr,  int& ni, int& nj);
template void linspace(int* arr, int& xi, int& xf, int& n);
template void linspace(double* arr, double& xi, double& xf, int& n);
template void write_array1d(int* arr,  int& n, std::string &path);
template void write_array1d(double* arr, int& n, std::string &path);
template void write_array1d(std::complex<double>* arr, int& n, std::string &path);
template void write_array2d(int** arr, int& ni,  int& nj, std::string &path);
template void write_array2d(double** arr, int& ni, int& nj, std::string &path);
template void write_array2d(std::complex<double>** arr, int& ni, int& nj, std::string &path);
