#include <iostream>
#include <sstream>
#include <algorithm>
#include "matrix2dv.h"

using namespace std;

template <class T>
Matrix<T>::Matrix() {
    resize(0, 0);
}

template <class T>
Matrix<T>::Matrix(typename vector< vector<T> >::size_type n) {
    resize(n, 1);
}

template <class T>
Matrix<T>::Matrix(typename vector< vector<T> >::size_type n, typename vector<T>::size_type m) {
    resize(n, m);
}

template <class T>
Matrix<T>::Matrix(vector< vector<T> > matData) {
    mat = matData;
    n = matData.size();
    m = 0;
    for (typename vector< vector<T> >::size_type i = 0; i < n; ++i) m = max(m, matData.at(i).size());
    for (typename vector< vector<T> >::size_type i = 0; i < n; ++i) matData.at(i).resize(m);
}

template <class T>
Matrix<T>::Matrix(const Matrix& p) {
    *this = p;
}

template <class T>
vector<T> Matrix<T>::at_row (typename vector< vector<T> >::size_type r) const {
    return mat.at(r);
}

template <class T>
vector<T> Matrix<T>::at_column(typename vector<T>::size_type c) const {
    vector<T> res(n);
    for (typename vector<T>::size_type i = 0; i < n; ++i) res.at(i) = mat.at(i).at(c);
    return res;
}

template <class T>
T& Matrix<T>::at(typename vector< vector<T> >::size_type r, typename vector<T>::size_type c) {
    return mat.at(r).at(c);
}

template <class T>
void Matrix<T>::resize(typename vector< vector<T> >::size_type nn, typename vector<T>::size_type mm) {
    mat = vector< vector<T> >(nn, vector<T>(mm));
    n = nn;
    m = mm;
}

template <class T>
Matrix<T> Matrix<T>::operator+() const {
    return *this;
}

template <class T>
Matrix<T> Matrix<T>::operator-() const {
    Matrix t(*this);
    t *= -1;
    return t;
}

template <class T>
Matrix<T>& Matrix<T>::operator=(Matrix<T> p) {
    mat = p.mat;
    n = p.n;
    m = p.m;
    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator+=(Matrix<T> p) {
    if (n != p.n || m != p.m) {
        ostringstream errMsg;
        errMsg << "Dimensions of matrices do not match (" << n << "x" << m << ", " << p.n << "x" << p.m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    for (typename vector< vector<T> >::size_type i = 0; i < n; ++i) for (typename vector< vector<T> >::size_type j = 0; j < m; ++j) (*this).at(i, j) += p.at(i, j);
    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator-=(Matrix<T> p) {
    if (n != p.n || m != p.m) {
        ostringstream errMsg;
        errMsg << "Dimensions of matrices do not match (" << n << "x" << m << ", " << p.n << "x" << p.m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    for (typename vector< vector<T> >::size_type i = 0; i < n; ++i) for (typename vector< vector<T> >::size_type j = 0; j < m; ++j) (*this).at(i, j) -= p.at(i, j);
    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator*=(Matrix<T> p) {
    cout << "Entered: *=\n";
    *this = *this * p;
    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator*=(T p) {
    for (typename vector< vector<T> >::size_type i = 0; i < n; ++i) for (typename vector<T>::size_type j = 0; j < m; ++j) mat.at(i).at(j) *= p;
    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator/=(T p) {
    for (typename vector< vector<T> >::size_type i = 0; i < n; ++i) for (typename vector<T>::size_type j = 0; j < m; ++j) mat.at(i).at(j) /= p;
    return *this;
}

template <class T>
Matrix<T> Matrix<T>::operator+(Matrix<T> p) const {
    Matrix<T> t(*this);
    t += p;
    return t;
}

template <class T>
Matrix<T> Matrix<T>::operator-(Matrix<T> p) const {
    Matrix<T> t(*this);
    t -= p;
    return t;
}

template <class T>
Matrix<T> Matrix<T>::operator*(Matrix<T> p) const {
    if (m != p.n) {
        ostringstream errMsg;
        errMsg << "Dimensions of matrices do not match (" << n << "x" << m << ", " << p.n << "x" << p.m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    Matrix<T> ans(n, p.m);
    for (typename vector< vector<T> >::size_type i = 0; i < n; ++i) for (typename vector< vector<T> >::size_type j = 0; j < m; ++j) for (typename vector<T>::size_type k = 0; k < p.m; ++k) ans.mat.at(i).at(k) += mat.at(i).at((typename vector<T>::size_type)j) * p.mat.at(j).at(k);
    return ans;
}

template <class T>
Matrix<T> Matrix<T>::operator*(T p) const {
    Matrix<T> t(*this);
    t *= p;
    return t;
}

template <class T>
Matrix<T> operator*(T n, const Matrix<T>& p) {
    return p * n;
}

template <class T>
Matrix<T> Matrix<T>::operator/(T p) const {
    Matrix<T> t(*this);
    t /= p;
    return t;
}

template <class T>
bool Matrix<T>::operator==(const Matrix<T>& p) const {
    if (n != p.n || m != p.m) return false;
    for (typename vector< vector<T> >::size_type i = 0; i < n; ++i) for (typename vector<T>::size_type j = 0; j < m; ++j) if (mat.at(i).at(j) != p.mat.at(i).at(j)) return false;
    return true;
}

template <class T>
bool Matrix<T>::operator!=(const Matrix<T>& p) const {
    return !(*this == p);
}

template class Matrix<short>;
template class Matrix<unsigned short>;
template class Matrix<int>;
template class Matrix<unsigned int>;
template class Matrix<long>;
template class Matrix<unsigned long>;
template class Matrix<long long>;
template class Matrix<unsigned long long>;
template class Matrix<float>;
template class Matrix<double>;
template class Matrix<long double>;
