#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <cmath>
#include "matrix.h"

using namespace std;

template <class T>
Matrix<T>::Matrix() {
    resize(0, 0);
}

template <class T>
Matrix<T>::Matrix(typename vector<T>::size_type n) {
    resize(n, 1);
}

template <class T>
Matrix<T>::Matrix(typename vector<T>::size_type n, typename vector<T>::size_type m) {
    resize(n, m);
}

template <class T>
Matrix<T>::Matrix(vector< vector<T> > matData) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: Constructor (vector< vector<int> >)\n";
    n = (typename vector<T>::size_type)matData.size();
    m = (typename vector<T>::size_type)0;
    for (typename vector< vector<T> >::size_type i = 0; i < n; ++i) m = max(m, (typename vector<T>::size_type)matData.at(i).size());
    mat = vector<T>(n * m);
    for (typename vector< vector<T> >::size_type i = 0; i < n; ++i) for (typename vector<T>::size_type j = 0; j < matData.at(i).size(); ++j) mat.at((typename vector<T>::size_type)i * m + j) = matData.at(i).at(j);
    if (___MATRIXINTARRAY_DEBUG_) cout << "End  : Constructor (vector< vector<int> >)\n";
}

template <class T>
Matrix<T>::Matrix(const Matrix<T>& p) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: Copy Constructor\n";
    mat = p.mat;
    n = p.n;
    m = p.m;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End  : Copy Constructor\n";
}

template <class T>
vector<T> Matrix<T>::at_row (typename vector<T>::size_type r) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: at_row(int): " << r << '\n';
    vector<T> res(m);
    for (typename vector<T>::size_type i = 0; i < m; ++i) res.at(i) = at(r, i);
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: at_row(int): " << r << '\n';
    return res;
}

template <class T>
vector<T> Matrix<T>::at_column(typename vector<T>::size_type c) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: at_column(int): " << c << '\n';
    vector<T> res(n);
    for (typename vector<T>::size_type i = 0; i < n; ++i) res.at(i) = at(i, c);
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: at_column(int): " << c << '\n';
    return res;
}


template <class T>
T& Matrix<T>::at(typename vector<T>::size_type r, typename vector<T>::size_type c) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: int& at(int, int): " << r << ", " << c << '\n';
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: int& at(int, int): " << r << ", " << c << '\n';
    return mat.at(r * m + c);
}

template <class T>
const T& Matrix<T>::at(typename vector<T>::size_type r, typename vector<T>::size_type c) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: const int& at(int, int): " << r << ", " << c << '\n';
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: const int& at(int, int): " << r << ", " << c << '\n';
    return mat.at(r * m + c);
}

template <class T>
void Matrix<T>::push_row(vector<T> p) {
    if (m != p.size()) {
        ostringstream errMsg;
        errMsg << "Dimensions of matrices do not match (" << n << "x" << m << ", " << 1 << "x" << p.size() << ")\n";
        throw invalid_argument(errMsg.str());
    }
    resize(n + 1);
    // Attention: the variable n is changed in "resize" method
    for (typename vector<T>::size_type j = 0; j < m; ++j) at(n - 1, j) = p.at(j);
}

template <class T>
void Matrix<T>::pop_row() {
    if (!n) {
        ostringstream errMsg;
        errMsg << "The matrix is empty\n";
        throw invalid_argument(errMsg.str());
    }
    resize(n - 1);
}

template <class T>
void Matrix<T>::push_column(vector<T> p) {
    if (n != p.size()) {
        ostringstream errMsg;
        errMsg << "Dimensions of matrices do not match (" << n << "x" << m << ", " << p.size() << "x" << 1 << ")\n";
        throw invalid_argument(errMsg.str());
    }
    resize(n, m + 1);
    // Attention: the variable m is changed in "resize" method
    for (typename vector<T>::size_type i = 0; i < n; ++i) at(i, m - 1) = p.at(i);
}

template <class T>
void Matrix<T>::pop_column() {
    if (!m) {
        ostringstream errMsg;
        errMsg << "The matrix is empty\n";
        throw invalid_argument(errMsg.str());
    }
    resize(n, m - 1);
}

template <class T>
void Matrix<T>::append_rows(Matrix<T> p) {
    if (m != p.m) {
        ostringstream errMsg;
        errMsg << "Dimensions of matrices do not match (" << n << "x" << m << ", " << p.n << "x" << p.m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    typename vector<T>::size_type nc = n;
    resize(n + p.n);
    // Attention: the variable n is changed in "resize" method
    for (typename vector<T>::size_type i = 0; i < p.n; ++i) for (typename vector<T>::size_type j = 0; j < m; ++j) at(i + nc, j) = p.at(i, j);
}

template <class T>
void Matrix<T>::append_columns(Matrix<T> p) {
    if (n != p.n) {
        ostringstream errMsg;
        errMsg << "Dimensions of matrices do not match (" << n << "x" << m << ", " << p.n << "x" << p.m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    typename vector<T>::size_type mc = m;
    resize(n, m + p.m);
    // Attention: the variable m is changed in "resize" method
    for (typename vector<T>::size_type i = 0; i < n; ++i) for (typename vector<T>::size_type j = 0; j < p.m; ++j) at(i, j + mc) = p.at(i, j);
}

template <class T>
void Matrix<T>::remove_rows(typename vector<T>::size_type p) {
    if (n < p) {
        ostringstream errMsg;
        errMsg << "The number of rows being removed is exceeding the matrix (" << p << " rows from " << n << "x" << m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    resize(n - p);
}

template <class T>
void Matrix<T>::remove_columns(typename vector<T>::size_type p) {
    if (m < p) {
        ostringstream errMsg;
        errMsg << "The number of columns being removed is exceeding the matrix (" << p << " columns from " << n << "x" << m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    resize(n, m - p);
}

template <class T>
void Matrix<T>::resize(typename vector<T>::size_type nn) {
    resize(nn, m);
}

template <class T>
void Matrix<T>::resize(typename vector<T>::size_type nn, typename vector<T>::size_type mm) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: resize(int, int): " << nn << ", " << mm << '\n';
    vector<T> temp(nn * mm);
    for (typename vector<T>::size_type i = 0; i < min(n, nn); ++i) for (typename vector<T>::size_type j = 0; j < min(m, mm); ++j) temp.at(i * mm + j) = mat.at(i * m + j);
    mat = temp;
    n = nn;
    m = mm;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End  : resize(int, int): " << nn << ", " << mm << '\n';
}

template <class T>
Matrix<T> Matrix<T>::operator+() const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: unary +\n";
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: unary +\n";
    return *this;
}

template <class T>
Matrix<T> Matrix<T>::operator-() const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: unary -\n";
    Matrix<T> t(*this);
    t *= -1;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: unary -\n";
    return t;
}

template <class T>
Matrix<T>& Matrix<T>::operator=(Matrix<T> p) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: assign =\n";
    mat = p.mat;
    n = p.n;
    m = p.m;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: assign =\n";
    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator+=(Matrix<T> p) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: assign +=\n";
    if (n != p.n || m != p.m) {
        ostringstream errMsg;
        errMsg << "Dimensions of matrices do not match (" << n << "x" << m << ", " << p.n << "x" << p.m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    for (typename vector<T>::size_type i = 0; i < n; ++i) for (typename vector<T>::size_type j = 0; j < m; ++j) at(i, j) += p.at(i, j);
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: assign +=\n";
    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator-=(Matrix<T> p) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: assign -=\n";
    if (n != p.n || m != p.m) {
        ostringstream errMsg;
        errMsg << "Dimensions of matrices do not match (" << n << "x" << m << ", " << p.n << "x" << p.m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    for (typename vector<T>::size_type i = 0; i < n; ++i) for (typename vector<T>::size_type j = 0; j < m; ++j) at(i, j) -= p.at(i, j);
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: assign -=\n";
    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator*=(Matrix<T> p) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: assign *= Matrix\n";
    *this = *this * p;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: assign *= Matrix\n";
    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator*=(T p) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: assign *= Scalar\n";
    for (typename vector<T>::size_type i = 0; i < n; ++i) for (typename vector<T>::size_type j = 0; j < m; ++j) at(i, j) *= p;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: assign *= Scalar\n";
    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator/=(Matrix<T> p) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: assign /= Matrix\n";
    *this = *this / p;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: assign /= Matrix\n";
    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator/=(T p) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: assign /= Scalar\n";
    for (typename vector<T>::size_type i = 0; i < n; ++i) for (typename vector<T>::size_type j = 0; j < m; ++j) at(i, j) /= p;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: assign /= Scalar\n";
    return *this;
}

template <class T>
Matrix<T> Matrix<T>::operator+(Matrix<T> p) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: binary +\n";
    Matrix<T> t(*this);
    t += p;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: binary +\n";
    return t;
}

template <class T>
Matrix<T> Matrix<T>::operator-(Matrix<T> p) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: binary -\n";
    Matrix<T> t(*this);
    t -= p;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: binary -\n";
    return t;
}

template <class T>
Matrix<T> Matrix<T>::operator*(Matrix<T> p) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: binary * Matrix\n";
    if (m != p.n) {
        ostringstream errMsg;
        errMsg << "Dimensions of matrices do not match (" << n << "x" << m << ", " << p.n << "x" << p.m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    Matrix<T> ans(n, p.m);
    for (typename vector<T>::size_type i = 0; i < n; ++i) for (typename vector<T>::size_type j = 0; j < m; ++j) for (typename vector<T>::size_type k = 0; k < p.m; ++k) ans.at(i, k) += at(i, j) * p.at(j, k);
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: binary * Matrix\n";
    return ans;
}

template <class T>
Matrix<T> Matrix<T>::operator*(T p) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: binary * Scalar\n";
    Matrix<T> t(*this);
    t *= p;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: binary * Scalar\n";
    return t;
}

// template <class T>
// Matrix<T> operator*(T n, const Matrix<T>& p) {
//     if (___MATRIXINTARRAY_DEBUG_) cout << "Start: friend binary * Scalar\n";
//     if (___MATRIXINTARRAY_DEBUG_) cout << "End r: friend binary * Scalar\n";
//     return p * n;
// }

template <class T>
Matrix<T> Matrix<T>::operator/(Matrix<T> p) const {
    return (*this) * p.inv();
}

template <class T>
Matrix<T> Matrix<T>::operator/(T p) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: binary / Scalar\n";
    Matrix<T> t(*this);
    t /= p;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: binary / Scalar\n";
    return t;
}

template <class T>
bool Matrix<T>::operator==(const Matrix<T>& p) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: comparison ==\n";
    if (n != p.n || m != p.m) return false;
    for (typename vector<T>::size_type i = 0; i < n; ++i) for (typename vector<T>::size_type j = 0; j < m; ++j) if (at(i, j) != p.at(i, j)) return false;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: comparison ==\n";
    return true;
}

template <class T>
bool Matrix<T>::operator!=(const Matrix<T>& p) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: comparison !=\n";
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: comparison !=\n";
    return !(*this == p);
}

template <class T>
bool Matrix<T>::is_zero() const {
    for (typename vector<T>::size_type i = 0; i < n; ++i) for (typename vector<T>::size_type j = 0; j < m; ++j) if (at(i, j)) return false;
    return true;
}

template <class T>
tuple<Matrix<T>, bool, typename vector<T>::size_type> Matrix<T>::gauss_count() const {
    Matrix<T> temp(*this);
    bool neg = false;
    typename vector<T>::size_type y = 0, x = 0;
    while (y < n && x < m) {
        typename vector<T>::size_type mPos = 0;
        T mVal = 0;
        for (typename vector<T>::size_type i = y; i < n; ++i) if (mVal < abs(temp.at(i, x))) mPos = i, mVal = abs(temp.at(i, x));
        if (!mVal) ++x;
        else {
            if (y != mPos) for (typename vector<T>::size_type j = 0; j < m; ++j) swap(temp.at(y, j), temp.at(mPos, j)), neg = !neg;
            for (typename vector<T>::size_type i = y + 1; i < n; ++i) {
                T f = temp.at(i, x) / temp.at(y, x);
                temp.at(i, x) = 0;
                for (typename vector<T>::size_type j = x + 1; j < m; ++j) temp.at(i, j) -= temp.at(y, j) * f;
            }
            ++y, ++x;
        }
    }
    return make_tuple(temp, neg, y);
}

template <class T>
Matrix<T> Matrix<T>::gauss() const {
    return get<0>(gauss_count());
}


template <class T>
Matrix<T> Matrix<T>::gauss_jordan() const {
    Matrix<T> temp = gauss();
    vector<typename vector<T>::size_type> lPos;
    lPos.reserve(m);
    typename vector<T>::size_type y = 0;
    for (typename vector<T>::size_type j = 0; y < n && j < m; ++j) if (temp.at(y, j)) lPos.push_back(j), ++y;
    --y;
    for (; y >= 0 && y < n; --y) {
        typename vector<T>::size_type x = lPos.at(y);
        for (typename vector<T>::size_type j = x + 1; j < m; ++j) temp.at(y, j) /= temp.at(y, x);
        temp.at(y, x) = 1;
        for (typename vector<T>::size_type i = 0; i < y; ++i) {
            for (typename vector<T>::size_type j = x + 1; j < m; ++j) temp.at(i, j) -= temp.at(y, j) * temp.at(i, x);
            temp.at(i, x) = 0;
        }
    }
    return temp;
}

template <class T>
T Matrix<T>::det() const {
    if (n != m) {
        ostringstream errMsg;
        errMsg << "Argument is not square matrix (" << n << "x" << m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    tuple<Matrix<T>, bool, typename vector<T>::size_type> temp = gauss_count();
    T ans = (get<1>(temp) ? -1 : 1);
    for (typename vector<T>::size_type i = 0; i < n; ++i) ans *= get<0>(temp).at(i, i);
    return ans;
}

template <class T>
Matrix<T> Matrix<T>::inv() const {
    if (n != m) {
        ostringstream errMsg;
        errMsg << "Argument is not square matrix (" << n << "x" << m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    AugmentedMatrix<T> temp1((*this), identity_matrix(n, (T)1));
    AugmentedMatrix<T> temp2 = temp1.gauss_jordan();
    if (temp2.left() != identity_matrix(n, (T)1)) {
        ostringstream errMsg;
        errMsg << "Argument is not invertible\n";
        throw invalid_argument(errMsg.str());
    }
    return temp2.right();
}

template <class T>
Matrix<T> Matrix<T>::pow(int r) const {
    if (n != m) {
        ostringstream errMsg;
        errMsg << "Argument is not square matrix (" << n << "x" << m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    if (r < 0) return inv().pow(-r);
    else if (!r) return identity_matrix(n, (T)1);
    else if (r % 2) return pow(r - 1) * (*this);
    else {
        Matrix<T> t = pow(r / 2);
        return t * t;
    }
}

template <class T>
T Matrix<T>::trace() const {
    if (n != m) {
        ostringstream errMsg;
        errMsg << "Argument is not square matrix (" << n << "x" << m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    T ans = 0;
    for (typename vector<T>::size_type i = 0; i < n; ++i) ans += at(i, i);
    return ans;
}

template <class T>
T Matrix<T>::rank() const {
    return get<2>(gauss_count());
}

template <class T>
T Matrix<T>::norm_squared() const {
    if (n != 1 && m != 1) {
        ostringstream errMsg;
        errMsg << "Arguments are neither row nor column vector (" << n << "x" << m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    if (m != 1) return transpose().norm_squared();
    T ans = 0;
    for (typename vector<T>::size_type i = 0; i < n; ++i) ans += at(i, 0) * at(i, 0);
    return ans;
}

template <class T>
double Matrix<T>::norm() const {
    return sqrt((double)norm_squared());
}

template <class T>
double Matrix<T>::distance(Matrix<T> p) const {
    return (*this - p).norm();
}

template <class T>
double Matrix<T>::angle(Matrix<T> p) const {
    return acos((double)dot(p) / norm() / p.norm());
}

template <class T>
T Matrix<T>::dot(Matrix<T> p) const {
    if (n != p.n || m != p.m) {
        ostringstream errMsg;
        errMsg << "Dimensions of matrices do not match (" << n << "x" << m << ", " << p.n << "x" << p.m << ")\n";
        throw invalid_argument(errMsg.str());
    } else if (n != 1 && m != 1) {
        ostringstream errMsg;
        errMsg << "Arguments are neither row nor column vector (" << n << "x" << m << ", " << p.n << "x" << p.m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    if (m != 1) return transpose().dot(p.transpose());
    T ans = 0;
    for (typename vector<T>::size_type i = 0; i < n; ++i) ans += at(i, 0) * p.at(i, 0);
    return ans;
}

template <class T>
bool Matrix<T>::orthogonal(Matrix<T> p) const {
    return (dot(p) == 0);
}

template <class T>
Matrix<T> Matrix<T>::projection(Matrix<T> p) const {
    return p * (dot(p) / p.norm_squared());
}

template <class T>
Matrix<T> Matrix<T>::basis() const {
    Matrix<T> temp1;
    typename vector<T>::size_type temp2;
    tie(temp1, ignore, temp2) = transpose().gauss_count();
    temp1.resize(temp2);
    return temp1;
}

// template <class T>
// Matrix<T> Matrix<T>::orthonormal() const {
//     vector< vector<T> > temp;
//     temp.reserve(n);
//     for (typename vector<T>::size_type i = 0; i < n; ++i) {
//         ;
//     }
//     return get<0>(temp1);
// }

template <class T>
Matrix<T> Matrix<T>::cross(Matrix<T> p) const {
    if (n != p.n || m != p.m) {
        ostringstream errMsg;
        errMsg << "Dimensions of matrices do not match (" << n << "x" << m << ", " << p.n << "x" << p.m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    if (n == 1 && m == 3) return transpose().cross(p.transpose());
    if (n == 3 && m == 1) {
        Matrix<T> ans(3);
        for (typename vector<T>::size_type i = 0; i < 3; ++i) ans.at(i, 0) = at((i + 1) % 3, 0) * p.at((i + 2) % 3, 0) - at((i + 2) % 3, 0) * p.at((i + 1) % 3, 0);
        return ans;
    }
    ostringstream errMsg;
    errMsg << "Arguments are neither 3x1 nor 1x3 vector (" << n << "x" << m << ", " << p.n << "x" << p.m << ")\n";
    throw invalid_argument(errMsg.str());
}

template <class T>
T Matrix<T>::minor(typename vector<T>::size_type y, typename vector<T>::size_type x) const {
    if (n != m) {
        ostringstream errMsg;
        errMsg << "Argument is not square matrix (" << n << "x" << m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    Matrix<T> temp(n - 1, m - 1);
    for (typename vector<T>::size_type i = 0; i < n; ++i) {
        if (y == i) continue;
        for (typename vector<T>::size_type j = 0; j < m; ++j) {
            if (x == j) continue;
            temp.at(i - (i > y), j - (j > x)) = at(i, j);
        }
    }
    return temp.det();
}

template <class T>
T Matrix<T>::cofactor(typename vector<T>::size_type y, typename vector<T>::size_type x) const {
    return minor(y, x) * ((y + x) % 2 ? -1 : 1);
}

template <class T>
Matrix<T> Matrix<T>::cofactor_matrix() const {
    Matrix<T> temp(n, m);
    for (typename vector<T>::size_type i = 0; i < n; ++i) for (typename vector<T>::size_type j = 0; j < m; ++j) temp.at(i, j) = cofactor(i, j);
    return temp;
}

template <class T>
Matrix<T> Matrix<T>::adjugate() const {
    return cofactor_matrix().transpose();
}

template <class T>
Matrix<T> Matrix<T>::transpose() const {
    Matrix<T> temp(m, n);
    for (typename vector<T>::size_type i = 0; i < n; ++i) for (typename vector<T>::size_type j = 0; j < m; ++j) temp.at(j, i) = at(i, j);
    return temp;
}

template <class T>
string Matrix<T>::to_string() const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: to_string\n";
    ostringstream s;
    for (typename vector<T>::size_type i = 0; i < n; ++i) for (typename vector<T>::size_type j = 0; j < m; ++j) s << at(i, j) << (j == m - 1 ? "\n" : ", ");
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: to_string\n";
    return s.str();
    // Count the number of digits of each component in a column, and forms matrix using ascii art
    // for (typename vector<T>::size_type j = 0; j < m; ++j) {
    //     string::size_type m = 0;
    //     for (typename vector<T>::size_type i = 0; i < n; ++i) m = max(m, to_string(at(i, j)).size())
    // }
}

template class Matrix<short>;
// template class Matrix<unsigned short>;
template class Matrix<int>;
// template class Matrix<unsigned int>;
template class Matrix<long>;
// template class Matrix<unsigned long>;
template class Matrix<long long>;
// template class Matrix<unsigned long long>;
template class Matrix<float>;
template class Matrix<double>;
template class Matrix<long double>;

template <class T>
AugmentedMatrix<T>::AugmentedMatrix(typename vector<T>::size_type n, typename vector<T>::size_type ml) {
    resize(n, ml, 1);
}

template <class T>
AugmentedMatrix<T>::AugmentedMatrix(typename vector<T>::size_type n, typename vector<T>::size_type ml, typename vector<T>::size_type mr) {
    resize(n, ml, mr);
}

template <class T>
AugmentedMatrix<T>::AugmentedMatrix(Matrix<T> mmatl) {
    matl = mmatl;
    matr = Matrix<T>(mmatl.row());
    n = mmatl.row();
    ml = mmatl.column();
    mr = 1;
}

template <class T>
AugmentedMatrix<T>::AugmentedMatrix(Matrix<T> mmatl, Matrix<T> mmatr) {
    if (mmatl.row() != mmatr.row()) {
        ostringstream errMsg;
        errMsg << "The number of rows does not match (" << mmatl.row() << "x" << mmatl.column() << ", " << mmatr.row() << "x" << mmatr.column() << ")\n";
        throw invalid_argument(errMsg.str());
    }
    matl = mmatl;
    matr = mmatr;
    n = mmatl.row();
    ml = mmatl.column();
    mr = mmatr.column();
}

template <class T>
AugmentedMatrix<T>::AugmentedMatrix(const AugmentedMatrix& p) {
    matl = p.matl;
    matr = p.matr;
    n = p.n;
    ml = p.ml;
    mr = p.mr;
}

template <class T>
pair< vector<T>, vector<T> > AugmentedMatrix<T>::at_row (typename vector<T>::size_type r) const {
    if (r < 0 || r >= n) {
        ostringstream errMsg;
        errMsg << "The argument is out of range (" << r << "-th row is called for [" << n << "x" << ml << "|" << n << "x" << mr << "] augmented matrix\n";
        throw out_of_range(errMsg.str());
    }
    return make_pair(matl.at_row(r), matr.at_row(r));
}

template <class T>
vector<T> AugmentedMatrix<T>::at_column (typename vector<T>::size_type c) const {
    if (c < 0 || c >= ml + mr) {
        ostringstream errMsg;
        errMsg << "The argument is out of range (" << c << "-th column is called for [" << n << "x" << ml << "|" << n << "x" << mr << "] augmented matrix\n";
        throw out_of_range(errMsg.str());
    }
    if (c < ml) return at_columnl(c);
    else return at_columnr(c - ml);
}

template <class T>
vector<T> AugmentedMatrix<T>::at_columnl(typename vector<T>::size_type c) const {
    return matl.at_column(c);
}

template <class T>
vector<T> AugmentedMatrix<T>::at_columnr(typename vector<T>::size_type c) const {
    return matr.at_column(c);
}

template <class T>
T& AugmentedMatrix<T>::at(typename vector<T>::size_type r, typename vector<T>::size_type c) {
    if (c < 0 || c >= ml + mr) {
        ostringstream errMsg;
        errMsg << "The argument is out of range ((" << r << ", " << c << ") is called for [" << n << "x" << ml << "|" << n << "x" << mr << "] augmented matrix\n";
        throw out_of_range(errMsg.str());
    }
    if (c < ml) return matl.at(r, c);
    else return matr.at(r, c - ml);
}

template <class T>
const T& AugmentedMatrix<T>::at(typename vector<T>::size_type r, typename vector<T>::size_type c) const {
    if (c < 0 || c >= ml + mr) {
        ostringstream errMsg;
        errMsg << "The argument is out of range ((" << r << ", " << c << ") is called for [" << n << "x" << ml << "|" << n << "x" << mr << "] augmented matrix\n";
        throw out_of_range(errMsg.str());
    }
    if (c < ml) return matl.at(r, c);
    else return matr.at(r, c - ml);
}

template <class T>
Matrix<T> AugmentedMatrix<T>::left() const {
    return matl;
}

template <class T>
Matrix<T> AugmentedMatrix<T>::right() const {
    return matr;
}

template <class T>
void AugmentedMatrix<T>::resize(typename vector<T>::size_type nn, typename vector<T>::size_type mml, typename vector<T>::size_type mmr) {
    matl = Matrix<T>(nn, mml);
    matr = Matrix<T>(nn, mmr);
    n = nn;
    ml = mml;
    mr = mmr;
}

template <class T>
bool AugmentedMatrix<T>::operator==(const AugmentedMatrix& p) const {
    return (matl == p.matl && matr == p.matr);
}

template <class T>
bool AugmentedMatrix<T>::operator!=(const AugmentedMatrix& p) const {
    return !(*this == p);
}

template <class T>
AugmentedMatrix<T> AugmentedMatrix<T>::gauss() const {
    AugmentedMatrix<T> temp(*this);
    typename vector<T>::size_type y = 0, x = 0;
    while (y < n && x < ml) {
        typename vector<T>::size_type mPos = 0;
        T mVal = 0;
        for (typename vector<T>::size_type i = y; i < n; ++i) if (mVal < abs(temp.at(i, x))) mPos = i, mVal = abs(temp.at(i, x));
        if (!mVal) ++x;
        else {
            for (typename vector<T>::size_type j = 0; j < ml + mr; ++j) swap(temp.at(y, j), temp.at(mPos, j));
            for (typename vector<T>::size_type i = y + 1; i < n; ++i) {
                T f = temp.at(i, x) / temp.at(y, x);
                temp.at(i, x) = 0;
                for (typename vector<T>::size_type j = x + 1; j < ml + mr; ++j) temp.at(i, j) -= temp.at(y, j) * f;
            }
            ++y, ++x;
        }
    }
    return temp;
}

template <class T>
AugmentedMatrix<T> AugmentedMatrix<T>::gauss_jordan() const {
    AugmentedMatrix<T> temp = gauss();
    vector<typename vector<T>::size_type> lPos;
    lPos.reserve(ml);
    typename vector<T>::size_type y = 0;
    for (typename vector<T>::size_type j = 0; y < n && j < ml; ++j) if (temp.at(y, j)) lPos.push_back(j), ++y;
    --y;
    for (; y >= 0 && y < n; --y) {
        typename vector<T>::size_type x = lPos.at(y);
        for (typename vector<T>::size_type j = x + 1; j < ml + mr; ++j) temp.at(y, j) /= temp.at(y, x);
        temp.at(y, x) = 1;
        for (typename vector<T>::size_type i = 0; i < y; ++i) {
            for (typename vector<T>::size_type j = x + 1; j < ml + mr; ++j) temp.at(i, j) -= temp.at(y, j) * temp.at(i, x);
            temp.at(i, x) = 0;
        }
    }
    return temp;
}

template class AugmentedMatrix<short>;
// template class AugmentedMatrix<unsigned short>;
template class AugmentedMatrix<int>;
// template class AugmentedMatrix<unsigned int>;
template class AugmentedMatrix<long>;
// template class AugmentedMatrix<unsigned long>;
template class AugmentedMatrix<long long>;
// template class AugmentedMatrix<unsigned long long>;
template class AugmentedMatrix<float>;
template class AugmentedMatrix<double>;
template class AugmentedMatrix<long double>;
