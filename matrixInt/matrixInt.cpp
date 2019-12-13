#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include "matrixInt.h"

using namespace std;

Matrix::Matrix() {
    resize(0, 0);
}

Matrix::Matrix(vector<int>::size_type n) {
    resize(n, 1);
}

Matrix::Matrix(vector<int>::size_type n, vector<int>::size_type m) {
    resize(n, m);
}

Matrix::Matrix(vector< vector<int> > matData) {
    n = (vector<int>::size_type)matData.size();
    m = (vector<int>::size_type)0;
    for (vector< vector<int> >::size_type i = 0; i < n; ++i) m = max(m, (vector<int>::size_type)matData.at(i).size());
    mat = vector<int>(n * m);
    for (vector< vector<int> >::size_type i = 0; i < n; ++i) for (vector<int>::size_type j = 0; j < matData.at(i).size(); ++j) mat.at((vector<int>::size_type)i * n + j) = matData.at(i).at(j);
}

Matrix::Matrix(const Matrix& p) {
    *this = p;
}

vector<int> Matrix::at_row(vector<int>::size_type r) const {
    vector<int> res(m);
    for (vector<int>::size_type i = 0; i < m; ++i) res.at(i) = at(r, i);
    return res;
}

vector<int> Matrix::at_column(vector<int>::size_type c) const {
    vector<int> res(n);
    for (vector<int>::size_type i = 0; i < n; ++i) res.at(i) = at(i, c);
    return res;
}

int Matrix::at(vector<int>::size_type r, vector<int>::size_type c) const {
    return mat.at(r * n + c);
}

int& Matrix::at(vector<int>::size_type r, vector<int>::size_type c) {
    return mat.at(r * n + c);
}

void Matrix::resize(vector<int>::size_type nn, vector<int>::size_type mm) {
    mat = vector<int>(nn * mm);
    n = nn;
    m = mm;
}

Matrix Matrix::operator+() const {
    return *this;
}

Matrix Matrix::operator-() const {
    Matrix t(*this);
    t *= -1;
    return t;
}

Matrix& Matrix::operator=(Matrix p) {
    mat = p.mat;
    n = p.n;
    m = p.m;
    return *this;
}

Matrix& Matrix::operator+=(Matrix p) {
    if (n != p.n || m != p.m) {
        ostringstream errMsg;
        errMsg << "Dimensions of matrices do not match (" << n << "x" << m << ", " << p.n << "x" << p.m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    for (vector<int>::size_type i = 0; i < n; ++i) for (vector<int>::size_type j = 0; j < m; ++j) (*this).at(i, j) += p.at(i, j);
    return *this;
}

Matrix& Matrix::operator-=(Matrix p) {
    if (n != p.n || m != p.m) {
        ostringstream errMsg;
        errMsg << "Dimensions of matrices do not match (" << n << "x" << m << ", " << p.n << "x" << p.m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    for (vector<int>::size_type i = 0; i < n; ++i) for (vector<int>::size_type j = 0; j < m; ++j) (*this).at(i, j) -= p.at(i, j);
    return *this;
}

Matrix& Matrix::operator*=(Matrix p) {
    cout << "Entered: *=\n";
    *this = *this * p;
    return *this;
}

Matrix& Matrix::operator*=(int p) {
    for (vector<int>::size_type i = 0; i < n; ++i) for (vector<int>::size_type j = 0; j < m; ++j) (*this).at(i, j) *= p;
    return *this;
}

Matrix& Matrix::operator/=(int p) {
    for (vector<int>::size_type i = 0; i < n; ++i) for (vector<int>::size_type j = 0; j < m; ++j) (*this).at(i, j) /= p;
    return *this;
}

Matrix Matrix::operator+(Matrix p) const {
    Matrix t(*this);
    t += p;
    return t;
}

Matrix Matrix::operator-(Matrix p) const {
    Matrix t(*this);
    t -= p;
    return t;
}

Matrix Matrix::operator*(Matrix p) const {
    if (m != p.n) {
        ostringstream errMsg;
        errMsg << "Dimensions of matrices do not match (" << n << "x" << m << ", " << p.n << "x" << p.m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    Matrix ans(n, p.m);
    for (vector<int>::size_type i = 0; i < n; ++i) for (vector<int>::size_type j = 0; j < m; ++j) for (vector<int>::size_type k = 0; k < p.m; ++k) ans.at(i, k) += (*this).at(i, j) * p.at(j, k);
    return ans;
}

Matrix Matrix::operator*(int p) const {
    Matrix t(*this);
    t *= p;
    return t;
}

Matrix operator*(int n, const Matrix& p) {
    return p * n;
}

Matrix Matrix::operator/(int p) const {
    Matrix t(*this);
    t /= p;
    return t;
}

bool Matrix::operator==(const Matrix& p) const {
    if (n != p.n || m != p.m) return false;
    for (vector<int>::size_type i = 0; i < n; ++i) for (vector<int>::size_type j = 0; j < m; ++j) if ((*this).at(i, j) != p.at(i, j)) return false;
    return true;
}

bool Matrix::operator!=(const Matrix& p) const {
    return !(*this == p);
}

string Matrix::to_string() const {
    ostringstream s;
    for (vector<int>::size_type i = 0; i < n; ++i) for (vector<int>::size_type j = 0; j < m; ++j) s << at(i, j) << (j == m - 1 ? "\n" : ", ");
    return s.str();
    // Count the number of digits of each component in a column, and forms matrix using ascii art
    // for (vector<int>::size_type j = 0; j < m; ++j) {
    //     string::size_type m = 0;
    //     for (vector<int>::size_type i = 0; i < n; ++i) m = max(m, to_string(at(i, j)).size())
    // }
}

ostream& operator<<(ostream& s, const Matrix& x) {
    return s << x.to_string();
}

istream& operator>>(istream& s,       Matrix& x) {
    for (vector<int>::size_type i = 0; i < x.row(); ++i) for (vector<int>::size_type j = 0; j < x.column(); ++j) s >> x.at(i, j);
    return s;
}
