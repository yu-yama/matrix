#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include "matrixIntArray.h"

using namespace std;

Matrix::Matrix(int n, int m = 1) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: Constructor (int, int): " << n << ", " << m << '\n';
    resize(n, m);
    if (___MATRIXINTARRAY_DEBUG_) cout << "End  : Constructor (int, int): " << n << ", " << m << '\n';
}

Matrix::Matrix(vector< vector<int> > matData) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: Constructor (vector< vector<int> >)\n";
    n = (int)matData.size();
    m = (int)0;
    for (vector< vector<int> >::size_type i = 0; i < (vector< vector<int> >::size_type)n; ++i) m = max(m, (int)matData.at(i).size());
    resize(n, m);
    for (vector< vector<int> >::size_type i = 0; i < (vector< vector<int> >::size_type)n; ++i) for (vector<int>::size_type j = 0; j < matData.at(i).size(); ++j) at((int)i, (int)j) = matData.at(i).at(j);
    if (___MATRIXINTARRAY_DEBUG_) cout << "End  : Constructor (vector< vector<int> >)\n";
}

Matrix::Matrix(const Matrix& p) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: Copy Constructor\n";
    resize(p.row(), p.column());
    for (int i = 0; i < p.row(); ++i) for (int j = 0; j < p.column(); ++j) at(i, j) = p.at(i, j);
    if (___MATRIXINTARRAY_DEBUG_) cout << "End  : Copy Constructor\n";
}

vector<int> Matrix::at_row(int r) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: at_row(int): " << r << '\n';
    vector<int> res(m);
    for (int i = 0; i < m; ++i) res.at(i) = at(r, i);
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: at_row(int): " << r << '\n';
    return res;
}

vector<int> Matrix::at_column(int c) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: at_column(int): " << c << '\n';
    vector<int> res(n);
    for (int i = 0; i < n; ++i) res.at(i) = at(i, c);
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: at_column(int): " << c << '\n';
    return res;
}

int& Matrix::at(int r, int c) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: int& at(int, int): " << r << ", " << c << '\n';
    if (r < 0 || c < 0 || r >= n || c >= m) throw out_of_range("Arguments given to at method is out of the range of matrix.");
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: int& at(int, int): " << r << ", " << c << '\n';
    return mat[r * n + c];
}

const int& Matrix::at(int r, int c) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: const int& at(int, int): " << r << ", " << c << '\n';
    if (r < 0 || c < 0 || r >= n || c >= m) throw out_of_range("Arguments given to at method is out of the range of matrix.");
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: const int& at(int, int): " << r << ", " << c << '\n';
    return mat[r * n + c];
}

void Matrix::resize(int nn, int mm) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: resize(int, int): " << nn << ", " << mm << '\n';
    if (nn <= 0 || mm <= 0) throw invalid_argument("The number of rows and columns of a matrix must be positive.");
    n = nn;
    m = mm;
    mat = new int[n * m];
    if (___MATRIXINTARRAY_DEBUG_) cout << "End  : resize(int, int): " << nn << ", " << mm << '\n';
}

Matrix Matrix::operator+() const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: unary +\n";
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: unary +\n";
    return *this;
}

Matrix Matrix::operator-() const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: unary -\n";
    Matrix t(*this);
    t *= -1;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: unary -\n";
    return t;
}

Matrix& Matrix::operator=(Matrix p) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: assign =\n";
    resize(p.row(), p.column());
    for (int i = 0; i < p.row(); ++i) for (int j = 0; j < p.column(); ++j) at(i, j) = p.at(i, j);
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: assign =\n";
    return *this;
}

Matrix& Matrix::operator+=(Matrix p) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: assign +=\n";
    if (n != p.n || m != p.m) {
        ostringstream errMsg;
        errMsg << "Dimensions of matrices do not match (" << n << "x" << m << ", " << p.n << "x" << p.m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    for (int i = 0; i < n; ++i) for (int j = 0; j < m; ++j) at(i, j) += p.at(i, j);
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: assign +=\n";
    return *this;
}

Matrix& Matrix::operator-=(Matrix p) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: assign -=\n";
    if (n != p.n || m != p.m) {
        ostringstream errMsg;
        errMsg << "Dimensions of matrices do not match (" << n << "x" << m << ", " << p.n << "x" << p.m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    for (int i = 0; i < n; ++i) for (int j = 0; j < m; ++j) at(i, j) -= p.at(i, j);
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: assign -=\n";
    return *this;
}

Matrix& Matrix::operator*=(Matrix p) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: assign *= Matrix\n";
    *this = *this * p;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: assign *= Matrix\n";
    return *this;
}

Matrix& Matrix::operator*=(int p) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: assign *= Scalar\n";
    for (int i = 0; i < n; ++i) for (int j = 0; j < m; ++j) at(i, j) *= p;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: assign *= Scalar\n";
    return *this;
}

Matrix& Matrix::operator/=(int p) {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: assign /= Scalar\n";
    for (int i = 0; i < n; ++i) for (int j = 0; j < m; ++j) at(i, j) /= p;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: assign /= Scalar\n";
    return *this;
}

Matrix Matrix::operator+(Matrix p) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: binary +\n";
    Matrix t(*this);
    t += p;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: binary +\n";
    return t;
}

Matrix Matrix::operator-(Matrix p) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: binary -\n";
    Matrix t(*this);
    t -= p;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: binary -\n";
    return t;
}

Matrix Matrix::operator*(Matrix p) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: binary * Matrix\n";
    if (m != p.n) {
        ostringstream errMsg;
        errMsg << "Dimensions of matrices do not match (" << n << "x" << m << ", " << p.n << "x" << p.m << ")\n";
        throw invalid_argument(errMsg.str());
    }
    Matrix ans(n, p.m);
    for (int i = 0; i < n; ++i) for (int j = 0; j < m; ++j) for (int k = 0; k < p.m; ++k) ans.at(i, k) += (*this).at(i, j) * p.at(j, k);
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: binary * Matrix\n";
    return ans;
}

Matrix Matrix::operator*(int p) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: binary * Scalar\n";
    Matrix t(*this);
    t *= p;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: binary * Scalar\n";
    return t;
}

Matrix operator*(int n, const Matrix& p) {
    // if (___MATRIXINTARRAY_DEBUG_) cout << "Start: friend binary * Scalar\n";
    // if (___MATRIXINTARRAY_DEBUG_) cout << "End r: friend binary * Scalar\n";
    return p * n;
}

Matrix Matrix::operator/(int p) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: binary / Scalar\n";
    Matrix t(*this);
    t /= p;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: binary / Scalar\n";
    return t;
}

bool Matrix::operator==(const Matrix& p) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: comparison ==\n";
    if (n != p.n || m != p.m) return false;
    for (int i = 0; i < n; ++i) for (int j = 0; j < m; ++j) if ((*this).at(i, j) != p.at(i, j)) return false;
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: comparison ==\n";
    return true;
}

bool Matrix::operator!=(const Matrix& p) const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: comparison !=\n";
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: comparison !=\n";
    return !(*this == p);
}

Matrix Matrix::row_echelon() const {
    Matrix temp(*this);
    int y = 0, x = 0;
    while (y < temp.n && x < temp.m) {
        int ymp = 0, ymv = 0;
        for (int i = y; i < temp.n; ++i) if (ymv < abs(temp.at(i, x))) ymp = i, ymv = abs(temp.at(i, x));
        if (ymv == 0) ++x;
        else {
            for (int j = 0; j < temp.m; ++j) swap(temp.at(y, j), temp.at(ymp, j));
            for (int i = y + 1; i < temp.n; ++i) {
                int f = temp.at(i, x) / temp.at(y, x);
                temp.at(i, x) = 0;
                for (int j = x + 1; j < temp.m; ++j) temp.at(i, j) -= temp.at(y, j) * f;
                ++y, ++x;
            }
        }
    }
    return temp;
}

string Matrix::to_string() const {
    if (___MATRIXINTARRAY_DEBUG_) cout << "Start: to_string\n";
    ostringstream s;
    for (int i = 0; i < n; ++i) for (int j = 0; j < m; ++j) s << at(i, j) << (j == m - 1 ? "\n" : ", ");
    if (___MATRIXINTARRAY_DEBUG_) cout << "End r: to_string\n";
    return s.str();
    // Count the number of digits of each component in a column, and forms matrix using ascii art
    // for (int j = 0; j < m; ++j) {
    //     string::size_type m = 0;
    //     for (int i = 0; i < n; ++i) m = max(m, to_string(at(i, j)).size())
    // }
}

ostream& operator<<(ostream& s, const Matrix& x) {
    // if (___MATRIXINTARRAY_DEBUG_) cout << "Start: ostream <<\n";
    // if (___MATRIXINTARRAY_DEBUG_) cout << "End r: ostream <<\n";
    return s << x.to_string();
}

istream& operator>>(istream& s,       Matrix& x) {
    // if (___MATRIXINTARRAY_DEBUG_) cout << "Start: istream >>\n";
    for (int i = 0; i < x.row(); ++i) for (int j = 0; j < x.column(); ++j) s >> x.at(i, j);
    // if (___MATRIXINTARRAY_DEBUG_) cout << "End r: istream >>\n";
    return s;
}
