#ifndef ___MATRIXINTARRAY_INCLUDED_
#define ___MATRIXINTARRAY_INCLUDED_

#define ___MATRIXINTARRAY_DEBUG_ false

#include <string>
#include <iostream>
#include <vector>

class Matrix {
private:
    int n, m;
    int* mat;

public:
    explicit Matrix(int n, int m);
    explicit Matrix(std::vector< std::vector<int> > matData);
    explicit Matrix(std::vector<int> matData, int n);
    Matrix(const Matrix& p);

    ~Matrix() {
        delete[] mat;
    }

    int row()    const { return n; }
    int column() const { return m; }

    std::vector<int> at_row   (int r) const;
    std::vector<int> at_column(int c) const;
    int& at(int r, int c);
    const int& at(int r, int c) const;

    void resize(int nn, int mm);
    Matrix  operator+() const;
    Matrix  operator-() const;

    Matrix& operator=(Matrix p);

    Matrix& operator+=(Matrix p);
    Matrix& operator-=(Matrix p);
    Matrix& operator*=(Matrix p);
    Matrix& operator*=(int p);
    // Matrix& operator/=(Matrix p);
    Matrix& operator/=(int p);

    Matrix  operator+(Matrix p) const;
    Matrix  operator-(Matrix p) const;
    Matrix  operator*(Matrix p) const;
    Matrix  operator*(int p) const;
    friend Matrix operator*(int n, const Matrix& p);
    // Matrix  operator/(Matrix p) const;
    Matrix  operator/(int p) const;
    // friend Matrix operator/(int n, const Matrix& p);

    bool operator==(const Matrix& p) const;
    bool operator!=(const Matrix& p) const;

    Matrix row_echelon() const;
    // Matrix reduced_row_echelon() const;
    // int det() const;
    // Matrix inv() const;
    // Matrix pow() const;
    //
    std::string to_string() const;
    // std::string to_latex() const;

    friend std::ostream& operator<<(std::ostream& s, const Matrix& x);

    friend std::istream& operator>>(std::istream& s,       Matrix& x);
};

#endif
