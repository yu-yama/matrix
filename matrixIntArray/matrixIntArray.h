#ifndef ___MATRIXINTARRAY_INCLUDED_
#define ___MATRIXINTARRAY_INCLUDED_

#define ___MATRIXINTARRAY_DEBUG_ false

#include <string>
#include <iostream>
#include <vector>

namespace std {
    class Matrix {
    private:
        int n, m;
        int* mat;

    public:
        explicit Matrix(int n, int m);
        explicit Matrix(vector< vector<int> > matData);
        explicit Matrix(vector<int> matData, int n);
        Matrix(const Matrix& p);

        ~Matrix() {
            delete[] mat;
        }

        int row()    const { return n; }
        int column() const { return m; }

        vector<int> at_row   (int r) const;
        vector<int> at_column(int c) const;
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
        string to_string() const;
        // string to_latex() const;

        friend ostream& operator<<(ostream& s, const Matrix& x);

        friend istream& operator>>(istream& s,       Matrix& x);
    };
}

#endif
