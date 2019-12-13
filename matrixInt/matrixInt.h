#ifndef ___MATRIXINT_INCLUDED_
#define ___MATRIXINT_INCLUDED_

#include <string>
#include <iostream>
#include <vector>

namespace std {
    class Matrix {
    private:
        vector<int> mat;
        vector<int>::size_type n, m;

    public:
        Matrix();
        Matrix(vector<int>::size_type n);
        Matrix(vector<int>::size_type n, vector<int>::size_type m);
        Matrix(vector< vector<int> > matData);
        Matrix(vector<int> matData, vector<int>::size_type n);
        Matrix(const Matrix& p);

        ~Matrix() {
            mat.clear();
        };

        vector<int>::size_type row()    const { return n; }
        vector<int>::size_type column() const { return m; }

        vector<int> at_row   (vector<int>::size_type r) const;
        vector<int> at_column(vector<int>::size_type c) const;
        int at(vector<int>::size_type r, vector<int>::size_type c) const;
        int& at(vector<int>::size_type r, vector<int>::size_type c);

        void resize(vector<int>::size_type nn, vector<int>::size_type mm);
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

        // Matrix row_echelon() const;
        // Matrix reduced_row_echelon() const;
        // int det() const;
        // Matrix inv() const;
        // Matrix pow() const;
        //
        string to_string() const;
        // string to_latex() const;
    };

    ostream& operator<<(ostream& s, const Matrix& x);

    istream& operator>>(istream& s,       Matrix& x);
}

#endif
