#ifndef ___MATRIX2DV_INCLUDED_
#define ___MATRIX2DV_INCLUDED_

#include <string>
#include <iostream>
#include <vector>

namespace std {
    template <typename T>
    class Matrix {
        // typedef typename vector< vector<T> >::size_type size_type;
    private:
        vector< vector<T> > mat;
        typename vector< vector<T> >::size_type n;
        typename vector<T>::size_type m;

    public:
        Matrix();
        Matrix(typename vector< vector<T> >::size_type n);
        Matrix(typename vector< vector<T> >::size_type n, typename vector<T>::size_type m);
        Matrix(vector< vector<T> > matData);
        Matrix(const Matrix& p);

        typename vector< vector<T> >::size_type row()      const { return n; }
        typename vector<T>::size_type column()   const { return m; }

        vector<T> at_row   (typename vector< vector<T> >::size_type r) const;
        vector<T> at_column(typename vector<T>::size_type c) const;
        T& at(typename vector< vector<T> >::size_type r, typename vector<T>::size_type c);

        void resize(typename vector< vector<T> >::size_type nn, typename vector<T>::size_type mm);
        Matrix  operator+() const;
        Matrix  operator-() const;

        Matrix& operator=(Matrix p);

        Matrix& operator+=(Matrix p);
        Matrix& operator-=(Matrix p);
        Matrix& operator*=(Matrix p);
        Matrix& operator*=(T p);
        // Matrix& operator/=(Matrix p);
        Matrix& operator/=(T p);

        Matrix  operator+(Matrix p) const;
        Matrix  operator-(Matrix p) const;
        Matrix  operator*(Matrix p) const;
        Matrix  operator*(T p) const;
        template <typename S>
        friend Matrix<S> operator*(S n, const Matrix<S>& p);
        // Matrix  operator/(Matrix p) const;
        Matrix  operator/(T p) const;
        // template <typename S>
        // friend Matrix<S> operator/(S n, const Matrix<S>& p);

        bool operator==(const Matrix& p) const;
        bool operator!=(const Matrix& p) const;

        // Matrix row_echelon() const;
        // Matrix reduced_row_echelon() const;
        // T det() const;
        // Matrix inv() const;
        // Matrix pow() const;
        //
        // string to_string() const;
        // string to_latex() const;
    };
}

// template <typename T>
// std::ostream& operator<<(std::ostream& s, const std::Matrix<T>& x);
// template <typename T>
// std::istream& operator>>(std::istream& s,       std::Matrix<T>& x);

#endif
