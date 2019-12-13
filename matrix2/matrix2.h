#ifndef ___MATRIX_INCLUDED_
#define ___MATRIX_INCLUDED_

#define ___MATRIXINTARRAY_DEBUG_ false

#include <string>
#include <iostream>
#include <vector>

template <typename T>
class Matrix;

template <typename T>
std::ostream& operator<<(std::ostream& s, const Matrix<T>& x);

template <typename T>
class Matrix {
    // typedef typename std::vector<T>::size_type size_type;
private:
    std::vector<T> mat;
    typename std::vector<T>::size_type n, m;

public:
    Matrix();
    Matrix(typename std::vector<T>::size_type n);
    Matrix(typename std::vector<T>::size_type n, typename std::vector<T>::size_type m);
    Matrix(std::vector< std::vector<T> > matData);
    Matrix(std::vector<T> matData, typename std::vector<T>::size_type n);
    Matrix(const Matrix& p);

    typename std::vector<T>::size_type row()    const { return n; }
    typename std::vector<T>::size_type column() const { return m; }

    std::vector<T> at_row   (typename std::vector<T>::size_type r) const;
    std::vector<T> at_column(typename std::vector<T>::size_type c) const;
    T& at(typename std::vector<T>::size_type r, typename std::vector<T>::size_type c);
    const T& at(typename std::vector<T>::size_type r, typename std::vector<T>::size_type c) const;

    void resize(typename std::vector<T>::size_type nn, typename std::vector<T>::size_type mm);
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
    std::string to_string() const;
    // std::string to_latex() const;

    friend std::ostream& operator<< <T>(std::ostream& s, const Matrix<T>& x);

    // friend std::istream& operator>>(std::istream& s,       Matrix& x) {
    //     if (___MATRIXINTARRAY_DEBUG_) cout << "Start: std::istream >>\n";
    //     for (typename std::vector<T>::size_type i = 0; i < x.row(); ++i) for (typename std::vector<T>::size_type j = 0; j < x.column(); ++j) s >> x.at(i, j);
    //     return s;
    //     if (___MATRIXINTARRAY_DEBUG_) cout << "End r: std::istream >>\n";
    // }
};

// template <typename T>
// std::std::ostream& operator<<(std::std::ostream& s, const std::Matrix<T>& x);
// template <typename T>
// std::std::istream& operator>>(std::std::istream& s,       std::Matrix<T>& x);

#endif
