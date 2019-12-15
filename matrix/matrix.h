#ifndef ___MATRIX_INCLUDED_
#define ___MATRIX_INCLUDED_

#define ___MATRIXINTARRAY_DEBUG_ false

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <tuple>

template <typename T>
class Matrix;

template <typename T>
std::ostream& operator<<(std::ostream&, const Matrix<T>&);

template <typename T>
std::istream& operator>>(std::istream&,       Matrix<T>&);

template <typename T>
class Matrix {
    // typedef typename std::vector<T>::size_type size_type;
private:
    std::vector<T> mat = {};
    typename std::vector<T>::size_type n = 0, m = 0;

    std::tuple<Matrix, bool, typename std::vector<T>::size_type> gauss_count() const;
    Matrix basis() const;

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

    void push_row(std::vector<T> p);
    void pop_row();
    void push_column(std::vector<T> p);
    void pop_column();

    void append_rows(Matrix p);
    void append_columns(Matrix p);

    void resize(typename std::vector<T>::size_type nn);
    void resize(typename std::vector<T>::size_type nn, typename std::vector<T>::size_type mm);
    Matrix  operator+() const;
    Matrix  operator-() const;

    Matrix& operator=(Matrix p);

    Matrix& operator+=(Matrix p);
    Matrix& operator-=(Matrix p);
    Matrix& operator*=(Matrix p);
    Matrix& operator*=(T p);
    Matrix& operator/=(Matrix p);
    Matrix& operator/=(T p);

    Matrix  operator+(Matrix p) const;
    Matrix  operator-(Matrix p) const;
    Matrix  operator*(Matrix p) const;
    Matrix  operator*(T p) const;
    // template <typename S>
    // friend Matrix<S> operator*(S n, const Matrix<S>& p);
    Matrix  operator/(Matrix p) const;
    Matrix  operator/(T p) const;
    // template <typename S>
    // friend Matrix<S> operator/(S n, const Matrix<S>& p);

    bool operator==(const Matrix& p) const;
    bool operator!=(const Matrix& p) const;

    bool is_zero() const;

    Matrix gauss() const;
    Matrix gauss_jordan() const;
    T det() const;
    Matrix inv() const;
    Matrix pow(int r) const;
    T trace() const;
    T rank() const;

    T norm_squared() const;
    double norm() const;
    double distance(Matrix p) const;
    double angle(Matrix p) const;

    T dot(Matrix p) const;
    bool orthogonal(Matrix p) const;
    Matrix projection(Matrix p) const;
    // Matrix orthonormal() const;
    Matrix cross(Matrix p) const;

    T minor(typename std::vector<T>::size_type y, typename std::vector<T>::size_type x) const;
    T cofactor(typename std::vector<T>::size_type y, typename std::vector<T>::size_type x) const;
    Matrix<T> cofactor_matrix() const;
    Matrix<T> adjugate() const;

    Matrix transpose() const;

    std::string to_string() const;
    // std::string to_latex() const;

    friend std::ostream& operator<< <> (std::ostream&, const Matrix&);
    friend std::istream& operator>> <> (std::istream&,       Matrix&);
};

template <typename T>
std::ostream& operator<<(std::ostream& s, const Matrix<T>& x) {
    if (___MATRIXINTARRAY_DEBUG_) std::cout << "Start: ostream <<\n";
    if (___MATRIXINTARRAY_DEBUG_) std::cout << "End r: ostream <<\n";
    return s << x.to_string();
}

template <typename T>
std::istream& operator>>(std::istream& s,       Matrix<T>& x) {
    if (___MATRIXINTARRAY_DEBUG_) std::cout << "Start: istream >>\n";
    for (typename std::vector<T>::size_type i = 0; i < x.row(); ++i) for (typename std::vector<T>::size_type j = 0; j < x.column(); ++j) s >> x.at(i, j);
    if (___MATRIXINTARRAY_DEBUG_) std::cout << "End r: istream >>\n";
    return s;
}

template <typename T>
Matrix<T> identity_matrix(typename std::vector<T>::size_type n, T val) {
    Matrix<T> t(n, n);
    for (typename std::vector<T>::size_type i = 0; i < n; ++i) t.at(i, i) = (T)1;
    return t;
}

template <typename T>
class AugmentedMatrix {
private:
    Matrix<T> matl, matr;
    typename std::vector<T>::size_type n, ml, mr;

public:
    AugmentedMatrix(typename std::vector<T>::size_type n, typename std::vector<T>::size_type ml);
    AugmentedMatrix(typename std::vector<T>::size_type n, typename std::vector<T>::size_type ml, typename std::vector<T>::size_type mr);
    AugmentedMatrix(Matrix<T> mmatl);
    AugmentedMatrix(Matrix<T> mmatl, Matrix<T> mmatr);
    AugmentedMatrix(const AugmentedMatrix& p);

    typename std::vector<T>::size_type row() const { return n; }
    typename std::pair<typename std::vector<T>::size_type, typename std::vector<T>::size_type> column() const { return {ml, mr}; }

    std::pair< std::vector<T>, std::vector<T> > at_row (typename std::vector<T>::size_type r) const;
    std::vector<T> at_column (typename std::vector<T>::size_type c) const;
    std::vector<T> at_columnl(typename std::vector<T>::size_type c) const;
    std::vector<T> at_columnr(typename std::vector<T>::size_type c) const;
    T& at(typename std::vector<T>::size_type r, typename std::vector<T>::size_type c);
    const T& at(typename std::vector<T>::size_type r, typename std::vector<T>::size_type c) const;

    Matrix<T> left() const;
    Matrix<T> right() const;

    void resize(typename std::vector<T>::size_type nn, typename std::vector<T>::size_type mml, typename std::vector<T>::size_type mmr);

    bool operator==(const AugmentedMatrix& p) const;
    bool operator!=(const AugmentedMatrix& p) const;

    AugmentedMatrix gauss() const;
    AugmentedMatrix gauss_jordan() const;
};

#endif
