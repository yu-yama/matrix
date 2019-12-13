#ifndef ___MATRIX_INCLUDED_
#define ___MATRIX_INCLUDED_

#include <string>
#include <sstream>
#include <iostream>
#include <vector>

namespace std {
    template <typename T>
    class Matrix {
        // typedef typename vector< vector<T> >::size_type size_type;
    private:
        vector< vector<T> > mat;
        typename vector< vector<T> >::size_type n, m;

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
        for (typename vector< vector<T> >::size_type i = 0; i < n; ++i) for (typename vector< vector<T> >::size_type j = 0; j < m; ++j) for (typename vector<T>::size_type k = 0; k < p.m; ++k) ans.mat.at(i).at(k) += mat.at(i).at(j) * p.mat.at(j).at(k);
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
        return !((*this) == p);
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

}

// template <typename T>
// std::ostream& operator<<(std::ostream& s, const std::Matrix<T>& x);
// template <typename T>
// std::istream& operator>>(std::istream& s,       std::Matrix<T>& x);

#endif
