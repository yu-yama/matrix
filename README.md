# Matrix Library for C++

![GitHub Action Build Status on Ubuntu](https://github.com/yu-yama/matrix/workflows/Matrix%20Library%20Build%20%26%20Run%20Check/badge.svg)

This library provides matrices computation features with a header file, `src/matrix.h`, and a cpp file `src/matrix.cpp`.

## License Information
The license for this repository is given in the LICENSE file.

## I. How to use
Add `#include "matrix.h"` in the preamble of your source code, then link your object file with the object file of `src/matrix.cpp`.

To compile the test code, `src/test.cpp` (which is not quite user-friendly yet), type `make test` in your terminal. Type `make run` or `./TEST_CPP` to run.

## II. Tested environment and compilation options
The sources in the `src` directory are tested using `g++-9 (Homebrew GCC 9.2.0_2) 9.2.0` on `macOS Catalina` with compilation options `-std=c++17 -O2`.

## III. Classes provided
Template classes `Matrix<T>` and `AugmentedMatrix<T>` are provided. In `src/test.cpp`, the classes are explicitly specialized for `short`, `int`, `long`, `long long`, `float`, `double`, and `long double`. Since a number of methods involve negative numbers, it is strongly discouraged to use `unsigned` types.

## IV. List of features
### 1. Matrix template class `Matrix<T>`
#### a. Member Variables
* `vector<T> mat` – stores the elements of a matrix
* `typename vector<T>::size_type n, m` – stores the number of rows, `n`, and columns, `m`

#### b. Member Functions
##### i. Constructors
* Default constructor `Matrix()` – generates a matrix with size `0 x 0`
* Constructor `Matrix(n)` – generates a column vector with the specified size, `n`
* Constructor `Matrix(n, m)` – generates an `n x m` matrix
* Constructor `Matrix(vector<T>)` – generates the specialized column vector
* Constructor `Matrix(vector< vector<T> >)` – generates a matrix that the specified vector is converted into. The number of columns will be equal to the greatest number of columns in the specified vector, and elements not contained in the specified vector will be filled with `0`.
* Constructor `Matrix(vector<T>, n)` – generates a matrix that the specified vector is converted into, by folding the vector by every n elements
* Copy constructor `Matrix(Matrix&)` – generates an identical matrix (deep copy)

##### ii. Accessors
###### ii-a. Getters
* `row()` – returns the number of rows
* `column()` – returns the number of columns
* `at_row(y)` – returns the `y`-th row of a matrix as a vector; throws an `out_of_range` exception when the `y`-th row does not exist
* `at_column(x)` – returns the `x`-th column of a matrix as a vector; throws an `out_of_range` exception when the `x`-th column does not exist
* `at(y, x)` – returns the element `(y, x)` of a matrix; throws an `out_of_range` exception when the element `(y, x)` does not exist

###### ii-b. Setters
* `push_row(vector<T>)` – adds the specified vector as a row at the bottom of a matrix; throws an `invalid_argument` exception when the number of columns does not match
* `pop_row()` – removes the bottom row of a matrix; throws `length_error` when the matrix does not have a row
* `push_column(vector<T>)` – adds the specified vector as a column at the rightmost of a matrix; throws an `invalid_argument` exception when the number of rows does not match
* `pop_column()` – removes the rightmost column of a matrix; throws `length_error` when the matrix does not have a column
* `append_rows(Matrix<T>)` – adds the matrix as rows at the bottom of a matrix; throws an `invalid_argument` exception when the number of columns does not match
* `remove_rows(p)` – removes `p` rows in the bottom of a matrix; throws `length_error` when the matrix does not have at least `p` rows
* `append_columns(Matrix<T>)` – adds the matrix as columns at the rightmost of a matrix; throws an `invalid_argument` exception when the number of rows does not match
* `remove_columns(p)` – removes `p` columns at the rightmost of a matrix; throws `length_error` when the matrix does not have at least `p` columns
* `insert_row(p, vector<T>)` – inserts the vector as the `p`-th row of a matrix, and consequently, increases the indexes of rows after the inserted row by `1`; throws an `invalid_argument` exception when the number of columns does not match
* `delete_row(p)` – deletes the `p`-th row of a matrix, and consequently, decreases the indexes of rows after the inserted row by `1`; throws `length_error` when the matrix does not have a row
* `insert_column(p, vector<T>)` – inserts the vector as the `p`-th column of a matrix, and consequently, increases the indexes of columns after the inserted column by `1`; throws an `invalid_argument` exception when the number of rows does not match
* `delete_column(p)` – deletes the `p`-th column of a matrix, and consequently, decreases the indexes of columns after the inserted column by `1`; throws a `length_error` exception when the matrix does not have a column
* `insert_rows(p, Matrix<T>)` – inserts the matrix as consecutive rows from the `p`-th row of a matrix, and consequently, increases the indexes of rows after the inserted row by the number of rows of the specified matrix; throws an `invalid_argument` exception when the number of columns does not match
* `delete_rows(p, q)` – deletes `q` consecutive rows from the `p`-th row of a matrix, and consequently, decreases the indexes of rows after the inserted row by `q`; throws a `length_error` exception when the matrix does not have at least `q` rows
* `insert_columns(p, Matrix<T>)` – inserts the matrix as consecutive columns from the `p`-th column of a matrix, and consequently, increases the indexes of columns after the inserted column by the number of rows of the specified matrix; throws an `invalid_argument` exception when the number of rows does not match
* `delete_columns(p, q)` – deletes `q` consecutive columns from the `p`-th column of a matrix, and consequently, decreases the indexes of columns after the inserted column by `q`; throws a `length_error` exception when the matrix does not have at least `q` columns
* `resize(nn)` – changes the number of rows of a matrix to `nn` and copies the elements contained in both matrices before and after this method
* `resize(nn, mm)` – changes the number of rows of a matrix to `nn` and the number of columns to `mm`, and copies the elements contained in both matrices before and after this method
* `resize(nn, mm, bool copy)` – equivalent to `resize(nn, mm)` if `copy = true`, resizes a matrix without copying otherwise

##### iii. Unary Operators
* `+` – returns the matrix itself
* `-` – returns the matrix times `-1`

##### iv. Binary Operators
###### iv-a. Assignment Operators
* `= Matrix<T>` – assigns the specified matrix
* `+= Matrix<T>` – adds the corresponding elements; throws an `invalid_argument` exception if the dimensions of two matrices do not match
* `-= Matrix<T>` – subtracts the corresponding elements; throws an `invalid_argument` exception if the dimensions of two matrices do not match
* `*= Matrix<T>` – equivalent to `= *this * Matrix<T>`
* `*= p` – computes a scalar multiplication
* `/= Matrix<T>` – equivalent to `= *this / Matrix<T>`
* `/= p` – computes a scalar division (multiplication of a reciprocal)

###### iv-b. Arithmetic Operators
* `+ Matrix<T>` – computes the sum of two matrices using `+= Matrix<T>`
* `- Matrix<T>` – computes the difference of two matrices using `-= Matrix<T>`
* `* Matrix<T>` – computes a matrix multiplication; throws an `invalid_argument` exception if the number of columns of one matrix and the number of rows of another matrix does not match
* `* p` – computes a scalar multiplication using `*= Matrix<T>`
* `/ Matrix<T>` – computes a matrix multiplication with an inverse (equivalent to `* Matrix<T>.inv()`)
* `/ p` – computes a scalar division (multiplication of a reciprocal)

###### iv-c. Comparison Operators
* `== Matrix<T>` – returns `true` if two matrices are equal (have the same number of rows and columns, and identical elements), and `false` otherwise
* `!= Matrix<T>` – equivalent to `!(== Matrix<T>)`

##### v. Arithmetic Methods
* `pow(int r)` – returns the matrix to the `r`-th power; throws an `invalid_argument` exception when the matrix is not a square matrix

##### vi. Comparison Methods
* `empty()` – returns `true` if the size of a matrix is zero, and `false` otherwise
* `is_zero()` – returns `true` if all the elements of a matrix is zero, and `false` otherwise
* `is_identity()` – returns `true` if a matrix is an identity matrix, and `false` otherwise; throws an `invalid_argument` exception when the matrix is not a square matrix

##### vii. Eliminations
* `gauss()` – returns a matrix gained after applying Gaussian elimination
* `gauss_jordan()` – returns a matrix gained after applying Gauss-Jordan elimination

##### viii. Attributes
* `det()` – returns the determinant of a matrix; throws an `invalid_argument` exception when the matrix is not a square matrix
* `inv()` – returns the inverse of a matrix; throws an `invalid_argument` exception when the matrix is not a square matrix
* `trace()` – returns the trace of a matrix; throws an `invalid_argument` exception when the matrix is not a square matrix
* `rank()` – returns the rank of a matrix
* `norm_squared()` – returns the square of the norm of a vector; throws an `invalid_argument` exception when the matrix is neither a column nor row vector
* `norm()` – returns the square root of `norm_squared()` in `double` type
* `distance(Matrix<T>)` – returns the norm of `*this - Matrix<T>` in `double` type
* `angle(Matrix<T>)` – returns the angle between two matrices in `double` type (in radian)
* `dot(Matrix<T>)` – returns a dot product of two matrices; throws an `invalid_argument` exception when two matrices are neither both row vectors nor both column vectors, or have different number of elements
* `orthogonal(Matrix<T>)` – returns whether two matrices are orthogonal or not
* `projection(Matrix<T>)` – returns a projection of a matrix on the specified matrix
* `orthonormal()` – returns an orthonormal basis of a matrix's column space
* `cross(Matrix<T>)` – returns a cross product of two matrices; throws an `invalid_argument` exception when two matrices are neither both `3 x 1` column vectors nor `1 x 3` row vectors
* `minor_at(y, x)` – returns a minor, _M<sub>yx</sub>_; throws an `invalid_argument` exception when the matrix is not a square matrix
* `cofactor_at(y, x)` – returns a cofactor, _C<sub>yx</sub>_
* `cofactor_matrix()` – returns the cofactor matrix
* `adjugate()` – returns the adjugate
* `transpose()` – returns the transpose

##### ix. Casting Methods
* `to_1d_vector()` – returns a one-dimensional vector representing a row or column vector; throws an `invalid_argument` exception when the matrix is neither row nor column vector
* `to_2d_vector()` – returns a two-dimensional vector representing the matrix
* `to_string()` – returns a string representing a matrix in CSV format

#### c. Non-Member Functions
##### i. Stream Operators
* `ostream << Matrix<T>` – writes `Matrix<T>.to_string()`
* `istream >> Matrix<T>` – reads consecutive numbers from the stream and write them in the matrix

##### ii. Identity Matrix
* `identity_matrix(n, val)` – returns the identity matrix of degree `n` in the type of `val`

##### iii. Vector Calculations
* `vdiff(vector<T> p, vector<T> q)` – returns `p - q`; throws an `invalid_argument` exception if `p` and `q` have different sizes
* `vscalar_mult(vector<T> p, T q)` – returns `q * p`
* `vnorm_squared(vector<T> p)` – returns the square of the norm of `p`
* `vnorm(vector<T> p)` – returns the norm of `p` in `double` type
* `vdistance(vector<T> p, vector<T> q)` – returns the norm of `p - q` in `double` type
* `vangle(vector<T> p, vector<T> q)` – returns the angle between `p` and `q` in `double` type (in radian)
* `vdot(vector<T> p, vector<T> q)` – returns the dot product of p and q; throws an `invalid_argument` exception if `p` and `q` have different sizes
* `vorthogonal(vector<T> p, vector<T> q)` – returns whether `p` and `q` are orthogonal or not
* `vprojection(vector<T> p, vector<T> q)` – returns the projection of `p` on `q`
* `vcross(vector<T> p, vector<T> q)` – returns the cross product of p and q; throws an `invalid_argument` exception if either `p` or `q` does not have exactly three elements

### 2. Augmented Matrix template class `AugmentedMatrix<T>`
#### a. Member Variables
* `Matrix<T> matl, matr` – stores the matrices on the left, `matl`, and the right, `matr`, side
* `typename vector<T>::size_type n, ml, mr` – stores the number of rows, `n`, columns of `matl`, `ml`, and columns of `matr`, `mr`

#### b. Member Functions
##### i. Constructors
* No default constructor provided
* Constructor `AugmentedMatrix(n, ml)` – generates an augmented matrix of an `n x ml` matrix and an `n x 1` column vector
* Constructor `AugmentedMatrix(n, ml, mr)` – generates an augmented matrix of an `n x ml` matrix and an `n x mr` matrix
* Constructor `AugmentedMatrix(Matrix<T>)` – generates an augmented matrix of the specified matrix and a column vector with the same number of rows
* Constructor `AugmentedMatrix(Matrix<T>, Matrix<T>)` – generates an augmented matrix of the two specified matrices; throws an `invalid_argument` exception when the number of rows of the two specified matrices does not match
* Copy Constructor `AugmentedMatrix(AugmentedMatrix&)` – generates an identical augmented matrix (deep copy)

##### ii. Accessors
###### ii-a. Getters
* `row()` – returns the number of rows
* `column()` – returns a pair of the number of columns of left and right matrices
* `at_row(y)` – returns a pair of the `y`-th row of the left and right matrices; throws `out_of_range` exception if the `y`-th row does not exist
* `at_column(x)` – returns `at_columnl(x)` if `0 <= x < ml` and `at_columnr(x - ml)` if `ml <= x < ml + mr`; throws an `out_of_range` exception otherwise
* `at_columnl(x)` – returns the `x`-th column of the left matrix
* `at_columnr(x)` – returns the `x`-th column of the right matrix
* `at(y, x)` – returns the element `(y, x)` of the left matrix if `0 <= x < ml` and `(y, x - ml)` of the right matrix if `ml <= x < ml + mr`; throws an `out_of_range` exception if `y`-th row does not exist, or `x` is not within `0` and `ml + mr`
* `left()` – returns the left matrix
* `right()` – returns the right matrix

###### ii-b. Setter
* `resize(nn, mml, mmr)` – resizes the left matrix to `nn x mml` and the right matrix to `nn x mmr`

##### iii. Comparison Operators
* `== AugmentedMatrix<T>` – returns `true` if both left and right matrices are equal, and `false` otherwise
* `!= AugmentedMatrix<T>` – equivalent to `!(== AugmentedMatrix<T>)`

##### iv. Eliminations
* `gauss()` – returns an augmented matrix after applying Gaussian elimination
* `gauss_jordan()` – returns an augmented matrix after applying Gauss-Jordan elimination

## V. Features to be added in the (near) future
* QR-decomposition, QR-algorithm
* Eigenvalue, eigenvector, eigenbasis, and eigenspace
* Diagonalization
