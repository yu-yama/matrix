#include <iostream>
#include "matrix2.h"

using namespace std;

int main() {
    Matrix<int> z();
    Matrix<int> v(2);
    Matrix<int> a(2, 2);
    cout << "000010\n" << a.to_string();
    cout << a;
    cin >> a;
    cout << "000020\n" << a;
    Matrix<int> b({{5, 6}, {7, 8}});
    cout << "000030\n";
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) cout << (a.at(i, j) = i + j) << '\n';
    a += b;
    cout << "000050\n" << a.to_string();
    // cout << "000070\n" << (2 * a).to_string();
    cout << "000070\n" << (a * 2);
    // cout << "000080\n" << (2 * a);
    cout << "000090\n" << (a + b);
    cout << "000110\n";
    cout << &a << '\n';
    Matrix<int> *c = &a;
    cout << "000130\n";
    cout << c << '\n';
    Matrix<int> d(a);
    cout << "000150\n" << d;
    cout << "000170\n" << &d << '\n';
    a = *(&a);
    cout << "000190\n" << a;
    a *= a;
    cout << "000210\n" << a;
    cout << "000230\n";
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) cout << a.at(i, j) << '\n';
    a *= b;
    cout << "000250\n" << a;
    cout << "000270\n";
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) cout << a.at(i, j) << '\n';
    return 0;
}
