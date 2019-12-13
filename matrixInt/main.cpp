#include <iostream>
#include "matrixInt.h"

using namespace std;

int main() {
    Matrix a(2, 2);
    cout << a.to_string();
    // cout << a;
    Matrix b({{5, 6}, {7, 8}});
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j)
    cout << (a.at(i, j) = i + j) << '\n';
    a += b;
    a + b;
    cout << "a\n";
    cout << &a << '\n';
    Matrix *c = &a;
    cout << c << '\n';
    // Matrix<int> c(b);
    cout << "a\n";
    a = *(&a);
    cout << "a\n";
    a *= a;
    cout << "a\n";
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) cout << a.at(i, j) << '\n';
    a *= b;
    cout << "a\n";
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) cout << a.at(i, j) << '\n';
    return 0;
}
