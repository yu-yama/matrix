#include <iostream>
#include "matrixIntArray.h"

using namespace std;

int main() {
    Matrix a(2, 2);
    cout << a.to_string();
    cin >> a;
    cout << a.to_string();
    cout << a << endl;
    Matrix b({{5, 6}, {7, 8}});
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) cout << (a.at(i, j) = i + j) << '\n';
    // a += b;
    // a + b;
    cout << "000010\n";
    cout << &a << '\n';
    Matrix *c = &a;
    cout << c << '\n';
    // Matrix<int> c(b);
    cout << "000030\n";
    cout << (*(&a)).to_string();
    cout << "000050\n";
    a *= a;
    cout << "000070\n";
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) cout << a.at(i, j) << '\n';
    a *= b;
    cout << "000090\n";
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) cout << a.at(i, j) << '\n';
    return 0;
}
