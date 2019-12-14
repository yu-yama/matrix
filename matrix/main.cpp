#include <iostream>
#include "matrix.h"

using namespace std;

int main() {
    Matrix<int> z();
    Matrix<int> v(2);
    Matrix<int> a(2, 2);
    cout << "000010\n" << a.to_string();
    cout << a;
    // cin >> a;
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
    cout << "000220\n" << a.pow(3) << a;
    cout << "000230\n";
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) cout << a.at(i, j) << '\n';
    a *= b;
    cout << "000250\n" << a;
    cout << "000270\n";
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) cout << a.at(i, j) << '\n';
    Matrix<double> e({{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}});
    cout << "000290\n" << e << e.gauss() << e.gauss_jordan();
    Matrix<double> f({{2, -1, 0}, {-1, 2, -1}, {3, 0, -1}});
    cout << "000310\n" << f << f.gauss() << f.gauss_jordan();
    Matrix<double> g({{3, 5}, {-2, 4}});
    cout << "000330\n" << g << g.det() << endl;
    Matrix<double> h({{-2, 1, 4}, {3, 5, -7}, {1, 6, 2}});
    cout << "000350\n" << h << h.det() << endl;

    Matrix<double> aa1({{1, 3, -2, 0, 2, 0}, {2, 6, -5, -2, 4, -3}, {0, 0, 5, 10, 0, 15}, {2, 6, 0, 8, 4, 18}});
    cout << "000360\n" << aa1.gauss() << aa1.gauss_jordan();
    Matrix<double> aa2({{0}, {-1}, {5}, {6}});
    AugmentedMatrix<double> aa(aa1, aa2);
    cout << "000370\n" << aa.left() << aa.right();
    cout << "000390\n" << aa.gauss().left() << aa.gauss().right();
    cout << "000410\n" << aa.gauss_jordan().left() << aa.gauss_jordan().right();

    Matrix<double> ab1({{1, 2, 3}, {2, 5, 3}, {1, 0, 8}});
    Matrix<double> ab2({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
    AugmentedMatrix<double> ab(ab1, ab2);
    cout << "000430\n" << ab.gauss_jordan().left() << ab.gauss_jordan().right();
    cout << "000450\n" << ab1.inv();

    // terminate called after throwing an instance of 'std::invalid_argument'
    // what():  Argument is not invertible
    // Matrix<double> ac1({{1, 6, 4}, {2, 4, -1}, {-1, 2, 5}});
    // cout << "000470\n" << ac1.inv();
    return 0;
}
