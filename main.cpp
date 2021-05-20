#include "matrix.hpp"
#include "vector.hpp"
#include "linear_system.hpp"
using std::cout;
int main()
{
    // int c[2] = {3, 3};
    // // Matrix<int> a = Matrix<int>(20, 20);

    // // std::cout << a;

    // auto b = Matrix<double>({{1, 2, 3, 4, 5}, {4, 5, 6, 7, 8}, {7, 8, 3, 10, 3}});
    // auto sim = Matrix<double>({{2, 3, 4}, {1, 1, 1}, {4, 2, 1}});
    // auto v = Vector<double>({5, 3, 7});
    // Vector<double> vcl({{2, 2, 2}});
    // auto ab = Vector<double>(vcl);
    // auto x = dot(ab, vcl);
    // Matrix<double> xx(ab);
    // auto cc = ab * (2); // vip compiler , from 2 it can create a Vector<double>2;
    // // cout << cc << cc.length();
    // std::cerr << v;
    // auto xy = sim.slice({0, 1}, {3, 1});
    // v = sim.slice({0, 1}, {3, 1});
    // cout << xy;
    // cout << "v is " << v;
    // // cout
    // //     << v;
    Vector<double> v({3, 3, 3});
    // auto ls = LinearSystem<double>(sim, v);
    Matrix<double> sym({{1, 2, 3}, {2, 4, 1}, {3, 1, 5}});
    Matrix<double> sym1({{-11, 2}, {2, 3}, {2, -1}});
    Matrix<double> sym2({{12, 10, 15}, {-4, 15, 9}, {-9, 13, 12}, {-8, 9, 10}});
    auto lmao = sym2.pseudoinverse();
    // PosSymLinSystem<double> psls(sym, v);
    // auto lmao = psls.solve();
    cout << lmao;
    // std::cout << sym;
    // auto solution = ls.solve();
    // cout << "solution is " << solution;

    // auto inv = sim.inverse();
    // cout << "inverse is " << inv;
    // cout << "b is " << b;
    // auto ge = gaussianElimination(b);
    // cout << "ge form is" << ge;
    // backwardSubstitution(ge);
    // cout << "rref form is " << ge;
    // auto n = b.slice({0, 0}, {3, 1});
    // cout << n;
    // auto m = b.slice({0, 3}, {3, 1});
    // m(2, 0) = 10;
    // auto oneDmat = Matrix<double>({
    //     {1},
    //     {2},
    //     {3},
    //     {4},
    // });
    // // m = oneDmat;
    // b.slice({0, 0}, {3, 1}) = b.slice({0, 1}, {3, 1});
    // cout << "new b is " << b;
    // auto haha1 = b.transpose();
    // cout << haha1;
    // auto haha = b.slice({0, 0}, {1, 4}) * oneDmat;
    // cout << haha;
    // auto gg = Matrix<int>({{1, 2, 3, 4}, {4, 5, 6, 7}, {7, 8, 9, 10}});
    // cout << "gg is" << gg;
    // auto hehe = gg / 2;
    // gg = hehe;

    // cout << gg;
    // cout << "checkpoint";
    // auto vec = new Vector<char>({'d', 'c', 'm'});
    // auto vec1 = vec->transpose();
    // cout << vec1;

    // infer_shape<int>({{1, 2, 3}, {4, 5, 6}}, c);
    // std::cout << c[0] << c[1];
    // return 0;
}