#ifndef _LS_H_
#define _LS_H_
#include <iostream>
#include "matrix.hpp"
#include "vector.hpp"
#include "matrix_utils.hpp"
template <typename T>
class LinearSystem
{
protected:
    Matrix<T> *A;
    Vector<T> *b;
    int mSize;

public:
    LinearSystem(Matrix<T> &mat, Vector<T> &vec)
    {
        check((mat.get_shape()[0] == vec.length()), "invalid Linear System");
        A = &mat;
        b = &vec;
        mSize = mat.get_shape()[1];
    };
    virtual Vector<T> solve()
    {
        check((A->get_shape()[0] >= A->get_shape()[1]), "this system has multiple solution");

        Matrix<T> augmented_mat = Matrix<T>(A->get_shape()[0], A->get_shape()[1] + 1);

        augmented_mat.slice({0, 0}, {A->get_shape()[0], A->get_shape()[1]}) = *A;
        // auto transpose_b = b->transpose();
        augmented_mat.slice({0, A->get_shape()[1]}, {A->get_shape()[0], 1}) = *b;

        Matrix<T> ge_mat = gaussianElimination(augmented_mat);
        bool has_solution = backwardSubstitution(ge_mat);

        check((has_solution), "this system is inconsistent");
        // check for multiple solutions here, we will also not deal with multiple solutions
        auto identity = Identity_matrix<T>(b->length());
        check(((ge_mat.slice({0, 0}, {A->get_shape()[1], A->get_shape()[1]}) == identity)), "this system has multiple solution");

        Vector<T> solution = Vector<T>(A->get_shape()[1]);
        solution = ge_mat.slice({0, ge_mat.get_shape()[1] - 1}, {A->get_shape()[1], 1}).transpose();
        return solution;
    };
};

template <typename T>
class PosSymLinSystem : public LinearSystem<T>
{
public:
    PosSymLinSystem(Matrix<T> &mat, Vector<T> &vec) : LinearSystem<T>(mat, vec)
    {
        check((mat.transpose() == mat), "the matrix is not symmetric");
    };
    Vector<T> solve()
    {
        auto A = *(this->A);
        auto b = *(this->b);
        Vector<T> solution(this->mSize);

        double TOLERANCE = 1.0e-10;
        int n = this->mSize;

        Vector<T> residual(b);
        Vector<T> search_direction(residual);
        debug("step");

        while (residual.norm() >= TOLERANCE)
        {
            debug("step");
            Vector<T> Rold = residual; // Store previous residual
            Vector<T> AP = A * search_direction;
            debug("step");
            // double y =
            //     std::cerr << x;
            T alpha = dot<T>(residual, residual) / dot<T>(search_direction, AP);
            solution = solution + (search_direction * alpha); // Next estimate of solution
            residual = b - A * solution;                      // Residual

            T beta = dot<T>(residual, residual) / dot<T>(Rold, Rold);
            search_direction = residual + search_direction * beta; // Next gradient
        }
        return solution;
    };
};
#endif