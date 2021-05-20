#ifndef _MATRIX_UTILS_H
#define _MATRIX_UTILS_H
#include "matrix.hpp"
#include "vector.hpp"
template <typename T>
Matrix<T> Identity_matrix(int size)
{
    Matrix<T> new_mat = Matrix<T>(size, size);
    for (int i = 0; i < size; i++)
    {
        new_mat(i, i) = 1;
    }
    return new_mat;
}
template <typename T>
Matrix<T> gaussianElimination(Matrix_base<T> &mat)
{
    // std::cout << "original mat " << mat;
    Matrix<T> placeholder = mat;
    // std::cout << "placeholder " << placeholder;
    int width = mat.get_shape()[1];
    int height = mat.get_shape()[0];
    int index_not_zero = 0;
    for (int j = 0; j < height; j++)
    {
        if (std::abs(mat(j, 0)) > 0)
        {
            index_not_zero = j;
            break;
        }
    }
    if (mat(index_not_zero, 0) != 0)
    {
        // swap the first row with the row having greatest 1st element
        View<T> first_row = placeholder.slice({0, 0}, {1, width});
        Matrix<T> temp = Matrix<T>(first_row);
        first_row = placeholder.slice({index_not_zero, 0}, {1, width});
        placeholder.slice({index_not_zero, 0}, {1, width}) = temp;
        // std::cout << "after swap:" << placeholder;
        // std::cout << "placeholder 1st element is " << placeholder(index_not_zero, 0) << "\n";
        first_row = first_row / placeholder(0, 0);
        // std::cout << "after divide:" << placeholder;

        for (int i = 1; i < height; i++)
        {
            T alpha = placeholder(i, 0);
            placeholder.slice({i, 0}, {1, width}) = placeholder.slice({i, 0}, {1, width}) - first_row * alpha;
        }
        // std::cout << "after normalize: " << placeholder;
        if (height == 1 || width == 1)
        {
            return placeholder;
        }
        View<T> recursive_GE_submat = placeholder.slice({1, 1}, {height - 1, width - 1});
        auto gaussian_eliminated_submat = gaussianElimination(recursive_GE_submat);

        recursive_GE_submat = gaussian_eliminated_submat;
        return placeholder;
    }
    if (height == 1 || width == 1)
    {
        return placeholder;
    }
    View<T> recursive_GE_submat = placeholder.slice({0, 1}, {height, width - 1});

    auto gaussian_eliminated_submat = gaussianElimination(recursive_GE_submat);
    recursive_GE_submat = gaussian_eliminated_submat;
    return placeholder;
}
template <typename T>
bool backwardSubstitution(Matrix<T> &mat)
{
    bool has_solution = true;
    int width = mat.get_shape()[1];
    int height = mat.get_shape()[0];
    for (int i = height - 1; i > 0; i--)
    {
        int leading_index = -1;
        // find the leading 1 in the matrix excluding the last column(you know why)
        for (int j = 0; j < width - 1; j++)
        {
            // std::cout << "i,j is" << i << j << "\n";
            if (mat(i, j) == 1)
            {
                leading_index = j;
                break;
            }
        }
        if (leading_index == -1)
        {
            if (mat(i, width - 1) != 0)
            {
                has_solution = false;
            }
            continue;
        } // no leading 1 in this row, skip this.
        for (int k = i - 1; k >= 0; k--)
        {
            mat.slice({k, 0}, {1, width}) = mat.slice({k, 0}, {1, width}) - mat.slice({i, 0}, {1, width}) * mat(k, leading_index);
        }
    }
    return has_solution;
}
template <typename T>
Matrix<T> extendedGaussianElimination(Matrix_base<T> &mat, T *det_quotient)
{
    Matrix<T> placeholder = mat;
    int width = mat.get_shape()[1];
    int height = mat.get_shape()[0];

    int index_not_zero = 0;
    for (int j = 0; j < height; j++)
    {
        if (std::abs(mat(j, 0)) > 0)
        {
            index_not_zero = j;
            break;
        }
    }
    if (mat(index_not_zero, 0) != 0)
    {
        // swap the first row with the row having greatest 1st element
        View<T> first_row = placeholder.slice({0, 0}, {1, width});
        Matrix<T> temp = Matrix<T>(first_row);
        first_row = placeholder.slice({index_not_zero, 0}, {1, width});
        placeholder.slice({index_not_zero, 0}, {1, width}) = temp;
        // std::cout << "after swap:" << placeholder;
        // std::cout << "placeholder 1st element is " << placeholder(index_not_zero, 0) << "\n";
        if (index_not_zero != 0)
            *det_quotient *= -1;
        *det_quotient *= placeholder(0, 0);
        first_row = first_row / placeholder(0, 0);
        // std::cout << "after divide:" << placeholder;

        for (int i = 1; i < height; i++)
        {
            T alpha = placeholder(i, 0);
            placeholder.slice({i, 0}, {1, width}) = placeholder.slice({i, 0}, {1, width}) - first_row * alpha;
        }
        // std::cout << "after normalzie: " << placeholder;
        if (height == 1 || width == 1)
        {
            return placeholder;
        }
        View<T> recursive_GE_submat = placeholder.slice({1, 1}, {height - 1, width - 1});
        auto gaussian_eliminated_submat = extendedGaussianElimination(recursive_GE_submat, det_quotient);
        recursive_GE_submat = gaussian_eliminated_submat;
        return placeholder;
    }
    if (height == 1 || width == 1)
    {
        return placeholder;
    }
    View<T> recursive_GE_submat = placeholder.slice({0, 1}, {height, width - 1});

    auto gaussian_eliminated_submat = extendedGaussianElimination(recursive_GE_submat, det_quotient);
    recursive_GE_submat = gaussian_eliminated_submat;
    return placeholder;
}
template <typename T>
T dot(Vector<T> v1, Vector<T> v2)
{
    check((v1.length() == v2.length()), "cannot dot product 2 different size vectors");
    Matrix<T> vec(v1);
    // std::cerr << "shape" << v.get_shape()[0] << v.get_shape()[1];
    Matrix<T> nice(v2.transpose() * vec); // why is v* trans double?
    return nice(0, 0);
}
#endif