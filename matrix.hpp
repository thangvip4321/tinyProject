#ifndef _MATRIX_H
#define _MATRIX_H
#include <iostream>
#include <cstddef>
#include "dbg.h"
template <typename T>
class Matrix;
template <typename T>
class View;

//=====code for matrix instantiation=================================
template <typename T, size_t N>
struct Matrix_init
{
    using type = std::initializer_list<typename Matrix_init<T, N - 1>::type>;
};

template <typename T>
struct Matrix_init<T, 1>
{
    using type = std::initializer_list<T>;
};

template <typename T, size_t N>
using Matrix_initializer = typename Matrix_init<T, N>::type; // what tf is typename here?
//
//
//
//
//
//
template <typename T>
class Matrix_base
{
protected:
    // height and width of the matrix
    int shape[2];
    T *data;

public:
    virtual T &operator()(int i, int j)
    {
        check((i < this->shape[0] && j < this->shape[1]), "index error,try to read value outside of matrix");
        return data[i * this->shape[1] + j];
    }; // for indexing inside the matrix
    Matrix<T> operator+(Matrix_base<T> &mat);
    Matrix<T> operator+(Matrix_base<T> &&mat) { return operator+(mat); };

    Matrix<T> operator*(Matrix_base<T> &mat);
    Matrix<T> operator*(Matrix_base<T> &&mat) { return operator*(mat); };

    Matrix<T> operator*(T scalar);
    Matrix<T> operator/(T scalar);
    Matrix<T> operator-(Matrix_base<T> &mat) { return *this + mat * (-1); };
    Matrix<T> operator-(Matrix_base<T> &&mat) { return operator-(mat); };

    Matrix<T> transpose();
    T determinant();
    Matrix<T> inverse();
    int rank();
    Matrix<T> pseudoinverse();
    bool operator==(Matrix_base<T> &mat);
    bool operator==(Matrix_base<T> &&mat) { return operator==(mat); };

    template <typename U>
    friend std::ostream &operator<<(std::ostream &s, Matrix_base<U> &mat);
    int *get_shape() { return this->shape; };
};

/* View
This class is for accessing the matrix easily. It does not contain any data by itself and thus should only be called by the original matrix.
Destruction of the original matrix will make the view unusable.
*/
template <typename T>
class View : public Matrix_base<T>
{
private:
    int start_pos[2];
    int orig_mat_size[2];
    View(int start_pos[2], int size[2], int orig_mat_size[2], T *data);

public:
    View &operator=(Matrix_base<T> &mat); //done
    View &operator=(Matrix_base<T> &&mat)
    {
        return operator=(mat);
    }
    View &operator=(View<T> view);
    T &operator()(int i, int j); //dones
    template <typename U>
    friend class Matrix;
};

// ====================================== CLASS FOR MATRIX========================================
// ===============================================================================================
template <typename T>
class Matrix : public Matrix_base<T>
{
public:
    Matrix(int h, int w);
    Matrix(Matrix_initializer<T, 2> data);
    Matrix(Matrix_base<T> &mat);
    Matrix(Matrix_base<T> &&mat) : Matrix(mat){};
    Matrix(Matrix<T> &mat);
    Matrix(Matrix<T> &&mat) : Matrix(mat){};
    ~Matrix()
    {
        // std::cerr << "calling matrix destructor\n";
        free(this->data);
    };
    Matrix<T> &operator=(Matrix_base<T> &mat);
    Matrix<T> &operator=(Matrix_base<T> &&mat) { return operator=(mat); };
    Matrix<T> &operator=(Matrix<T> &mat);
    // Matrix<T> &operator=(Matrix<T> &&mat) { return operator=(mat); };

    View<T> slice(std::initializer_list<int> start_pos, std::initializer_list<int> size);

private:
    void _insert_flat(Matrix_initializer<T, 2> data);
};

#include "matrix.tpp"
#endif
