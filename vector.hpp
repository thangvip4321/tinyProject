#ifndef _VECTOR_H_
#define _VECTOR_H_
#include "matrix.hpp"
#include "dbg.h"
#include "math.h"
template <typename T>
class Vector : public Matrix<T>
{
public:
    Vector(int length) : Matrix<T>(length, 1){};
    Vector(std::initializer_list<T> vec_init) : Matrix<T>({vec_init})
    {
        int temp = this->shape[0];
        this->shape[0] = this->shape[1];
        this->shape[1] = temp;
    };
    Vector(Matrix_base<T> &mat) : Matrix<T>(mat)
    {

        if (mat.get_shape()[1] != 1)
        {
            throw std::invalid_argument("bạn ơi vector chỉ có 1 chiều thôi");
        }
    };
    Vector(Matrix_base<T> &&mat) : Vector(mat){};
    // Vector(Vector<T> &vec)
    Vector<T> &operator=(Matrix_base<T> &mat)
    {
        debug("hmm");
        check((mat.get_shape()[1] == 1 && mat.get_shape()[0] == this->length()), "mismatch dimension mat: (%d,%d),vector(%d,%d)", mat.get_shape()[0], mat.get_shape()[1], this->shape[0], this->shape[1]);

        for (int i = 0; i < mat.get_shape()[0]; i++)
        {
            this->operator[](i) = mat(i, 0);
        }
        return *this;
    };
    Vector<T> &operator=(Matrix_base<T> &&mat)
    {
        return this->operator=(mat);
    }

    int length()
    {
        return this->shape[0];
    }
    T &operator[](int i)
    {
        check((i < this->length()), "the index lies outside of array range");
        return this->data[i];
    }
    T norm()
    {
        T sum = 0;
        for (int i = 0; i < length(); i++)
        {
            sum += pow(this->operator[](i), 2);
        }
        return sqrt(sum);
    }

    Vector<T> operator*(T v)
    {
        return Matrix<T>(*this) * v;
    }
};
#endif