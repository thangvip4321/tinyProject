#include "dbg.h"
#include "matrix_utils.hpp"
// ================================================================================================================================
// ====================================================     MATRIX_BASE     =======================================================
// ================================================================================================================================
// ================================================================================================================================
// ================================================================================================================================
// ================================================================================================================================
template <typename T>
Matrix<T> Matrix_base<T>::operator+(Matrix_base<T> &adder)
{
    check((adder.shape[0] == this->shape[0] && adder.shape[1] == this->shape[1]), "mismatch dimension");
    Matrix<T> new_mat = Matrix<T>(this->shape[0], this->shape[1]);
    for (int i = 0; i < this->shape[0]; i++)
    {
        for (int j = 0; j < this->shape[1]; j++)
        {
            new_mat(i, j) = this->operator()(i, j) + adder(i, j);
        }
    }
    return new_mat;
}

template <typename T>
Matrix<T> Matrix_base<T>::operator*(T scalar)
{
    Matrix<T> new_mat = Matrix<T>(*this);
    for (int i = 0; i < this->shape[0]; i++)
    {
        for (int j = 0; j < this->shape[1]; j++)
        {
            new_mat(i, j) *= scalar;
        }
    }
    return new_mat;
}
template <typename T>
Matrix<T> Matrix_base<T>::operator*(Matrix_base<T> &mat)
{
    check((mat.shape[0] == this->shape[1]), "mismatch dimension (%d,%d) and (%d,%d)", this->shape[0], this->shape[1], mat.shape[0], mat.shape[1]);
    Matrix<T> new_mat = Matrix<T>(this->shape[0], mat.shape[1]);

    for (int i = 0; i < this->shape[0]; i++)

    {
        for (int j = 0; j < mat.shape[1]; j++)
        {
            T sum = this->operator()(i, 0) * mat(0, j);
            for (int k = 1; k < this->shape[1]; k++)
            {
                sum += this->operator()(i, k) * mat(k, j);
            }
            new_mat(i, j) = sum;
        }
    }

    return new_mat;
}
template <typename T>
Matrix<T> Matrix_base<T>::operator/(T scalar)
{
    Matrix<T> new_mat = Matrix<T>(*this);
    for (int i = 0; i < this->shape[0]; i++)
    {
        for (int j = 0; j < this->shape[1]; j++)
        {
            new_mat(i, j) /= scalar;
        }
    }
    return new_mat;
}
template <typename T>
bool Matrix_base<T>::operator==(Matrix_base<T> &mat)
{
    if (this->get_shape()[0] != mat.get_shape()[0] || this->get_shape()[1] != mat.get_shape()[1])
        return false;
    for (int i = 0; i < this->shape[0]; i++)
    {
        for (int j = 0; j < this->shape[1]; j++)
        {
            if (this->operator()(i, j) != mat(i, j))
            {
                std::cout << "gotcha!: "
                          << this->operator()(i, j) << " and " << mat(i, j) << "\n";
                return false;
            }
        }
    }
    return true;
}

template <typename T>
std::ostream &operator<<(std::ostream &s, Matrix_base<T> &mat)
{
    s << "[ ";
    for (int i = 0; i < mat.shape[0]; i++)
    {
        for (int j = 0; j < mat.shape[1]; j++)
        {
            s << mat(i, j) << " ";
        }
        if (i == mat.shape[0] - 1)
            s << "]";
        s << "\n";
    }
    return s;
};

template <typename T>
T Matrix_base<T>::determinant()
{
    int width = this->get_shape()[1];
    int height = this->get_shape()[0];
    if (width != height)
    {
        throw std::runtime_error("cannot calculate determinant of non-square matrix");
    }
    T det_quotient = 1;
    T det = 1;
    auto ge_mat = extendedGaussianElimination(*this, &det_quotient);
    for (int i = 0; i < det_quotient; i++)
    {
        det *= ge_mat(i, i);
    }
    det /= det_quotient;
    return det;
}

template <typename T>
Matrix<T> Matrix_base<T>::inverse()
{
    int width = this->get_shape()[1];
    int height = this->get_shape()[0];
    if (width != height)
    {
        throw std::runtime_error("cannot find inverse of non-square matrix");
    }
    Matrix<T> identity_mat = Identity_matrix<T>(width);
    Matrix<T> placeholder_mat = Matrix<T>(width, width * 2);
    placeholder_mat.slice({0, 0}, {width, width}) = *this;
    placeholder_mat.slice({0, width}, {width, width}) = identity_mat;
    auto ge_mat = gaussianElimination(placeholder_mat);
    backwardSubstitution(ge_mat);

    if (!(ge_mat.slice({0, 0}, {width, width}) == identity_mat))
    {
        throw std::runtime_error("cannot find inverse of this");
    }
    auto view = ge_mat.slice({0, width}, {width, width});
    auto inverse = Matrix<T>(view);
    return inverse;
}
template <typename T>
int Matrix_base<T>::rank()
{
    int width = this->get_shape()[1];
    int height = this->get_shape()[0];
    auto ge_mat = gaussianElimination(*this);
    int num_leading_ones = 0;
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            std::cout << this->operator()(i, j) << " ";
            if (ge_mat(i, j) == 1)
            {
                num_leading_ones += 1;
                break;
            } //found a leading one
        }
    }
    return num_leading_ones;
}
template <typename T>
Matrix<T> Matrix_base<T>::pseudoinverse()
{
    int width = this->get_shape()[1];
    int height = this->get_shape()[0];
    int rank_mat = this->rank();
    std::cerr << "rank: " << rank_mat << "\n";
    if (rank_mat < std::min(width, height))
    {
        throw std::runtime_error("matrix not full rank so we cannot do moore-penrose pseudoinverse");
    }
    auto new_mat = Matrix<T>(width, height); //matrix automatically call destructor?
    if (width > height)
    {
        Matrix<T> res = (this->operator*(this->transpose()));
        debug("x");
        new_mat = this->transpose() * res.inverse();
        debug("y");
    }
    else if (width < height)
    {
        Matrix<T> res = (this->transpose() * (*this));
        new_mat = res.inverse() * this->transpose();
    }
    else
    {
        new_mat = this->inverse();
    }
    return new_mat;
}

// ================================================================================================================================
// ====================================================     MATRIX PART     =======================================================
// ================================================================================================================================
// ================================================================================================================================
// ================================================================================================================================
// ================================================================================================================================

template <typename T>
Matrix<T> &Matrix<T>::operator=(Matrix_base<T> &mat)
{
    int *shape = this->get_shape();
    check((shape[0] == mat.get_shape()[0] && shape[1] == mat.get_shape()[1]), "mismatch dimension");
    for (int i = 0; i < shape[0]; i++)
    {
        for (int j = 0; j < shape[1]; j++)
        {
            this->operator()(i, j) = mat(i, j);
        }
    }
    return *this;
}
template <typename T>
Matrix<T> &Matrix<T>::operator=(Matrix<T> &mat)
{

    int *shape = this->get_shape();
    check((shape[0] == mat.get_shape()[0] && shape[1] == mat.get_shape()[1]), "mismatch dimension");

    for (int i = 0; i < shape[0]; i++)
    {
        for (int j = 0; j < shape[1]; j++)
        {
            this->operator()(i, j) = mat(i, j);
        }
    }
    return *this;
}

template <typename T>
void infer_shape(Matrix_initializer<T, 2> data, int shape[2])
{
    auto start_row = data.begin();
    // std::cout << std::type_info<start_row>() << "\n";
    int height = data.size();
    int width = start_row[0].size();
    for (int i = 1; i < data.size(); i++)
    {
        if (width != start_row[i].size())
        {
            throw std::out_of_range("bạn ơi xin vui lòng gõ đúng số mảng`");
        }
    }
    shape[0] = height;
    shape[1] = width;
}
template <typename T>
void Matrix<T>::_insert_flat(Matrix_initializer<T, 2> data)
{
    auto start_mat = data.begin();
    for (int i = 0; i < this->shape[0]; i++)
    {
        auto ath_row = start_mat[i];
        auto start_row = ath_row.begin();
        for (int j = 0; j < this->shape[1]; j++)
        {
            this->data[i * this->shape[1] + j] = start_row[j];
        }
    }
}
template <typename T>
Matrix<T>::Matrix(int h, int w)
{
    this->shape[0] = h;
    this->shape[1] = w;
    this->data = (T *)calloc(this->shape[1] * this->shape[0], sizeof(T));
};
template <typename T>
Matrix<T>::Matrix(Matrix_initializer<T, 2> data)
{
    infer_shape<T>(data, this->shape);
    int size = this->shape[0] * this->shape[1];
    this->data = (T *)calloc(size, sizeof(T));
    _insert_flat(data);
};
template <typename T>
Matrix<T>::Matrix(Matrix_base<T> &mat)
{

    this->shape[0] = mat.get_shape()[0];
    this->shape[1] = mat.get_shape()[1];
    this->data = (T *)calloc(this->shape[0] * this->shape[1], sizeof(T));
    for (int i = 0; i < this->shape[0]; i++)
    {
        for (int j = 0; j < this->shape[1]; j++)
        {
            this->operator()(i, j) = mat(i, j);
        }
    }
};
template <typename T>
Matrix<T>::Matrix(Matrix<T> &mat)
{

    this->shape[0] = mat.get_shape()[0];
    this->shape[1] = mat.get_shape()[1];
    this->data = (T *)calloc(this->shape[0] * this->shape[1], sizeof(T));
    for (int i = 0; i < this->shape[0]; i++)
    {
        for (int j = 0; j < this->shape[1]; j++)
        {
            this->operator()(i, j) = mat(i, j);
        }
    }
};
template <typename T>
View<T> Matrix<T>::slice(std::initializer_list<int> start_pos, std::initializer_list<int> size)
{
    if (start_pos.size() != size.size() || start_pos.size() != 2)
    {
        throw std::logic_error("we only support 2d matrix for now");
    }
    int dim = start_pos.size();
    int start_view[2];
    int size_view[2];
    int orig_mat_size_view[2];
    const int *pos_start = start_pos.begin();
    const int *size_start = size.begin();
    for (int i = 0; i < dim; i++)
    {
        if (pos_start[i] < 0 || pos_start[i] + size_start[i] > this->shape[i] || size_start[i] > this->shape[i])
        {
            throw std::out_of_range("invalid size for setting a view");
        }
        start_view[i] = pos_start[i];
        size_view[i] = size_start[i];
        orig_mat_size_view[i] = this->shape[i];
    }
    // std::cout << "view is " << start_view[0] << start_view[1] << "\n";
    // std::cout << "size is " << size_view[0] << size_view[1] << "\n";

    View<T> new_view = View<T>(start_view, size_view, orig_mat_size_view, this->data);
    return new_view;
};
template <typename T>
Matrix<T> Matrix_base<T>::transpose()
{
    int new_height = this->shape[1];
    int new_width = this->shape[0];
    auto new_mat = Matrix<T>(new_height, new_width);
    for (int i = 0; i < this->shape[0]; i++)
    {
        for (int j = 0; j < this->shape[1]; j++)
        {
            new_mat(j, i) = this->operator()(i, j);
        }
    }
    return new_mat;
}

// ================================================================================================================================
// ================================================================================================================================
// ======================================================== VIEW PART==============================================================
// ================================================================================================================================
// ================================================================================================================================
// ================================================================================================================================
template <typename T>
View<T>::View(int start_pos[2], int size[2], int orig_mat_size[2], T *data)
{
    for (int i = 0; i < 2; i++)
    {
        this->start_pos[i] = start_pos[i];
        this->shape[i] = size[i];
        this->orig_mat_size[i] = orig_mat_size[i];
    }
    this->data = data;
}

template <typename T>
View<T> &View<T>::operator=(Matrix_base<T> &mat)
{
    int *shape = this->get_shape();
    check((shape[0] == mat.get_shape()[0] && shape[1] == mat.get_shape()[1]), "mismatch dimension");

    for (int i = 0; i < shape[0]; i++)
    {
        for (int j = 0; j < shape[1]; j++)
        {
            this->operator()(i, j) = mat(i, j);
        }
    }
    return *this;
}
template <typename T>
View<T> &View<T>::operator=(View<T> mat)
{
    int *shape = this->get_shape();
    check((shape[0] == mat.get_shape()[0] && shape[1] == mat.get_shape()[1]), "mismatch dimension");
    for (int i = 0; i < shape[0]; i++)
    {
        for (int j = 0; j < shape[1]; j++)
        {
            this->operator()(i, j) = mat(i, j);
        }
    }
    return *this;
}

template <typename T>
T &View<T>::operator()(int i, int j)
{

    if (j >= this->shape[1] || i >= this->shape[0] || j < 0 || i < 0)
    {
        throw std::out_of_range("out of range for this view");
    }
    return this->data[(i + this->start_pos[0]) * this->orig_mat_size[1] + (this->start_pos[1] + j)];
}
