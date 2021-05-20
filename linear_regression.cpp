#include <iostream>
#include <fstream>
#include "linear_system.hpp"
#include "vector.hpp"
#include "matrix.hpp"
using namespace std;

int main()
{
    ifstream fin;
    fin.open("machine.data", ifstream::in);

    // Just read 170 instances
    Matrix<double> A(170, 7);
    Vector<double> b(170);
    char cNum[10];

    if (!fin.is_open())
    {
        cout << "Error to open file";
        exit(-1);
    }

    // Read the training instances
    for (int i = 0; i < 170; i++)
    {
        fin.ignore(100, ',')
            .ignore(100, ',');
        for (int j = 0; j < 7; j++)
        {
            if (j != 6)
            {
                fin.getline(cNum, 20, ',');
                A(i, j) = atof(cNum);
            }
            else
                A(i, j) = 1;
        }
        fin.getline(cNum, 20, ',');
        b[i] = atof(cNum);
        fin.ignore(20, '\n');
    }

    //  Read the testing instances
    Matrix<double> A_test(39, 7);
    Vector<double> b_test(39);
    for (int i = 0; i < 39; i++)
    {
        fin.ignore(100, ',').ignore(100, ',');
        for (int j = 0; j < 7; j++)
        {
            if (j != 6)
            {
                fin.getline(cNum, 20, ',');
                A_test(i, j) = atof(cNum);
            }
            else
                A_test(i, j) = 1;
        }
        fin.getline(cNum, 20, ',');
        b_test[i] = atof(cNum);
        fin.ignore(20, '\n');
    }
    fin.close();

    Vector<double> x = A.pseudoinverse() * b;
    Vector<double> y_hat = A_test * x;
    // Calculate RMSE and deviation
    double sum = 0;
    double sum_2 = 0;
    for (int i = 0; i < 39; i++)
    {
        // cout << b_test(i) << " " << c(i) << endl;
        sum += (y_hat[i] - b_test[i]) * (y_hat[i] - b_test[i]);
        sum_2 += fabs(y_hat[i] - b_test[i]) / b_test[i] * 100;
    }
    // cout << sum;
    double RMSE = sqrt(sum / 39);
    double average_deviation = sum_2 / 39;
    cout << RMSE << " " << average_deviation << " " << endl;
}
// return 0;
