#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <chrono>
using namespace std::chrono;

double frand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void FillUpper(std::vector<std::vector<double>> &v)
{
    for (int i = 0; i < v.size(); i++)
    {
        for (int j = i; j < v[i].size(); j++)
        {
            v[i][j] = frand(-4, 4); // [-4, 4]
        }
    }
}

void FillLower(std::vector<std::vector<double>> &v)
{
    for (int i = 0; i < v.size(); i++)
    {
        for (int j = 0; j < i; j++)
        {
            v[i][j] = v[j][i];
        }
    }
}

void FillDiagonal(std::vector<std::vector<double>> &v)
{
    for (int i = 0; i < v.size(); i++)
    {
        double sum = 0;
        for (int j = 0; j < v[i].size(); j++)
        {
            if (i != j)
            {
                sum += fabs(v[i][j]);
            }
        }
        v[i][i] = sum + 1;
    }
}

void FillMatrix(std::vector<std::vector<double>> &v)
{
    FillUpper(v);
    FillLower(v);
    FillDiagonal(v);
}

void PrintMatrix(std::vector<std::vector<double>> &v)
{
    for (int i = 0; i < v.size(); i++)
    {
        for (int j = 0; j < v[i].size(); j++)
        {
            std::cout << v[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// print matrix to ofstream
void PrintMatrix(std::vector<std::vector<double>> &v, std::ofstream &out)
{
    for (int i = 0; i < v.size(); i++)
    {
        for (int j = 0; j < v[i].size(); j++)
        {
            out << v[i][j] << " ";
        }
        out << std::endl;
    }
    out << std::endl;
}

void FillVector(std::vector<double> &v)
{
    for (int i = 0; i < v.size(); i++)
    {
        v[i] = frand(-4, 4); // [-4, 4]
    }
}

void PrintVector(const std::vector<double> &v)
{
    for (int i = 0; i < v.size(); i++)
    {
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
}

void MatrixMultiply(const std::vector<std::vector<double>> &a,
                    const std::vector<double> &b,
                    std::vector<double> &c)
{
    for (int i = 0; i < a.size(); i++)
    {
        for (int j = 0; j < b.size(); j++)
        {
            double sum = 0;
            for (int k = 0; k < a[0].size(); k++)
            {
                sum += a[i][k] * b[k];
            }
            c[i] = sum;
        }
    }
}

// vector subtraction and retirn the result
std::vector<double> VectorSubtract(const std::vector<double> &a, const std::vector<double> &b)
{
    std::vector<double> c(a.size());
    for (int i = 0; i < a.size(); i++)
    {
        c[i] = a[i] - b[i];
    }
    return c;
}

void MatrixMultiply(const std::vector<std::vector<double>> &a,
                    const std::vector<std::vector<double>> &b,
                    std::vector<std::vector<double>> &c)
{
    for (int i = 0; i < a.size(); i++)
    {
        for (int j = 0; j < b[0].size(); j++)
        {
            double sum = 0;
            for (int k = 0; k < a[0].size(); k++)
            {
                sum += a[i][k] * b[k][j];
            }
            c[i][j] = sum;
        }
    }
}

// find inverse matrix using Gauss-Jordan elimination
std::vector<std::vector<double>> InverseMatrixGaussianElimination(const std::vector<std::vector<double>> &a)
{
    std::vector<std::vector<double>> b(a.size(), std::vector<double>(a.size()));
    std::vector<std::vector<double>> a1 = a;

    for (int i = 0; i < a.size(); ++i)
    {
        b[i][i] = 1;
    }

    for (int i = 0; i < a.size() - 1; ++i)
    {
        double fact = a1[i][i];

        for (int j = i + 1; j < a.size(); ++j)
        {
            double coff = a1[j][i] / fact;

            for (int k = 0; k < a.size(); ++k)
            {
                b[j][k] -= coff * b[i][k];
                a1[j][k] -= coff * a1[i][k];
                // if (a1[j][k] < 0.00001) a1[j][k] = 0;
            }
        }
    }

    for (int j = a.size() - 1; j > 0; --j)
    {
        double fact = a1[j][j];

        for (int i = j - 1; i >= 0; --i)
        {
            // std::cout << -a[i][j] << ' ' << fact << std::endl;
            double coff = a1[i][j] / fact;

            for (int k = 0; k < a.size(); ++k)
            {
                b[i][k] -= coff * b[j][k];
                a1[i][k] -= coff * a1[j][k];
                // if (a1[i][k] < 0.00001) a1[i][k] = 0;
            }
        }
    }

    for (int i = 0; i < a.size(); ++i)
    {
        for (int j = 0; j < a.size(); ++j)
        {
            b[i][j] /= a1[i][i];
        }
    }

    return b;
}

// get norm of matrix
double FindNorm(const std::vector<std::vector<double>> &a)
{
    double max = -1;
    double sum = 0;

    for (int i = 0; i < a.size(); i++)
    {
        sum = 0;
        for (int j = 0; j < a[i].size(); j++)
        {
            if (a[i][j] >= 0)
                sum += a[i][j];
            else
                sum -= a[i][j];
        }

        if (sum > max)
        {
            max = sum;
        }
    }

    return max;
}

// get condition number of matrix
double ConditionNumber(const std::vector<std::vector<double>> &a)
{
    std::vector<std::vector<double>> inverse_matrix = InverseMatrixGaussianElimination(a);
    return FindNorm(a) * FindNorm(inverse_matrix);
}

// solve Ax = b using Gaussian elimination and return x
std::vector<double> SolveMatrixGaussianElimination(const std::vector<std::vector<double>> &a,
                                                   const std::vector<double> &b)
{
    std::vector<std::vector<double>> a1 = a;
    std::vector<double> b1 = b;
    std::vector<double> x(a.size());

    int n = a.size();
    int k = 0, index;
    const double eps = 0.00001;
    double max = 0;

    while (k < n)
    {
        // find max a[i][k]
        max = fabs(a1[k][k]);
        index = k;

        for (int i = k + 1; i < n; ++i)
        {
            if (fabs(a1[i][k]) > max)
            {
                max = fabs(a1[i][k]);
                index = i;
            }
        }
        // ???????????????????????? ??????????
        if (max < eps)
        {
            std::cout << "Matrix is singular" << std::endl; // no not-zero diagonal elements
            return x;
        }

        for (int j = 0; j < n; ++j)
        {
            double temp = a1[k][j];
            a1[k][j] = a1[index][j];
            a1[index][j] = temp;
        }

        double temp = b1[k];
        b1[k] = b1[index];
        b1[index] = temp;

        // Normalize equation
        for (int i = k; i < n; ++i)
        {
            double temp = a1[i][k];
            // std::cout << temp << ' ' << abs(temp) << std::endl;
            if (fabs(temp) < eps) // skip zero elements
                continue;
            for (int j = k; j < n; j++)
                a1[i][j] = a1[i][j] / temp;
            b1[i] = b1[i] / temp;
            if (i == k)
                continue; // don't self-minus
            for (int j = k; j < n; j++)
                a1[i][j] = a1[i][j] - a1[k][j];
            b1[i] = b1[i] - b1[k];
        }
        ++k;
    }

    // reverse
    for (k = n - 1; k >= 0; k--)
    {
        x[k] = b1[k];
        for (int i = 0; i < k; i++)
            b1[i] = b1[i] - a1[i][k] * x[k];
    }

    return x;
}

// get LUP decomposition of matrix
void LUP(std::vector<std::vector<double>> a, std::vector<std::vector<double>> &l,
         std::vector<std::vector<double>> &u, std::vector<double> &p)
{
    u = a;
    int n = a.size();

    for (int i = 0; i < n; ++i)
        p[i] = i;

    for (int i = 0; i < n; ++i)
        l[i][i] = 1;

    for (int i = 0; i < n; i++)
        for (int j = i; j < n; j++)
            l[j][i] = u[j][i] / u[i][i];

    for (int k = 1; k < n; k++)
    {
        double max = -5;
        int index = 0;

        for (int i = k; i < n; ++i)
        {
            if (fabs(u[i][k]) > max)
            {
                max = fabs(u[i][k]);
                index = i;
            }
        }

        // ???????????????????????? ??????????
        for (int j = 0; j < n; ++j)
        {
            double temp = u[k][j];
            u[k][j] = u[index][j];
            u[index][j] = temp;
        }
        double temp = p[k];
        p[k] = p[index];
        p[index] = temp;

        for (int j = 0; j < k; ++j)
        {
            double temp = l[k][j];
            l[k][j] = l[index][j];
            l[index][j] = temp;
        }

        for (int i = k - 1; i < n; i++)
            for (int j = i; j < n; j++)
                l[j][i] = u[j][i] / u[i][i];

        for (int i = k; i < n; i++)
            for (int j = k - 1; j < n; j++)
                u[i][j] = u[i][j] - l[i][k - 1] * u[k - 1][j];
    }
}

// solve Ax = b using LUP decomposition
std::vector<double> LUP_solve(const std::vector<std::vector<double>> &a, const std::vector<double> &b)
{
    std::vector<std::vector<double>> l(a.size(), std::vector<double>(a.size()));
    std::vector<std::vector<double>> u(a.size(), std::vector<double>(a.size()));
    std::vector<double> p(a.size());

    LUP(a, l, u, p);

    std::vector<double> x(a.size());
    std::vector<double> y(a.size());

    for (int i = 0; i < a.size(); i++)
        y[i] = b[p[i]];

    for (int i = 0; i < a.size(); i++)
    {
        double sum = 0;
        for (int j = 0; j < i; j++)
            sum += l[i][j] * y[j];
        y[i] = y[i] - sum;
    }

    for (int i = a.size() - 1; i >= 0; i--)
    {
        double sum = 0;
        for (int j = i + 1; j < a.size(); j++)
            sum += u[i][j] * x[j];
        x[i] = (y[i] - sum) / u[i][i];
    }

    return x;
}

// get transpose of matrix
std::vector<std::vector<double>> Transpose(const std::vector<std::vector<double>> &a)
{
    int n = a.size();
    std::vector<std::vector<double>> b(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            b[i][j] = a[j][i];
    return b;
}

// slove Ax = b by sqare root method and return x
std::vector<double> Sqrt_solve(const std::vector<std::vector<double>> &a, const std::vector<double> &b)
{
    int n = a.size();
    // std::vector<std::vector<double>> a1 = a;
    // std::vector<double> b1 = b;
    std::vector<std::vector<double>> s(n, std::vector<double>(n));

    // fill first element s[1][1]
    s[0][0] = sqrt(a[0][0]);

    // fill first column s[1][j]
    for (int j = 1; j < n; ++j)
    {
        s[0][j] = a[0][j] / s[0][0];
    }

    // fill main diagonal s[i][i]
    // PrintMatrix(s);
    for (int i = 1; i < n; ++i)
    {
        double cof = 0;
        for (int k = 0; k < i; ++k)
        {
            cof += pow((double)s[k][i], 2.0);
        }
        // std::cout << cof << std::endl;
        s[i][i] = sqrt(a[i][i] - cof);
    }

    // PrintMatrix(s);
    //  fill upper triangle, lower triangle = 0
    for (int i = 1; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            double cof = 0;
            for (int k = 0; k < i; ++k)
            {
                cof += s[k][i] * s[k][j];
            }
            s[i][j] = (a[i][j] - cof) / s[i][i];
        }
    }

    // PrintMatrix(s);
    std::vector<double> y(n);
    // fill first y[1][1]
    // y[0] = b[0] / s[0][0];

    // fill y
    std::vector<std::vector<double>> s1 = Transpose(s);
    // PrintMatrix(s1);
    for (int i = 0; i < n; ++i)
    {
        double cof = 0;
        for (int k = 0; k < i; ++k)
        {
            cof += s1[k][i] * y[k];
        }
        y[i] = (b[i] - cof) / s1[i][i];
    }
    // PrintVector(y);

    std::vector<double> x(n);
    // fill last x
    // x[n - 1] = y[n - 1] / s[n - 1][n - 1];

    // fill x
    for (int i = n - 1; i >= 0; --i)
    {
        double cof = 0;
        for (int k = n - 1; k > i; --k)
        {
            cof += s[i][k] * x[k];
        }
        x[i] = (y[i] - cof) / s[i][i];
    }
    // PrintVector(x);
    // std::vector<std::vector<double>> wtf(n, std::vector<double>(n));
    // MatrixMultiply(s, s1, wtf);
    // PrintMatrix(wtf);
    return x;
}

// get LDLT decomposition of a matrix
bool LDLT(const std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &l, std::vector<double> &d)
{
    int n = a.size();
    double sum = 0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = i; j < n; ++j)
        {
            sum = a[j][i];
            for (int k = 0; k < i; ++k)
            {
                sum -= l[i][k] * d[k] * l[j][k];
            }
            if (i == j) // diagonal element
            {
                if (sum <= 0)
                {
                    std::cout << "Error: matrix is not positive definite" << std::endl;
                    return false;
                }
                d[i] = sum;
                l[i][i] = 1;
            }
            else
            {
                l[j][i] = sum / d[i]; // not diagonal element
            }
        }
    }
    return true;
}

// multiply matrix by vector
std::vector<double> MatrixMultiply(const std::vector<std::vector<double>> &a, const std::vector<double> &b)
{
    int n = a.size();
    std::vector<double> c(n);
    for (int i = 0; i < n; ++i)
    {
        double sum = 0;
        for (int j = 0; j < n; ++j)
        {
            sum += a[i][j] * b[j];
        }
        c[i] = sum;
    }
    return c;
}

// multiply matrix by matrix
std::vector<std::vector<double>> MatrixMultiply(const std::vector<std::vector<double>> &a, const std::vector<std::vector<double>> &b)
{
    int n = a.size();
    std::vector<std::vector<double>> c(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            double sum = 0;
            for (int k = 0; k < n; ++k)
            {
                sum += a[i][k] * b[k][j];
            }
            c[i][j] = sum;
        }
    }
    return c;
}

// find norm of two vectors
double Norm(const std::vector<double> &a, const std::vector<double> &b)
{
    int n = a.size();
    double sum = 0;
    for (int i = 0; i < n; ++i)
    {
        sum += pow(a[i] - b[i], 2.0);
    }
    return sqrt(sum);
}

// solve matrix by relaxation method with parameter alpha
std::vector<double> Relaxation_solve(const std::vector<std::vector<double>> &a, const std::vector<double> &b, double alpha, int &count)
{
    int n = a.size();
    std::vector<double> x(n);
    std::vector<double> x_old(n);
    std::vector<double> x_new(n);
    for (int i = 0; i < n; ++i)
    {
        x[i] = b[i];
    }
    do
    {
        for (int i = 0; i < n; ++i)
        {
            x_old[i] = x[i];
        }
        x_new = x_old;
        for (int i = 0; i < n; ++i)
        {
            double sum = 0;
            for (int j = 0; j < n; ++j)
            {
                if (i != j)
                {
                    sum += a[i][j] * x_old[j];
                }
            }
            x_new[i] = (1 - alpha) * x_old[i] + alpha * (b[i] - sum) / a[i][i];
        }
        x = x_new;
        ++count;
    } while (Norm(x_new, x_old) > 1e-6);
    return x;
}

// find cube norm of a vector
double Norm(const std::vector<double> &a)
{
    int n = a.size();
    double max = 0;
    for (int i = 0; i < n; ++i)
    {
        if (fabs(a[i]) > max)
        {
            max = fabs(a[i]);
        }
    }
    return max;
}

int main()
{
    const int size = 3;
    std::vector<std::vector<double>> matrix(size, std::vector<double>(size));
    std::vector<double> vec(size);
    std::vector<double> mult(size);
    FillMatrix(matrix);
    std::cout << "\tFirst task:\n";
    std::cout << "Matrix A:" << std::endl;
    PrintMatrix(matrix);

    FillVector(vec);
    std::cout << "\tSecond task:\n";
    std::cout << "Vector y:" << std::endl;
    PrintVector(vec);
    MatrixMultiply(matrix, vec, mult);
    std::cout << "Vector b (multiplication result):\n";
    PrintVector(mult);

    std::vector<std::vector<double>> inverse_matrix = InverseMatrixGaussianElimination(matrix);
    std::cout << "\tThird task:\n";
    std::cout << "Inverse matrix:" << std::endl;
    PrintMatrix(inverse_matrix);

    double conditionNumber = ConditionNumber(matrix);
    std::cout << "Condition number: " << conditionNumber << std::endl
              << std::endl;

    std::vector<double> x = SolveMatrixGaussianElimination(matrix, mult);
    std::cout << "\tFourth task:\n";
    std::cout << "Vector x (solution to Ax = b):\n";
    PrintVector(x);

    std::cout << "\tFifth task:\n";
    std::vector<std::vector<double>> l(size, std::vector<double>(size)), u(size, std::vector<double>(size));
    std::vector<double> p(size);

    LUP(matrix, l, u, p);
    std::cout << "LUP decomposition:\n";
    std::cout << "Matrix L:" << std::endl;
    PrintMatrix(l);
    std::cout << "Matrix U:" << std::endl;
    PrintMatrix(u);
    std::cout << "Vector P:" << std::endl;
    PrintVector(p);
    std::vector<std::vector<double>> P(size, std::vector<double>(size));

    for (int i = 0; i < size; ++i)
        P[i][p[i]] = 1;

    std::cout << "Matrix P:" << std::endl;
    PrintMatrix(P);
    std::cout << "Solve LUP:\n";
    PrintVector(LUP_solve(matrix, mult));

    std::cout << "\tSixth task:\n";
    std::vector<double> b = Sqrt_solve(matrix, mult);
    std::cout << "Vector b (solution to Ax = b):\n";
    PrintVector(b);

    std::cout << "LDLT decomposition:\n";
    std::vector<std::vector<double>> l1(size, std::vector<double>(size));
    std::vector<double> d(size);
    LDLT(matrix, l1, d);
    std::cout << "Matrix L:" << std::endl;
    PrintMatrix(l1);
    std::cout << "Vector D (diagonal):" << std::endl;
    PrintVector(d);

    std::cout << "\tSeventh task:\n";
    std::cout << "Solution to Ax=b by relaxation method with param 1 - 8/40\n";
    int count = 0;
    std::vector<double> xx = Relaxation_solve(matrix, mult, 1 - 8 / 40, count);
    PrintVector(xx);

    std::cout << "\tEight task\n";
    double max = 0, sum = 0, min = 1000000000, time = 0, time1 = 0, time2 = 0, time3 = 0, time4 = 0, time5 = 0;
    double max_r = 0, min_r = 10000000, sum_r = 0;
    double max_gauss = 0, min_gauss = 10000000, sum_gauss = 0;
    double max_sqrt = 0, min_sqrt = 10000000, sum_sqrt = 0;
    double max_lup = 0, min_lup = 10000000, sum_lup = 0;
    double max_relax = 0, min_relax = 10000000, sum_relax = 0;

    std::ofstream out("matrix_max_condition_number.txt");
    std::vector<std::vector<double>> matrix_max_cond_num(size, std::vector<double>(size));

    for (int i = 0; i < 100; ++i)
    {
        FillMatrix(matrix);
        double cn = ConditionNumber(matrix);
        sum += cn;
        if (cn > max)
        {
            max = cn;
            out.close();
            out.open("matrix_max_condition_number.txt", std::ios::out | std::ios::in);
            matrix_max_cond_num = matrix;
            PrintMatrix(matrix, out);
        }
        if (cn < min)
        {
            min = cn;
        }
        auto start = high_resolution_clock::now();
        InverseMatrixGaussianElimination(matrix);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        time += duration.count();

        auto start1 = high_resolution_clock::now();
        std::vector<double> gauss_solution = SolveMatrixGaussianElimination(matrix, mult);
        auto stop1 = high_resolution_clock::now();
        auto duration1 = duration_cast<microseconds>(stop1 - start1);
        time1 += duration1.count();

        auto start2 = high_resolution_clock::now();
        LUP(matrix, l, u, p);
        auto stop2 = high_resolution_clock::now();
        auto duration2 = duration_cast<microseconds>(stop2 - start2);
        time2 += duration2.count();

        auto start3 = high_resolution_clock::now();
        std::vector<double> lup_solution = LUP_solve(matrix, mult);
        auto stop3 = high_resolution_clock::now();
        auto duration3 = duration_cast<microseconds>(stop3 - start3);
        time3 += duration3.count();

        auto start4 = high_resolution_clock::now();
        std::vector<double> sqrt_solve = Sqrt_solve(matrix, mult);
        auto stop4 = high_resolution_clock::now();
        auto duration4 = duration_cast<microseconds>(stop4 - start4);
        time4 += duration4.count();

        int count = 0;
        auto start5 = high_resolution_clock::now();
        xx = Relaxation_solve(matrix, mult, 1 - 8 / 40, count);
        auto stop5 = high_resolution_clock::now();
        auto duration5 = duration_cast<microseconds>(stop5 - start5);
        time5 += duration5.count();
        sum_r += count;

        if (count > max_r)
            max_r = count;
        if (count < min_r)
            min_r = count;


        double gauss_norm = Norm(VectorSubtract(vec, gauss_solution)); // find norm of x_gauss - y

        if (gauss_norm > max_gauss)
            max_gauss = gauss_norm; // find max norm of difference
        if (gauss_norm < min_gauss)
            min_gauss = gauss_norm; // find min norm of difference
        sum_gauss += gauss_norm;    // sum norm of difference

        double lup_norm = Norm(VectorSubtract(vec, lup_solution)); // find norm of x_lup - y

        if (lup_norm > max_lup)
            max_lup = lup_norm; // find max norm of difference
        if (lup_norm < min_lup)
            min_lup = lup_norm; // find min norm of difference
        sum_lup += lup_norm;    // sum norm of difference

        double sqrt_norm = Norm(VectorSubtract(vec, sqrt_solve)); // find norm of x_sqrt - y

        if (sqrt_norm > max_sqrt)
            max_sqrt = sqrt_norm; // find max norm of difference
        if (sqrt_norm < min_sqrt)
            min_sqrt = sqrt_norm; // find min norm of difference
        sum_sqrt += sqrt_norm;    // sum norm of difference

        double relax_norm = Norm(VectorSubtract(vec, xx)); // find norm of x_relax - y

        if (relax_norm > max_relax)
            max_relax = relax_norm; // find max norm of difference
        if (relax_norm < min_relax)
            min_relax = relax_norm; // find min norm of difference
        sum_relax += relax_norm;    // sum norm of difference
    }
    out.close();

    std::ofstream out1("task8.txt");
    out1 << "Average condition number: " << sum / 100 << std::endl;
    out1 << "Max condition number: " << max << std::endl;
    out1 << "Min condition number: " << min << std::endl;
    out1 << std::endl;

    out1 << "Average time for inverse matrix: " << time / 100 << " microseconds" << std::endl
         << std::endl;

    out1 << "Max norm of difference for Gauss method: " << max_gauss << std::endl;
    out1 << "Min norm of difference for Gauss method: " << min_gauss << std::endl;
    out1 << "Average norm of difference for Gauss method: " << sum_gauss / 100 << std::endl;
    out1 << "Max norm of difference for LUP method: " << max_lup << std::endl;
    out1 << "Min norm of difference for LUP method: " << min_lup << std::endl;
    out1 << "Average norm of difference for LUP method: " << sum_lup / 100 << std::endl;
    out1 << "Max norm of difference for Sqrt method: " << max_sqrt << std::endl;
    out1 << "Min norm of difference for Sqrt method: " << min_sqrt << std::endl;
    out1 << "Average norm of difference for Sqrt method: " << sum_sqrt / 100 << std::endl;
    out1 << "Max norm of difference for Relaxation method: " << max_relax << std::endl;
    out1 << "Min norm of difference for Relaxation method: " << min_relax << std::endl;
    out1 << "Average norm of difference for Relaxation method: " << sum_relax / 100 << std::endl;
    out1 << std::endl;

    out1 << "Average time for Gauss method: " << time1 / 100 << " microseconds" << std::endl
         << std::endl;

    out1 << "Average time for LUP method: " << time2 / 100 << " microseconds" << std::endl
         << std::endl;

    out1 << "Average time for LUP_solve method: " << time3 / 100 << " microseconds" << std::endl
         << std::endl;

    out1 << "Average time for Sqrt_solve method: " << time4 / 100 << " microseconds" << std::endl
         << std::endl;

    out1 << "Average time for Relaxation_solve method: " << time5 / 100 << " microseconds" << std::endl
         << std::endl;

    out1 << "Max iterations for Relaxation_solve method: " << max_r << std::endl;
    out1 << "Min iterations for Relaxation_solve method: " << min_r << std::endl;
    out1 << "Average iterations for Relaxation_solve method: " << sum_r / 100 << std::endl
         << std::endl;
    out1.close();

    std::cout << "\tTask 9" << std::endl;
    double N = 8;
    std::vector<std::vector<double>> a1 = {{pow(N, 2) + 15, N - 1, -1, -2}, {N - 1, -15 - pow(N, 2), -N + 4, -4}, {-1, -N + 4, pow(N, 2) + 8, -N}, {-2, -4, -N, pow(N, 2) + 10}};

    std::vector<std::vector<double>> a2_t = {{1, 1 + N, 2 + N, 3 + N, 4 + N, 5 + N, 6 + N, 7 + N},
                                             {100 * N, 1000 * N, 10000 * N, 100000 * N, -1000 * N, -10000 * N, -100000 * N, 1},
                                             {N, -1 + N, -2 + N, -3 + N, -4 + N, -5 + N, -6 + N, -7 + N},
                                             {N - 1000, 10 * N - 1000, 100 * N - 1000, 1000 * N - 1000, 10000 * N - 1000, -N, -N + 1, -N + 2},
                                             {-2 * N, 0, -1, -2, -3, -4, -5, -6},
                                             {N - 2019, -N + 2020, N - 2021, -N + 2022, N - 2023, -N + 2024, N - 2025, -N + 2026},
                                             {2 * N - 2000, 4 * N - 2005, 8 * N - 2010, 16 * N - 2015, 32 * N - 2020, 2019 * N, -2020 * N, 2021 * N},
                                             {1020 - 2 * N, -2924 + 896 * N, 1212 + 9808 * N, -2736 + 98918 * N, 1014 - 11068 * N, -1523 - 8078 * N, 2625 - 102119 * N, -1327 + 1924 * N}};

    std::vector<std::vector<double>> a2 = MatrixMultiply(Transpose(a2_t), a2_t);

    std::cout << "For matrix A1:" << std::endl;
    PrintMatrix(a1);
    int size_a1 = a1.size();
    std::vector<double> y_a1(size_a1);
    FillVector(y_a1);
    std::cout << "Vector y:" << std::endl;
    PrintVector(y_a1);

    std::vector<double> b_a1 = MatrixMultiply(a1, y_a1);
    std::cout << "Vector b:" << std::endl;
    PrintVector(b_a1);

    std::cout << "Condition number for A1: " << ConditionNumber(a1) << std::endl;

    std::cout << "Gauss for matrix A1:" << std::endl;
    PrintVector(SolveMatrixGaussianElimination(a1, b_a1));

    std::cout << "LUP for matrix A1:" << std::endl;
    std::vector<std::vector<double>> l_a1(size_a1, std::vector<double>(size_a1)), u_a1(size_a1, std::vector<double>(size_a1));
    std::vector<double> p_a1(size_a1);

    LUP(a1, l_a1, u_a1, p_a1);
    std::cout << "LUP decomposition:\n";
    std::cout << "Matrix L:" << std::endl;
    PrintMatrix(l_a1);
    std::cout << "Matrix U:" << std::endl;
    PrintMatrix(u_a1);
    std::cout << "Vector P:" << std::endl;
    PrintVector(p_a1);
    std::cout << "Vector b (solution using LUP):" << std::endl;
    PrintVector(LUP_solve(a1, b_a1));

    std::cout << "Solution using square root method:" << std::endl;
    PrintVector(Sqrt_solve(a1, b_a1));

    std::cout << "LDLT decomposition:\n";
    std::vector<std::vector<double>> l1_a1(size_a1, std::vector<double>(size_a1));
    std::vector<double> d_a1(size_a1);
    if (LDLT(a1, l1_a1, d_a1))
    {
        std::cout << "Matrix L:" << std::endl;
        PrintMatrix(l1_a1);
        std::cout << "Vector D (diagonal):" << std::endl;
        PrintVector(d_a1);
    }
    else
    {
        std::cout << "Matrix is not positive definite" << std::endl
                  << std::endl;
    }

    std::cout << "Solution using Relaxation method:" << std::endl;
    PrintVector(Relaxation_solve(a1, b_a1, 1 - (8 / 40), count));

    std::cout << "For matrix A2:" << std::endl;
    PrintMatrix(a2);
    int size_a2 = a2.size();
    std::vector<double> y_a2(size_a2);
    FillVector(y_a2);
    std::cout << "Vector y:" << std::endl;
    PrintVector(y_a2);

    std::vector<double> b_a2 = MatrixMultiply(a2, y_a2);
    std::cout << "Vector b:" << std::endl;
    PrintVector(b_a2);

    std::cout << "Condition number for A2: " << ConditionNumber(a2) << std::endl;

    std::cout << "Gauss for matrix A2:" << std::endl;
    PrintVector(SolveMatrixGaussianElimination(a2, b_a2));

    std::cout << "LUP for matrix A2:" << std::endl;
    std::vector<std::vector<double>> l_a2(size_a2, std::vector<double>(size_a2)), u_a2(size_a2, std::vector<double>(size_a2));
    std::vector<double> p_a2(size_a2);

    LUP(a2, l_a2, u_a2, p_a2);
    std::cout << "LUP decomposition:\n";
    std::cout << "Matrix L:" << std::endl;
    PrintMatrix(l_a2);
    std::cout << "Matrix U:" << std::endl;
    PrintMatrix(u_a2);
    std::cout << "Vector P:" << std::endl;
    PrintVector(p_a2);
    std::cout << "Vector b (solution using LUP):" << std::endl;
    PrintVector(LUP_solve(a2, b_a2));

    std::cout << "Solution using square root method:" << std::endl;
    PrintVector(Sqrt_solve(a2, b_a2));

    std::cout << "LDLT decomposition:\n";
    std::vector<std::vector<double>> l1_a2(size_a2, std::vector<double>(size_a2));
    std::vector<double> d_a2(size_a2);
    if (LDLT(a2, l1_a2, d_a2))
    {
        std::cout << "Matrix L:" << std::endl;
        PrintMatrix(l1_a2);
        std::cout << "Vector D (diagonal):" << std::endl;
        PrintVector(d_a2);
    }
    else
    {
        std::cout << "Matrix is not positive definite" << std::endl
                  << std::endl;
    }

    std::cout << "Solution using Relaxation method:" << std::endl;
    PrintVector(Relaxation_solve(a2, b_a2, 1 - (8 / 40), count));

    std::cout << "\tTask 10:" << std::endl;
    std::cout << "Investigate cases for matrix with max condition number:\n"
              << std::endl;

    std::vector<double> vec_for_cm(matrix_max_cond_num.size());
    FillVector(vec_for_cm);

    std::vector<double> solve_r = Relaxation_solve(matrix_max_cond_num, vec_for_cm, 1-(8/40), count);
    for (int i = 0; i < 5; ++i)
    {
        for (int k = 0; k < matrix_max_cond_num.size(); ++k)
        {
            vec_for_cm[i] += 0.1;
        }

        std::vector<double> solution = Relaxation_solve(matrix_max_cond_num, vec_for_cm,  1-(8/40), count);
        double diff_norm = Norm(VectorSubtract(solution, solve_r));
        std::cout << "Norm for diff vector = " << diff_norm << std::endl;
    }

    std::cout << "\nFind errors for Relaxation method:\n" << std::endl;

    std::cout << "Param = 0.8:\n";
    std::vector<double> solve_1 = Relaxation_solve(matrix, vec_for_cm, 0.8, count);
    for (int i = 0; i < 5; ++i)
    {
        FillVector(vec_for_cm);

        std::vector<double> solution = Relaxation_solve(matrix_max_cond_num, vec_for_cm,  0.8, count);
        double diff_norm = Norm(VectorSubtract(solution, solve_1));
        std::cout << "Norm for diff vector = " << diff_norm << std::endl;
    }

    std::cout << "\nParam = 1.0\n";
    solve_1 = Relaxation_solve(matrix, vec_for_cm, 1.0, count);
    for (int i = 0; i < 5; ++i)
    {
        FillVector(vec_for_cm);

        std::vector<double> solution = Relaxation_solve(matrix_max_cond_num, vec_for_cm,  1.0, count);
        double diff_norm = Norm(VectorSubtract(solution, solve_1));
        std::cout << "Norm for diff vector = " << diff_norm << std::endl;
    }

    std::cout << "\nParam = 1.2\n";
    solve_1 = Relaxation_solve(matrix, vec_for_cm, 1.2, count);
    for (int i = 0; i < 5; ++i)
    {
        FillVector(vec_for_cm);

        std::vector<double> solution = Relaxation_solve(matrix_max_cond_num, vec_for_cm,  1.2, count);
        double diff_norm = Norm(VectorSubtract(solution, solve_1));
        std::cout << "Norm for diff vector = " << diff_norm << std::endl;
    }

    return 0;
}