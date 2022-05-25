#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <complex>
#include <string>
#include <cmath>
#include <algorithm>

using namespace std::chrono;
using namespace std;

#define sqr(x) x*x
#define PI 3.14159265358979323846
#define E 2.71828182845904523536

struct QR_Matrix{
    QR_Matrix()= default;
    vector<vector<double>> R;
    vector<vector<double>> T;
    vector<vector<double>> RQ;
};

struct Power_Matrix{
    Power_Matrix()= default;
    vector<complex<double>> prev_y1;
    vector<complex<double>> prev_y2;
    vector<complex<double>> prev_y3;
    vector<complex<double>> prev_y4;
    vector<complex<double>> y;
    vector<complex<double>> u;
    vector<complex<double>> prev_u1;
    vector<complex<double>> prev_u2;
    vector<complex<double>> prev_u3;
    vector<complex<double>> prev_u4;
    int iter_y1;
    int iter_y2;
    int iter_y3;
    int iter_y4;
};

// func that returns a random matrix
vector<vector<double>> GetRandomMatrix(int n, int m) {
    vector<vector<double>> matrix(n, vector<double>(m));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            matrix[i][j] = rand() % 10;
        }
    }
    return matrix;
}

// func that prints a matrix
void PrintMatrix(const vector<vector<double>>& matrix, ofstream& out) {
    for (auto & i : matrix) {
        for (auto j : i) {
            if (imag(j) == 0) {
                out << left << setw(15) << setprecision(6) << fixed << real(j) << " ";
            }
            else {
                out << left << setw(15) << setprecision(6) << fixed << real(j) << " + " << imag(j) << "i ";
            }
            //out << left << setw(15) << setprecision(6) << fixed << j << " ";
        }
        out << endl;
    }
}

// print complex vector
void PrintEigenvalues(const vector<complex<double>>& vector, ofstream& out) {
    for (auto i : vector) {
        if (imag(i) == 0) {
            out << left << setw(15) << setprecision(6) << fixed << real(i) << " ";
        }
        else {
            string s = to_string(real(i)) + " + " + to_string(imag(i)) + "i       ";
            out << s;
        }
    }
    out << endl;
}

// print complex vector
void PrintVector(const vector<complex<double>>& vector, ofstream& out) {
    for (auto i : vector) {
        if (imag(i) == 0) {
            out << left << setw(15) << setprecision(6) << fixed << real(i) << " ";
        }
        else {
            out << left << setw(15) << setprecision(6) << fixed << real(i) << " + " << imag(i) << "i ";
        }
    }
    out << endl;
}

// func that returns a random vector
vector<complex<double>> GetRandomVector(int n) {
    vector<complex<double>> vector(n);
    for (int i = 0; i < n; i++) {
        vector[i] = rand() % 10;
    }
    return vector;
}

// func that return square vector norm
double VectorNorm(const vector<double>& vector) {
    double norm = 0;
    for (double i : vector) {
        norm += i * i;
    }
    return sqrt(norm);
}

complex<double> VectorNorm(const vector<complex<double>>& vector) {
    complex<double> norm = 0;
    for (complex<double> i : vector) {
        norm += i * conj(i);
    }
    return sqrt(norm);
}

complex<double> VectorNormInf(const vector<complex<double>>& vector) {
    complex<double> norm = 0;
    for (complex<double> i : vector) {
        norm += abs(i);
    }
    return norm;
}

// func that divides a vector by divider
vector<complex<double>> VectorDivide(vector<complex<double>> vec, complex<double> divider){
    for (auto & i : vec) {
        i /= divider;
    }
    return vec;
}

// multiply matrix by matrix
vector<vector<double>> MatrixMultiplication(const vector<vector<double>>& matrix1, const vector<vector<double>>& matrix2) {
    vector<vector<double>> result(matrix1.size(), vector<double>(matrix2[0].size()));
    for (int i = 0; i < matrix1.size(); i++) {
        for (int j = 0; j < matrix2[0].size(); j++) {
            for (int k = 0; k < matrix1[0].size(); k++) {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
    return result;
}

// multiply matrix by vector
vector<double> MatrixVectorMultiplication(const vector<double>& vec, const vector<vector<double>>& matrix) {
    vector<double> result(matrix.size());
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            result[i] += vec[j] * matrix[i][j];
        }
    }
    return result;
}

vector<complex<double>> MatrixVectorMultiplication(const vector<complex<double>>& vec, const vector<vector<double>>& matrix) {
    vector<complex<double>> result(matrix.size());
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            result[i] += real(vec[j]) * matrix[i][j];
        }
    }
    return result;
}

vector<complex<double>> MatrixVectorMultiplication(const vector<complex<double>>& vec, const vector<vector<complex<double>>>& matrix) {
    vector<complex<double>> result(matrix.size());
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            result[i] += vec[j] * matrix[i][j];
        }
    }
    return result;
}

// divide two vectors
vector<double> VectorVectorDivision(const vector<double>& vec1, const vector<double>& vec2) {
    vector<double> result(vec1.size());
    for (int i = 0; i < vec1.size(); i++) {
        result[i] = vec1[i] / vec2[i];
    }
    return result;
}

vector<complex<double>> VectorVectorDivision(const vector<complex<double>>& vec1, const vector<complex<double>>& vec2) {
    vector<complex<double>> result(vec1.size());
    for (int i = 0; i < vec1.size(); i++) {
        result[i] = vec1[i] / vec2[i];
    }
    return result;
}

// multiply two vectors
vector<double> VectorVectorMultiplication(const vector<double>& vec1, const vector<double>& vec2) {
    vector<double> result(vec1.size());
    for (int i = 0; i < vec1.size(); i++) {
        result[i] = vec1[i] * vec2[i];
    }
    return result;
}

// subtract two vectors
vector<double> VectorVectorSubtraction(const vector<double>& vec1, const vector<double>& vec2) {
    vector<double> result(vec1.size());
    for (int i = 0; i < vec1.size(); i++) {
        result[i] = vec1[i] - vec2[i];
    }
    return result;
}

vector<complex<double>> VectorVectorSubtraction(const vector<complex<double>>& vec1, const vector<complex<double>>& vec2) {
    vector<complex<double>> result(vec1.size());
    for (int i = 0; i < vec1.size(); i++) {
        result[i] = vec1[i] - vec2[i];
    }
    return result;
}

// add two vectors
vector<double> VectorVectorAddition(const vector<double>& vec1, const vector<double>& vec2) {
    vector<double> result(vec1.size());
    for (int i = 0; i < vec1.size(); i++) {
        result[i] = vec1[i] + vec2[i];
    }
    return result;
}

vector<complex<double>> VectorVectorAddition(const vector<complex<double>>& vec1, const vector<complex<double>>& vec2) {
    vector<complex<double>> result(vec1.size());
    for (int i = 0; i < vec1.size(); i++) {
        result[i] = vec1[i] + vec2[i];
    }
    return result;
}

// scalar multiplication of two vectors
double ScalarVectorMultiplication(const vector<double>& vec1, const vector<double>& vec2) {
    double result = 0;
    for (int i = 0; i < vec1.size(); i++) {
        result += vec1[i] * vec2[i];
    }
    return result;
}

complex<double> ScalarVectorMultiplication(const vector<complex<double>>& vec1, const vector<complex<double>>& vec2) {
    complex<double> result = 0;
    for (int i = 0; i < vec1.size(); i++) {
        result += vec1[i] * vec2[i];
    }
    return result;
}

// multiply vector by scalar
vector<double> NumberVectorMultiplication(double number, const vector<double>& vec) {
    vector<double> result(vec.size());
    for (int i = 0; i < vec.size(); i++) {
        result[i] = number * vec[i];
    }
    return result;
}

vector<complex<double>> ComplexVectorMultiplication(complex<double> number, const vector<complex<double>>& vec) {
    vector<complex<double>> result(vec.size());
    for (int i = 0; i < vec.size(); i++) {
        result[i] = number * vec[i];
    }
    return result;
}

// add scalar to a vector
vector<double> NumberVectorAddition(double number, const vector<double>& vec) {
    vector<double> result(vec.size());
    for (int i = 0; i < vec.size(); i++) {
        result[i] = number + vec[i];
    }
    return result;
}

void PrintForPowerIteration(complex<double> lambda, const vector<complex<double>>& eigenvector, int count, auto time, ofstream& out){
    out << "Dominant eigenvalue: " << lambda << endl;
    out << "\n -------------------------------------------------- \n\n";
    out << "Eigenvector:" << endl;
    PrintVector(eigenvector, out);
    out << "\n -------------------------------------------------- \n\n";
    out << "Number of iterations: " << count << endl;
    out << "\n -------------------------------------------------- \n\n";
    out << "Time elapsed: " << time.count() << "s" << endl;
    out << "\n -------------------------------------------------- \n\n";
}

void PrintForPowerIteration(complex<double> lambda1, complex<double> lambda2, const vector<complex<double>>& eigenvector1,
                            const vector<complex<double>>& eigenvector2, int count, auto time, ofstream& out){
    out << "Dominant eigenvalue1: " << lambda1 << endl;
    out << "\n -------------------------------------------------- \n\n";
    out << "Dominant eigenvalue2: " << lambda2 << endl;
    out << "\n -------------------------------------------------- \n\n";
    out << "Eigenvector1:" << endl;
    PrintVector(eigenvector1, out);
    out << "\n -------------------------------------------------- \n\n";
    out << "Eigenvector2:" << endl;
    PrintVector(eigenvector2, out);
    out << "\n -------------------------------------------------- \n\n";
    out << "Number of iterations: " << count << endl;
    out << "\n -------------------------------------------------- \n\n";
    out << "Time elapsed: " << time.count() << "s" << endl;
    out << "\n -------------------------------------------------- \n\n";
}

int GetMaxIndexNorm(const vector<complex<double>> &vec) {
    int index = 0;
    complex<double> max = (0, 0);
    for (int i = 0; i < vec.size(); i++) {
        if (real(vec[i]) > real(max)) {
            max = abs(vec[i]);
            index = i;
        }
    }
    return index;
}

// not working...
// func that uses power iteration to find eigenvector and max eigenvalue of a matrix
complex<double> PowerIteration(vector<vector<double>> matrix, double epsilon, ofstream& out) {
    // cast matrix to complex matrix
    vector<vector<complex<double>>> complex_matrix(matrix.size());
    for (int i = 0; i < matrix.size(); i++) {
        complex_matrix[i].resize(matrix.size());
        for (int j = 0; j < matrix.size(); j++) {
            complex_matrix[i][j] = complex<double>(matrix[i][j], 0);
        }
    }

    Power_Matrix power_matrix{};

    power_matrix.y = GetRandomVector(matrix.size());
    power_matrix.u = VectorDivide(power_matrix.y, VectorNorm(power_matrix.y));
    power_matrix.prev_y1 = power_matrix.y;
    power_matrix.prev_y2 = power_matrix.y;
    power_matrix.prev_y3 = power_matrix.y;
    power_matrix.prev_u1 = power_matrix.u;
    power_matrix.prev_u2 = power_matrix.u;
    power_matrix.prev_u3 = power_matrix.u;
    power_matrix.iter_y1 = 0;
    power_matrix.iter_y2 = 0;
    power_matrix.iter_y3 = 0;

    vector<complex<double>> u1;

    complex<double> lambda = ScalarVectorMultiplication(power_matrix.u, MatrixVectorMultiplication(power_matrix.u, matrix));
    int count = 0;

    auto start = chrono::steady_clock::now();

    while (true) {
        // fill previous steps
        if (power_matrix.iter_y4 == count - 4) {
            power_matrix.prev_y4 = power_matrix.prev_y3;
            power_matrix.prev_u4 = power_matrix.prev_u3;
            power_matrix.iter_y4 = count;
        }
        if (power_matrix.iter_y3 == count - 3) {
            power_matrix.prev_y3 = power_matrix.prev_y2;
            power_matrix.prev_u3 = power_matrix.prev_u2;
            power_matrix.iter_y3 = count;
        }
        if (power_matrix.iter_y2 == count - 2) {
            power_matrix.prev_y2 = power_matrix.prev_y1;
            power_matrix.prev_u2 = power_matrix.prev_u1;
            power_matrix.iter_y2 = count;
        }
        if (power_matrix.iter_y1 == count - 1) {
            power_matrix.prev_y1 = power_matrix.y;
            power_matrix.prev_u1 = power_matrix.u;
            power_matrix.iter_y1 = count;
        }

        power_matrix.y = MatrixVectorMultiplication(power_matrix.u, complex_matrix);
        power_matrix.u = VectorDivide(power_matrix.y, VectorNorm(power_matrix.y));
        lambda = ScalarVectorMultiplication(power_matrix.u, MatrixVectorMultiplication(power_matrix.u, complex_matrix));
        ++count;

        if (real(VectorNorm(VectorVectorSubtraction(MatrixVectorMultiplication(power_matrix.u, complex_matrix),
                                               ComplexVectorMultiplication(lambda, power_matrix.u)))) < epsilon) {
            PrintForPowerIteration(lambda, power_matrix.u, count, chrono::steady_clock::now() - start, out);
            return lambda;
        }

        // case two is the same as case one (lambda is the same)

        // check if there are two equal but opposite eigenvalues (case 3)
        if (count > 3) {
            // k starts from 0
            vector<complex<double>> check_max = VectorVectorDivision(
                    ComplexVectorMultiplication(VectorNormInf(power_matrix.prev_y2), power_matrix.prev_y1),
                    power_matrix.prev_u3);
            lambda = sqrt(VectorNormInf(check_max));

            power_matrix.u = VectorVectorAddition(power_matrix.prev_y2, ComplexVectorMultiplication(lambda, power_matrix.prev_u3));
            u1 = VectorVectorSubtraction(power_matrix.prev_y2,
                                         ComplexVectorMultiplication(lambda, power_matrix.prev_u3));

            if (real(VectorNorm(VectorVectorSubtraction(MatrixVectorMultiplication(power_matrix.u, complex_matrix),
                                                   ComplexVectorMultiplication(lambda, power_matrix.u)))) < epsilon) {
                PrintForPowerIteration(lambda, -lambda, power_matrix.u, u1, count, chrono::steady_clock::now() - start, out);
                return lambda;
            }
        }

        // check for complex pair (case 4)
        if (count > 5) {
            int index = GetMaxIndexNorm(power_matrix.prev_y1);

            double yk = real(power_matrix.prev_y3[index]);
            double yk_2 = real(power_matrix.prev_y1[index]);
            double yk_1 = real(power_matrix.prev_y2[index]);
            double uk = real(power_matrix.prev_u3[index]);
            double uk1 = real(power_matrix.prev_u4[index]);

            double r = sqrt(abs(((yk * yk_2 * real(VectorNormInf(power_matrix.prev_y2))) - (yk_1 * yk_1 * real(VectorNormInf(power_matrix.prev_y3)))) /
                    ((uk1 * yk_1) - uk * uk * real(VectorNormInf(power_matrix.prev_y3)))));

            //double cos_ = max(-1, min((yk_2 + r * r * yk) / (2 * r * yk_1), 1));
            double cos_ = clamp((yk_2 * real(VectorNormInf(power_matrix.prev_y2)) + r * r * uk) / (2 * r * yk_1), -1.0, 1.0);
            double sin_ = sqrt(abs(1 - cos_ * cos_));

            lambda = complex<double>(r * cos_, r * sin_);
            complex<double> lambda_2 = complex<double>(r * cos_, -r * sin_);

            power_matrix.u = VectorVectorSubtraction(power_matrix.prev_y2,
                                          ComplexVectorMultiplication(lambda_2, power_matrix.prev_u3));
            u1 = VectorVectorSubtraction(power_matrix.prev_y2,
                                          ComplexVectorMultiplication(lambda, power_matrix.prev_u3));

            if (real(VectorNorm(VectorVectorSubtraction(MatrixVectorMultiplication(power_matrix.u, complex_matrix),
                                                        ComplexVectorMultiplication(lambda, power_matrix.u)))) < epsilon) {
                auto stop = chrono::steady_clock::now();
                std::chrono::duration<double> elapsed = stop - start;
                PrintForPowerIteration(lambda, lambda_2, power_matrix.u, u1, count, elapsed, out);
                return 0;
            }
        }

        if (count > 1500) {
            out << "Too many iterations" << endl;
            return -1;
        }
    }

    out << "WTF case." << endl;
    return 0;
}

// used to find Hessenberg form of a matrix
vector<vector<double>> HessenbergForm(const vector<vector<double>>& matrix, ofstream &out) {
    double c, s;
    vector<vector<double>> matrix_new(matrix.size(), vector<double>(matrix.size()));
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix.size(); j++) {
            matrix_new[i][j] = matrix[i][j];
        }
    }

    for (int i = 1; i < matrix_new[0].size() - 1; i++) {
        for (int j = 1 + i; j < matrix_new.size(); j++) {
            if (fabs(matrix_new[i][i - 1]) + fabs(matrix_new[j][i - 1]) != 0) {
                c = matrix_new[i][i - 1] / sqrt(sqr(matrix_new[i][i - 1]) + sqr(matrix_new[j][i - 1]));
                s = matrix_new[j][i - 1] / sqrt(sqr(matrix_new[i][i - 1]) + sqr(matrix_new[j][i - 1]));

                // find matrix R
                for (int m = i - 1; m < matrix_new[i].size(); m++) { ;
                    double temp = matrix_new[i][m];
                    matrix_new[i][m] = c * temp + s * matrix_new[j][m];
                    matrix_new[j][m] = -s * temp + c * matrix_new[j][m];

                }

                for (int m = 0; m < matrix_new[i].size(); m++) {// multiply by transpose
                    double temp = matrix_new[m][i];
                    matrix_new[m][i] = c * temp + s * matrix_new[m][j];
                    matrix_new[m][j] = -s * temp + c * matrix_new[m][j];

                }

                // create new matrix
                vector<vector<complex<double>>> MatrixNew(matrix.size(), vector<complex<double>>(matrix[0].size()));
                for (int i = 0; i < matrix.size(); i++) {
                    for (int j = 0; j < matrix[0].size(); j++) {
                        if (j == i) {
                            MatrixNew[i][j] = 1;
                        } else MatrixNew[i][j] = 0;
                    }
                }

                MatrixNew[i][i] = c;
                MatrixNew[j][j] = c;
                MatrixNew[i][j] = s;
                MatrixNew[j][i] = -s;

            }
        }

    }

    out << "Hessenberg Matrix:" << endl;
    PrintMatrix(matrix_new, out);
    std::cout << std::endl;
    return matrix_new;
}

vector<vector<double>> product_rotation_matrix(const vector<vector<double>> &matrix1, const vector<vector<double>> &matrix2, int str, int row){
    vector<vector<double>> product_rotation_matrix(matrix2.size(), vector<double>(matrix2[0]. size()));
    for (int i = 0; i < matrix2.size(); i++){
        for (int j = 0; j < matrix2[0].size(); j++){
            product_rotation_matrix[i][j] = matrix2[i][j];
        }
    }
    for (int k = 0; k < matrix2[0]. size(); k ++){
        product_rotation_matrix[str][k] = matrix2[str][k] * matrix1[str][str] + matrix2[row][k] * matrix1[str][row];
        product_rotation_matrix[row][k] = matrix2[row][k] * matrix1[str][str] - matrix2[str][k] * matrix1[str][row];
    }

    return product_rotation_matrix;
}

// multiply R(upper-triangle) by T(transpose)
vector<vector<double>> product_matrix(const vector<vector<double>> &matrix1, const vector<vector<double>> &matrix2){
    int new_row = matrix2.size();
    int new_col = matrix2[0].size();
    vector<vector<double>> matrix_product(new_row, vector<double>(new_col));
    for (int i = 0; i < new_row; i++){
        for (int j = 0; j < new_col; j++){
            matrix_product[i][j] = 0;
            for (int k = 0; k < new_col; k++){
                matrix_product[i][j] = matrix_product[i][j] + (matrix1[i][k]*matrix2[k][j]);
            }
        }

    }
    return matrix_product;
}

// implementation of QR algorithm
QR_Matrix qr_method(const vector<vector<double>> &matrix) {
    double c, s;
    QR_Matrix struct_Matrix;
    struct_Matrix.T = vector<vector<double>>(matrix.size(),vector<double>(matrix.size()));
    struct_Matrix.R = vector<vector<double>>(matrix.size(),vector<double>(matrix[0].size()));
    struct_Matrix.RQ = vector<vector<double>>(matrix.size(),vector<double>(matrix.size()));


    vector<vector<double>> matrix_r_new(matrix.size(), vector<double>(matrix[0].size()));
    // write input data
    for (int i = 0; i < matrix.size(); i ++){
        for (int j = 0; j < matrix[0].size(); j ++){
            matrix_r_new[i][j] = matrix[i][j];
        }
    }
    // create matrix of ones (square)
    vector<vector<double>> MatrixInitial(matrix.size(), vector<double>(matrix.size()));
    for (int i = 0; i < matrix.size(); i++){
        for (int j = 0; j < matrix.size(); j++){
            if (i == j){
                MatrixInitial[i][j] = 1;
            }
            else{
                MatrixInitial[i][j] = 0;
            }
        }
    }
    for (int i = 0; i < matrix_r_new[0].size() - 1; i++){
        for(int j = 1 + i; j < matrix_r_new.size(); j++){
            if (abs(matrix_r_new[i][i]) + abs(matrix_r_new[j][i]) != 0){
                c = matrix_r_new[i][i] / sqrt(sqr(matrix_r_new[i][i]) + sqr(matrix_r_new[j][i]));
                s = matrix_r_new[j][i] / sqrt(sqr(matrix_r_new[i][i]) + sqr(matrix_r_new[j][i]));
                // find R
                for (int m = i; m < matrix_r_new[i].size(); m++){
                    double temp = matrix_r_new[i][m];
                    matrix_r_new[i][m] = c * temp + s * matrix_r_new[j][m];
                    matrix_r_new[j][m] = -s * temp + c * matrix_r_new[j][m];
                }

                // get new matrix
                vector<vector<double>> MatrixNew(matrix_r_new.size(), vector<double>(matrix_r_new[0].size()));
                for (int i = 0; i < matrix_r_new.size(); i++){
                    for (int j = 0; j < matrix_r_new[0].size(); j++){
                        if (j == i){
                            MatrixNew[i][j] = 1;
                        }
                        else MatrixNew[i][j] = 0;
                    }
                }

                MatrixNew[i][i] = c;
                MatrixNew[j][j] = c;
                MatrixNew[i][j] = s;
                MatrixNew[j][i] = -s;

                //multiply rotate matrix
                MatrixInitial = product_rotation_matrix(MatrixNew, MatrixInitial, i , j); // get matrix T
            }
        }
    }

    // transpose T => get Q
    vector<vector<double>> transpose_matrix(MatrixInitial.size(),vector<double>(MatrixInitial.size()));
    for(int i = 0; i < MatrixInitial.size(); i++){
        for (int j = 0; j < MatrixInitial.size(); j++){
            transpose_matrix[i][j] = MatrixInitial[j][i];

        }
    }

    // multiply R = TA by T(transpose) = RQ
    vector<vector<double>> matrix_new;
    matrix_new = product_matrix(matrix_r_new, transpose_matrix);

    for(int i = 0; i < matrix.size(); i++){
        for (int j = 0; j < matrix[0].size(); j++){
            struct_Matrix.R[i][j] = matrix_r_new[i][j];
        }
    }

    for(int i = 0; i < matrix.size(); i++){
        for (int j = 0; j < matrix.size(); j++){
            struct_Matrix.T[i][j] = MatrixInitial[i][j];
        }
    }


    for(int i = 0; i < matrix.size(); i++){
        for (int j = 0; j < matrix.size(); j++){
            struct_Matrix.RQ[i][j] = matrix_new[i][j];
        }
    }

    return struct_Matrix;
}

// solve square equation to get a pair of complex eigenvalues
void solve_square_equation(complex<double> cof1, complex<double> cof2, complex<double>& x1, complex<double> &x2){
    complex<double> d = sqr(cof1) - 4 * real(cof2);
    x1 = (-cof1 + sqrt(d));
    x2 = (-cof1 - sqrt(d));
    (x1).real((x1).real() / 2);
    (x2).real((x2).real() / 2);
    (x1).imag((x1).imag() / 2);
    (x2).imag((x2).imag() / 2);
}

// function that returns eigenvalues of a matrix found by QR algorithm
vector<complex<double>> get_eigenvalues(vector<vector<double>> matrix){
    vector<complex<double>> eigenvalues;
    int i;
    for (i = 0; i < matrix.size()-1; i++){
        if (abs(real(matrix[i+1][i])) < 0.001){
            eigenvalues.emplace_back(matrix[i][i]);
        }
        else{
            complex<double> cof1 = matrix[i][i+1] * matrix[i+1][i];
            complex<double> cof2 = matrix[i][i] * matrix[i+1][i+1];
            complex<double> cof3 = cof1 - cof2;
            complex<double> cof4 = matrix[i][i] + matrix[i+1][i+1];
            complex<double> x1, x2;
            solve_square_equation(abs(cof4), abs(cof3), x1, x2);
            eigenvalues.push_back(x1);
            eigenvalues.push_back(x2);
            i += 1;
        }
    }
    if (i == matrix.size()-1){
        eigenvalues.emplace_back(matrix[i][i]);
    }
    return eigenvalues;
}

// nonlinear function given in task
double nonlinear_func(double x){
    return ((pow(x, 9) + PI) * cos(log(x * x + 1))) / (pow(E, x*x)) - (x / 2022);
}

double nonlinear_func_derivative(double x){
    return (-(2*pow(E,(-x*x))*(pow(x,9) + PI)*x*sin(log(x*x + 1)))/(x*x + 1)) - (2*pow(E,(-x*x))*(pow(x,9) + PI)*x*cos(log(x*x + 1))) +
    9*pow(E,(-x*x))*pow(x,8)*cos(log(x*x + 1)) - 1/2022;
}

// function to narrow down the interval that contains roots of nonlinear function
void bisection_method(double& a, double& b, double eps, ofstream& out){
    auto start = std::chrono::high_resolution_clock::now();
    int count = 0;
    while (fabs(b - a) > 2 * eps){
        double x = (a + b) / 2;
        if (nonlinear_func(x) * nonlinear_func(a) < 0){
            b = x;
        }
        else{
            a = x;
        }
        ++count;
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    out << "Number of iterations for bisection method: " << count << endl;
    out << "Time elapsed for bisection method: " << elapsed.count() << " s" << endl;
}

double Newton_method(double x, double eps, ofstream& out){
    int count = 0;
    auto start = std::chrono::high_resolution_clock::now();
    double x_new = x - nonlinear_func(x) / nonlinear_func_derivative(x);
    while (fabs(x_new - x) > eps) {
        x = x_new;
        x_new = x - nonlinear_func(x) / nonlinear_func_derivative(x);
        ++count;
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    out << "Number of iterations for Newton method: " << count << endl;
    out << "Time elapsed for Newton method: " << elapsed.count() << " s" << endl;
    return x_new;
}

vector<vector<double>> get_a1(){
    return {{1,-2,1,0,-1,1,-2,2,0,-2},
            {0,2,0,0,2,1,-1,-1,-1,-2},
            {0,1,0,-1,1,-1,0,-1,1,-1},
            {-2,-1,2,-1,0,0,0,0,1,0},
            {1,-2,0,1,0,-2,-1,0,2,2},
            {-2,-2,0,-2,0,1,1,-2,1,1},
            {-1,-2,-1,-1,-2,-1,-2,1,-1,2},
            {-2,1,2,-2,0,2,1,-1,-2,2},
            {0,1,0,1,1,-2,2,0,1,1},
            {0,0,2,-1,-1,0,-2,2,-1,-1}
    };
}

vector<vector<double>> get_a2(){
    return {{-1,1,-1,0,-1,0,-1,1,1,-1,0,-1,-1,1,0,0,1,1,1,1},
            {-1,0,-1,1,-1,0,0,0,0,-1,0,0,-1,1,0,-1,1,-1,-1,0},
            {1,0,-1,1,0,1,-1,-1,-1,0,-1,-1,1,-1,1,1,-1,1,-1,0},
            {-1,1,0,0,-1,0,0,-1,0,-1,1,1,-1,-1,1,1,-1,1,-1,0},
            {1,0,-1,0,0,-1,1,1,0,0,0,1,1,1,0,0,-1,0,0,1},
            {0,0,0,0,-1,1,1,0,0,1,1,0,-1,0,1,1,0,1,0,0},
            {-1,0,1,1,1,-1,-1,0,-1,1,-1,-1,-1,0,-1,0,0,0,-1,1},
            {0,0,-1,-1,0,1,1,1,1,-1,0,0,-1,1,1,1,1,0,0,-1},
            {0,0,1,1,0,1,1,0,1,-1,1,0,0,0,1,1,0,0,0,1},
            {0,-1,0,0,1,0,-1,0,-1,0,-1,0,-1,0,1,-1,0,0,1,1},
            {1,-1,1,-1,-1,-1,1,0,-1,0,1,1,-1,0,1,1,1,0,0,0},
            {0,1,0,0,-1,0,1,0,1,0,0,1,1,-1,-1,0,-1,1,1,-1},
            {-1,-1,-1,-1,0,1,-1,0,0,-1,0,0,0,1,1,0,0,0,-1,0},
            {-1,0,1,0,-1,0,0,1,-1,1,1,-1,1,1,1,-1,1,-1,-1,0},
            {1,-1,0,-1,-1,0,-1,-1,0,0,1,0,1,1,-1,1,0,0,-1,0},
            {-1,-1,1,0,-1,1,1,-1,1,0,0,-1,1,-1,-1,0,0,1,1,1},
            {0,0,-1,0,0,0,0,-1,1,1,0,-1,1,-1,0,0,0,-1,-1,1},
            {-1,0,-1,-1,-1,1,1,-1,1,-1,1,-1,1,-1,1,1,0,-1,0,-1},
            {-1,0,1,0,0,0,0,-1,1,-1,1,-1,0,-1,-1,1,0,1,0,0},
            {0,-1,-1,1,-1,1,-1,-1,-1,1,1,-1,0,-1,-1,0,1,0,-1,-1}
    };
}

// find dominant eigenvalue and eigenvector using Power Iteration method
// test speed on random matrix
bool task1(){
    // Power Method for matrix a1
    ofstream out("power_iterations.txt");

    PowerIteration(get_a1(), 1e-6, out);
    // end for a1

    // Power Method for matrix a2
    PowerIteration(get_a2(), 1e-6, out);
    // end for a2

    vector<vector<double>> test1 = GetRandomMatrix(10, 10);
    PowerIteration(test1, 1e-6, out);

    vector<vector<double>> test2 = GetRandomMatrix(30, 30);
    PowerIteration(test1, 1e-6, out);

    out.close();

    return true;
}

// realize QR algorithm for matrix a1, a2 and test speed on random matrix
bool task2(){
    ofstream out("qr_algorithm.txt");
    ofstream out1 ("hessenberg_forms.txt");

    out << "QR algorithm:" << endl;

    // hessenberg form for a1
    out1 << "Hessenberg form of matrix a1:" << endl;
    vector<vector<double>> hessenberg_a1 = HessenbergForm(get_a1(), out1);
    out1 << "\n -------------------------------------------------- \n\n";

    // hessenberg form for a2
    out1 << "Hessenberg form of matrix a2:" << endl;
    vector<vector<double>> hessenberg_a2 = HessenbergForm(get_a2(), out1);
    out1 << "\n -------------------------------------------------- \n\n";

    /*
    //test
    vector<vector<complex<double>>> test_matrix = {
            { 5, 3, -5, -3},
            { 1, -24.4, 32.4, 23.4},
            { 4, -38.4, 49.4, 38.4},
            { -8, -8, 0, 7},
    };

    QR_Matrix test = qr_method(test_matrix);
    PrintMatrix(test.RQ, out);

    //end test
     */

    // check time for QR method for matrix a1
    auto start1 = std::chrono::high_resolution_clock::now();
    QR_Matrix a1_ = qr_method(hessenberg_a1);

    for (int i = 0; i < 1000; ++i) {
        a1_ = qr_method(a1_.RQ);
    }

    vector<complex<double>> eigenvector1 = get_eigenvalues(a1_.RQ);
    auto finish1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed1 = finish1 - start1;
    // end for a1

    // check time for QR method a2
    auto start2 = std::chrono::high_resolution_clock::now();
    QR_Matrix a2_ = qr_method(hessenberg_a2);

    for (int i = 0; i < 7000; ++i) {
        a2_ = qr_method(a2_.RQ);
    }

    vector<complex<double>> eigenvector2 = get_eigenvalues(a2_.RQ);
    auto finish2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed2 = finish2 - start2;
    // end for a2

    // output for a1
    out << "Time for QR algorithm (matrix a1): " << elapsed1.count() << " s" << endl;

    PrintMatrix(a1_.RQ, out);
    out << "\n -------------------------------------------------- \n\n";

    out << "Eigenvalues for matrix a1:" << endl;
    PrintEigenvalues(eigenvector1, out);
    out << "\n -------------------------------------------------- \n\n";
    // end for a1

    // output for a2
    out << "Time for QR algorithm (matrix a2): " << elapsed2.count() << " s" << endl;
    PrintMatrix(a2_.RQ, out);
    out << "\n -------------------------------------------------- \n\n";

    out << "Eigenvalues for matrix a2:" << endl;
    PrintEigenvalues(eigenvector2, out);
    out << "\n -------------------------------------------------- \n\n";
    // end for a2

    out1.close();
    out.close();

    // time tests
    ofstream out3 ("time_tests.txt");
    ofstream out4 ("test_hessenberg.txt");
    out3 << "Time tests:" << endl;

    vector<vector<double>> test1 = GetRandomMatrix(10, 10);
    vector<vector<double>> test2 = GetRandomMatrix(20, 20);
    vector<vector<double>> test3 = GetRandomMatrix(30, 30);
    vector<vector<double>> test4 = GetRandomMatrix(40, 40);

    out4 << "Hessenberg form of random matrix 10x10:" << endl;
    vector<vector<double>> hessenberg_test1 = HessenbergForm(test1, out4);
    out4 << "\n -------------------------------------------------- \n\n";
    out4 << "Hessenberg form of random matrix 20x20:" << endl;
    vector<vector<double>> hessenberg_test2 = HessenbergForm(test2, out4);
    out4 << "\n -------------------------------------------------- \n\n";
    out4 << "Hessenberg form of random matrix 30x30:" << endl;
    vector<vector<double>> hessenberg_test3 = HessenbergForm(test3, out4);
    out4 << "\n -------------------------------------------------- \n\n";
    out4 << "Hessenberg form of random matrix 40x40:" << endl;
    vector<vector<double>> hessenberg_test4 = HessenbergForm(test4, out4);
    out4 << "\n -------------------------------------------------- \n\n";

    auto start3 = std::chrono::high_resolution_clock::now();
    QR_Matrix test1_ = qr_method(hessenberg_test1);

    for (int i = 0; i < 100; ++i) {
        test1_ = qr_method(test1_.RQ);
    }
    vector<complex<double>> eigenvector_test1 = get_eigenvalues(test1_.RQ);

    auto finish3 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed3 = finish3 - start3;

    auto start4 = std::chrono::high_resolution_clock::now();
    QR_Matrix test2_ = qr_method(hessenberg_test2);

    for (int i = 0; i < 100; ++i) {
        test2_ = qr_method(test2_.RQ);
    }
    vector<complex<double>> eigenvector_test2 = get_eigenvalues(test2_.RQ);

    auto finish4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed4 = finish4 - start4;

    auto start5 = std::chrono::high_resolution_clock::now();
    QR_Matrix test3_ = qr_method(hessenberg_test3);

    for (int i = 0; i < 100; ++i) {
        test3_ = qr_method(test3_.RQ);
    }
    vector<complex<double>> eigenvector_test3 = get_eigenvalues(test3_.RQ);

    auto finish5 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed5 = finish5 - start5;

    auto start6 = std::chrono::high_resolution_clock::now();
    QR_Matrix test4_ = qr_method(hessenberg_test4);

    for (int i = 0; i < 100; ++i) {
        test4_ = qr_method(test4_.RQ);
    }
    vector<complex<double>> eigenvector_test4 = get_eigenvalues(test4_.RQ);

    auto finish6 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed6 = finish6 - start6;

    out3 << "Time for QR algorithm (10x10): " << elapsed3.count() << " s" << endl;
    out3 << "Time for QR algorithm (20x20): " << elapsed4.count() << " s" << endl;
    out3 << "Time for QR algorithm (30x30): " << elapsed5.count() << " s" << endl;
    out3 << "Time for QR algorithm (40x40): " << elapsed6.count() << " s" << endl;

    return true;
}

bool task3(){
    ofstream nonlinear_solve("nonlinear_solve.txt");
    nonlinear_solve << "Nonlinear solve:" << endl;
    // first interval: (-2, -1.5)
    // second interval: (-1.5, -1)
    // third interval: (1, 2)
    double a1_ = -2, b1 = -1.5, a2_ = -1.5, b2 = -1, a3 = 1, b3 = 2;
    double eps = 1e-4;

    bisection_method(a1_, b1, eps, nonlinear_solve);
    nonlinear_solve << "Interval found by bisection method: [" << setprecision(15) << fixed << a1_ << ", " << setprecision(15) << fixed << b1 << "]" << endl;
    double x1 = Newton_method(a1_, 1e-16, nonlinear_solve);
    nonlinear_solve << "x1 = " << setprecision(15) << fixed << x1 << endl;
    nonlinear_solve << "\n -------------------------------------------------- \n\n";

    bisection_method(a2_, b2, eps, nonlinear_solve);
    nonlinear_solve << "Interval found by bisection method: [" << setprecision(15) << fixed << a2_ << ", " << setprecision(15) << fixed << b2 << "]" << endl;
    double x2 = Newton_method(a2_, 1e-16, nonlinear_solve);
    nonlinear_solve << "x2 = " << setprecision(15) << fixed << x2 << endl;
    nonlinear_solve << "\n -------------------------------------------------- \n\n";

    bisection_method(a3, b3, eps, nonlinear_solve);
    nonlinear_solve << "Interval found by bisection method: [" << setprecision(15) << fixed << a3 << ", " << setprecision(15) << fixed << b3 << "]" << endl;
    double x3 = Newton_method(a3, 1e-16, nonlinear_solve);
    nonlinear_solve << "x3 = " << setprecision(15) << fixed << x3 << endl;
    nonlinear_solve.close();
    return true;
}

int main() {

    // task 1
    if (task1()){
        cout << "Task 1 is done" << endl;
    } else {
        cout << "Task 1 is failed" << endl;
    }

    // task 2
    if (task2()){
        cout << "Task 2 is done" << endl;
    } else {
        cout << "Task 2 is failed" << endl;
    }

    // task 3
    if (task3()){
        cout << "Task 3 is done" << endl;
    } else {
        cout << "Task 3 is failed" << endl;
    }

    return 0;
}
