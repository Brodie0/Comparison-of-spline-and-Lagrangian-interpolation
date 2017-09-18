#include <iostream>
#include <iomanip>
#include <fstream>
#include "Matrix.h"

	Matrix::Matrix(const int r, const int c)
		:rows(r), cols(c) {
		matrix = new double *[r];
		for (int y = 0; y < r; y++) {
			matrix[y] = new double[c];
			for (int x = 0; x < c; x++) {
				matrix[y][x] = 0;
			}
		}
	}

	Matrix::Matrix(const Matrix& m)
		:rows(m.rows), cols(m.cols) {
		matrix = new double *[rows];
		for (int y = 0; y < rows; y++) {
			matrix[y] = new double[cols];
			for (int x = 0; x < cols; x++) {
				matrix[y][x] = m.matrix[y][x];
			}
		}
	}

	Matrix::~Matrix() {
		if (matrix != nullptr) {
			for (int i = 0; i < rows; i++) {
				if (matrix[i] != nullptr) {
					delete[] matrix[i];
					matrix[i] = nullptr;
				}
			}
		}
		delete[] matrix;
		matrix = nullptr;
	}

	void Matrix::print() {
		for (int y = 0; y < rows; y++) {
			for (int x = 0; x < cols; x++) {
				std::cout << std::setprecision(30) << matrix[y][x]<<" ";
			}
			std::cout << "\n";
		}
		std::cout << "\n";
	}

	void Matrix::switchRows(int r1, int r2) {
		if (r1 > rows || r2 > rows) {
			std::cout << "Blad zamiany wierszy";
			return;
		}
		double tmp;
		for (int i = 0; i < cols; i++) {
			tmp = matrix[r1][i];
			matrix[r1][i] = matrix[r2][i];
			matrix[r2][i] = tmp;
		}
	}

	int Matrix::findMaxInCol(int row, int from, int to) {
		int currRow = from;
		double cur = matrix[row][from];
		for (int i = from; i <= to; i++) {
			double tmp = matrix[row][i];
			if (cur < tmp) {
				currRow = i;
				cur = tmp;
			}
		}
		return currRow;
	}

	Matrix& Matrix::operator=(const Matrix& m) {
		if (matrix != nullptr) {
			for (int i = 0; i < rows; i++) {
				if (matrix[i] != nullptr) {
					delete[] matrix[i];
					matrix[i] = nullptr;
				}
			}
		}
		delete[] matrix;
		rows = m.rows;
		cols = m.cols;
		matrix = new double *[rows];
		for (int y = 0; y < rows; y++) {
			matrix[y] = new double[cols];
			for (int x = 0; x < cols; x++) {
				matrix[y][x] = m.matrix[y][x];
			}
		}
		return *this;
	}

	Matrix Matrix::operator+(const Matrix& m) {
		Matrix ret(rows, cols);
		if (m.rows != rows && m.cols != cols) {
			std::cout << "Blad wymiarow macierzy";
			return ret;
		}
		for (int y = 0; y < rows; y++) {
			for (int x = 0; x < cols; x++) {
				ret.matrix[y][x] = matrix[y][x] + m.matrix[y][x];
			}
		}
		return ret;
	}

	Matrix Matrix::operator-(const Matrix& m) {
		Matrix ret(rows, cols);
		if (m.rows != rows && m.cols != cols) {
			std::cout << "Blad wymiarow macierzy";
			return ret;
		}

		for (int y = 0; y < rows; y++) {
			for (int x = 0; x < cols; x++) {
				ret.matrix[y][x] = matrix[y][x] - m.matrix[y][x];
			}
		}
		return ret;
	}

	Matrix Matrix::operator*(const Matrix& right) {
		Matrix ret(rows, right.cols);
		if (cols != right.rows) {
			std::cout << "Blad rozmiarow mnozenie";
			return ret;
		}
		for (int y = 0; y < rows; y++) {
			for (int x = 0; x < right.cols; x++) {
				for (int r = 0; r < right.rows; r++) {
					ret.matrix[y][x] += matrix[y][r] * right.matrix[r][x];
				}
			}
		}
		return ret;
	}

	Matrix Matrix::JacobiMethod::solve(Matrix A, Matrix b, int iterations = 5000) {
		Matrix ret(A.rows, 1);
		Matrix retNew(A.rows, 1);
		Matrix residuum = Matrix(A.cols, 1);
		double norm = 0.0;
		for (int iter = 0; iter < iterations; iter++) {
			for (int y = 0; y < A.rows; y++) {
				double actual = 0;
				for (int j = 0; j < y; j++) {		//dolny trojkat
					actual -= (A.matrix[y][j] * ret.matrix[j][0]);
				}
				for (int j = y + 1; j < A.rows; j++) {	//gorny, wszystko przez x(k-1)
					actual -= (A.matrix[y][j] * ret.matrix[j][0]);
				}
				actual += b.matrix[y][0];
				actual /= A.matrix[y][y];	//dzielimy przez wartosci z diagonali
				retNew.matrix[y][0] = actual;
			}
			ret = retNew;
			residuum = A*ret - b;
			double sum = 0.0;
			for (int j = 0; j < residuum.cols; j++)
				sum += powl(residuum.matrix[j][0], 2.0);
			norm = sqrtl(sum);
			if (norm <= 1e-9) {
				std::cout << std::endl << "Osiagnieto norme w Jacobim, wynosi ona: " << iter << std::endl;
				return ret;
			}
		}
		return ret;
	}

	Matrix Matrix::GaussSiedelMethod::solve(Matrix A, Matrix b, int iterations = 5000) {
		Matrix ret(A.rows, 1);
		Matrix retNew(A.rows, 1);
		Matrix residuum = Matrix(A.cols, 1);
		double norm = 0.0;
		//
		for (int iter = 0; iter < iterations; iter++) {
			for (int y = 0; y < A.rows; y++) {
				double actual = 0;
				for (int j = 0; j < y; j++) {		//dla dolnego trojkata odejmujemy wartosci przemnozone przez x(k)
					actual -= (A.matrix[y][j] * retNew.matrix[j][0]);
				}
				for (int j = y + 1; j < A.rows; j++) {	//dla gornego x(k-1)
					actual -= (A.matrix[y][j] * ret.matrix[j][0]);
				}
				actual += b.matrix[y][0];
				actual /= A.matrix[y][y];
				retNew.matrix[y][0] = actual;
			}
			ret = retNew;
			//
			residuum = A*ret - b;
			double sum = 0.0;
			for (int j = 0; j < residuum.cols; j++)
				sum += powl(residuum.matrix[j][0], 2.0);
			norm = sqrtl(sum);
			if (norm <= 1e-9) {
				std::cout << std::endl << "Osiagnieto norme w Gauss-Siedlu, wynosi ona: " << iter << std::endl;
				return ret;
			}
		}
		return ret;
	}

	Matrix Matrix::GaussMethod::reduce(Matrix A, Matrix& b) {
		int rMax;
		double reducer, toReduct;
		for (int r = 0; r < A.rows; r++) {
			if (A.matrix[r][r] == 0) {
				rMax = A.findMaxInCol(r, 0, A.rows - 1);
				if (rMax == r) {
					continue;
				}
				else {
					A.switchRows(r, rMax);
					b.switchRows(r, rMax);
				}
			}
			reducer = A.matrix[r][r];
			for (int y = r + 1; y < A.rows; y++) {
				if (A.matrix[y][r] == 0) {
					continue;
				}
				toReduct = A.matrix[y][r];
				for (int x = r; x < A.cols; x++) {
					A.matrix[y][x] = A.matrix[y][x] - A.matrix[r][x] * toReduct / (reducer);
				}
				b.matrix[y][0] = b.matrix[y][0] - b.matrix[r][0] * toReduct / (reducer);
			}
		}
		return A;
	}

	Matrix Matrix::GaussMethod::solve(Matrix A, Matrix b) {
		A = reduce(A, b);
		int dim = A.rows;
		Matrix results(dim, 1);
		results.matrix[dim - 1][0] = b.matrix[dim - 1][0] / A.matrix[dim - 1][dim - 1];		//musi byc kwadratowa
		for (int y = dim - 2; y >= 0; y--) {
			results.matrix[y][0] = b.matrix[y][0];
			for (int x = dim - 1; x >= y + 1; x--) {
				results.matrix[y][0] -= (A.matrix[y][x] * results.matrix[x][0]);
			}
			results.matrix[y][0] /= A.matrix[y][y];
		}
		double norm = 0.0;
		Matrix residuum = A*results - b;
		double sum = 0.0;
		for (int j = 0; j < residuum.cols; j++)
			sum += powl(residuum.matrix[j][0], 2.0);
		std::cout << "sum: " << sum << std::endl;
		norm = sqrtl(sum);
		std::cout << "Norma w Gaussie wynosi: " << std::setprecision(30) << norm;
		return results;
	}