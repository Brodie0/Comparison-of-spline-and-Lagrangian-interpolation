#pragma once

class Matrix {
public:
	double ** matrix;
	int rows;
	int cols;

	Matrix(int, int);
	Matrix(const Matrix&);
	~Matrix();

	void print();
	void switchRows(int, int);
	int findMaxInCol(int, int, int);
	Matrix& operator=(const Matrix&);
	Matrix operator+(const Matrix&);
	Matrix operator-(const Matrix&);
	Matrix operator*(const Matrix&);

	class JacobiMethod {
	public:
		static Matrix solve(Matrix, Matrix, int);
	};

	class GaussSiedelMethod {
	public:
		static Matrix solve(Matrix, Matrix, int);
	};

	class GaussMethod {
	public:
		static Matrix reduce(Matrix, Matrix&);
		static Matrix solve(Matrix, Matrix);
	};
};