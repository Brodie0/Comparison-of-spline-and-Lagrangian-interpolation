// Projekt3MN.cpp: Okreœla punkt wejœcia dla aplikacji konsoli.
//
#include <fstream>
#include <iomanip>
#include "Matrix.h"

struct Factors {
	double a;
	double b;
	double c;
	double d;
};

double calculateS(Factors *f, Matrix nodes, int j, double x) {
	double xj = nodes.matrix[0][j];
	return f[j].a + f[j].b*(x - xj) + f[j].c*(x - xj)*(x - xj) + f[j].d*(x - xj)*(x - xj)*(x - xj);
}

double h(Matrix nodes, int j) {
	return nodes.matrix[0][j] - nodes.matrix[0][j - 1];
}

double mi(Matrix nodes, int j) {
	if (j == 0 || j == nodes.cols - 1)
	{
		return 0.0;
	}
	return h(nodes, j) / (h(nodes, j) + h(nodes, j + 1));
}

double lambda(Matrix nodes, int j) {
	if (j == 0 || j == nodes.cols - 1)
	{
		return 0.0;
	}
	return h(nodes, j + 1) / (h(nodes, j) + h(nodes, j + 1));
}

double delta(Matrix nodes, Matrix values, int j) {
	if (j == 0 || j == nodes.cols - 1)
	{
		return 0.0;
	}
	return (6.0 / (h(nodes, j) + h(nodes, j + 1)))
		*(((values.matrix[0][j + 1] - values.matrix[0][j]) / h(nodes, j + 1))
			- ((values.matrix[0][j] - values.matrix[0][j - 1]) / h(nodes, j)));
}

Matrix createMmatrix(Matrix nodes) {
	Matrix ret(nodes.cols, nodes.cols);
	ret.matrix[0][0] = 2.0;
	ret.matrix[0][1] = lambda(nodes, 0);
	for (int y = 1; y < nodes.cols - 1; y++) {
		ret.matrix[y][y - 1] = mi(nodes, y);
		ret.matrix[y][y] = 2.0;
		ret.matrix[y][y + 1] = lambda(nodes, y);
	}
	ret.matrix[nodes.cols - 1][nodes.cols - 1] = 2.0;
	ret.matrix[nodes.cols - 1][nodes.cols - 2] = mi(nodes, nodes.cols - 1);
	return ret;
}

double li(Matrix X, double xS, int i) {
	double li = 1.0;
	for (unsigned int j = 0; j < X.cols; j++) {
		if (j != i)
			li *= (xS - X.matrix[0][j]) / (X.matrix[0][i] - X.matrix[0][j]);
	}
	return li;
}

double lagrangeInterpolation(Matrix X, Matrix Y, double xS) {

	double yS = 0.0;
	for (unsigned int i = 0; i < X.cols; i++) {
		yS += Y.matrix[0][i] * li(X, xS, i);
	}
	return yS;
}

void loadFile(char* filename, Matrix &x, Matrix &fx) {
	std::fstream file;
	file.open(filename, std::ios::in);

	for (int r = 0; r < x.rows; r++) {
		if (file.eof()) {
			break;
		}
		for (int c = 0; c < x.cols; c++) {
			char d = 0;
			file >> x.matrix[r][c];
			file >> d;	//przecinek
			file >> fx.matrix[r][c];
		}
	}
	file.close();
}

int main() {
	char *path = "someNodes\\middleSahara20evenly.txt";
	//char *path = "full\\middleSahara.txt";
	const int nodesCount = 20;
	Matrix nodes(1, nodesCount), functionValues(1, nodesCount);
	loadFile(path, nodes, functionValues);
	Matrix left = createMmatrix(nodes);
	Matrix right(nodesCount, 1);
	for (int y = 0; y < right.rows; y++) {
		right.matrix[y][0] = delta(nodes, functionValues, y);
	}
	Matrix M = Matrix::GaussMethod::solve(left, right);
	Factors factors[nodesCount];

	for (int i = 0; i < nodesCount - 1; i++) {
		factors[i].a = functionValues.matrix[0][i];
		factors[i].b = ((functionValues.matrix[0][i + 1]) - (functionValues.matrix[0][i])) / (h(nodes, i + 1))
			- (((2 * M.matrix[i][0]) + (M.matrix[i + 1][0])) / 6.0) * h(nodes, i + 1);
		factors[i].c = M.matrix[i][0] / 2;
		factors[i].d = (M.matrix[i + 1][0] - M.matrix[i][0]) / (6.0 * h(nodes, i + 1));
	}

	const int count = 512;
	double x[count];
	double lagrangeResults[count];
	double splineResults[count];
	double start = nodes.matrix[0][0];
	double end = nodes.matrix[0][nodesCount-1];
	double step = (end - start) / (double)count;
	for (int i = 0; i < count; i++) {
		x[i] = start + i*step;
	}

	for (int i = 0; i < count; i++) {
		lagrangeResults[i] = lagrangeInterpolation(nodes, functionValues, x[i]);
	}

	for (int i = 0; i < count; i++) {
		for (int j = 1; j < nodesCount; j++) {
			if (x[i] < nodes.matrix[0][j]) {
				splineResults[i] = calculateS(factors, nodes, j - 1, x[i]);
				break;
			}
		}
	}

	std::fstream file;
	file.open("splineInterp.txt", std::ios::out);
	for (int i = 0; i < count; i++) {
		file << splineResults[i] << "\n";
	}
	//for (int i = 0; i < count; i++) {
	//	file << std::setprecision(10)<<nodes.matrix[0][i] << "\n";
	//}
	file.close();

	file.open("lagrangeInterp.txt", std::ios::out);
	for (int i = 0; i < count; i++) {
		file << lagrangeResults[i] << "\n";
	}
	//for (int i = 0; i < count; i++) {
	//	file << std::setprecision(10) << functionValues.matrix[0][i] << "\n";
	//}
	file.close();

	return 0;
}
