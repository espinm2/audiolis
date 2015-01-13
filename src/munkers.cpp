#include "munkers.h"

#include <iostream>
#include <limits.h>
using namespace std;

void step_one(int** matrix, int n, int& step) {
	int min_in_row;

	for (int r = 0; r < n; r++) {
		min_in_row = matrix[r][0];
		for (int c = 0; c < n; c++)
			if (matrix[r][c] < min_in_row)
				min_in_row = matrix[r][c];
		for (int c = 0; c < n; c++)
			matrix[r][c] -= min_in_row;
	}
	step = 2;
}

void step_two(int** matrix, int n, int* rowCover, int* colCover, int** m,
		int& step) {
	for (int r = 0; r < n; r++)
		for (int c = 0; c < n; c++) {
			if (matrix[r][c] == 0 && rowCover[r] == 0 && colCover[c] == 0) {
				m[r][c] = 1;
				rowCover[r] = 1;
				colCover[c] = 1;
			}
		}
	for (int r = 0; r < n; r++)
		rowCover[r] = 0;
	for (int c = 0; c < n; c++)
		colCover[c] = 0;
	step = 3;
}

void step_three(int n, int* colCover, int** m, int& step) {
	int colcount;
	for (int r = 0; r < n; r++)
		for (int c = 0; c < n; c++)
			if (m[r][c] == 1)
				colCover[c] = 1;

	colcount = 0;
	for (int c = 0; c < n; c++)
		if (colCover[c] == 1)
			colcount += 1;
	if (colcount >= n || colcount >= n)
		step = 7;
	else
		step = 4;
}

void find_a_zero(int** matrix, int n, int* rowCover, int* colCover, int& row,
		int& col) {
	int r = 0;
	int c;
	bool done;
	row = -1;
	col = -1;
	done = false;
	while (!done) {
		c = 0;
		while (true) {
			if (matrix[r][c] == 0 && rowCover[r] == 0 && colCover[c] == 0) {
				row = r;
				col = c;
				done = true;
			}
			c += 1;
			if (c >= n || done)
				break;
		}
		r += 1;
		if (r >= n)
			done = true;
	}
}

bool star_in_row(int** m, int n, int row) {
	bool tmp = false;
	for (int c = 0; c < n; c++)
		if (m[row][c] == 1)
			tmp = true;
	return tmp;
}

void find_star_in_row(int** m, int n, int row, int& col) {
	col = -1;
	for (int c = 0; c < n; c++)
		if (m[row][c] == 1)
			col = c;
}

void step_four(int** matrix, int n, int* rowCover, int* colCover, int** m,
		int& path_row_0, int& path_col_0, int& step) {
	int row = -1;
	int col = -1;
	bool done;

	done = false;
	while (!done) {
		find_a_zero(matrix, n, rowCover, colCover, row, col);
		if (row == -1) {
			done = true;
			step = 6;
		} else {
			m[row][col] = 2;
			if (star_in_row(m, n, row)) {
				find_star_in_row(m, n, row, col);
				rowCover[row] = 1;
				colCover[col] = 0;
			} else {
				done = true;
				step = 5;
				path_row_0 = row;
				path_col_0 = col;
			}
		}
	}
}

void find_star_in_col(int n, int** m, int c, int& r) {
	r = -1;
	for (int i = 0; i < n; i++)
		if (m[i][c] == 1)
			r = i;
}

void find_prime_in_row(int n, int** m, int r, int& c) {
	for (int j = 0; j < n; j++)
		if (m[r][j] == 2)
			c = j;
}

void augment_path(int** m, int path_count, int** path) {
	for (int p = 0; p < path_count; p++)
		if (m[path[p][0]][path[p][1]] == 1)
			m[path[p][0]][path[p][1]] = 0;
		else
			m[path[p][0]][path[p][1]] = 1;
}

void clear_covers(int n, int* rowCover, int* colCover) {
	for (int r = 0; r < n; r++)
		rowCover[r] = 0;
	for (int c = 0; c < n; c++)
		colCover[c] = 0;
}

void erase_primes(int n, int** m) {
	for (int r = 0; r < n; r++)
		for (int c = 0; c < n; c++)
			if (m[r][c] == 2)
				m[r][c] = 0;
}

void step_five(int n, int* rowCover, int* colCover, int** m, int& path_row_0,
		int& path_col_0, int& path_count, int** path, int& step) {
	bool done;
	int r = -1;
	int c = -1;

	path_count = 1;
	path[path_count - 1][0] = path_row_0;
	path[path_count - 1][1] = path_col_0;
	done = false;
	while (!done) {
		find_star_in_col(n, m, path[path_count - 1][1], r);
		if (r > -1) {
			path_count += 1;
			path[path_count - 1][0] = r;
			path[path_count - 1][1] = path[path_count - 2][1];
		} else
			done = true;
		if (!done) {
			find_prime_in_row(n, m, path[path_count - 1][0], c);
			path_count += 1;
			path[path_count - 1][0] = path[path_count - 2][0];
			path[path_count - 1][1] = c;
		}
	}
	augment_path(m, path_count, path);
	clear_covers(n, rowCover, colCover);
	erase_primes(n, m);
	step = 3;
}

void find_smallest(int** matrix, int n, int* rowCover, int* colCover,
		int& minval) {
	for (int r = 0; r < n; r++)
		for (int c = 0; c < n; c++)
			if (rowCover[r] == 0 && colCover[c] == 0)
				if (minval > matrix[r][c])
					minval = matrix[r][c];
}

void step_six(int** matrix, int n, int* rowCover, int* colCover, int& step) {
	int minval = INT_MAX;
	find_smallest(matrix, n, rowCover, colCover, minval);
	for (int r = 0; r < n; r++)
		for (int c = 0; c < n; c++) {
			if (rowCover[r] == 1)
				matrix[r][c] += minval;
			if (colCover[c] == 0)
				matrix[r][c] -= minval;
		}
	step = 4;
}

int** runMunkers(int** matrix, int n, bool max) {

	if (max == true) {
		int maxValue = matrix[0][0];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (matrix[i][j] > maxValue) {
					maxValue = matrix[i][j];
				}
			}
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				matrix[i][j] = maxValue - matrix[i][j];
			}
		}
	}

	bool done = false;
	int step = 1;
	int* rowCover = new int[n];
	int* colCover = new int[n];
	int** m = new int*[n];
	int path_row_0 = 0;
	int path_col_0 = 0;
	int path_count = 0;
	for (int i = 0; i < n; i++) {
		rowCover[i] = colCover[i] = 0;
		m[i] = new int[n];
		for (int j = 0; j < n; j++) {
			m[i][j] = 0;
		}
	}
	int** path = new int*[n * 2];
	for (int i = 0; i < n * 2; i++) {
		path[i] = new int[2];
	}
	while (!done) {
    // cout << "Matrix: " << endl;
    // for (int i = 0; i < n; i++) {
    //   for (int j = 0; j < n; j++) {
    //     cout << matrix[i][j] << " ";
    //   }
    //   cout << endl;
    // }
    // cout << "Mask: " << endl;
    // for (int i = 0; i < n; i++) {
    //   for (int j = 0; j < n; j++) {
    //     cout << m[i][j] << " ";
    //   }
    //   cout << endl;
    // }
		switch (step) {
		case 1:
			step_one(matrix, n, step);
			break;
		case 2:
			step_two(matrix, n, rowCover, colCover, m, step);
			break;
		case 3:
			step_three(n, colCover, m, step);
			break;
		case 4:
			step_four(matrix, n, rowCover, colCover, m, path_row_0, path_col_0,
					step);
			break;
		case 5:
			step_five(n, rowCover, colCover, m, path_row_0, path_col_0,
					path_count, path, step);
			break;
		case 6:
			step_six(matrix, n, rowCover, colCover, step);
			break;
		case 7:
			done = true;
			break;
		}
	}

	return m;
}

// int main() {
// 
//   const int row = 4;
//   const int col = 4;
// 
//   int ** c = new int*[row];
// 
//   for(int i = 0; i < row; i++){
//     c[i] = new int[col];
//   }
// 
//   c[0][0]= 3;
//   c[0][1]= 3;
//   c[0][2]= 3;
//   c[0][3]= 100;
// 
//   c[1][0]= 3;
//   c[1][1]= 2;
//   c[1][2]= 3;
//   c[1][3]= 100;
// 
//   c[2][0]= 3;
//   c[2][1]= 3;
//   c[2][2]= 2;
//   c[2][3]= 100;
// 
//   c[3][0]= 1;
//   c[3][1]= 1;
//   c[3][2]= 2;
//   c[3][3]= 1;
// 
// 	c = runMunkers(c, 4, false);
// 
// 	for (int i = 0; i < 4; i++) {
// 		for (int j = 0; j < 4; j++) {
// 			cout << c[i][j] << " ";
// 		}
// 		cout << endl;
// 	}
// 	return 0;
// }

