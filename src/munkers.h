#ifndef MUNKERS_H
#define MUNKERS_H

void step_one(int** matrix, int n, int& step);
void step_two(int** matrix, int n, int* rowCover, int* colCover, int** m, int& step);
void step_three(int n, int* colCover, int** m, int& step);
void find_a_zero(int** matrix, int n, int* rowCover, int* colCover, int& row, int& col);
bool star_in_row(int** m, int n, int row);
void find_star_in_row(int** m, int n, int row, int& col);
void step_four(int** matrix, int n, int* rowCover, int* colCover, int** m, int& path_row_0, int& path_col_0, int& step);
void find_star_in_col(int n, int** m, int c, int& r);
void find_prime_in_row(int n, int** m, int r, int& c);
void augment_path(int** m, int path_count, int** path);
void clear_covers(int n, int* rowCover, int* colCover);
void erase_primes(int n, int** m);
void step_five(int n, int* rowCover, int* colCover, int** m, int& path_row_0, int& path_col_0, int& path_count, int** path, int& step);
void find_smallest(int** matrix, int n, int* rowCover, int* colCover, int& minval);
void step_six(int** matrix, int n, int* rowCover, int* colCover, int& step);
int** runMunkers(int** matrix, int n, bool max);

#endif


