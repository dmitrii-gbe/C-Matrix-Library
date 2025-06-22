#pragma once
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define EXIT_CODE_SUCCESS 0
#define EXIT_CODE_INCORRECT_MATRIX 1
#define EXIT_CODE_CALCULATION_ERROR 2
#define EPSILON 1e-7
#define COMP_RESULT_FAILURE 0
#define COMP_RESULT_SUCCESS 1

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);
void print_row(FILE *f, const double *const row, int columns);
void dump_matrix(const char *const message, FILE *f, const matrix_t *const m);
void read_matrix(matrix_t *const m, FILE *f, int rows, int cols);
int find_proper_row(const matrix_t *const m, const int col);
int copy_matrix(const matrix_t *const from, matrix_t *const to);
void swap_two_rows(matrix_t *const m, const int row1, const int row2);
int make_row_echelon_form(matrix_t *A, int *const swap_count);
void make_identity_matrix(const int rows, const int cols,
                          matrix_t *const result);
void make_identity_diag(matrix_t *const m, matrix_t *const result);
void multiply_row_by_number(double *const row, const double number,
                            const int columns);
void forward_elimination(matrix_t *const m, matrix_t *const result);
void reverse_elimination(matrix_t *const m, matrix_t *const result);
void fill_minor(const matrix_t *const m, const int row_to_avoid,
                const int col_to_avoid, matrix_t *const minor);
int s21_eq_internal(const matrix_t *const m1, const matrix_t *const m2);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
double abs_double(const double val);
double determinant_by_definition(const matrix_t *const m);
void make_row_echelon_form_two_matrices(matrix_t *const m,
                                        matrix_t *const result);
int s21_eq_matrix_relative(const matrix_t *const A, const matrix_t *const B);
int s21_eq_relative_internal(const matrix_t *const m1,
                             const matrix_t *const m2);
void add_num_to_matrix(matrix_t *const m, const double num);
matrix_t create_random_matrix(const int rows, const int columns);
void make_random_negative(double *const val);
