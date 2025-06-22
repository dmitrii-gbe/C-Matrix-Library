#pragma once
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_double.h>
#include <time.h>

#include "../s21_matrix.h"

gsl_matrix *from_s21_matrix_to_gsl(const matrix_t *const s21_m);
matrix_t from_gsl_to_s21_matrix(const gsl_matrix *const gsl_m);
void print_gsl_row(const gsl_matrix *const m, const int row);
void dump_gsl_matrix(const char *const message, const gsl_matrix *const m);
double gsl_determinant(const gsl_matrix *m);
int gsl_invert_matrix(const gsl_matrix *A, gsl_matrix *inverse);
void compute_cofactor_matrix(const gsl_matrix *A, gsl_matrix *cofactor);
gsl_matrix *gsl_transpose(const gsl_matrix *m);
