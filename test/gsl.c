#include "gsl.h"

gsl_matrix *from_s21_matrix_to_gsl(const matrix_t *const s21_m) {
  gsl_matrix *gsl_m = gsl_matrix_calloc(s21_m->rows, s21_m->columns);
  for (int row = 0; row < s21_m->rows; ++row) {
    for (int col = 0; col < s21_m->columns; ++col) {
      gsl_matrix_set(gsl_m, row, col, s21_m->matrix[row][col]);
    }
  }
  return gsl_m;
}

matrix_t from_gsl_to_s21_matrix(const gsl_matrix *const gsl_m) {
  matrix_t s21_m;
  s21_create_matrix(gsl_m->size1, gsl_m->size2, &s21_m);
  for (int row = 0; row < gsl_m->size1; ++row) {
    for (int col = 0; col < gsl_m->size2; ++col) {
      s21_m.matrix[row][col] = gsl_matrix_get(gsl_m, row, col);
    }
  }
  return s21_m;
}

void print_gsl_row(const gsl_matrix *const m, const int row) {
  char first = 1;
  for (int col = 0; col < m->size2; ++col) {
    if (first) {
      printf("%8.2f", gsl_matrix_get(m, row, col));
      first = 0;
    } else {
      printf(" %8.2f", gsl_matrix_get(m, row, col));
    }
  }
  printf("\n");
}

void dump_gsl_matrix(const char *const message, const gsl_matrix *const m) {
  printf("%s\n", message);
  for (int row = 0; row < m->size1; ++row) {
    print_gsl_row(m, row);
  }
  printf("\n");
}

double gsl_determinant(const gsl_matrix *m) {
  size_t n = m->size1;
  if (n != m->size2) {
    return 0.0;
  }

  gsl_matrix *tmp = gsl_matrix_alloc(n, n);
  gsl_matrix_memcpy(tmp, m);

  gsl_permutation *p = gsl_permutation_alloc(n);
  int signum;

  gsl_linalg_LU_decomp(tmp, p, &signum);
  double det = gsl_linalg_LU_det(tmp, signum);

  gsl_permutation_free(p);
  gsl_matrix_free(tmp);

  return det;
}

int gsl_invert_matrix(const gsl_matrix *A, gsl_matrix *inverse) {
  int n = A->size1;

  if (A->size1 != A->size2 || inverse->size1 != n || inverse->size2 != n) {
    return -1;
  }

  gsl_matrix *A_copy = gsl_matrix_alloc(n, n);
  gsl_matrix_memcpy(A_copy, A);

  gsl_permutation *perm = gsl_permutation_alloc(n);
  int signum;

  int status = gsl_linalg_LU_decomp(A_copy, perm, &signum);
  if (status != 0) {
    gsl_permutation_free(perm);
    gsl_matrix_free(A_copy);
    return status;
  }

  status = gsl_linalg_LU_invert(A_copy, perm, inverse);

  gsl_permutation_free(perm);
  gsl_matrix_free(A_copy);
  return status;
}

void get_minor(const gsl_matrix *src, gsl_matrix *minor, size_t row,
               size_t col) {
  size_t m_row = 0, m_col = 0;
  for (size_t i = 0; i < src->size1; ++i) {
    if (i == row) continue;
    m_col = 0;
    for (size_t j = 0; j < src->size2; ++j) {
      if (j == col) continue;
      gsl_matrix_set(minor, m_row, m_col, gsl_matrix_get(src, i, j));
      m_col++;
    }
    m_row++;
  }
}

void compute_cofactor_matrix(const gsl_matrix *A, gsl_matrix *cofactor) {
  size_t n = A->size1;
  gsl_matrix *minor = gsl_matrix_alloc(n - 1, n - 1);

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      get_minor(A, minor, i, j);

      gsl_permutation *p = gsl_permutation_alloc(n - 1);
      int signum;
      gsl_matrix *tmp = gsl_matrix_alloc(n - 1, n - 1);
      gsl_matrix_memcpy(tmp, minor);

      gsl_linalg_LU_decomp(tmp, p, &signum);
      double det = gsl_linalg_LU_det(tmp, signum);

      double cofactor_ij = ((i + j) % 2 == 0 ? 1 : -1) * det;
      gsl_matrix_set(cofactor, i, j, cofactor_ij);

      gsl_permutation_free(p);
      gsl_matrix_free(tmp);
    }
  }

  gsl_matrix_free(minor);
}

gsl_matrix *gsl_transpose(const gsl_matrix *m) {
  size_t rows = m->size1;
  size_t cols = m->size2;
  gsl_matrix *result = gsl_matrix_alloc(cols, rows);
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      double val = gsl_matrix_get(m, i, j);
      gsl_matrix_set(result, j, i, val);
    }
  }
  return result;
}
