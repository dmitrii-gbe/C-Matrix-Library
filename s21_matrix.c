#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int status = rows > 0 && columns > 0 && result ? EXIT_CODE_SUCCESS
                                                 : EXIT_CODE_INCORRECT_MATRIX;
  if (status == EXIT_CODE_SUCCESS) {
    result->columns = columns;
    result->rows = rows;
    result->matrix = (double **)malloc(sizeof(double *) * rows +
                                       sizeof(double) * columns * rows);
    status = result->matrix
                 ? status
                 : EXIT_CODE_INCORRECT_MATRIX;  // if malloc returns null
    if (status == EXIT_CODE_SUCCESS) {
      for (int i = 0; i < rows; ++i) {
        result->matrix[i] = (double *)(result->matrix + rows) + columns * i;
      }
      for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
          result->matrix[i][j] = 0.0;
        }
      }
    }
  }
  return status;
}

void s21_remove_matrix(matrix_t *A) {
  free(A->matrix);
  A->matrix = NULL;
  A->rows = 0;
  A->columns = 0;
}

void read_matrix(matrix_t *const m, FILE *f, const int rows, const int cols) {
  s21_create_matrix(rows, cols, m);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      fscanf(f, "%lf", &m->matrix[i][j]);
    }
  }
}

void add_rows(double *to, const double *from, const int cols,
              const double multiplyer) {
  for (int i = 0; i < cols; ++i) {
    *to++ += multiplyer * *from++;
  }
}

int copy_matrix(const matrix_t *const from, matrix_t *const to) {
  int status = s21_create_matrix(from->rows, from->columns, to);
  if (!status) {
    for (int i = 0; i < from->rows; ++i) {
      for (int j = 0; j < from->columns; ++j) {
        to->matrix[i][j] = from->matrix[i][j];
      }
    }
  }
  return status;
}

int find_proper_row(const matrix_t *const m, const int col) {
  int result = -1;
  double max = 0.0;
  for (int row = col; row < m->rows; ++row) {
    double val = abs_double(m->matrix[row][col]);
    if (val > max) {
      max = val;
      result = row;
    }
  }
  return max == 0.0 ? -1 : result;
}

void swap_two_rows(matrix_t *const m, const int row1, const int row2) {
  double *tmp = m->matrix[row1];
  m->matrix[row1] = m->matrix[row2];
  m->matrix[row2] = tmp;
}

int make_row_echelon_form(matrix_t *const A, int *const swap_count) {
  int error = 0;
  for (int col = 0; col < A->columns && !error; ++col) {
    int idx;
    if ((idx = find_proper_row(A, col)) != -1) {
      *swap_count += idx != col;
      swap_two_rows(A, col, idx);
    }

    if (A->matrix[col][col] != 0.0) {
      for (int j = col + 1; j < A->rows; ++j) {  // j - row
        add_rows(A->matrix[j], A->matrix[col], A->columns,
                 -A->matrix[j][col] / A->matrix[col][col]);
      }
    } else {
      error = 1;
    }
  }
  return error;
}

double accumulate_determinant(const matrix_t *const m) {
  double result = 1.0;
  for (int i = 0; i < m->columns; ++i) {  // i - column
    result *= m->matrix[i][i];
  }
  return result;
}

int s21_determinant(matrix_t *A, double *result) {
  int status = result && A && A->columns == A->rows && A->columns > 0
                   ? EXIT_CODE_SUCCESS
                   : EXIT_CODE_INCORRECT_MATRIX;
  if (status == EXIT_CODE_SUCCESS) {
    matrix_t m;
    copy_matrix(A, &m);
    int swap_count = 0;
    if (!make_row_echelon_form(&m, &swap_count)) {
      *result = accumulate_determinant(&m);
      *result = swap_count % 2 == 0 ? *result : -(*result);
    } else {
      *result = 0.0;
    }
    s21_remove_matrix(&m);
  }
  return status;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int status = (A && A->matrix && A->columns > 0 && A->rows > 0 && result)
                   ? EXIT_CODE_SUCCESS
                   : EXIT_CODE_INCORRECT_MATRIX;
  if (status == EXIT_CODE_SUCCESS) {
    status = s21_create_matrix(A->rows, A->columns, result);
    for (int row = 0; row < A->rows && !status; ++row) {
      for (int col = 0; col < A->columns; ++col) {
        result->matrix[row][col] = A->matrix[row][col] * number;
      }
    }
  }
  return status;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = (A && B && A->columns > 0 && A->rows > 0 && A->matrix &&
                B->matrix && result)
                   ? EXIT_CODE_SUCCESS
                   : EXIT_CODE_INCORRECT_MATRIX;
  if (!status && (A->rows != B->rows || A->columns != B->columns))
    status = EXIT_CODE_CALCULATION_ERROR;
  if (status == EXIT_CODE_SUCCESS) {
    status = s21_create_matrix(A->rows, A->columns, result);
    for (int row = 0; row < A->rows && !status; ++row) {
      for (int col = 0; col < A->columns; ++col) {
        result->matrix[row][col] = A->matrix[row][col] + B->matrix[row][col];
      }
    }
  }
  return status;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = (A && B && A->columns > 0 && A->rows > 0 && A->matrix &&
                B->matrix && result)
                   ? EXIT_CODE_SUCCESS
                   : EXIT_CODE_INCORRECT_MATRIX;
  if (!status && (A->rows != B->rows || A->columns != B->columns))
    status = EXIT_CODE_CALCULATION_ERROR;
  if (status == EXIT_CODE_SUCCESS) {
    status = s21_create_matrix(A->rows, A->columns, result);
    for (int row = 0; row < A->rows && !status; ++row) {
      for (int col = 0; col < A->columns; ++col) {
        result->matrix[row][col] = A->matrix[row][col] - B->matrix[row][col];
      }
    }
  }
  return status;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = (A && B && A->columns > 0 && A->rows > 0 && A->matrix &&
                B->matrix && result)
                   ? EXIT_CODE_SUCCESS
                   : EXIT_CODE_INCORRECT_MATRIX;
  if (!status && (A->columns != B->rows)) status = EXIT_CODE_CALCULATION_ERROR;
  if (status == EXIT_CODE_SUCCESS) {
    status = s21_create_matrix(A->rows, B->columns, result);
    for (int i = 0; i < A->rows && !status; ++i) {
      for (int j = 0; j < B->columns; ++j) {
        for (int k = 0; k < A->columns; ++k) {
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      }
    }
  }
  return status;
}

void make_identity_matrix(const int rows, const int cols,
                          matrix_t *const result) {
  s21_create_matrix(rows, cols, result);
  for (int i = 0; i < cols; ++i) {
    result->matrix[i][i] = 1.0;
  }
}

void multiply_row_by_number(double *const row, const double number,
                            const int columns) {
  for (int i = 0; i < columns; ++i) {
    row[i] *= number;
  }
}

void forward_elimination(matrix_t *const m, matrix_t *const result) {
  for (int col = 0; col < m->columns; ++col) {
    int idx;
    if ((idx = find_proper_row(m, col)) != -1) {
      swap_two_rows(m, col, idx);
      swap_two_rows(result, col, idx);
    }
    double multiplier = 1.0 / m->matrix[col][col];
    multiply_row_by_number(m->matrix[col] + col, multiplier, m->columns - col);
    multiply_row_by_number(result->matrix[col], multiplier, result->columns);
    for (int row = col + 1; row < m->rows; ++row) {
      multiplier = m->matrix[row][col];
      add_rows(m->matrix[row] + col, m->matrix[col] + col, m->columns - col,
               -multiplier);
      add_rows(result->matrix[row], result->matrix[col], result->columns,
               -multiplier);
    }
  }
}

void reverse_elimination(matrix_t *const m, matrix_t *const result) {
  for (int col = m->columns - 1; col > 0; --col) {
    for (int row = col - 1; row >= 0; --row) {
      double multiplier = m->matrix[row][col];
      add_rows(m->matrix[row] + col, m->matrix[col] + col, m->columns - col,
               -multiplier);
      add_rows(result->matrix[row], result->matrix[col], result->columns,
               -multiplier);
    }
  }
}

void make_row_echelon_form_two_matrices(matrix_t *const m,
                                        matrix_t *const result) {
  forward_elimination(m, result);
  reverse_elimination(m, result);
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int status = (A && A->columns > 0 && A->rows > 0 && A->matrix && result)
                   ? EXIT_CODE_SUCCESS
                   : EXIT_CODE_INCORRECT_MATRIX;
  double det = 0.0;
  if (!status && (s21_determinant(A, &det) || det == 0.0))
    status = EXIT_CODE_CALCULATION_ERROR;
  if (status == EXIT_CODE_SUCCESS) {
    make_identity_matrix(A->rows, A->columns, result);
    matrix_t m;
    status = copy_matrix(A, &m);
    make_row_echelon_form_two_matrices(&m, result);
    s21_remove_matrix(&m);
  }
  return status;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int status = (A && A->columns > 0 && A->rows > 0 && A->matrix && result)
                   ? EXIT_CODE_SUCCESS
                   : EXIT_CODE_INCORRECT_MATRIX;
  if (status == EXIT_CODE_SUCCESS) {
    s21_create_matrix(A->columns, A->rows, result);
    for (int row = 0; row < A->rows; ++row) {
      for (int col = 0; col < A->columns; ++col) {
        result->matrix[col][row] = A->matrix[row][col];
      }
    }
  }
  return status;
}

void fill_minor(const matrix_t *const m, const int row_to_avoid,
                const int col_to_avoid, matrix_t *const minor) {
  int minor_row = 0, minor_col = 0;
  for (int row = 0; row < m->rows; ++row) {
    for (int col = 0; col < m->columns; ++col) {
      if (row != row_to_avoid && col != col_to_avoid) {
        minor->matrix[minor_row][minor_col++] = m->matrix[row][col];
      }
    }
    minor_col = 0;
    minor_row += row != row_to_avoid;
  }
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int status = (A && A->columns > 0 && A->rows > 0 && A->matrix && result)
                   ? EXIT_CODE_SUCCESS
                   : EXIT_CODE_INCORRECT_MATRIX;
  if (!status && (A->rows != A->columns || A->rows < 2))
    status = EXIT_CODE_CALCULATION_ERROR;
  if (status == EXIT_CODE_SUCCESS) {
    s21_create_matrix(A->rows, A->columns, result);
    matrix_t minor;
    s21_create_matrix(A->rows - 1, A->columns - 1, &minor);
    for (int i = 0; i < A->rows; ++i) {
      for (int j = 0; j < A->columns; ++j) {
        double det = 0.0;
        fill_minor(A, i, j, &minor);
        s21_determinant(&minor, &det);
        result->matrix[i][j] = (i + j) % 2 == 0 ? det : -det;
      }
    }
    s21_remove_matrix(&minor);
  }
  return status;
}

double abs_double(const double val) { return val < 0.0 ? -val : val; }

int s21_eq_internal(const matrix_t *const m1, const matrix_t *const m2) {
  int result = 1;
  for (int row = 0; row < m1->rows && result; ++row) {
    for (int col = 0; col < m1->columns && result; ++col) {
      if (abs_double(m1->matrix[row][col] - m2->matrix[row][col]) >= EPSILON)
        result = 0;
    }
  }
  return result;
}

int s21_eq_relative_internal(const matrix_t *const m1,
                             const matrix_t *const m2) {
  int result = 1;
  for (int row = 0; row < m1->rows && result; ++row) {
    for (int col = 0; col < m1->columns && result; ++col) {
      if ((m2->matrix[row][col] != 0.0 &&
           abs_double(1.0 - m1->matrix[row][col] / m2->matrix[row][col]) >=
               EPSILON) ||
          (m2->matrix[row][col] == 0.0 &&
           abs_double(m1->matrix[row][col]) >= EPSILON))
        result = 0;
    }
  }
  return result;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  return A && B && A->rows == B->rows && A->columns == B->columns &&
         A->matrix && B->matrix && s21_eq_internal(A, B);
}

int s21_eq_matrix_relative(const matrix_t *const A, const matrix_t *const B) {
  return A && B && A->rows == B->rows && A->columns == B->columns &&
         s21_eq_relative_internal(A, B);
}

double determinant_by_definition(const matrix_t *const m) {
  if (m->rows == 2 && m->columns == 2) {
    return m->matrix[0][0] * m->matrix[1][1] -
           m->matrix[0][1] * m->matrix[1][0];
  }
  if (m->rows == 1 && m->columns == 1) {
    return m->matrix[0][0];
  }
  double determinant = 0.0;
  double sign = 1.0;
  for (int col = 0; col < m->columns; ++col) {
    matrix_t minor;
    s21_create_matrix(m->rows - 1, m->columns - 1, &minor);
    fill_minor(m, 0, col, &minor);
    determinant += sign * m->matrix[0][col] * determinant_by_definition(&minor);
    s21_remove_matrix(&minor);
    sign = -sign;
  }
  return determinant;
}

void print_row(FILE *f, const double *const row, const int columns) {
  char first = 1;
  for (int i = 0; i < columns; ++i) {
    if (first) {
      fprintf(f, "%19.16f", row[i]);
      first = 0;
    } else {
      fprintf(f, " %19.16f", row[i]);
    }
  }
  fprintf(f, "\n");
}

void dump_matrix(const char *const message, FILE *f, const matrix_t *const m) {
  fprintf(f, "%s\n", message);
  for (int i = 0; i < m->rows; ++i) {
    print_row(f, m->matrix[i], m->columns);
  }
  fprintf(f, "\n");
}

void add_num_to_matrix(matrix_t *const m, const double num) {
  for (int row = 0; row < m->rows; ++row) {
    for (int col = 0; col < m->columns; ++col) {
      m->matrix[row][col] += num;
    }
  }
}

void make_random_negative(double *const val) {
  int num = rand();
  if (num % 3 == 0) {
    *val = -(*val);
  }
  if (num % 4 == 0) {
    *val = 0.0;
  }
}

matrix_t create_random_matrix(const int rows, const int columns) {
  matrix_t result;
  s21_create_matrix(rows, columns, &result);
  for (int row = 0; row < rows; ++row) {
    for (int col = 0; col < columns; ++col) {
      double a = rand() % 10000;
      double b;
      while ((b = rand() % 100) == 0.0) {
      }
      make_random_negative(&a);
      result.matrix[row][col] = a / b;
    }
  }
  return result;
}
