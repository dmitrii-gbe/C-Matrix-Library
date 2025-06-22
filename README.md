# s21_matrix — Matrix Library in C

[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Overview

`s21_matrix` is a lightweight C library for basic matrix operations, implemented for educational purposes. It supports creation, manipulation, and arithmetic operations on matrices with double-precision floating-point numbers.

---

## Features

- Create and free matrices
- Copy matrices
- Compare matrices for equality
- Add, subtract, and multiply matrices
- Multiply a matrix by a scalar
- Transpose matrices
- Calculate the determinant of square matrices
- Calculate the inverse of invertible square matrices
- Calculate the matrix of algebraic complements (cofactors)

---

## Functions Overview

Below is a summary of the main functions in the library:

| Function                         | Description                                                      |
|---------------------------------|------------------------------------------------------------------|
| `int s21_create_matrix(int rows, int cols, matrix_t *result)` | Creates a matrix with given dimensions. Initializes memory. Returns error code if invalid. |
| `void s21_remove_matrix(matrix_t *m)`                        | Frees matrix memory and resets fields.                          |
| `int s21_eq_matrix(const matrix_t *m1, const matrix_t *m2)`  | Compares two matrices for equality within a small epsilon.      |
| `int s21_sum_matrix(const matrix_t *m1, const matrix_t *m2, matrix_t *result)` | Adds two matrices of the same size.                             |
| `int s21_sub_matrix(const matrix_t *m1, const matrix_t *m2, matrix_t *result)` | Subtracts second matrix from first, same size required.        |
| `int s21_mult_matrix(const matrix_t *m1, const matrix_t *m2, matrix_t *result)` | Multiplies two matrices if dimensions allow.                    |
| `int s21_mult_number(const matrix_t *m, double num, matrix_t *result)` | Multiplies all elements of a matrix by a scalar.                |
| `int s21_transpose(const matrix_t *m, matrix_t *result)`      | Transposes the input matrix.                                    |
| `int s21_determinant(const matrix_t *m, double *result)`      | Calculates determinant of a square matrix.                      |
| `int s21_inverse_matrix(const matrix_t *m, matrix_t *result)` | Calculates inverse matrix if determinant is not zero.          |
| `int s21_calc_complements(const matrix_t *m, matrix_t *result)` | Calculates matrix of algebraic complements (cofactors).       |

---

## Usage Examples

```c
#include "s21_matrix.h"
#include <stdio.h>

int main() {
    matrix_t A, B, C;
    int status;

    // Create 3x3 matrix A
    status = s21_create_matrix(3, 3, &A);
    if (status != 0) return status;

    // Initialize A
    for (int i = 0; i < A.rows; ++i) {
        for (int j = 0; j < A.columns; ++j) {
            A.matrix[i][j] = i * A.columns + j + 1;
        }
    }

    // Copy A into B
    status = s21_create_matrix(3, 3, &B);
    if (status != 0) return status;

    for (int i = 0; i < B.rows; ++i) {
        for (int j = 0; j < B.columns; ++j) {
            B.matrix[i][j] = A.matrix[i][j];
        }
    }

    // Add A and B into C
    status = s21_sum_matrix(&A, &B, &C);
    if (status == 0) {
        printf("Sum matrix C:\n");
        for (int i = 0; i < C.rows; ++i) {
            for (int j = 0; j < C.columns; ++j) {
                printf("%8.2f ", C.matrix[i][j]);
            }
            printf("\n");
        }
        s21_remove_matrix(&C);
    }

    // Calculate determinant of A
    double det = 0;
    status = s21_determinant(&A, &det);
    if (status == 0) {
        printf("Determinant of A: %lf\n", det);
    }

    // Clean up
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);

    return 0;
}
```
---
# Build and test
- make all      # build library, run tests, generate coverage report
- make test     # run tests only
- make clean    # clean build artifacts
- make clang-format  # format source files
- make gcov-report     # run tests and see code coverage
---
# Contact
For questions or suggestions, feel free to open an issue or contact the maintainer.
---
# License

This project is licensed under the MIT License.

---
