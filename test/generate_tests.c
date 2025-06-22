#include <gsl/gsl_matrix_double.h>
#include <stdio.h>

#include "../s21_matrix.h"
#include "gsl.h"

void gen_s21_inverse_tests(int test_count) {
  FILE* f = fopen("test_files/s21_inverse_test", "w");
  if (!f) {
    fprintf(stderr, "Output file is not found\n");
    return;
  }
  fprintf(f, "%d\n", test_count);  // test count
  while (test_count != 0) {
    matrix_t m;
    gsl_matrix* gsl_m;
    int dimention = 1 + rand() % 30;
    m = create_random_matrix(dimention, dimention);
    gsl_m = from_s21_matrix_to_gsl(&m);
    gsl_matrix* gsl_inverse = gsl_matrix_calloc(m.rows, m.columns);
    if (gsl_invert_matrix(gsl_m, gsl_inverse) == 0) {
      fprintf(f, "%d\n", dimention);  // dimention
      dump_matrix("", f, &m);
      matrix_t inverted_conv_matrix_t = from_gsl_to_s21_matrix(gsl_inverse);
      dump_matrix("", f, &inverted_conv_matrix_t);
      s21_remove_matrix(&inverted_conv_matrix_t);
      --test_count;
    }
    gsl_matrix_free(gsl_inverse);
    s21_remove_matrix(&m);
    gsl_matrix_free(gsl_m);
  }
  fclose(f);
}

void gen_s21_determinant_tests(int test_count) {
  FILE* f = fopen("test_files/s21_determinant_test", "w");
  if (!f) {
    fprintf(stderr, "Output file is not found\n");
    return;
  }
  fprintf(f, "%d\n", test_count);  // test count
  matrix_t m;
  gsl_matrix* gsl_m;
  while (test_count != 0) {
    int dimention = 1 + rand() % 30;
    m = create_random_matrix(dimention, dimention);
    gsl_m = from_s21_matrix_to_gsl(&m);
    fprintf(f, "%d\n", dimention);  // dimention
    dump_matrix("", f, &m);
    double det = gsl_determinant(gsl_m);
    fprintf(f, "%.16f\n", det);
    s21_remove_matrix(&m);
    gsl_matrix_free(gsl_m);
    --test_count;
  }
  fclose(f);
}

void gen_s21_sum_tests(int test_count) {
  FILE* f = fopen("test_files/s21_sum_test", "w");
  if (!f) {
    fprintf(stderr, "Output file is not found\n");
    return;
  }
  fprintf(f, "%d\n", test_count);  // test count
  matrix_t m, m1;
  gsl_matrix *gsl_m, *gsl_m1;
  while (test_count != 0) {
    int row = 1 + rand() % 30;
    int col = 1 + rand() % 30;
    m = create_random_matrix(row, col);
    m1 = create_random_matrix(row, col);
    gsl_m = from_s21_matrix_to_gsl(&m);
    gsl_m1 = from_s21_matrix_to_gsl(&m1);
    fprintf(f, "%d %d\n", row, col);
    dump_matrix("", f, &m);
    dump_matrix("", f, &m1);
    gsl_matrix_add(gsl_m, gsl_m1);
    matrix_t result = from_gsl_to_s21_matrix(gsl_m);
    dump_matrix("", f, &result);
    s21_remove_matrix(&m);
    s21_remove_matrix(&m1);
    s21_remove_matrix(&result);
    gsl_matrix_free(gsl_m);
    gsl_matrix_free(gsl_m1);
    --test_count;
  }
  fclose(f);
}

void gen_s21_sub_tests(int test_count) {
  FILE* f = fopen("test_files/s21_sub_test", "w");
  if (!f) {
    fprintf(stderr, "Output file is not found\n");
    return;
  }
  fprintf(f, "%d\n", test_count);  // test count
  matrix_t m, m1;
  gsl_matrix *gsl_m, *gsl_m1;
  while (test_count != 0) {
    int row = 1 + rand() % 30;
    int col = 1 + rand() % 30;
    m = create_random_matrix(row, col);
    m1 = create_random_matrix(row, col);
    gsl_m = from_s21_matrix_to_gsl(&m);
    gsl_m1 = from_s21_matrix_to_gsl(&m1);
    fprintf(f, "%d %d\n", row, col);
    dump_matrix("", f, &m);
    dump_matrix("", f, &m1);
    gsl_matrix_sub(gsl_m, gsl_m1);
    matrix_t result = from_gsl_to_s21_matrix(gsl_m);
    dump_matrix("", f, &result);
    s21_remove_matrix(&m);
    s21_remove_matrix(&m1);
    s21_remove_matrix(&result);
    gsl_matrix_free(gsl_m);
    gsl_matrix_free(gsl_m1);
    --test_count;
  }
  fclose(f);
}

void gen_s21_mul_tests(int test_count) {
  FILE* f = fopen("test_files/s21_mul_test", "w");
  if (!f) {
    fprintf(stderr, "Output file is not found\n");
    return;
  }
  fprintf(f, "%d\n", test_count);  // test count
  matrix_t m, m1;
  gsl_matrix *gsl_m, *gsl_m1, *gsl_res;
  while (test_count != 0) {
    int row1 = 1 + rand() % 30;
    int col = 1 + rand() % 30;
    int col2 = 1 + rand() % 30;
    m = create_random_matrix(row1, col);
    m1 = create_random_matrix(col, col2);
    gsl_m = from_s21_matrix_to_gsl(&m);
    gsl_m1 = from_s21_matrix_to_gsl(&m1);
    gsl_res = gsl_matrix_calloc(row1, col2);
    fprintf(f, "%d %d %d\n", row1, col, col2);
    dump_matrix("", f, &m);
    dump_matrix("", f, &m1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, gsl_m, gsl_m1, 0.0,
                   gsl_res);
    matrix_t result = from_gsl_to_s21_matrix(gsl_res);
    dump_matrix("", f, &result);
    s21_remove_matrix(&m);
    s21_remove_matrix(&m1);
    s21_remove_matrix(&result);
    gsl_matrix_free(gsl_m);
    gsl_matrix_free(gsl_m1);
    gsl_matrix_free(gsl_res);
    --test_count;
  }
  fclose(f);
}

void gen_s21_scale_tests(int test_count) {
  FILE* f = fopen("test_files/s21_scale_test", "w");
  if (!f) {
    fprintf(stderr, "Output file is not found\n");
    return;
  }
  fprintf(f, "%d\n", test_count);  // test count
  matrix_t m;
  gsl_matrix* gsl_m;
  while (test_count != 0) {
    int row = 1 + rand() % 30;
    int col = 1 + rand() % 30;
    double multiplyer = ((double)rand() / (double)rand());
    m = create_random_matrix(row, col);
    gsl_m = from_s21_matrix_to_gsl(&m);
    gsl_matrix_scale(gsl_m, multiplyer);
    fprintf(f, "%d %d %.20lf\n", row, col, multiplyer);
    dump_matrix("", f, &m);
    matrix_t result = from_gsl_to_s21_matrix(gsl_m);
    dump_matrix("", f, &result);
    s21_remove_matrix(&m);
    s21_remove_matrix(&result);
    gsl_matrix_free(gsl_m);
    --test_count;
  }
  fclose(f);
}

void gen_s21_complements_tests(int test_count) {
  FILE* f = fopen("test_files/s21_complements_test", "w");
  if (!f) {
    fprintf(stderr, "Output file is not found\n");
    return;
  }
  fprintf(f, "%d\n", test_count);
  matrix_t m;
  gsl_matrix *gsl_m, *gsl_result;
  while (test_count != 0) {
    int dimention = 2 + rand() % 30;
    m = create_random_matrix(dimention, dimention);
    gsl_m = from_s21_matrix_to_gsl(&m);
    fprintf(f, "%d\n", dimention);
    dump_matrix("", f, &m);
    gsl_result = gsl_matrix_calloc(dimention, dimention);
    compute_cofactor_matrix(gsl_m, gsl_result);
    matrix_t result = from_gsl_to_s21_matrix(gsl_result);
    dump_matrix("", f, &result);
    s21_remove_matrix(&m);
    s21_remove_matrix(&result);
    gsl_matrix_free(gsl_m);
    gsl_matrix_free(gsl_result);
    --test_count;
  }
  fclose(f);
}

void gen_s21_transpose_tests(int test_count) {
  FILE* f = fopen("test_files/s21_transpose_test", "w");
  if (!f) {
    fprintf(stderr, "Output file is not found\n");
    return;
  }
  fprintf(f, "%d\n", test_count);
  matrix_t m;
  gsl_matrix *gsl_m, *gsl_result;
  while (test_count != 0) {
    int row = 1 + rand() % 30;
    int col = 1 + rand() % 30;
    m = create_random_matrix(row, col);
    gsl_m = from_s21_matrix_to_gsl(&m);
    fprintf(f, "%d %d\n", row, col);
    dump_matrix("", f, &m);
    gsl_result = gsl_transpose(gsl_m);
    matrix_t result = from_gsl_to_s21_matrix(gsl_result);
    dump_matrix("", f, &result);
    s21_remove_matrix(&m);
    s21_remove_matrix(&result);
    gsl_matrix_free(gsl_m);
    gsl_matrix_free(gsl_result);
    --test_count;
  }
  fclose(f);
}

int main() {
  srand(time(NULL));
  gen_s21_transpose_tests(1000);
  gen_s21_complements_tests(300);
  gen_s21_scale_tests(1000);
  gen_s21_mul_tests(1000);
  gen_s21_sum_tests(1000);
  gen_s21_sub_tests(1000);
  gen_s21_determinant_tests(1000);
  gen_s21_inverse_tests(1000);
  return 0;
}
