#include <check.h>
/* #include <gsl/gsl_matrix_double.h> */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../s21_matrix.h"
/* #include "gsl.h" */

/* void CalculateDet() { */
/*   printf(__func__); */
/*   printf("\n"); */
/*   matrix_t m1; */
/*   int rows, cols; */
/*   scanf("%d %d", &rows, &cols); */
/*   read_matrix(&m1, stdin, rows, cols); */
/*   dump_matrix("", stdout, &m1); */
/*   double det; */
/*   s21_determinant(&m1, &det); */
/*   printf("%.19lf\n", det); */
/*   s21_remove_matrix(&m1); */
/* } */
/**/
/* void InverseMatrix() { */
/*   printf(__func__); */
/*   printf("\n"); */
/*   matrix_t m1, m2; */
/*   int rows, cols; */
/*   scanf("%d %d", &rows, &cols); */
/*   read_matrix(&m1, stdin, rows, cols); */
/*   dump_matrix("Initial matrix", stdout, &m1); */
/*   s21_inverse_matrix(&m1, &m2); */
/*   dump_matrix("Inversed matrix", stdout, &m2); */
/*   matrix_t m3; */
/*   s21_mult_matrix(&m1, &m2, &m3); */
/*   dump_matrix("Initial matrix multiplied by inverse matrix", stdout, &m3); */
/*   s21_remove_matrix(&m1); */
/*   s21_remove_matrix(&m2); */
/*   s21_remove_matrix(&m3); */
/* } */
/**/
/* void MultiplyMatrix() { */
/*   printf(__func__); */
/*   printf("\n"); */
/*   matrix_t m1, m2, res; */
/*   int rows, cols; */
/*   scanf("%d %d", &rows, &cols); */
/*   read_matrix(&m1, stdin, rows, cols); */
/*   scanf("%d %d", &rows, &cols); */
/*   read_matrix(&m2, stdin, rows, cols); */
/*   s21_mult_matrix(&m2, &m1, &res); */
/*   dump_matrix("", stdout, &m2); */
/*   dump_matrix("", stdout, &m1); */
/*   dump_matrix("", stdout, &res); */
/*   s21_remove_matrix(&m1); */
/*   s21_remove_matrix(&m2); */
/*   s21_remove_matrix(&res); */
/* } */
/**/
/* void SumMatrix() { */
/*   printf(__func__); */
/*   printf("\n"); */
/*   matrix_t m1, m2, res; */
/*   int rows, cols; */
/*   scanf("%d %d", &rows, &cols); */
/*   read_matrix(&m1, stdin, rows, cols); */
/*   scanf("%d %d", &rows, &cols); */
/*   read_matrix(&m2, stdin, rows, cols); */
/*   dump_matrix("", stdout, &m2); */
/*   dump_matrix("", stdout, &m1); */
/*   s21_sum_matrix(&m2, &m1, &res); */
/*   dump_matrix("", stdout, &res); */
/*   s21_remove_matrix(&m1); */
/*   s21_remove_matrix(&m2); */
/*   s21_remove_matrix(&res); */
/* } */
/**/
/* void SubMatrix() { */
/*   printf(__func__); */
/*   printf("\n"); */
/*   matrix_t m1, m2, res; */
/*   int rows, cols; */
/*   scanf("%d %d", &rows, &cols); */
/*   read_matrix(&m1, stdin, rows, cols); */
/*   scanf("%d %d", &rows, &cols); */
/*   read_matrix(&m2, stdin, rows, cols); */
/*   dump_matrix("", stdout, &m2); */
/*   dump_matrix("", stdout, &m1); */
/*   s21_sub_matrix(&m2, &m1, &res); */
/*   dump_matrix("", stdout, &res); */
/*   s21_remove_matrix(&m1); */
/*   s21_remove_matrix(&m2); */
/*   s21_remove_matrix(&res); */
/* } */
/**/
/* void TransposeTest() { */
/*   printf(__func__); */
/*   printf("\n"); */
/*   matrix_t m1, m2; */
/*   int rows, cols; */
/*   scanf("%d %d", &rows, &cols); */
/*   read_matrix(&m1, stdin, rows, cols); */
/*   dump_matrix("Original matrix", stdout, &m1); */
/*   s21_transpose(&m1, &m2); */
/*   dump_matrix("Transposed matrix", stdout, &m2); */
/*   s21_remove_matrix(&m1); */
/*   s21_remove_matrix(&m2); */
/* } */
/**/
/* void CalcComplementsTest() { */
/*   printf(__func__); */
/*   printf("\n"); */
/*   matrix_t m1, m2; */
/*   int rows, cols; */
/*   scanf("%d %d", &rows, &cols); */
/*   read_matrix(&m1, stdin, rows, cols); */
/*   dump_matrix("Original matrix", stdout, &m1); */
/*   s21_calc_complements(&m1, &m2); */
/*   dump_matrix("Complements matrix", stdout, &m2); */
/*   s21_remove_matrix(&m1); */
/*   s21_remove_matrix(&m2); */
/* } */

/* void test_determinant_gsl() { */
/*   matrix_t m; */
/*   gsl_matrix *gsl_m; */
/*   int counter = 0; */
/*   while (counter++ < 100) { */
/*     int dimention = 1 + rand() % 30; */
/*     m = create_random_matrix(dimention, dimention); */
/*     gsl_m = from_s21_matrix_to_gsl(&m); */
/*     double mine_det = determinant_by_definition(&m); */
/*     double det_gauss, gsl_det = gsl_determinant(gsl_m); */
/*     s21_determinant(&m, &det_gauss); */
/*     if (!(mine_det == 0.0 || mine_det == -0.0 || */
/*           (abs_double(det_gauss / gsl_det - 1.0) < EPSILON && */
/*            abs_double(mine_det / gsl_det - 1.0) < EPSILON))) { */
/*       dump_matrix("Random Matrix", stdout, &m); */
/*       dump_gsl_matrix("GSL Matrix", gsl_m); */
/*       printf("%6.10f %6.10f %6.10f\n", mine_det, det_gauss, gsl_det); */
/*       break; */
/*     } */
/*     printf("%6.10f %6.10f %6.10f\n", mine_det, det_gauss, gsl_det); */
/*     s21_remove_matrix(&m); */
/*     gsl_matrix_free(gsl_m); */
/*   } */
/* } */
/**/
/* void test_inverse_gsl() { */
/*   matrix_t m; */
/*   gsl_matrix *gsl_m; */
/*   int counter = 0; */
/*   while ((counter++ < 100)) { */
/*     int dimention = 1 + rand() % 30; */
/*     m = create_random_matrix(dimention, dimention); */
/*     dump_matrix("Original matrix", stdout, &m); */
/*     gsl_m = from_s21_matrix_to_gsl(&m); */
/*     matrix_t mine_inverse; */
/*     int status = s21_inverse_matrix(&m, &mine_inverse); */
/*     double det; */
/*     s21_determinant(&m, &det); */
/*     printf("Determinant: %6.15f\n", det); */
/*     if (!status && abs_double(det) > EPSILON) { */
/*       printf("Status: %d\n", status); */
/*       dump_matrix("Mine inversed", stdout, &mine_inverse); */
/*       gsl_matrix *gsl_inverse = gsl_matrix_calloc(m.rows, m.columns); */
/*       gsl_invert_matrix(gsl_m, gsl_inverse); */
/*       matrix_t inverted_conv_matrix_t = from_gsl_to_s21_matrix(gsl_inverse);
 */
/*       dump_matrix("GSL inversed", stdout, &inverted_conv_matrix_t); */
/*       assert(s21_eq_matrix(&mine_inverse, &inverted_conv_matrix_t)); */
/*       s21_remove_matrix(&inverted_conv_matrix_t); */
/*       gsl_matrix_free(gsl_inverse); */
/*       s21_remove_matrix(&mine_inverse); */
/*     } */
/*     s21_remove_matrix(&m); */
/*     gsl_matrix_free(gsl_m); */
/*   } */
/* } */

START_TEST(test_precalculated_inverse) {
  fprintf(stdout, __func__);
  fprintf(stdout, "\n");
  FILE *f = fopen("test_files/s21_inverse_test", "r");
  if (!f) {
    fprintf(stderr, "Input file is not found\n");
    return;
  }
  int test_count, test_count_bak;
  fscanf(f, "%d\n", &test_count);
  test_count_bak = test_count;
  matrix_t m, inversed_m, precalc_answer;
  int successful_tests = 0;
  while (test_count != 0) {
    int dimention;
    fscanf(f, "%d\n", &dimention);
    read_matrix(&m, f, dimention, dimention);
    read_matrix(&precalc_answer, f, dimention, dimention);
    int status = s21_inverse_matrix(&m, &inversed_m);
    successful_tests +=
        status == 0 && s21_eq_matrix(&inversed_m, &precalc_answer);
    s21_remove_matrix(&m);
    s21_remove_matrix(&inversed_m);
    s21_remove_matrix(&precalc_answer);
    --test_count;
  }
  fclose(f);
  fprintf(stdout, "Tests performed: %d Successfull tests: %d\n", test_count_bak,
          successful_tests);
  ck_assert_int_eq(test_count_bak, successful_tests);
}
END_TEST

START_TEST(test_precalculated_determinant) {
  fprintf(stdout, __func__);
  fprintf(stdout, "\n");
  FILE *f = fopen("test_files/s21_determinant_test", "r");
  if (!f) {
    fprintf(stderr, "Input file is not found\n");
    return;
  }
  int test_count, test_count_bak;
  fscanf(f, "%d\n", &test_count);
  test_count_bak = test_count;
  matrix_t m;
  int successful_tests = 0;
  while (test_count != 0) {
    int dimention;
    fscanf(f, "%d\n", &dimention);
    read_matrix(&m, f, dimention, dimention);
    double precalculated_determinant;
    fscanf(f, "%lf\n", &precalculated_determinant);
    double determinant;
    int status = s21_determinant(&m, &determinant);
    successful_tests +=
        status == 0 &&
        ((determinant == 0.0 && precalculated_determinant == 0.0) ||
         abs_double(determinant / precalculated_determinant - 1.0) < EPSILON);
    --test_count;
    s21_remove_matrix(&m);
  }
  fclose(f);
  fprintf(stdout, "Tests performed: %d Successfull tests: %d\n", test_count_bak,
          successful_tests);
  ck_assert_int_eq(test_count_bak, successful_tests);
}
END_TEST

START_TEST(test_precalculated_sum) {
  fprintf(stdout, __func__);
  fprintf(stdout, "\n");
  FILE *f = fopen("test_files/s21_sum_test", "r");
  if (!f) {
    fprintf(stderr, "Input file is not found\n");
    return;
  }
  int test_count, test_count_bak;
  fscanf(f, "%d\n", &test_count);
  test_count_bak = test_count;
  matrix_t lhs, rhs, precalculated_result, result;
  int successful_tests = 0;
  while (test_count != 0) {
    int row, col;
    fscanf(f, "%d %d\n", &row, &col);
    read_matrix(&lhs, f, row, col);
    read_matrix(&rhs, f, row, col);
    read_matrix(&precalculated_result, f, row, col);
    successful_tests += !s21_sum_matrix(&lhs, &rhs, &result) &&
                        s21_eq_matrix(&result, &precalculated_result);
    --test_count;
    s21_remove_matrix(&lhs);
    s21_remove_matrix(&rhs);
    s21_remove_matrix(&precalculated_result);
    s21_remove_matrix(&result);
  }
  fclose(f);
  fprintf(stdout, "Tests performed: %d Successfull tests: %d\n", test_count_bak,
          successful_tests);
  ck_assert_int_eq(test_count_bak, successful_tests);
}
END_TEST

START_TEST(test_precalculated_sub) {
  fprintf(stdout, __func__);
  fprintf(stdout, "\n");
  FILE *f = fopen("test_files/s21_sub_test", "r");
  if (!f) {
    fprintf(stderr, "Input file is not found\n");
    return;
  }
  int test_count, test_count_bak;
  fscanf(f, "%d\n", &test_count);
  test_count_bak = test_count;
  matrix_t lhs, rhs, precalculated_result, result;
  int successful_tests = 0;
  while (test_count != 0) {
    int row, col;
    fscanf(f, "%d %d\n", &row, &col);
    read_matrix(&lhs, f, row, col);
    read_matrix(&rhs, f, row, col);
    read_matrix(&precalculated_result, f, row, col);
    successful_tests += !s21_sub_matrix(&lhs, &rhs, &result) &&
                        s21_eq_matrix(&result, &precalculated_result);
    --test_count;
    s21_remove_matrix(&lhs);
    s21_remove_matrix(&rhs);
    s21_remove_matrix(&precalculated_result);
    s21_remove_matrix(&result);
  }
  fclose(f);
  fprintf(stdout, "Tests performed: %d Successfull tests: %d\n", test_count_bak,
          successful_tests);
  ck_assert_int_eq(test_count_bak, successful_tests);
}
END_TEST

START_TEST(test_precalculated_mul) {
  fprintf(stdout, __func__);
  fprintf(stdout, "\n");
  FILE *f = fopen("test_files/s21_mul_test", "r");
  if (!f) {
    fprintf(stderr, "Input file is not found\n");
    return;
  }
  int test_count, test_count_bak;
  fscanf(f, "%d\n", &test_count);
  test_count_bak = test_count;
  matrix_t lhs, rhs, precalculated_result, result;
  int successful_tests = 0;
  while (test_count != 0) {
    int row1, col, col2;
    fscanf(f, "%d %d %d\n", &row1, &col, &col2);
    read_matrix(&lhs, f, row1, col);
    read_matrix(&rhs, f, col, col2);
    read_matrix(&precalculated_result, f, row1, col2);
    successful_tests += !s21_mult_matrix(&lhs, &rhs, &result) &&
                        s21_eq_matrix(&result, &precalculated_result);
    --test_count;
    s21_remove_matrix(&lhs);
    s21_remove_matrix(&rhs);
    s21_remove_matrix(&precalculated_result);
    s21_remove_matrix(&result);
  }
  fclose(f);
  fprintf(stdout, "Tests performed: %d Successfull tests: %d\n", test_count_bak,
          successful_tests);
  ck_assert_int_eq(test_count_bak, successful_tests);
}
END_TEST

START_TEST(test_precalculated_scale) {
  fprintf(stdout, __func__);
  fprintf(stdout, "\n");
  FILE *f = fopen("test_files/s21_scale_test", "r");
  if (!f) {
    fprintf(stderr, "Input file is not found\n");
    return;
  }
  int test_count, test_count_bak;
  fscanf(f, "%d\n", &test_count);
  test_count_bak = test_count;
  matrix_t m, precalculated_result, result;
  int successful_tests = 0;
  while (test_count != 0) {
    int row, col;
    double multiplyer;
    fscanf(f, "%d %d %lf\n", &row, &col, &multiplyer);
    read_matrix(&m, f, row, col);
    read_matrix(&precalculated_result, f, row, col);
    successful_tests += !s21_mult_number(&m, multiplyer, &result) &&
                        s21_eq_matrix(&result, &precalculated_result);
    --test_count;
    s21_remove_matrix(&m);
    s21_remove_matrix(&precalculated_result);
    s21_remove_matrix(&result);
  }
  fclose(f);
  fprintf(stdout, "Tests performed: %d Successfull tests: %d\n", test_count_bak,
          successful_tests);
  ck_assert_int_eq(test_count_bak, successful_tests);
}
END_TEST

START_TEST(test_precalculated_complements) {
  fprintf(stdout, __func__);
  fprintf(stdout, "\n");
  FILE *f = fopen("test_files/s21_complements_test", "r");
  if (!f) {
    fprintf(stderr, "Input file is not found\n");
    return;
  }
  int test_count, test_count_bak;
  fscanf(f, "%d\n", &test_count);
  test_count_bak = test_count;
  matrix_t m, precalculated_result, result;
  int successful_tests = 0;
  while (test_count != 0) {
    int dimention;
    fscanf(f, "%d\n", &dimention);
    read_matrix(&m, f, dimention, dimention);
    read_matrix(&precalculated_result, f, dimention, dimention);
    /* fprintf(stdout, "Status: %d\n", s21_calc_complements(&m, &result)); */
    /* fprintf(stdout, "Comparison: %d\n", */
    /*         s21_eq_matrix(&result, &precalculated_result)); */
    successful_tests += !s21_calc_complements(&m, &result) &&
                        s21_eq_matrix_relative(&result, &precalculated_result);
    if (!s21_eq_matrix_relative(&precalculated_result, &result)) {
      dump_matrix("Result: ", stdout, &result);
      dump_matrix("Precalculated result: ", stdout, &precalculated_result);
    }
    --test_count;
    s21_remove_matrix(&m);
    s21_remove_matrix(&precalculated_result);
    s21_remove_matrix(&result);
  }
  fclose(f);
  fprintf(stdout, "Tests performed: %d Successfull tests: %d\n", test_count_bak,
          successful_tests);
  ck_assert_int_eq(test_count_bak, successful_tests);
}
END_TEST

START_TEST(test_precalculated_transpose) {
  fprintf(stdout, __func__);
  fprintf(stdout, "\n");
  FILE *f = fopen("test_files/s21_transpose_test", "r");
  if (!f) {
    fprintf(stderr, "Input file is not found\n");
    return;
  }
  int test_count, test_count_bak;
  fscanf(f, "%d\n", &test_count);
  test_count_bak = test_count;
  matrix_t m, precalculated_result, result;
  int successful_tests = 0;
  while (test_count != 0) {
    int row, col;
    fscanf(f, "%d %d\n", &row, &col);
    read_matrix(&m, f, row, col);
    read_matrix(&precalculated_result, f, col, row);
    successful_tests += !s21_transpose(&m, &result) &&
                        s21_eq_matrix(&result, &precalculated_result);
    --test_count;
    s21_remove_matrix(&m);
    s21_remove_matrix(&precalculated_result);
    s21_remove_matrix(&result);
  }
  fclose(f);
  ck_assert_int_eq(test_count_bak, successful_tests);
}
END_TEST

START_TEST(test_eq) {
  int test_count = 100000;
  fprintf(stdout, __func__);
  fprintf(stdout, "\n");
  int test_count_bak;
  test_count_bak = test_count;
  matrix_t m, m1, m2, m3, m4, m5, m6, m7;
  int successful_tests = 0;
  while (test_count != 0) {
    int row = 1 + rand() % 30, col = 1 + rand() % 30;
    m = create_random_matrix(row, col);
    copy_matrix(&m, &m1);
    copy_matrix(&m, &m2);
    copy_matrix(&m, &m3);
    copy_matrix(&m, &m4);
    copy_matrix(&m, &m5);
    copy_matrix(&m, &m6);
    copy_matrix(&m, &m7);
    add_num_to_matrix(&m2, EPSILON * 1.1);
    add_num_to_matrix(&m3, EPSILON / 2.0);
    add_num_to_matrix(&m4, EPSILON * 2.0);
    m5.rows += 1;
    m6.columns += 1;
    free(m7.matrix);
    m7.matrix = NULL;
    successful_tests += s21_eq_matrix(&m, &m1) && !s21_eq_matrix(&m, &m2) &&
                        s21_eq_matrix(&m, &m3) && !s21_eq_matrix(&m, &m4) &&
                        !s21_eq_matrix(&m, &m5) && !s21_eq_matrix(&m, &m6) &&
                        !s21_eq_matrix(&m, NULL) && !s21_eq_matrix(NULL, &m) &&
                        !s21_eq_matrix(NULL, NULL) && !s21_eq_matrix(&m, &m7);
    --test_count;
    s21_remove_matrix(&m);
    s21_remove_matrix(&m1);
    s21_remove_matrix(&m2);
    s21_remove_matrix(&m3);
    s21_remove_matrix(&m4);
    s21_remove_matrix(&m5);
    s21_remove_matrix(&m6);
    s21_remove_matrix(&m7);
  }
  fprintf(stdout, "Tests performed: %d Successfull tests: %d\n", test_count_bak,
          successful_tests);
  ck_assert_int_eq(test_count_bak, successful_tests);
}
END_TEST

START_TEST(test_error_codes_create) {
  matrix_t m;
  ck_assert_int_eq(s21_create_matrix(0, 1, &m), EXIT_CODE_INCORRECT_MATRIX);
  ck_assert_int_eq(s21_create_matrix(1, 0, &m), EXIT_CODE_INCORRECT_MATRIX);
  ck_assert_int_eq(s21_create_matrix(0, 0, &m), EXIT_CODE_INCORRECT_MATRIX);
  ck_assert_int_eq(s21_create_matrix(2, 2, NULL), EXIT_CODE_INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_error_codes_sum_and_sub) {
  matrix_t a = {}, b = {}, r = {};
  s21_create_matrix(2, 2, &a);
  s21_create_matrix(3, 2, &b);
  ck_assert_int_eq(s21_sum_matrix(&a, &b, &r), EXIT_CODE_CALCULATION_ERROR);
  ck_assert_int_eq(s21_sub_matrix(&a, &b, &r), EXIT_CODE_CALCULATION_ERROR);
  s21_remove_matrix(&a);
  s21_remove_matrix(&b);
}
END_TEST

START_TEST(test_error_codes_mult_matrix) {
  matrix_t a = {}, b = {}, r = {};
  s21_create_matrix(2, 3, &a);
  s21_create_matrix(4, 2, &b);
  ck_assert_int_eq(s21_mult_matrix(&a, &b, &r), EXIT_CODE_CALCULATION_ERROR);
  s21_remove_matrix(&a);
  s21_remove_matrix(&b);
}
END_TEST

START_TEST(test_error_codes_determinant) {
  matrix_t a = {};
  double result;
  ck_assert_int_eq(s21_determinant(NULL, &result), EXIT_CODE_INCORRECT_MATRIX);
  s21_create_matrix(2, 3, &a);
  ck_assert_int_eq(s21_determinant(&a, &result), EXIT_CODE_INCORRECT_MATRIX);
  s21_remove_matrix(&a);
}
END_TEST

START_TEST(test_error_codes_inverse) {
  matrix_t a = {}, r = {};
  s21_create_matrix(2, 2, &a);
  a.matrix[0][0] = 1;
  a.matrix[0][1] = 2;
  a.matrix[1][0] = 2;
  a.matrix[1][1] = 4;
  ck_assert_int_eq(s21_inverse_matrix(&a, &r), EXIT_CODE_CALCULATION_ERROR);
  s21_remove_matrix(&a);
}
END_TEST

START_TEST(test_error_codes_complements) {
  matrix_t a = {}, r = {};
  s21_create_matrix(2, 3, &a);
  ck_assert_int_eq(s21_calc_complements(&a, &r), EXIT_CODE_CALCULATION_ERROR);
  s21_remove_matrix(&a);

  s21_create_matrix(1, 1, &a);
  ck_assert_int_eq(s21_calc_complements(&a, &r), EXIT_CODE_CALCULATION_ERROR);
  s21_remove_matrix(&a);
}
END_TEST

START_TEST(test_error_codes_transpose) {
  ck_assert_int_eq(s21_transpose(NULL, NULL), EXIT_CODE_INCORRECT_MATRIX);
}
END_TEST

START_TEST(test_success_create_and_remove) {
  matrix_t m;
  ck_assert_int_eq(s21_create_matrix(3, 3, &m), EXIT_CODE_SUCCESS);
  ck_assert_ptr_nonnull(m.matrix);
  ck_assert_int_eq(m.rows, 3);
  ck_assert_int_eq(m.columns, 3);
  s21_remove_matrix(&m);
  ck_assert_ptr_null(m.matrix);
}
END_TEST

START_TEST(test_success_sum_sub) {
  matrix_t a = {}, b = {}, result = {};
  s21_create_matrix(2, 2, &a);
  s21_create_matrix(2, 2, &b);
  a.matrix[0][0] = 1;
  b.matrix[0][0] = 2;
  a.matrix[1][1] = 3;
  b.matrix[1][1] = 4;
  ck_assert_int_eq(s21_sum_matrix(&a, &b, &result), EXIT_CODE_SUCCESS);
  ck_assert_double_eq(result.matrix[0][0], 3);
  ck_assert_double_eq(result.matrix[1][1], 7);
  s21_remove_matrix(&result);
  ck_assert_int_eq(s21_sub_matrix(&a, &b, &result), EXIT_CODE_SUCCESS);
  ck_assert_double_eq(result.matrix[0][0], -1);
  ck_assert_double_eq(result.matrix[1][1], -1);
  s21_remove_matrix(&a);
  s21_remove_matrix(&b);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_success_mult_matrix) {
  matrix_t a = {}, b = {}, result = {};
  s21_create_matrix(2, 3, &a);
  s21_create_matrix(3, 2, &b);

  int val = 1;
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j) a.matrix[i][j] = val++;

  val = 1;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 2; ++j) b.matrix[i][j] = val++;

  ck_assert_int_eq(s21_mult_matrix(&a, &b, &result), EXIT_CODE_SUCCESS);
  ck_assert_int_eq(result.rows, 2);
  ck_assert_int_eq(result.columns, 2);
  s21_remove_matrix(&a);
  s21_remove_matrix(&b);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_success_determinant) {
  matrix_t m = {};
  s21_create_matrix(2, 2, &m);
  m.matrix[0][0] = 4;
  m.matrix[0][1] = 6;
  m.matrix[1][0] = 3;
  m.matrix[1][1] = 8;
  double det = 0.0;
  ck_assert_int_eq(s21_determinant(&m, &det), EXIT_CODE_SUCCESS);
  ck_assert_double_eq_tol(det, 14.0, EPSILON);
  s21_remove_matrix(&m);
}
END_TEST

START_TEST(test_success_inverse) {
  matrix_t m = {}, inv = {};
  s21_create_matrix(2, 2, &m);
  m.matrix[0][0] = 4;
  m.matrix[0][1] = 7;
  m.matrix[1][0] = 2;
  m.matrix[1][1] = 6;
  ck_assert_int_eq(s21_inverse_matrix(&m, &inv), EXIT_CODE_SUCCESS);
  double expected[2][2] = {{0.6, -0.7}, {-0.2, 0.4}};
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      ck_assert_double_eq_tol(inv.matrix[i][j], expected[i][j], EPSILON);
  s21_remove_matrix(&m);
  s21_remove_matrix(&inv);
}
END_TEST

START_TEST(test_success_transpose) {
  matrix_t m = {}, tr = {};
  s21_create_matrix(2, 3, &m);
  int val = 1;
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j) m.matrix[i][j] = val++;
  ck_assert_int_eq(s21_transpose(&m, &tr), EXIT_CODE_SUCCESS);
  ck_assert_int_eq(tr.rows, 3);
  ck_assert_int_eq(tr.columns, 2);
  ck_assert_double_eq(tr.matrix[0][0], 1);
  ck_assert_double_eq(tr.matrix[1][0], 2);
  ck_assert_double_eq(tr.matrix[2][1], 6);
  s21_remove_matrix(&m);
  s21_remove_matrix(&tr);
}
END_TEST

START_TEST(test_success_complements) {
  matrix_t m = {}, comp = {};
  s21_create_matrix(3, 3, &m);
  int vals[3][3] = {{1, 2, 3}, {0, 4, 2}, {5, 2, 1}};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) m.matrix[i][j] = vals[i][j];
  ck_assert_int_eq(s21_calc_complements(&m, &comp), EXIT_CODE_SUCCESS);
  s21_remove_matrix(&m);
  s21_remove_matrix(&comp);
}
END_TEST

Suite *s21_test_suite(void) {
  Suite *s;
  TCase *tc_core;
  s = suite_create("s21_test_suite");
  tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_precalculated_inverse);
  tcase_add_test(tc_core, test_precalculated_scale);
  tcase_add_test(tc_core, test_eq);
  tcase_add_test(tc_core, test_precalculated_transpose);
  tcase_add_test(tc_core, test_precalculated_mul);
  tcase_add_test(tc_core, test_precalculated_inverse);
  tcase_add_test(tc_core, test_precalculated_complements);
  tcase_add_test(tc_core, test_precalculated_determinant);
  tcase_add_test(tc_core, test_precalculated_sum);
  tcase_add_test(tc_core, test_precalculated_sub);
  tcase_add_test(tc_core, test_error_codes_create);
  tcase_add_test(tc_core, test_error_codes_sum_and_sub);
  tcase_add_test(tc_core, test_error_codes_mult_matrix);
  tcase_add_test(tc_core, test_error_codes_determinant);
  tcase_add_test(tc_core, test_error_codes_inverse);
  tcase_add_test(tc_core, test_error_codes_complements);
  tcase_add_test(tc_core, test_error_codes_transpose);
  tcase_add_test(tc_core, test_success_create_and_remove);
  tcase_add_test(tc_core, test_success_sum_sub);
  tcase_add_test(tc_core, test_success_mult_matrix);
  tcase_add_test(tc_core, test_success_determinant);
  tcase_add_test(tc_core, test_success_inverse);
  tcase_add_test(tc_core, test_success_transpose);
  tcase_add_test(tc_core, test_success_complements);
  suite_add_tcase(s, tc_core);
  return s;
}

int main() {
  srand(time(NULL));
  Suite *s;
  SRunner *sr;
  s = s21_test_suite();
  sr = srunner_create(s);
  srunner_set_fork_status(sr, CK_NOFORK);
  srunner_run_all(sr, CK_NORMAL);
  int number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
