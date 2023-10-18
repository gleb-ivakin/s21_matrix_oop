#include "s21_matrix_oop.h"

int main(int argc, char const* argv[]) {
  S21Matrix m1;
  S21Matrix m2(3, 5);

  m1.FillMatrix();
  m2.FillMatrix();

  m1.PrintMatrix();
  m2.PrintMatrix();

  S21Matrix m3 = m1;
  m3.PrintMatrix();

  bool res1, res2;
  res1 = m1.EqMatrix(m2);
  res2 = m1.EqMatrix(m3);

  printf("res1: %s, res2: %s\n", (res1 ? "true" : "false"),
         (res2 ? "true" : "false"));

  m1.SumMatrix(m3);
  m1.PrintMatrix();

  m3.SubMatrix(m1);
  m3.PrintMatrix();

  m1.MulNumber(3.50);
  m1.PrintMatrix();

  S21Matrix m4 = m2.Transpose();
  m4.PrintMatrix();

  S21Matrix m5(2, 3);
  m5.FillMatrix();
  m5.PrintMatrix();
  S21Matrix m6(3, 1);
  m6.FillMatrix();
  m6.PrintMatrix();
  m5.MulMatrix(m6);

  m3.FillMatrixRandom();
  m3.PrintMatrix();
  double det = m3.Determinant();
  printf("det: %lf\n", det);

  double* tesr_arr[0];

  S21Matrix m7 = m3.CalcComplements();
  m7.PrintMatrix();

  printf("----------- Перегрузка операторов ------------\n");
  S21Matrix m8(2, 3);
  S21Matrix m9(2, 3);
  m8.FillMatrixRandom();
  m9.FillMatrixRandom();
  m8.PrintMatrix();
  m9.PrintMatrix();
  S21Matrix m10 = m8 + m9;
  m10.PrintMatrix();

  S21Matrix m11 = m8 - m9;
  m11.PrintMatrix();

  S21Matrix m13(3, 2);
  m13.FillMatrixRandom();
  m13.PrintMatrix();
  S21Matrix m12 = m8 * m13;

  m12.PrintMatrix();
  S21Matrix m14 = m12 * 3.50;
  m14.PrintMatrix();

  printf("----------- += ------------\n");
  S21Matrix m15;
  S21Matrix m16;
  m15.FillMatrixRandom();
  m16.FillMatrixRandom();
  m15.PrintMatrix();
  m16.PrintMatrix();
  m15 += m16;
  m15.PrintMatrix();

  m15 -= m16;
  m15.PrintMatrix();

  m15 *= m16;

  m15 *= 3.50;
  m15.PrintMatrix();

  printf("operator(): %lf\n", m15(2, 1));

  printf("----------- m15 ------------\n");
  S21Matrix m17 = std::move(m15);
  m17.PrintMatrix();
  m15.PrintMatrix();

  printf("----------- Mul ------------\n");
  S21Matrix m18(4, 4);
  S21Matrix m19(4, 4);
  m18.FillMatrixRandom();
  m19.FillMatrixRandom();
  m18.PrintMatrix();
  m19.PrintMatrix();
  m18.MulMatrix(m19);
  m18.PrintMatrix();

  // int* test = (int*) malloc(350);

  // m19.Free();
  m19 = m18;
  printf("----------- m19 ------------\n");
  m19.PrintMatrix();

  // S21Matrix m20(0,0);
  // m20.Determinant();

  S21Matrix m20(3, 3);
  m19 = m20;
  m19.PrintMatrix();

  S21Matrix M1(1, 2), M2(2, 3);
  M1.FillMatrix();
  M1.PrintMatrix();
  M2.FillMatrix();
  M2.PrintMatrix();
  M1.MulMatrix(M2);
  M1.PrintMatrix();

  S21Matrix M3;
  M3.FillMatrix();
  M3.PrintMatrix();
  M3.EditSize(3, 3);
  M3.PrintMatrix();

  return 0;
}
