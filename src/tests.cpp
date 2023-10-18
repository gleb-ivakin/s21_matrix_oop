#include <iostream>

#include "gtest/gtest.h"
#include "s21_matrix_oop.h"

// ---------- Конструкторы ----------

TEST(Constructors, Base) {
  S21Matrix M;
  EXPECT_EQ(M.GetRows(), 3);
  EXPECT_EQ(M.GetCols(), 3);
}

TEST(Constructors, Parameterized) {
  S21Matrix M1(10, 15);
  EXPECT_EQ(M1.GetRows(), 10);
  EXPECT_EQ(M1.GetCols(), 15);

  // пробуем создать матрицу (0, 0)
  bool exception = false;
  try {
    S21Matrix M2(0, 0);
  } catch (const std::out_of_range &e) {
    exception = true;
  }
  EXPECT_TRUE(exception);
}

TEST(Constructors, Copy) {
  S21Matrix M1;
  M1.FillMatrixRandom();
  S21Matrix M2 = M1;
  for (int y = 0; y < M1.GetRows(); y++)
    for (int x = 0; x < M1.GetRows(); x++) EXPECT_EQ(M1(y, x), M2(y, x));
}

TEST(Constructors, Move) {
  S21Matrix M1;
  M1.FillMatrixRandom();
  S21Matrix M2 = M1;
  S21Matrix M3 = std::move(M1);

  for (int y = 0; y < M3.GetRows(); y++)
    for (int x = 0; x < M3.GetRows(); x++) EXPECT_EQ(M3(y, x), M2(y, x));

  EXPECT_EQ(M1.GetRows(), 0);
  EXPECT_EQ(M1.GetCols(), 0);
}

TEST(Constructors, Destructor) {
  S21Matrix M1;
  M1.~S21Matrix();

  EXPECT_EQ(M1.GetRows(), 0);
  EXPECT_EQ(M1.GetCols(), 0);
}

// ----------- Методы -----------

TEST(Methods, EqMatrix) {
  S21Matrix M1;
  M1.FillMatrixRandom();
  S21Matrix M2 = M1;
  EXPECT_TRUE(M1.EqMatrix(M2));

  S21Matrix M3;
  M3.FillMatrix();
  EXPECT_FALSE(M3.EqMatrix(M2));
}

TEST(Methods, SumMatrix) {
  S21Matrix M1, M2;
  M1.FillMatrixRandom();
  M2.FillMatrixRandom();
  S21Matrix M3 = M1;
  M1.SumMatrix(M2);

  for (int y = 0; y < M1.GetRows(); y++)
    for (int x = 0; x < M1.GetRows(); x++)
      EXPECT_EQ(M1(y, x), M2(y, x) + M3(y, x));

  // матрицы с разными размерами
  bool exception = false;
  try {
    S21Matrix M4(5, 5);
    S21Matrix M5(2, 8);
    M4.SumMatrix(M5);
  } catch (const std::out_of_range &e) {
    exception = true;
  }
  EXPECT_TRUE(exception);
}

TEST(Methods, SubMatrix) {
  S21Matrix M1, M2;
  M1.FillMatrixRandom();
  M2.FillMatrixRandom();
  S21Matrix M3 = M1;
  M1.SubMatrix(M2);

  for (int y = 0; y < M1.GetRows(); y++)
    for (int x = 0; x < M1.GetRows(); x++)
      EXPECT_EQ(M1(y, x), M3(y, x) - M2(y, x));

  // матрицы с разными размерами
  bool exception = false;
  try {
    S21Matrix M4(5, 5);
    S21Matrix M5(2, 8);
    M4.SubMatrix(M5);
  } catch (const std::out_of_range &e) {
    exception = true;
  }
  EXPECT_TRUE(exception);
}

TEST(Methods, MulNumber) {
  S21Matrix M1;
  M1.FillMatrixRandom();
  S21Matrix M2 = M1;
  double number = 3.50;
  M1.MulNumber(number);

  for (int y = 0; y < M1.GetRows(); y++)
    for (int x = 0; x < M1.GetRows(); x++)
      EXPECT_EQ(M1(y, x), M2(y, x) * number);
}

TEST(Methods, MulMatrix) {
  // [0, 1] * [[0, 1, 2] [3, 4, 5]] == [3, 4, 5]
  S21Matrix M1(1, 2), M2(2, 3);
  M1.FillMatrix();
  M2.FillMatrix();
  M1.MulMatrix(M2);

  EXPECT_EQ(M1(0, 0), 3);
  EXPECT_EQ(M1(0, 1), 4);
  EXPECT_EQ(M1(0, 2), 5);

  // M3.cols != M4.rows
  bool exception = false;
  try {
    S21Matrix M3(5, 5);
    S21Matrix M4(2, 7);
    M3.MulMatrix(M4);
  } catch (const std::out_of_range &e) {
    exception = true;
  }
  EXPECT_TRUE(exception);
}

TEST(Methods, Transpose) {
  // [0, 1, 2] == [[0], [1], [2]]
  S21Matrix M1(1, 3);
  M1.FillMatrix();
  S21Matrix M2 = M1.Transpose();

  EXPECT_EQ(M2(0, 0), 0);
  EXPECT_EQ(M2(1, 0), 1);
  EXPECT_EQ(M2(2, 0), 2);
}

TEST(Methods, CalcComplements) {
  S21Matrix M1;
  M1.SetValue(0, 0, 1);
  M1.SetValue(0, 1, 2);
  M1.SetValue(0, 2, 3);
  M1.SetValue(1, 0, 0);
  M1.SetValue(1, 1, 4);
  M1.SetValue(1, 2, 2);
  M1.SetValue(2, 0, 5);
  M1.SetValue(2, 1, 2);
  M1.SetValue(2, 2, 1);

  double complement[9] = {0, 10, -20, 4, -14, 8, -8, -2, 4};
  S21Matrix M2 = M1.CalcComplements();

  for (int y = 0, n = 0; y < M1.GetRows(); y++)
    for (int x = 0; x < M1.GetCols(); x++, n++)
      EXPECT_EQ(M2(y, x), complement[n]);

  // Определитель равен нулю
  bool exception = false;
  try {
    S21Matrix M1;
    M1.FillMatrix();
    M1.CalcComplements();
  } catch (const std::out_of_range &e) {
    exception = true;
  }
  EXPECT_TRUE(exception);

  // Неквадратная матрица
  try {
    S21Matrix M1(4, 5);
    M1.FillMatrix();
    M1.CalcComplements();
  } catch (const std::out_of_range &e) {
    exception = true;
  }
  EXPECT_TRUE(exception);
}

TEST(Methods, Determinant) {
  double det;
  S21Matrix M1(3, 3);
  det = M1.Determinant();
  EXPECT_EQ(det, 0);

  S21Matrix M2(2, 2);
  M2.SetValue(0, 0, 5);
  M2.SetValue(0, 1, 2);
  M2.SetValue(1, 0, 3);
  M2.SetValue(1, 1, 4);
  det = M2.Determinant();
  EXPECT_EQ(det, 14);

  // Неквадратная матрица
  bool exception = false;
  try {
    S21Matrix M1(4, 5);
    M1.FillMatrix();
    M1.Determinant();
  } catch (const std::out_of_range &e) {
    exception = true;
  }
  EXPECT_TRUE(exception);
}

TEST(Methods, InverseMatrix) {
  S21Matrix M1(2, 2);
  M1.SetValue(0, 0, 0);
  M1.SetValue(0, 1, 1);
  M1.SetValue(1, 0, 4);
  M1.SetValue(1, 1, 4);
  S21Matrix Inv = M1.InverseMatrix();
  EXPECT_EQ(Inv(0, 0), -1);
  EXPECT_EQ(Inv(0, 1), 0.25);
  EXPECT_EQ(Inv(1, 0), 1);
  EXPECT_EQ(Inv(1, 1), 0);

  // Определитель равен нулю
  bool exception = false;
  try {
    S21Matrix M1;
    M1.FillMatrix();
    M1.CalcComplements();
  } catch (const std::out_of_range &e) {
    exception = true;
  }
  EXPECT_TRUE(exception);

  // Неквадратная матрица
  try {
    S21Matrix M1(4, 5);
    M1.FillMatrix();
    M1.CalcComplements();
  } catch (const std::out_of_range &e) {
    exception = true;
  }
  EXPECT_TRUE(exception);
}

// ----------- Перегрузка операторов ------------

TEST(Operators, Plus) {
  S21Matrix M1, M2;
  M1.FillMatrixRandom();
  M2.FillMatrixRandom();
  S21Matrix M3 = M1;
  S21Matrix M4 = M1 + M2;

  for (int y = 0; y < M4.GetRows(); y++)
    for (int x = 0; x < M4.GetRows(); x++)
      EXPECT_EQ(M4(y, x), M2(y, x) + M3(y, x));

  // матрицы с разными размерами
  bool exception = false;
  try {
    S21Matrix M5(5, 5);
    S21Matrix M6(2, 7);
    S21Matrix M7 = M5 + M6;
  } catch (const std::out_of_range &e) {
    exception = true;
  }
  EXPECT_TRUE(exception);
}

TEST(Operators, Minus) {
  S21Matrix M1, M2;
  M1.FillMatrixRandom();
  M2.FillMatrixRandom();
  S21Matrix M3 = M1;
  S21Matrix M4 = M1 - M2;

  for (int y = 0; y < M4.GetRows(); y++)
    for (int x = 0; x < M4.GetRows(); x++)
      EXPECT_EQ(M4(y, x), M3(y, x) - M2(y, x));

  // матрицы с разными размерами
  bool exception = false;
  try {
    S21Matrix M5(5, 5);
    S21Matrix M6(2, 7);
    S21Matrix M7 = M5 - M6;
  } catch (const std::out_of_range &e) {
    exception = true;
  }
  EXPECT_TRUE(exception);
}

TEST(Operators, MulNumber) {
  S21Matrix M1;
  M1.FillMatrixRandom();
  double number = 3.50;
  S21Matrix M2 = M1;
  S21Matrix M3 = M1 * number;

  for (int y = 0; y < M3.GetRows(); y++)
    for (int x = 0; x < M3.GetRows(); x++)
      EXPECT_EQ(M3(y, x), M2(y, x) * number);
}

TEST(Operators, MulMatrix) {
  // [0, 1] * [[0, 1, 2] [3, 4, 5]] == [3, 4, 5]
  S21Matrix M1(1, 2), M2(2, 3);
  M1.FillMatrix();
  M2.FillMatrix();
  S21Matrix M3 = M1 * M2;

  EXPECT_EQ(M3(0, 0), 3);
  EXPECT_EQ(M3(0, 1), 4);
  EXPECT_EQ(M3(0, 2), 5);

  // M3.cols != M4.rows
  bool exception = false;
  try {
    S21Matrix M3(5, 5);
    S21Matrix M4(2, 7);
    S21Matrix M5 = M3 * M4;
  } catch (const std::out_of_range &e) {
    exception = true;
  }
  EXPECT_TRUE(exception);
}

TEST(Operators, Equal) {
  S21Matrix M1;
  M1.FillMatrixRandom();
  S21Matrix M2 = M1;
  EXPECT_TRUE(M1 == M2);

  S21Matrix M3;
  M3.FillMatrix();
  EXPECT_FALSE(M3 == M2);
}

TEST(Operators, Assignment) {
  S21Matrix M1;
  M1.FillMatrixRandom();
  S21Matrix M2 = M1;
  EXPECT_TRUE(M1 == M2);
}

TEST(Operators, PlusEquals) {
  S21Matrix M1, M2;
  M1.FillMatrixRandom();
  M2.FillMatrixRandom();
  S21Matrix M3 = M1;
  M3 += M2;

  for (int y = 0; y < M3.GetRows(); y++)
    for (int x = 0; x < M3.GetRows(); x++)
      EXPECT_EQ(M3(y, x), M1(y, x) + M2(y, x));

  // матрицы с разными размерами
  bool exception = false;
  try {
    S21Matrix M5(5, 5);
    S21Matrix M6(2, 7);
    M5 += M6;
  } catch (const std::out_of_range &e) {
    exception = true;
  }
  EXPECT_TRUE(exception);
}

TEST(Operators, MinusEquals) {
  S21Matrix M1, M2;
  M1.FillMatrixRandom();
  M2.FillMatrixRandom();
  S21Matrix M3 = M1;
  M3 -= M2;

  for (int y = 0; y < M3.GetRows(); y++)
    for (int x = 0; x < M3.GetRows(); x++)
      EXPECT_EQ(M3(y, x), M1(y, x) - M2(y, x));

  // матрицы с разными размерами
  bool exception = false;
  try {
    S21Matrix M5(5, 5);
    S21Matrix M6(2, 7);
    M5 -= M6;
  } catch (const std::out_of_range &e) {
    exception = true;
  }
  EXPECT_TRUE(exception);
}

TEST(Operators, MulNumEquals) {
  S21Matrix M1;
  M1.FillMatrixRandom();
  double number = 3.50;
  S21Matrix M2 = M1;
  M1 *= number;

  for (int y = 0; y < M1.GetRows(); y++)
    for (int x = 0; x < M1.GetRows(); x++)
      EXPECT_EQ(M1(y, x), M2(y, x) * number);
}

TEST(Operators, MulMatrixEquals) {
  // [0, 1] * [[0, 1, 2] [3, 4, 5]] == [3, 4, 5]
  S21Matrix M1(1, 2), M2(2, 3);
  M1.FillMatrix();
  M2.FillMatrix();
  M1 *= M2;

  EXPECT_EQ(M1(0, 0), 3);
  EXPECT_EQ(M1(0, 1), 4);
  EXPECT_EQ(M1(0, 2), 5);

  // M3.cols != M4.rows
  bool exception = false;
  try {
    S21Matrix M3(5, 5);
    S21Matrix M4(2, 7);
    M3 *= M4;
  } catch (const std::out_of_range &e) {
    exception = true;
  }
  EXPECT_TRUE(exception);
}

TEST(Operators, Indexation) {
  S21Matrix M1(1, 3);
  M1.FillMatrix();

  EXPECT_EQ(M1(0, 0), 0);
  EXPECT_EQ(M1(0, 1), 1);
  EXPECT_EQ(M1(0, 2), 2);
}

// --------- Accessors и mutators ---------

TEST(Accessors_and_mutators, GetRows) {
  S21Matrix M1;
  int rows = M1.GetRows();
  EXPECT_EQ(rows, 3);
}

TEST(Accessors_and_mutators, GetCols) {
  S21Matrix M1;
  int cols = M1.GetCols();
  EXPECT_EQ(cols, 3);
}

TEST(Accessors_and_mutators, SetValue) {
  S21Matrix M1;
  M1.SetValue(2, 2, 3.50);
  EXPECT_EQ(M1(2, 2), 3.50);
}

TEST(Accessors_and_mutators, EditSize) {
  S21Matrix M1;
  M1.FillMatrix();
  M1.EditSize(2, 2);
  double m1[4] = {0, 1, 3, 4};
  for (int y = 0, n = 0; y < M1.GetRows(); y++)
    for (int x = 0; x < M1.GetCols(); x++, n++) EXPECT_EQ(M1(y, x), m1[n]);

  M1.EditSize(3, 3);
  double m2[9] = {0, 1, 0, 3, 4, 0, 0, 0, 0};
  for (int y = 0, n = 0; y < M1.GetRows(); y++)
    for (int x = 0; x < M1.GetCols(); x++, n++) EXPECT_EQ(M1(y, x), m2[n]);

  M1.EditSize(2, 5);
  double m3[10] = {0, 1, 0, 0, 0, 3, 4, 0, 0, 0};
  for (int y = 0, n = 0; y < M1.GetRows(); y++)
    for (int x = 0; x < M1.GetCols(); x++, n++) EXPECT_EQ(M1(y, x), m3[n]);

  M1.EditSize(5, 2);
  double m4[10] = {0, 1, 3, 4, 0, 0, 0, 0, 0, 0};
  for (int y = 0, n = 0; y < M1.GetRows(); y++)
    for (int x = 0; x < M1.GetCols(); x++, n++) EXPECT_EQ(M1(y, x), m4[n]);

  // Тот же размер
  M1.EditSize(5, 2);
  for (int y = 0, n = 0; y < M1.GetRows(); y++)
    for (int x = 0; x < M1.GetCols(); x++, n++) EXPECT_EQ(M1(y, x), m4[n]);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  std::cout << "Running tests:" << std::endl;
  int ret{RUN_ALL_TESTS()};
  if (!ret)
    std::cout << "Well done." << std::endl;
  else
    std::cout << "Fail." << std::endl;
    
  return 0;
}