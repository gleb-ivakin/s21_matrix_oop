#ifndef __S21MATRIX_H__
#define __S21MATRIX_H__

#include <cmath>
#include <cstring>
#include <iostream>

class S21Matrix {
 private:
  // Атрибуты
  int rows_, cols_;
  double** matrix_;

 public:
  // Конструкторы
  S21Matrix();  // Базовый конструктор, инициализирующий матрицу некоторой
                // заранее заданной размерностью
  S21Matrix(int rows, int cols);  // Параметризированный конструктор с
                                  // количеством строк и столбцов
  S21Matrix(const S21Matrix& other);  // Конструктор копирования
  S21Matrix(S21Matrix&& other) noexcept;  // Конструктор переноса
  ~S21Matrix();                           // Деструктор

  // Методы
  bool EqMatrix(
      const S21Matrix& other);  // Проверяет матрицы на равенство между собой
  void SumMatrix(
      const S21Matrix& other);  // Прибавляет вторую матрицу к текущей
  void SubMatrix(const S21Matrix& other);  // Вычитает из текущей матрицы другую
  void MulNumber(const double num);  // Умножает текущую матрицу на число
  void MulMatrix(const S21Matrix& other);  // Умножает текущую матрицу на вторую
  S21Matrix Transpose();  // Создает новую транспонированную матрицу из текущей
                          // и возвращает ее
  S21Matrix CalcComplements();  // Вычисляет матрицу алгебраических дополнений
                                // текущей матрицы и возвращает ее
  double Determinant();  // Вычисляет и возвращает определитель текущей матрицы
  S21Matrix InverseMatrix();  // Вычисляет и возвращает обратную матрицу

  // Перегрузка операторов
  S21Matrix operator+(const S21Matrix& other);  // Сложение двух матриц
  S21Matrix operator-(
      const S21Matrix& other);  // Вычитание одной матрицы из другой
  S21Matrix operator*(const S21Matrix& other);  // Умножение матриц
  S21Matrix operator*(const double num);  // Умножение матрицы на число
  bool operator==(
      const S21Matrix& other);  // Проверка матриц на равенство (EqMatrix)
  S21Matrix& operator=(
      const S21Matrix& other);  // Присвоение матрице значений другой матрицы
  S21Matrix& operator+=(
      const S21Matrix& other);  // Присвоение сложения (SumMatrix)
  S21Matrix& operator-=(
      const S21Matrix& other);  // Присвоение разности (SubMatrix)
  S21Matrix& operator*=(
      const S21Matrix& other);  // Присвоение умножения (MulMatrix)
  S21Matrix& operator*=(const double num);  // Присвоение умножения (MulNumber)
  double operator()(
      int i, int j);  // Индексация по элементам матрицы (строка, колонка)

  // Accessors и mutators
  int GetRows();
  int GetCols();
  void SetValue(int row, int col, double value);
  void EditSize(int rows, int cols);

  // Вспомогательные
  void
  FillMatrix();  // Заполняет матрицу значениями от 0 до (rows_ * cols_ - 1)
  void FillMatrixRandom();  // Заполняет матрицу случайными значениями
  void PrintMatrix();  // Выводит матрицу в консоль

 private:
  // Вспомогательные
  void FillMatrixByZero();  // Заполняет матрицу нулями
  void Free();  // Освобождает память, выделенную под матрицы
  S21Matrix GetMinor(int oy, int ox);  // Возвращает минор матрицы
  void CopyMatrix(
      const S21Matrix& other);  // Копирует содержимое одной матрицы в другую
};

#endif
