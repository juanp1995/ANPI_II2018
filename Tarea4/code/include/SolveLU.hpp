#include <cmath>
#include <limits>
#include <functional>

#include <math.h>

#include "Exception.hpp"
#include "Matrix.hpp"
#include "LU.hpp"

#define TINY 1e-15

#ifndef ANPI_SOLVE_LU_HPP
#define ANPI_SOLVE_LU_HPP

namespace anpi {

  /**
   * Swaps the solution vector using the same permutation vector of the decomposition LU
   *
   * The solution vector has to swap because the result of the decomposition
   * is not A = LU, instead the result is PA = LU; where the vector P contains
   * the permutations done as result of the partial pivoting
   */
  template<typename T>
  void permuteVector(std::vector<T>& vector, const std::vector<size_t>& permut) {

    if (vector.size() != permut.size())
      throw anpi::Exception("Vectors has to be of the same size");

    std::vector<size_t> permutTemp = permut;

    for (size_t i = 0; i < vector.size(); ++i) {
      size_t indx = permutTemp[i];
      if (indx != i) {
        T temp = vector[i];
        vector[i] = vector[indx];
        vector[indx] = temp;

        permutTemp[i] = i;
        permutTemp[indx] = indx;
      }
    }
  } // permutVector

  /**
   * Finds a vector "y" such that satisfies " Ly=b "
   * using forward substitution.
   *
   * @param[in] m Lower triangular matrix of coefficients
   * @param[out] y Vector to be calculated
   * @param[in] b Vector with the result values
   *
   * @throws anpi::Exception if the size of the vector with the
   *         result values is different than the number rows
   *         of the matrix
   */
  template<typename T>
  void forwardSubs(const anpi::Matrix<T>& m, std::vector<T>& y, const std::vector<T>& b) {

    if (m.rows() != b.size()) throw anpi::Exception("Forward substitution can't be done");

    size_t n = m.rows();
    T sum, diag = T(0.0);
    y.resize(n);

    for (size_t i = 0; i < n; ++i) {

      if (i == 0) y[i] = b[i] / m[i][i];
      else {
        diag = m[i][i];
        if (diag == T(0)) diag = TINY; //Avoid division by 0
        sum = b[i];
        for (size_t j = 0; j < i; ++j)
          sum -= m[i][j] * y[j];
        y[i] = sum / diag;
      }
    }
  }

  /**
   * Finds a vector "x" such that satisfies " Ux=y "
   * using backward substitution.
   *
   * @param[in] m Upper triangular matrix of coefficients
   * @param[out] x Vector to be calculated
   * @param[in] y Vector with the result values
   *
   * @throws anpi::Exception if the size of the vector with the
   *         result values is different than the number rows
   *         of the matrix
   */
  template<typename T>
  void backwardSubs(const anpi::Matrix<T>& m, std::vector<T>& x, const std::vector<T>& y) {

    if (m.rows() != y.size()) throw anpi::Exception("Backward substitution can't be done");

    size_t n = m.rows();
    T sum, diag = T(0.0);
    x.resize(n);
    size_t i = n - 1;

    for (auto it = x.rbegin(); it != x.rend(); ++it) {
      if (i == (n - 1))
        x[i] = y[i] / m[i][i];
      else {
        diag = m[i][i];
        if (diag == T(0)) diag = TINY; //Avoid division by 0
        sum = y[i];
        for (size_t j = (i + 1); j < n; ++j)
          sum -= m[i][j] * x[j];
        x[i] = sum / diag;
      }
      --i;
    }
  }


  /**
   * Solves a linear equation system using LU decomposition
   *
   * @param[in] A Square matrix of coefficients
   * @param[out] x Solution vector of the system
   * @param[in] b Vector with the result values of the equations
   *
   * @throws anpi::Exception if the matrix is not square or if
   *         there are less result values in vector b than equations
   *         in the system
   */
  template<typename T>
  bool solveLU(const anpi::Matrix<T>& A, std::vector<T>& x, const std::vector<T>& b) {

    if (A.rows() != A.cols()) throw anpi::Exception("Matrix has to be square");
    if (b.size() != A.cols()) throw anpi::Exception("Equations system can't be solve");

    anpi::Matrix<T> LU;
    std::vector<size_t> permut;
    std::vector<T> bPermut = b;
    anpi::Matrix<T> L, U;
    
    //Try to decompose the matrix
    try {
      lu(A, LU, permut);
    } catch (anpi::Exception& exc) {
      return false;
    }

    unpackDoolittle(LU, L, U);
    permuteVector(bPermut, permut);

    std::vector<T> y;
    //Solves Ly=b
    forwardSubs<T>(L, y, bPermut);

    //Solves Ux=y
    backwardSubs<T>(U, x, y);
    return true;
  } // solveLU

}

#endif
