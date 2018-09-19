/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: 
 * @Date  : 03.03.2018
 */

#include <cmath>
#include <limits>
#include <functional>

#include "Exception.hpp"
#include "Matrix.hpp"

#ifndef ANPI_LU_CROUT_HPP
#define ANPI_LU_CROUT_HPP

namespace anpi {

  /**
   * Auxiliary method used to debug LU decomposition.
   *
   * It separates a packed LU matrix into the lower triangular matrix
   * L and the upper triangular matrix U, such that the diagonal of U
   * is composed by 1's.
   */
  template<typename T>
  void unpackCrout(const Matrix<T>& LU,
                   Matrix<T>& L,
                   Matrix<T>& U) {

    if (LU.rows() != LU.cols()) throw anpi::Exception("Matrix is not a LU decomposition");

    for (size_t i = 0; i < LU.rows(); ++i) {
      for (size_t j = 0; j < LU.cols(); ++j) {
        if (i == j) {
          L[i][j] = LU[i][j];
          U[i][j] = T(1);
        }
        else if (j > i) {
          L[i][j] = T(0);
          U[i][j] = LU[i][j];
        }
        else {
          L[i][j] = LU[i][j];
          U[i][j] = T(0);
        }
      } // for j
    } // for i

  }
  
  /**
   * Decompose the matrix A into a lower triangular matrix L and an
   * upper triangular matrix U.  The matrices L and U are packed into
   * a single matrix LU.  
   *
   * Crout's way of packing assumes a diagonal of
   * 1's in the U matrix.
   *
   * @param[in] A a square matrix 
   * @param[out] LU matrix encoding the L and U matrices
   * @param[out] permut permutation vector, holding the indices of the
   *             original matrix falling into the corresponding element.
   *             For example if permut[5]==3 holds, then the fifth row
   *             of the LU decomposition in fact is dealing with the third
   *             row of the original matrix.
   *
   * @throws anpi::Exception if matrix cannot be decomposed, or input
   *         matrix is not square.
   */
  template<typename T>
  void luCrout(const Matrix<T>& A,
               Matrix<T>& LU,
               std::vector<size_t>& permut) {

    if (A.rows() != A.cols()) throw anpi::Exception("Matrix has to be square");

    size_t n = A.rows();
    size_t i, j, k, imax = 0;
    T pivot, temp, sum = T(0.0);
    LU = A;

    for (size_t m = 0; m < n; ++m) {
      permut[m] = m;
    }

    for (j = 0; j < n; ++j) {
      pivot = T(0.0);
      //Search for the pivot
      for (i = j; i < n; ++i) {
        if ((temp = std::fabs(LU[i][j])) >= pivot) {
          pivot = temp;
          imax = i;
        }
      }
      if (pivot == T(0.0)) throw anpi::Exception("Singular matrix");

      // Interchange rows
      if (j != imax) {
        for (k = 0; k < n; ++k) {
          temp = LU[imax][k];
          LU[imax][k] = LU[j][k];
          LU[j][k] = temp;
        }
        permut[imax] = j;
        permut[j] = imax;
      }


      if (j == 0) {
        for (k = j + 1; k < n; ++k)
          LU[j][k] = LU[j][k] / pivot;
      }

      else if (j > 0 && j < n - 1) {
        for (i = j; i < n; ++i) {
          sum = LU[i][j];
          for (k = 0; k < j; ++k)
            sum -= LU[i][k] * LU[k][j];
          LU[i][j] = sum;
        }

        for (k = j + 1; k < n; ++k) {
          sum = LU[j][k];
          for (i = 0; i < j; ++i)
            sum -= LU[j][i] * LU[i][k];
          LU[j][k] = sum / pivot;
        }
      }

      else {
        sum = LU[n - 1][n - 1];
        for (k = 0; k < n - 1; ++k)
          sum -= LU[n - 1][k] * LU[k][n - 1];
        LU[n - 1][n - 1] = sum;
      }
    }
  }

}
  
#endif

