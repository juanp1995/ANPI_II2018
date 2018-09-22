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
#include <algorithm>

#include "Exception.hpp"
#include "Matrix.hpp"

#ifndef ANPI_LU_DOOLITTLE_HPP
#define ANPI_LU_DOOLITTLE_HPP

namespace anpi {


  /**
   * Auxiliary method used to debug LU decomposition.
   *
   * It separates a packed LU matrix into the lower triangular matrix
   * L and the upper triangular matrix U, such that the diagonal of L
   * is composed by 1's.
   */
  template<typename T>
  void unpackDoolittle(const Matrix<T>& LU,
                       Matrix<T>& L,
                       Matrix<T>& U) {
	if(LU.rows() != LU.cols())
		throw anpi::Exception("The matrix LU is not a square matrix");

	L.allocate(LU.rows(), LU.rows());
   	U.allocate(LU.rows(), LU.rows());

	std::size_t matrixSize = LU.rows();
	for(std::size_t i = 0 ; i < matrixSize; i++){
		for (std::size_t j = 0 ; j < matrixSize; j++){
			if(i > j){
				L[i][j] = LU[i][j];
				U[i][j] = T(0);
			}
			else if (i == j){
				L[i][j] = T(1);
				U[i][j] = LU[i][j];
			}
			else{
				U[i][j] = LU[i][j];
				L[i][j] = T(0);
			}
		}
	}
  }
  
  /**
   * Decompose the matrix A into a lower triangular matrix L and an
   * upper triangular matrix U.  The matrices L and U are packed into
   * a single matrix LU. 
   *
   * The L matrix will have in the Doolittle's LU decomposition a
   * diagonal of 1's
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
  void luDoolittle(const Matrix<T>& A,
                   Matrix<T>& LU,
                   std::vector<size_t>& permut) {

	if(A.rows() != A.cols())					//First corroboration.
		throw anpi::Exception("Matrix not square");
	const std::size_t matrixSize = A.rows();

	permut.resize(matrixSize);
	LU = A;

	for (std::size_t i = 0 ; i < matrixSize ; ++i) 			//Filling the permutation vector.
		permut[i] = T(i);
	
									//Copy of the original matrix to one that can be worked on.
	for (std::size_t i = 0 ; i < matrixSize-1; ++i){ 		//Iterates through rows. Main for.
		std::size_t j = i; 					//It works diagonally.
		T pivot = LU[i][j];					//The starting pivot always starts in an element of the diagonal.
		std::size_t row_pivot = i;				//It's necessary to know in which row is the pivot.

		for (std::size_t i_temp = i ; i_temp < matrixSize; ++i_temp){
			/*Finds the biggest absolute value in the column
		 		* and sets it as the pivot.*/		
			if(std::abs(LU[i_temp][j]) > std::abs(pivot)){ 	//A biggest value is found.
				pivot = LU[i_temp][j];			//Changes the pivot.
				row_pivot = i_temp;			//Sets the row of the element found as the new row of the pivot.
				permut[i] = permut[row_pivot];  				//Changes the permutation vector.
				permut[row_pivot] = i;
			}
		}

		if(pivot == T(0))
			throw anpi::Exception("Singular matrix");
		swapRows(LU, i,row_pivot,j);				//Swaps the row with the old pivot with the row of the new one.
		for (std::size_t i2 = i ; i2 < matrixSize -1 ; ++i2){ 	//Calculates the factors according to the row.
		T factor = LU[i2+1][j]/LU[i][j];						//Calculates the factor.
			for (std::size_t j2 = i ; j2 < matrixSize ; ++j2) 	//Executes the changes in the U matrix with the right factor.
				LU[i2+1][j2] = LU[i2+1][j2]-LU[i][j2]*factor;	//Executes the changes.
			LU[i2+1][j] = factor;								//Sets the factor calculated in the L matrix.
		}
	}
  }
	template<class T> 
	void swapRows(Matrix<T>& matriz, std::size_t row1, std::size_t row2, std::size_t column){ //Swap two rows given in the given matrix.

		for (std::size_t j = column; j < matriz.rows(); j++){
			T temporal = matriz[row1][j];
			matriz[row1][j] = matriz[row2][j];
			matriz[row2][j] = temporal;
		}
	}

}
  
#endif

