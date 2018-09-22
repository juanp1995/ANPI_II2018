/**
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 */

#include <boost/test/unit_test.hpp>

#include "LUCrout.hpp"
#include "LUDoolittle.hpp"
#include "SolveLU.hpp"

#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>

#include <functional>

#include <cmath>

namespace anpi {
  namespace test {
  
  	enum UnpackEnum {
  		unpackCrout, 
  		unpackDoolittle
  		};
    
    /// Test the given closed root finder
    template<typename T>
    void luTest(const std::function<void(const Matrix<T>&,
                                         Matrix<T>&,
                                         std::vector<size_t>&)>& decomp,
                const std::function<void(const Matrix<T>&,
                                         Matrix<T>&,
                                         Matrix<T>&)>& unpack) {

      // The result
      Matrix<T> LU;

      // Test if a non-square matrix is successfully detected
      {
        Matrix<T> A = {{1,7,6,4},{2,17,27,17}};
        std::vector<size_t> p;
        try {
          decomp(A,LU,p);
          BOOST_CHECK_MESSAGE(false,"Rectangular matrix not properly catched");
        }
        catch(anpi::Exception& exc) {
          BOOST_CHECK_MESSAGE(true,"Rectangular matrix properly detected");
        }
      }

      // Test pivoting
      {
        anpi::Matrix<T> A = { {-1,-2,1,2},{ 2, 0,1,2},{-1,-1,0,1},{ 1, 1,1,1} };
        std::vector<size_t> p;
        decomp(A,LU,p);

        std::vector<size_t> gp= {1,0,3,2};
        BOOST_CHECK(gp==p);
      }
      
      // Test decomposition
      {
        // same matrix as before, but already permuted to force a
        // clean decomposition
        anpi::Matrix<T> A = { { 2, 0,1,2},{-1,-2,1,2},{ 1, 1,1,1},{-1,-1,0,1} };
        std::vector<size_t> p;
        decomp(A,LU,p);
        Matrix<T> L,U;
        unpack(LU,L,U);
        Matrix<T> Ar=L*U;

        const T eps = std::numeric_limits<T>::epsilon();

        BOOST_CHECK(Ar.rows()==A.rows());
        BOOST_CHECK(Ar.cols()==A.cols());

        for (size_t i=0;i<Ar.rows();++i) {
          for (size_t j=0;j<Ar.cols();++j) {
            BOOST_CHECK(std::abs(Ar(i,j)-A(i,j)) < eps);
          }
        }
      }
    }//luTest
    
    
    /*
     * Test the given unpack method
     */
    template<typename T>
    void unpackTest(const std::function<void(const Matrix<T>&,
    																				Matrix<T>&,
    																				Matrix<T>&)>& unpack,
										const UnpackEnum method) {
			
			//Result							
			anpi::Matrix<T> L, U;	
			
			//Test if a non-square matrix is successfully detected 
			{
				anpi::Matrix<T> LU = { {1,2,3,4}, {5,6,7,8} };
				try{
					unpack(LU, L, U);
					BOOST_CHECK_MESSAGE(false, "Rectangular matrix not properly catched");
				}
				catch(anpi::Exception& exc){
					BOOST_CHECK_MESSAGE(true, "Rectangular matrix properly detected");
				}
			}
			
			//Test lower and upper triangular matrix
			{
				anpi::Matrix<T> LU = { {1,2,3,4}, {5,6,7,8}, {9,10,11,12}, {13,14,15,16} };
				unpack(LU,L,U);
								
				for(size_t i=0; i<LU.rows(); ++i){
					for(size_t j=0; j<LU.cols(); ++j){
						//Test unitary diagonal
						if(i==j){
							if(method == UnpackEnum::unpackCrout)
								BOOST_CHECK(U[i][j] == T(1));
							else
								BOOST_CHECK(L[i][j] == T(1));
						}
						
						else if(j>i)//Test lower triangular matrix
							BOOST_CHECK( L[i][j] == T(0) );
						
						else //Test upper triangular matrix
							BOOST_CHECK( U[i][j] == T(0) );
					}//for j
				}// for i
			}
										
		}//unpackTest

  } // test
}  // anpi

BOOST_AUTO_TEST_SUITE( LU )

BOOST_AUTO_TEST_CASE(Doolittle) 
{
  anpi::test::luTest<float>(anpi::luDoolittle<float>,
                            anpi::unpackDoolittle<float>);
  anpi::test::luTest<double>(anpi::luDoolittle<double>,
                             anpi::unpackDoolittle<double>);
}

BOOST_AUTO_TEST_CASE(Crout) 
{
  anpi::test::luTest<float>(anpi::luCrout<float>,anpi::unpackCrout<float>);
  anpi::test::luTest<double>(anpi::luCrout<double>,anpi::unpackCrout<double>);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( Unpack )

BOOST_AUTO_TEST_CASE(UnpackCrout)
{
	anpi::test::unpackTest<float>(anpi::unpackCrout<float>, anpi::test::UnpackEnum::unpackCrout);
	anpi::test::unpackTest<double>(anpi::unpackCrout<double>, anpi::test::UnpackEnum::unpackCrout);
}

BOOST_AUTO_TEST_CASE(UnpackDoolittle)
{
	anpi::test::unpackTest<float>(anpi::unpackDoolittle<float>, anpi::test::UnpackEnum::unpackDoolittle);
	anpi::test::unpackTest<double>(anpi::unpackDoolittle<double>, anpi::test::UnpackEnum::unpackDoolittle);
}

BOOST_AUTO_TEST_SUITE_END()

