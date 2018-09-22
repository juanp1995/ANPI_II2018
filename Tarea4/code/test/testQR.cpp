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
#include "QR.hpp"
#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>
#include <functional>
#include <cmath>

namespace anpi {
  namespace test {


        template<typename T>
        void qrTest(){
	anpi::Matrix<double> a;
	anpi::Matrix<double> q;
	anpi::Matrix<double> r;
	anpi::Matrix<double> qt;
	anpi::Matrix<double> qqt;
	anpi::Matrix<double> qr;
	a={{12,-51,4},{6,167,-68},{-4,24,-41}};
	anpi::Matrix<T> I={{1,0,0},{0,1,0},{0,0,1}};
	anpi::qr(a,q,r);
	T eps = std::numeric_limits<T>::epsilon();
	
	//test que Q*Qt=I
	{
		anpi::get_transpuesta(q,qt);
		
		qqt=anpi::operator*(q,qt);
		anpi::normalizar(qqt);
	
		for (size_t i=0;i<qqt.rows();++i) {
			for (size_t j=0;j<qqt.cols();++j) {
	            		BOOST_CHECK(qqt(i,j)-I(i,j) < eps);
			}
		}  
	}

	//test que R es triangular superior
	{
		for (size_t i=0;i<r.rows();++i) {
			for (size_t j=0;j<r.cols();++j) {
	            		if(i>j){
					BOOST_CHECK(r[i][j] <eps);
				}
			}
		}  
	}

	//test que A=QR
	{
		qr=anpi::operator*(q,r);
		for (size_t i=0;i<qr.rows();++i) {
			for (size_t j=0;j<qr.cols();++j) {
	            		BOOST_CHECK((a(i,j)-qr(i,j)) < 0.0001);
			}
		}  
	}


        }//qrTest

       
    }//test
}//anpi

BOOST_AUTO_TEST_SUITE( QR )

    BOOST_AUTO_TEST_CASE(QR)
    {
        anpi::test::qrTest<float>();
        anpi::test::qrTest<double>();
    }

  
BOOST_AUTO_TEST_SUITE_END()



