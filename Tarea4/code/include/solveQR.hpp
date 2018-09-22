/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Emilly Sancho Murillo
 * @Date  : 15.09.2018
 */


#include <iostream>
#include "Matrix.hpp"
#include <math.h>
#include <cmath>
#include <limits>
#include "QR.hpp"

using namespace std;
namespace anpi{


template<typename T>
void solveQR2(const anpi::Matrix<T>& A, std::vector<T>& x, const std::vector<T>& b){
	//Variables necesarias para el computo
	anpi::Matrix<T> Q;
	anpi::Matrix<T> Qt;
	anpi::Matrix<T> R;
	std::vector<T> bprima;
	size_t jr;//para almacenar el # de columnas de R
	T acumulado;

	anpi::qr(A,Q,R);//Descomposicion de A
	anpi::get_transpuesta(Q,Qt);//Calculo de la trasnpuesta de Q

	bprima=anpi::operator *(Qt,b);// bpivote=Qt*b
	jr=R.cols();
	x.resize(bprima.size());
	x[jr-1]=(bprima[jr-1]/R[jr-1][jr-1]);//sustitucion hacia atras x=b/R, para x[n-1]

	//cout<<"antes del for de solveQR2 todo bien "<<endl;
	for(size_t j=jr-2;j>0;--j){
		//cout<<"j "<<j<<endl;
		acumulado=bprima[j];
		for(size_t jj=j+1;jj<R.cols();++jj){
			acumulado-=(R[j][jj]*x[jj]);
			//cout<<" jj"<<jj<<"\n"<<endl;
		}
		x[j]=(T(1)/R[j][j])*acumulado;
	}
}


/*
*Dadas const anpi::Matrix<T>& A, std::vector<T>& x, const std::vector<T>& b
*toma A y los descompone usando qr, utiliza Qt para calcular un b prima y este lo 
*sustituye hacia atras, que sea "igual" al b y cuando esto se cumple el x econtrado en solveQR2 es la
*respuesta del sistema.
*/
template<typename T>
bool solveQR(const anpi::Matrix<T>& A,
		std::vector<T>& x, const std::vector<T>& b){
	T epsilon = std::numeric_limits<T>::epsilon();

	//Variables necesarias para el computo
	std::vector<T> bpivote;
	std::vector<T> bprima;
	std::vector<T> berror;
	std::vector<T> xerror;
	bpivote=b;
	int ultimo_intento=1000;

	float bmagnitud=anpi::get_magnitud(bpivote);

	//cout<<"aqui casual ";
	anpi::solveQR2(A,x,b); //primer calculo de x
	//cout<<"after casual ";
	while(bmagnitud>epsilon){
		bprima=anpi::operator *(A,x);

		for(size_t i=0;i<berror.size();++i){
			berror[i]=bprima[i]-b[i];
		}


		bmagnitud=anpi::get_magnitud(berror);
		if(bmagnitud <=epsilon){
			return true;
		}anpi::solveQR2(A,xerror,berror);

		ultimo_intento++;
		if(ultimo_intento<=0){
			return false;
		}
	}return false;

}


}//fin namesapace anpi

