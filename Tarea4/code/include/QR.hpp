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
#include "Exception.hpp"
#include <math.h>
#include <cmath>
#include <limits>

using namespace std;
namespace anpi{

/*
 * Imprime una matriz en consola
 */
template<typename T>
void printMatrix(anpi::Matrix<T>& m){
	for(size_t i=0;i<m.rows();++i){
		for(size_t j=0;j<m.cols();++j){
			std::cout<<m[i][j]<<"\t";
		}std::cout<<"\n";
	}
}

/*
 * Imprime un vector en consola
 */
template<typename T>
void printVector(std::vector<T>& v){
	for(size_t i=0;i<v.size();++i){
			printf("%g  ",v[i]);
		}printf("\n");
}
/*
 *Dado std::vector<T>& u retorna la
 *magnitud (tipo float) de u
 */
template<typename T>
float get_magnitud(std::vector<T>& u){
	float m=0;
	for(size_t i=0;i<u.size();++i){
				m=m+u[i]*u[i];
		} return sqrtf(m);
}



/*
 * Dados std::vector<T>& u, const size_t e
 * aplica la formula de Householder y calcula
 * ||u||e, siendo e un vector {1,0,...,0}
 */
template<typename T>
std::vector<T> get_ae(std::vector<T>& u, const size_t e){
	std::vector<T> me;
	me.resize(u.size());


	for(size_t i=0;i<u.size();++i){
		if(i>=e){
			me[i]=0;
		}else{
			me[i]=get_magnitud(u);
		}
	}return me;

}

/*
 * Dados std::vector<T>& a,std::vector<T>& b
 * utiliza la formula de Householder para calcular
 * u=a-||a||e y retorna el vector u
 */
template<typename T>
std::vector<T> get_u(std::vector<T>& a,std::vector<T>& b){
	std::vector<T> c;
	c.resize(a.size());
	for(size_t i=0; i< a.size();++i){
		c[i]=a[i]-b[i];
	}return c;
}

/*
 * Dados const std::vector<T>& u,anpi::Matrix<T>& uut
 * multiplica u por su tranpuesta y retorna la matriz resultante
 */
template<typename T>
void get_uut(const std::vector<T>& u,anpi::Matrix<T>& uut){
	uut.allocate(u.size(),u.size());
	for(size_t i=0; i<u.size();++i){
		for(size_t j=0; j<u.size();++j){
			uut[i][j]=u[i]*u[j];
		}
	}
}

/*
 * Dadas anpi::Matrix<T>& q,anpi::Matrix<T>& uut,const float mag, size_t e
 * calcula Q usando transfomaciones de HouseHolder.
 * q es la matriz a calcular, uut es el vectoru por su transpuesta(matriz)
 * mag es la magnitud de u, e corresponde al numero de iteracion actual
 */
template<typename T>
void get_q(anpi::Matrix<T>& q,anpi::Matrix<T>& uut,const float mag,size_t e){

	for(size_t i=0;i<uut.rows();++i){
		for(size_t j=0;j<uut.cols();++j){
			if(i!=j ){
				q[i+e][j+e]=(-2/(mag*mag))*uut[i][j];

				}else if(i==j){
					q[i+e][j+e]=1+(-2/(mag*mag))*uut[i][j];

				}
			}
	}
	if(e>0){
			for (size_t i=0;i<uut.rows()+e;++i){
				for(size_t j=0;j<uut.cols()+e;++j){
					if(i<e || j<e){
						q[i][j]=0;
						}if(i==j && i<e){
							q[i][j]=1;
						}
				}
			}
		}
}

/*
 * Dadas anpi::Matrix<T>& q,const size_t e
 * inserta 0 en una columna apartir de la fila e
 */
template<typename T>
void normalizar(anpi::Matrix<T>& q,const size_t e){
	for(size_t i=0;i<q.rows();++i){
		for(size_t j=0;j<q.rows();++j){
			if(i>j && j<e){
				q[i][j]=0;
			}
		}
	}
}

/*
*Dada anpi::Matrix<T>& q
*la recore y coloca ceros en los lugares donde 
*el numero es menor que 0.0000001
*/
template<typename T>
void normalizar(anpi::Matrix<T>& q){

	for(size_t i=0;i<q.rows();++i){
		for(size_t j=0;j<q.rows();++j){
			if((i!=j) &&(q[i][j]<0.0000001)){
				q[i][j]=0;
			}
		}
	}
}

/*
 * Toma anpi::Matrix<T>& qa y le reduce e filas y e columnas
 * la nueva matriz se almacena en anpi::Matrix<T>& ap
 */
template<typename T>
void reducir(anpi::Matrix<T>& qa,anpi::Matrix<T>& ap,size_t e){
	ap.allocate(qa.rows()-(e+1),qa.cols()-(e+1));
	for(size_t i=0;i<ap.rows();++i){
		for(size_t j=0; j<ap.cols();++j){
			ap[i][j]=qa[i+e+1][j+e+1];
		}

	}
}

/*
 * Dadas const anpi::Matrix<T>& m,anpi::Matrix<T>& trans
 * toma la matriz m la transpone y guarda el resultado en trans
 *
 */
template<typename T>
void get_transpuesta(const anpi::Matrix<T>& m,anpi::Matrix<T>& trans){
	trans.allocate(m.cols(),m.rows());
	for(size_t i=0;i<m.rows();++i){
		for(size_t j=0; j<m.cols();++j){
			trans[j][i]=m[i][j];
		}
	}
}



/* 
*Dadas const anpi::Matrix<T>& A,anpi::Matrix<T>& Q,anpi::Matrix<T>& R
*toma cada columna de Q y usando las tranformaciones de Householder calcula un Qn, 
*luego lo transpone ylo multiplica por a, ese resultado lo inserta  en R la cual es una matriz 
*triangular superior, estas transformaciones se hacen hasta Qn-1, y la multiplicacion de todos los Q
*da como resultado R
*/

template<typename T>
void qr ( const anpi::Matrix<T>& A,
		anpi::Matrix<T>& Q,
		anpi::Matrix<T>& R ) {

	//Variables necesarias para el computo
	const size_t i = A.rows();
	std::vector<T> u;
	std::vector<T> a;
	std::vector<T> b;
	anpi::Matrix<T> uut;
	anpi::Matrix<T> QA;
	anpi::Matrix<T> a_pivot = A;
	anpi::Matrix<T> q_pivot;
	anpi::Matrix<T> Qt;
	Q.allocate(A.rows(),A.cols());
	q_pivot.allocate(A.rows(),A.cols());
	QA.allocate(A.rows(),A.cols());

	for(size_t ii=0;ii<i-1;++ii){
		// calculo del vector a:n
		a.resize(a_pivot.rows());
		for(size_t x=0;x<i;++x){
			a[x]=a_pivot[x][0];
		}

		//calculo el vector ||a||e0
		b.resize(a_pivot.rows());
		b=anpi::get_ae(a,1);

		//calculo el vector u
		u=anpi::get_u(a,b);

		//Calculo de la matriz uut
		anpi::get_uut(u,uut);

		//Calculo de la matriz Q
		anpi::get_q(q_pivot,uut,get_magnitud(u),ii);
			if(ii>0){
				Q=anpi::operator *(Q,q_pivot);
			}else{
				Q=q_pivot;
			}

		//calculo de QA
		QA=anpi::operator *(q_pivot,A);

		//coloco 0 en las filas inferiores ii
		anpi::normalizar(QA,ii+1);

		//reduzco QA para usarla como nueva matriz
		anpi::reducir(QA,a_pivot,ii);
	}// fin del for

	anpi::get_transpuesta(Q,Qt);//Calculo la transpuesta de Q
	R=anpi::operator *(Qt,A); //Calculo Qt*A
	anpi::normalizar(R,i); //Coloco ceros en R

	//cout<<"\n"<<"Beautiful R "<<endl;
	//anpi::printMatrix(R);
	//cout<<"\n"<<"Acabo de hacer una descomposicion QR"<<endl;

	}//fin de qr


}//fin namesapace anpi
