/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 04.08.2018
 */

#include <cmath>
#include <limits>
#include <functional>

#include "Exception.hpp"

#ifndef ANPI_ROOT_RIDDER_HPP
#define ANPI_ROOT_RIDDER_HPP

namespace anpi {
  
  /**
   * Find a root of the function funct looking for it starting at xi
   * by means of the secant method.
   *
   * @param funct a functor of the form "T funct(T x)"
   * @param xi initial position
   * @param xii second initial position 
   *
   * @return root found, or NaN if no root could be found
   */
  template<typename T>
   	T rootRidder(const std::function<T(T)>& funct,T xi,T xii,const T eps) {

   		const int MaxIteraciones = numeric_limits<T>::digits;
   		T raiz = 0;
   		T funcion_x3=0;
   		T funcion_pibote=0;
   		T denominador =0;
   		T numerador = 0;
   		T x3=0;
   		T x4=0;
   		T funcion_x1=funct(xi); //valor de la funcion en el extremo inferior
   		T funcion_x2=funct(xii); //valor de la funcion en el extremo superior

   		//comprobar que los signos sean opuestos, para que encierrela raiz
   		if((funcion_x1*funcion_x2<0)){
			for(int i=0;i<MaxIteraciones;++i){

				x3=(xi+xii)/2; //calculo del punto medio numerical
				funcion_x3=funct(x3); //valor de la funcion evaluada en el punto medio
				numerador=(x3-xi)*funcion_x3*(funcion_x1>=funcion_x2?1:-1);
				denominador=sqrt(funcion_x3*funcion_x3-funcion_x1*funcion_x2);
				x4=x3+(numerador/denominador);//vuelvo a calcular x4

				if(abs(x4-raiz)<=eps){
					return raiz;
				}
				raiz=x4;
				funcion_pibote=funct(raiz);

				//comprobar que los signos seas distintos y dependendiendo del caso, cambiar los parametros
				if(((funcion_x3*funcion_pibote)>0 && funcion_x3<0) ||((funcion_x3*funcion_pibote)<0 && funcion_x3>0)){
					xi=x3;
					funcion_x1=funcion_x3;
					xii=raiz;
					funcion_x2=funcion_pibote;
				}

				if(((funcion_x1*funcion_pibote)>0 && funcion_x1<0) ||((funcion_x1*funcion_pibote)<0 && funcion_x1>0)){
					xii=raiz;
					funcion_x2=funcion_pibote;
				}

				if(((funcion_x2*funcion_pibote)>0 && funcion_x2<0) ||((funcion_x2*funcion_pibote)<0 && funcion_x2>0)){
					xi=raiz;
					funcion_x1=funcion_pibote;
				}

				if(abs(xii-xi)<=eps){
					return raiz;
				}
			}//Cierra el for
   		}
			else {
   			//comprobar si alguno de los extremos es raiz
	   		if(funcion_x1==0) return xi;
	   		if(funcion_x2==0) return xii;
   		}
	    // Return NaN if no root was found
	    return std::numeric_limits<T>::quiet_NaN();
    }

}
  
#endif

