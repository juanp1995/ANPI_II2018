/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 */

#include <cmath>
#include <limits>
#include <functional>

#include "Exception.hpp"

#ifndef ANPI_ROOT_BRENT_HPP
#define ANPI_ROOT_BRENT_HPP
using namespace std;
namespace anpi {
  
  /**
   * Find the roots of the function funct looking for it in the
   * interval [xl,xu], using the Brent's method.
   *
   * @param funct a std::function of the form "T funct(T x)"
   * @param xl lower interval limit
   * @param xu upper interval limit
   *
   * @return root found, or NaN if none could be found.
   *
   * @throws anpi::Exception if inteval is reversed or both extremes
   *         have same sign.
   */
	template<typename T>
	T rootBrent(const std::function<T(T)>& funct,T xl,T xu,const T eps) {

		//-------------AQUI empieza mi codigo ----------------------
		const int MaxIteraciones = numeric_limits<T>::digits;
		T raiz = T(0);//El valor que ando buscando
		T funcion_inferior = funct(xl); //La funcion evaluada en el limite inferior
		T funcion_superior = funct (xu); //La funcion evaluada en el limite superior

		//reviso si los terminos del intervalo estan ordenados
		if(xu<xl){
			throw anpi::Exception("Inverta los valores del intervalo");
		}

		//primero hay que ver si el intervalo encierra una raiz, usando la logica de biseccion (signos)
		if(funcion_superior*funcion_inferior>0){
			throw Exception("No hay raiz en el intervalo");
		}

		//revisar si alguno de los limites es una raiz
		if(abs(funcion_inferior)<eps){
			return xl;
		}
		if(abs(funcion_superior)<eps){
			return xu;
		}

		/*Primero intento interpolacion inversa cuadratica
		 *para hacer la interpolacion necesito tres puntos y aqui solo tengo 2,
		 *me voy a sacar uno de la manga.*/
		T x_medio = (xu+xl)/2;
		T funcion_media = funct(xu); 
		T funcion_raiz;

		for(int i=0;i<MaxIteraciones;++i){
			if(abs(xu-xl)<(eps)){//condicion de parada para la biseccion, lecc 05 pag 19
				return raiz;
			}

			//intento calcular la raiz por interpolacion
			if(funcion_inferior!=funcion_media && funcion_media!=funcion_superior){//no puede tener 2 o mas y iguales (y=f(x))
			//Formula lecc 06, pag 28
			raiz = (xl*funcion_superior*funcion_media)/((funcion_inferior-funcion_superior)*(funcion_inferior-funcion_media))
						 	+ (x_medio*funcion_inferior*funcion_superior)/((funcion_media-funcion_inferior)*(funcion_media-funcion_superior))
							+(xu*funcion_inferior*funcion_media)/((funcion_superior-funcion_inferior)*(funcion_superior-funcion_media));

			}
			else { //si no cumple para interpolacion se calcula por secante
				raiz= xu-(funcion_superior*(xl-xu))/(funcion_superior-funcion_inferior);//Formula lecc 06, pag 20
			}
			if(funcion_inferior==funcion_media && funcion_media==funcion_superior
					&& abs(funct(raiz))<(eps-numeric_limits<T>::epsilon()) ){//si no se puede con interpolacion o secante, se calcula por biseccion
						raiz = (xl+xu)/T(2);
			}

			funcion_raiz = funct(raiz);
			x_medio=xu;
			funcion_media=funcion_superior;

			if(abs(funcion_inferior)<abs(funcion_superior)){//Cuando hay que invertir los limites
				swap(xl,xu);
				swap(funcion_inferior,funcion_superior);
			}

			if(funcion_inferior*funcion_raiz<0){//reacomodo de signos
				xu=raiz;
				funcion_superior=funcion_raiz;
			}
			else {
				xl=raiz;
				funcion_inferior=funcion_raiz;
			}

		}//Cierra el for
		return numeric_limits<T>::quiet_NaN();
	}//Cierra el rootBrent
}
  
#endif

