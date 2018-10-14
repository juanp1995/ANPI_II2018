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

#ifndef ANPI_ROOT_BISECTION_HPP
#define ANPI_ROOT_BISECTION_HPP

namespace anpi {
  
  /**
   * Find the roots of the function funct looking for it in the
   * interval [xl,xu], using the bisection method.
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
T rootBisection(const std::function<T(T)>& funct,T xl,T xu,const T eps) {

  const int MaxIterationValue = std::numeric_limits<T>::digits;                   //Valor máximo de iteraciones permitidas.

  if(xl > xu){
        throw anpi::Exception("El intervalo se encuentra invertido");       
  }    

  T distancia,imagen_1, imagen_2, x_Medio, raiz;

  imagen_1 = funct(xl);
  imagen_2 = funct(xu);
  
  if(imagen_1 * imagen_2 > T(0)){
  	throw anpi::Exception("Las preimágenes no tienen imágenes con signos opuestos");       
  }
  raiz = imagen_1 < T(0) ? (distancia = xu - xl,xl) : (distancia = xl - xu,xu);   //Si la imagen de xl es negativa, la distancia
                                                                                  //será positiva para ir sumando y la raíz se ubica
                                                                           //en xl, para ir acercándose al 0. Igual al revés.
  for (int i = 1 ; i <= MaxIterationValue ; i++){

    imagen_2 = funct(x_Medio = raiz+(distancia *= static_cast<T>(0.5)));        //Se divide la distancia a la mitad. Se le suma a la raíz.
                                                                                  //Se le asigna el resultado a x_Medio y después se evalúa.
    if (imagen_2 <= T(0))
      raiz = x_Medio;
    if (std::abs(distancia) < eps || imagen_2 == T(0))
      return raiz; 
  }
  if(std::abs(funct(x_Medio)) < eps) return x_Medio;                                        //En caso de que se alcance el máximo de iteraciones y.
                                                                                  //el valor de la función es menor al épsilon.
  return std::numeric_limits<T>::quiet_NaN();                                     //En caso de que no se encuentra, retorna NaN.
  }

}
  
#endif

