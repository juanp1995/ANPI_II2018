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

#ifndef ANPI_ROOT_INTERPOLATION_HPP
#define ANPI_ROOT_INTERPOLATION_HPP

namespace anpi {
  
  /**
   * Find the roots of the function funct looking for it in the
   * interval [xl,xu], by means of the interpolation method.
   *
   * @param funct a functor of the form "T funct(T x)"
   * @param xl lower interval limit
   * @param xu upper interval limit
   *
   * @return root found, or NaN if none could be found.
   *
   * @throws anpi::Exception if inteval is reversed or both extremes
   *         have same sign.
   */
  template<typename T>
  T rootInterpolation(const std::function<T(T)>& funct,T xl,T xu,const T eps) {

    if(xl > xu){
        throw anpi::Exception("El intervalo se encuentra invertido");       
    }
    
    const int MaxIterationValue = std::numeric_limits<T>::digits;       //Valor máximo de iteraciones que puede haber.       			

		T x1, x2, y1, y2, distancia_X, del, raiz, y_raiz;
    
    x1 = xl;
    x2 = xu;
    y1 = funct(x1);
    y2 = funct(x2);

    if(y1 * y2 > 0){
        throw anpi::Exception("Las preimágenes no tienen imágenes con signos opuestos");       
    }
    
    distancia_X = x2 - x1;

    for (int i = 1 ; i <= MaxIterationValue ; i++){

      raiz = x1+(distancia_X*y1)/(y1-y2);
      y_raiz = funct(raiz);

      if (y_raiz < T(0)){
        del = x1 - raiz;
        x1 = raiz;
        y1 = y_raiz;
      }

      else{
        del = x2 - raiz;
        x2 = raiz;
        y2 = y_raiz;
      }

      distancia_X = x2 - x1;
      if (std::abs(del) < eps || y_raiz == T(0)) return raiz;
    }
 
    if(std::abs(funct(raiz)) < eps) return raiz;                            //Si el valor de la imagen es menor a épsilon, la retorna.    
    // Return NaN if no root was found
    return std::numeric_limits<T>::quiet_NaN();                         //Si no se encuentra el número retorna un NaN.
  }

}
  
#endif

