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

#ifndef ANPI_ROOT_SECANT_HPP
#define ANPI_ROOT_SECANT_HPP

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
  T rootSecant(const std::function<T(T)>& funct,T xi,T xii,const T eps) {

    T xPrev = xi;
    T x0 = xii;
    T error = 0;
    T xNext = 0;
    T maxIterations = std::numeric_limits<T>::digits;

    for (int i = maxIterations; i > 0; --i) {

      // Approximation of the derivative using backward difference
      T fdx = (funct(xPrev) - funct(x0)) / (xPrev - x0);

      // Calculates a new approximation of the root
      xNext = x0 - funct(x0) / fdx;
      error = xNext - x0;

      // Return the approximated value of the root
      if(std::abs(error) < eps){
        return xNext;
      }

      xPrev = x0;
      x0 = xNext;
    }

    // Return NaN if no root was found
    return std::numeric_limits<T>::quiet_NaN();
  }

}
  
#endif

