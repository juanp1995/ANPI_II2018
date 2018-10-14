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

#ifndef ANPI_NEWTON_RAPHSON_HPP
#define ANPI_NEWTON_RAPHSON_HPP

namespace anpi {
  
  /**
   * Find the roots of the function funct looking by means of the
   * Newton-Raphson method
   *
   * @param funct a functor of the form "T funct(T x)"
   * @param xi initial root guess
   * 
   * @return root found, or NaN if none could be found.
   *
   * @throws anpi::Exception if inteval is reversed or both extremes
   *         have same sign.
   */
  template<typename T>
  T rootNewtonRaphson(const std::function<T(T)>& funct,T xi,const T eps) {

    // Step size in terms of eps
    T h = eps * static_cast<T>(0.5);
    T maxIterations = std::numeric_limits<T>::digits;
    T xNext = 0;
    T error = 0;

    for (int i = maxIterations; i > 0; --i) {

      // Approximation of the derivative using central difference
      T derivative = (funct(xi + h) - funct(xi - h)) / (T(2) * h);

      // Calculate the next approximation of the root
      xNext = xi - (funct(xi) / derivative);
      error = xNext - xi;

      if (std::abs(error) < eps) {
        break;
      }

      xi = xNext;
    }

    // Check if the method converged to a root, according to the
    // tolerance specified
    if (std::abs(funct(xNext)) < eps) {
      return xNext;
    }

    // Return NaN if no root was found
    return std::numeric_limits<T>::quiet_NaN();
    }

}
  
#endif
