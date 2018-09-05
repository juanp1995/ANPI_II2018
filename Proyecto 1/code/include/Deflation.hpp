/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @author Pablo Alvarado
 * @date   18.08.2018
 */

#ifndef ANPI_DEFLATION_HPP
#define ANPI_DEFLATION_HPP

#include <vector>
#include <type_traits>

#include <boost/type_traits/is_complex.hpp>
#include <boost/math/tools/polynomial.hpp>


namespace anpi {
  namespace bmt=boost::math::tools; // bt as alias for boost::math::tools
  
  /**
   * Deflate polynomial
   *
   * @param[in] poly Input polynomial to be deflated with the provided root
   * @param[in] root Root of poly to be deflated with.
   * @param[out] residuo Residual of the polynomial deflation
   * @return deflated polynomial
   */
template<class T>
bmt::polynomial<T> deflate(const bmt::polynomial<T>& poly, const T& root,
    T& residuo, T tolerance = std::numeric_limits<T>::epsilon()) {

  residuo = poly[poly.degree()];
  bmt::polynomial<T> quotient = poly;
  quotient[poly.degree()] = 0;


  T coef;

  for (int i = poly.degree() - 1; i >= 0; --i) {
    coef = quotient[i];
    quotient[i] = residuo;
    residuo = coef + residuo * root;
  }
  quotient.normalize();
  return quotient;

}

template<class T> //Delfaciona una raiz compleja
void aux_deflate(std::vector<std::complex<T>>& coeficientes, int n,
    std::complex<T>& residuo, const std::complex<T>& root) {

  std::complex<T> pivote;
  for (int i = n - 1; i >= 0; --i) {
    pivote = coeficientes[i];
    coeficientes[i] = residuo;
    residuo = pivote + residuo * root;
  }
}

/**
 * Deflate polynomial with a second order polynomial.
 *
 * The second order polynomial equals x^2 -2 Re(root)x + |root|^2.
 *
 * @param[in] poly Input polynomial to be deflated with the provided root
 * @param[in] root Root of poly to be deflated with.
 * @param[out] residuo Residual of the polynomial deflation
 * @return deflated polynomial
 */
template<class T>
bmt::polynomial<T> deflate2(const bmt::polynomial<T>& poly,
    const std::complex<T>& root, bmt::polynomial<T>& residuo, T tolerance =
        std::numeric_limits<T>::epsilon()) {

  int n = int(poly.size() - 1);
  bmt::polynomial<T> coeficiente(poly);
  std::vector<std::complex<T>> coeficientes_complejos;
  std::complex<T> residuo_complejo;
  std::complex<T> raiz_conjugada;

  //pasar los coeficientes del polinomio a complejos para usarlos en aux_deflate
  for (int i = 0; i < int(poly.size()); ++i) {
    std::complex<T> complejo(poly[i], T(0));
    coeficientes_complejos.push_back(complejo);
  }

  //deflacionar la raiz por 1era vez
  residuo_complejo = coeficientes_complejos[n];
  coeficientes_complejos[n].real(T(0));
  coeficientes_complejos[n].imag(T(0));
  anpi::aux_deflate(coeficientes_complejos, n, residuo_complejo, root);

  //deflacionar la raiz por 2da vez
  raiz_conjugada = std::conj(root);
  residuo_complejo = coeficientes_complejos[n];
  coeficientes_complejos[n].real(T(0));
  coeficientes_complejos[n].imag(T(0));
  anpi::aux_deflate(coeficientes_complejos, n, residuo_complejo,
      raiz_conjugada);

  //se devuelven los coeficientes a reales para retornarlos
  for (int j = 0; j < int(coeficiente.size()); ++j) {
    coeficiente[j] = coeficientes_complejos[j].real();
  }

  residuo = bmt::polynomial<T>(residuo_complejo.real());
  coeficiente.normalize();
  return coeficiente;

}

}



#endif
