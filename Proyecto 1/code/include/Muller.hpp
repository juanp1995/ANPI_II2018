/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @author Pablo Alvarado
 * @date   18.08.2018
 */

#ifndef ANPI_MULLER_HPP
#define ANPI_MULLER_HPP

#include <vector>
#include <type_traits>

#include <boost/type_traits/is_complex.hpp>
#include <boost/math/tools/polynomial.hpp>

#include <Deflation.hpp>
#include <NewtonRaphson.hpp>



namespace anpi {
  
  /// Enumerator makes explicit polishing roots
  enum PolishEnum {
    DoNotPolish,
    PolishRoots
  };

  namespace bmt=boost::math::tools; // for polynomial

template<typename T>
struct extractType;

/**
 * Extracts the inner type of a template
 * Example:
 *    C = std::complex<double> //Outer type
 *    D = double //Inner type
 */
template<template<typename ...> class C, typename D>
struct extractType<C<D>>
{ using subType = D;};

/**
 * Approximates one real root of a polynomial
 * Return in complex format in case of finding a complex root, which needs
 * needs to be deflated
 *
 * @param[in] poly polynomial to be analyzed for roots
 * @param[in] x0 first root approximation
 * @param[in] x1 second root approximation
 * @param[in] x2 third root approximation
 *
 * return root found
 */
template<typename T>
std::complex<T> calculateRealRoot(const bmt::polynomial<T>& poly, T& x0, T& x1,
    T& x2) {

  T root = T(0);

  for (size_t i = 0; i < 10000; ++i) {

    T h1 = x1 - x0;
    T h2 = x2 - x1;

    T f0 = poly.evaluate(x0);
    T f1 = poly.evaluate(x1);
    T f2 = poly.evaluate(x2);

    T d1 = (f1 - f0) / h1;
    T d2 = (f2 - f1) / h2;

    // Coefficients to approximate the root
    T A = (d2 - d1) / (h2 + h1);
    T B = d2 + h2 * A;
    T C = f2;

    T discriminant = B * B - 4 * A * C;
    if (discriminant < 0) { //Root needs to be removed
      std::complex<T> compSqrt = std::complex<T>(0,
          std::sqrt(std::abs(discriminant)));
      std::complex<T> compDenominator = B - compSqrt;
      if (std::abs(B) > 0) {
        compDenominator = B + compSqrt;
      }
      std::complex<T> compFormula = (T(-2) * C) / compDenominator;
      std::complex<T> compRoot = x2 + compFormula;
      return compRoot;
    }

    T sqrtDiscriminant = std::sqrt(discriminant);
    T denominator = B - sqrtDiscriminant;

    if (B > 0) {
      denominator = B + sqrtDiscriminant;
    }

    T quadFormula = (T(-2) * C) / denominator;
    T root = x2 + quadFormula; // New approximation of the root

    bool convergenceCond1 = std::abs(quadFormula)
        < std::abs(std::numeric_limits<T>::epsilon());
    bool convergenceCond2 = std::abs(poly.evaluate(root))
        < std::abs(std::numeric_limits<T>::epsilon());
    if (convergenceCond1 || convergenceCond2) {
      return root;
    }
    x0 = x1;
    x1 = x2;
    x2 = root;

    f0 = poly.evaluate(x0);
    f1 = poly.evaluate(x1);
    f2 = poly.evaluate(x2);

    // Perturb the points if they have the same value
    if ((f0 == f1) && (f1 == f2)) {
      x0 *= T(0.85);
      x1 *= T(0.85);
      x2 *= T(0.85);
    }
  }
  return std::complex<T>(root);
} //end calculateRealRoot

/**
 * Approximates the value of a complex root
 *
 * @param[in] poly polynomial to be analyzed for roots
 * @param[in] x0 first root approximation
 * @param[in] x1 second root approximation
 * @param[in] x2 third root approximation
 *
 * @return Complex root found
 */
template<typename T, typename U>
std::complex<T> calculateComplexRoot(
    const bmt::polynomial<std::complex<T>>& poly,
    std::complex<U>& x0, std::complex<U>& x1, std::complex<U>& x2) {

  std::complex<T> root(0, 0);
  std::complex<T> x_0 = static_cast<std::complex<T>>(x0);
  std::complex<T> x_1 = static_cast<std::complex<T>>(x1);
  std::complex<T> x_2 = static_cast<std::complex<T>>(x2);

  for (size_t i = 0; i < 10000; ++i) {
    std::complex<T> h1 = x_1 - x_0;
    std::complex<T> h2 = x_2 - x_1;

    std::complex<T> f0 = poly.evaluate(x_0);
    std::complex<T> f1 = poly.evaluate(x_1);
    std::complex<T> f2 = poly.evaluate(x_2);

    std::complex<T> d1 = (f1 - f0) / h1;
    std::complex<T> d2 = (f2 - f1) / h2;

    // Coefficients to approximate the root
    std::complex<T> A = (d2 - d1) / (h2 + h1);
    std::complex<T> B = d2 + h2 * A;
    std::complex<T> C = f2;

    std::complex<T> discriminant = B * B
        - static_cast<std::complex<T>>(4) * A * C;

    std::complex<T> sqrtDiscriminant = std::sqrt(discriminant);
    std::complex<T> denominator = B - sqrtDiscriminant;

    if (std::abs(B - sqrtDiscriminant) < std::abs(B + A)) {
      denominator = B + sqrtDiscriminant;
    }

    std::complex<T> quadFormula = (static_cast<std::complex<T>>(-2) * C)
        / denominator;
    std::complex<T> root = x_2 + quadFormula; // New approximation of the root

    bool convergenceCond1 = std::abs(quadFormula)
        < std::abs(std::numeric_limits<float>::epsilon());
    bool convergenceCond2 = std::abs(poly.evaluate(root))
        < std::abs(std::numeric_limits<float>::epsilon());

    if (typeid(U) == typeid(double)) {
      convergenceCond1 = std::abs(quadFormula)
          < std::abs(std::numeric_limits<double>::epsilon());
      convergenceCond2 = std::abs(poly.evaluate(root))
          < std::abs(std::numeric_limits<double>::epsilon());
    }

    if (convergenceCond1 || convergenceCond2) {
      return root;
    }
    x_0 = x_1;
    x_1 = x_2;
    x_2 = root;

    f0 = f1;
    f1 = f2;
    f2 = poly.evaluate(x_2);

    // Perturb the points if they have the same value
    if ((f0 == f1) && (f1 == f2)) {
      x_0 *= T(0.85);
      x_1 *= T(0.85);
      x_2 *= T(0.85);
    }
  }
  return root;
} //end calculateComplexRoot

/**
 * Finds all the roots of polynomials with real coefficients and
 * real roots only (complex roots -if any- are deflate)
 *
 * @param[in] poly Polynomial to be analyzed for roots
 * @param[in] x0 First root approximation
 * @param[in] x1 Second root approximation
 * @param[in] x2 Third root approximation
 * @param[out] roots Vector where found roots are stored
 * @param[in] polish Flag to enable/disable roots polishing
 */
template<typename T>
void findAllReal(const bmt::polynomial<T>& poly, T& x0, T& x1, T& x2,
    std::vector<T>& roots, const PolishEnum polish) {

  bmt::polynomial<T> quotient = poly;
  std::vector<std::complex<T>> allRoots(poly.degree());
  bmt::polynomial<T> polyResiduo = poly;
  T residuo = 0;
  for (size_t i = 0; i < poly.degree(); ++i) {
    std::complex<T> root = calculateRealRoot(quotient, x0, x1, x2);

    if (root.imag() != 0) {
      quotient = deflate2(quotient, root, polyResiduo);
    }
    else {
      T realRoot = root.real();
      if (polish == PolishEnum::PolishRoots) {
        realRoot = newtonRaphson(poly, realRoot);
      }
			roots.push_back(root.real());
      quotient = deflate(quotient, root.real(), residuo);
    }
		if(quotient.degree()==0) return;
  }
} // findAllReal


/**
 * Finds all the complex roots of a polynomial with complex
 * coefficients
 *
 * @param[in] poly Polynomial to be analyzed for roots
 * @param[in] x0 First root approximation
 * @param[in] x1 Second root approximation
 * @param[in] x2 Third root approximation
 * @param[out] roots Vector where found roots are stored
 * @param[in] polish Flag to enable/disable roots polishing
 */
template<typename T, typename U>
void findAllComplex(const bmt::polynomial<std::complex<T>>& poly,
    std::complex<U>& x0, std::complex<U>& x1, std::complex<U>& x2,
    std::vector<std::complex<U>>& roots, const PolishEnum polish) {

  bmt::polynomial<std::complex<T>> quotient = poly;
  bmt::polynomial<std::complex<T>> compResiduo = poly;
  std::complex<T> residuo = 0;

  for (size_t i = 0; i < poly.degree(); ++i) {
    std::complex<T> root = calculateComplexRoot<T, U>(quotient, x0, x1, x2);
    roots.push_back(static_cast<std::complex<U>>(root));
    quotient = deflate(quotient, root, residuo);
  }
} //findAllComplex

/**
 * Enable only when: Coefficients -> float/double
 *                   Roots        -> float/double
 *
 * Compute the roots of the given polynomial using the Muller method.
 * @param[in] poly polynomial to be analyzed for roots
 * @param[out] roots all roots found
 * @param[in] start initial point for finding the roots
 * @param[in] polish indicate if polishing is needed or not.
 *
 * @return the number of roots found
 */
template<class T, class U>
typename std::enable_if<
    (std::is_floating_point<T>::value && std::is_floating_point<U>::value), void>::type muller(
    const bmt::polynomial<T>& poly, std::vector<U>& roots,
    const PolishEnum polish, const U start) {

  static_assert(std::is_floating_point<T>::value ||
      boost::is_complex<T>::value,
      "T must be floating point or complex");
  static_assert(std::is_floating_point<U>::value ||
      boost::is_complex<U>::value,
      "U must be floating point or complex");

  U dx = U(1);
  U x0 = start;
  U x1 = start - dx;
  U x2 = start - U(2) * dx;


  findAllReal<U>(poly, x0, x1, x2, roots, polish);
}


/**
 * Enable only when: Coefficients -> dComplex/fComplex
 *                   Roots        -> dComplex/fComplex
 *
 * Compute the roots of the given polynomial using the Muller method.
 * @param[in] poly polynomial to be analyzed for roots
 * @param[out] roots all roots found
 * @param[in] start initial point for finding the roots
 * @param[in] polish indicate if polishing is needed or not.
 *
 * @return the number of roots found
 */
template<class T, class U>
typename std::enable_if<
    (boost::is_complex<T>::value && boost::is_complex<U>::value), void>::type muller(
    const bmt::polynomial<T>& poly, std::vector<U>& roots,
    const PolishEnum polish, const U start) {

  static_assert(std::is_floating_point<T>::value ||
      boost::is_complex<T>::value,
      "T must be floating point or complex");
  static_assert(std::is_floating_point<U>::value ||
      boost::is_complex<U>::value,
      "U must be floating point or complex");

  typedef typename extractType<T>::subType TSub; //Inner type of T
  typedef typename extractType<U>::subType USub; //Inner type of U

  USub dx = USub(1);
  std::complex<USub> x0 = start;
  std::complex<USub> x1 = start - dx;
  std::complex<USub> x2 = start - USub(2) * dx;

  findAllComplex<TSub, USub>(poly, x0, x1, x2, roots, polish);
}



/**
 * Enable only when: Coefficients -> dComplex/fComplex
 *                   Roots        -> float/double
 *
 * Compute the roots of the given polynomial using the Muller method.
 * @param[in] poly polynomial to be analyzed for roots
 * @param[out] roots all roots found
 * @param[in] start initial point for finding the roots
 * @param[in] polish indicate if polishing is needed or not.
 *
 * @return the number of roots found
 */
template<class T, class U>
typename std::enable_if<
    (boost::is_complex<T>::value && !boost::is_complex<U>::value), void>::type muller(
    const bmt::polynomial<T>& poly, std::vector<U>& roots,
    const PolishEnum polish, const U start) {

  static_assert(std::is_floating_point<T>::value ||
      boost::is_complex<T>::value,
      "T must be floating point or complex");
  static_assert(std::is_floating_point<U>::value ||
      boost::is_complex<U>::value,
      "U must be floating point or complex");

  typedef typename extractType<T>::subType TSub; //Inner type of T

  std::complex<U> dx = static_cast<std::complex<U>>(1.0);
  std::complex<U> x0 = { start, 0 };
  std::complex<U> x1 = x0 - dx;
  std::complex<U> x2 = x0 - static_cast<std::complex<U>>(2.0) * dx;

  bmt::polynomial<T> quotient = poly;
  T residuo = 0;

  for (size_t i = 0; i < poly.degree(); ++i) {
    T root = calculateComplexRoot<TSub, U>(quotient, x0, x1, x2);
    //If imaginary part is negligible
    if (std::abs(root.imag()) < std::abs(std::numeric_limits<U>::epsilon() * root.real())) {
      roots.push_back(static_cast<U>(root.real()));
    }
    quotient = deflate(quotient, root, residuo);
  }
}




/**
 * Enable only when: Coefficients -> float/double
 *                   Roots        -> dComplex/fComplex
 *
 * Compute the roots of the given polynomial using the Muller method.
 * @param[in] poly polynomial to be analyzed for roots
 * @param[out] roots all roots found
 * @param[in] start initial point for finding the roots
 * @param[in] polish indicate if polishing is needed or not.
 *
 * @return the number of roots found
 */
template<class T, class U>
typename std::enable_if<(!boost::is_complex<T>::value && boost::is_complex<U>::value), void>::type muller(
    const bmt::polynomial<T>& poly, std::vector<U>& roots, const PolishEnum polish, const U start) {

  static_assert(std::is_floating_point<T>::value ||
      boost::is_complex<T>::value,
      "T must be floating point or complex");
  static_assert(std::is_floating_point<U>::value ||
      boost::is_complex<U>::value,
      "U must be floating point or complex");

  if (start.imag() != 0) {
    throw anpi::Exception("Polynomial with real coefficients can't have a complex initial point");
  }

  typedef typename extractType<U>::subType USub; //Inner type of U

  bmt::polynomial<T> quotient = poly;
  bmt::polynomial<T> polyResiduo = poly;
  T residuo = 0;

  T dx = T(1.0);
  T x0 = static_cast<T>(start.real());
  T x1 = x0 - dx;
  T x2 = x0 - T(2.0) * dx;
  size_t degree = poly.degree();

  for (size_t i = 0; i < degree; ++i) {
    std::complex<T> root = calculateRealRoot(quotient, x0, x1, x2);
    if (root.imag() == 0) {
      quotient = deflate(quotient, static_cast<T>(root.real()), residuo);
      roots.push_back(std::complex<USub>(static_cast<USub>(root.real()), USub(0)));
    }
    else {
      roots.push_back(static_cast<U>(root));
      std::complex<USub> conjugate = static_cast<U>(std::conj(root));
      roots.push_back(conjugate);
      quotient = deflate2(quotient, root, polyResiduo);
    }
    if (quotient.degree() == 0) return;
  }
}
 

}


#endif
