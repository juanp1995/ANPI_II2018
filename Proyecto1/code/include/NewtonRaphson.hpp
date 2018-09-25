#include <limits>

#include <boost/type_traits/is_complex.hpp>
#include <boost/math/tools/polynomial.hpp>

namespace anpi {

namespace bmt = boost::math::tools;

  /**
   * Calculate the analytic derivative of a polynomial with
   * real coefficients
   *
   * @param[in] realPoly Polynomial to be derived
   *
   * @return Derivative of the polynomial
   */
  template<class T>
  bmt::polynomial<T> realDerivative(const bmt::polynomial<T>& realPoly) {

    bmt::polynomial<T> HPolinomio0(realPoly);
    for (int unsigned (i) = 1; i <= realPoly.degree(); i++) {
      HPolinomio0[i - 1] = HPolinomio0[i] * i;
    }
    HPolinomio0[HPolinomio0.degree()] = 0;
    HPolinomio0.normalize();
    return HPolinomio0;
  }

  /**
   * Polishes a found root of a polynomial with real
   * coefficients, using the Newton-Raphson method
   *
   * @param[in] poly Polynomial whose root was found
   * @param[in] root Root found
   *
   * @return Root polished if the method converges
   */
  template<typename T>
  T newtonRaphson(const bmt::polynomial<T>& poly, T& root) {

    T maxIterations = std::numeric_limits<T>::digits;
    bmt::polynomial<T> derivative = realDerivative(poly);
    T xPrev = root;
    T error = 1;
    T xNext;

    for (int i = 0; i < maxIterations; ++i) {
      xNext = xPrev - (poly.evaluate(xPrev) / derivative.evaluate(xPrev));
      error = xNext - xPrev;

      if (std::abs(error) < std::numeric_limits<T>::epsilon()) {
        break;
      }
      xPrev = xNext;
    }

    if (std::abs(poly.evaluate(xNext)) < std::numeric_limits<T>::epsilon()) {
      return xNext;
    }
    return root;
  }

}

