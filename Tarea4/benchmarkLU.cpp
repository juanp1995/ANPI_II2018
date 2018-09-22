/**
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @author Juan Pablo Brenes
 * @date   19.9.2018
 */


#include <boost/test/unit_test.hpp>


#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>

/**
 * Unit tests for the matrix class
 */
#include "benchmarkFramework.hpp"
#include "Allocator.hpp"
#include "LUCrout.hpp"
#include "LUDoolittle.hpp"


/// Benchmark for addition operations
template<typename T>
class benchLU {
protected:
  /// Maximum allowed size for the square matrices
  const size_t _maxSize;

  /// A large matrix holding 
  anpi::Matrix<T> _data;

  /// State of the benchmarked evaluation
  anpi::Matrix<T> _a;
  anpi::Matrix<T> _b;
  anpi::Matrix<T> _c;
public:
  /// Construct
  benchLU(const size_t maxSize)
    : _maxSize(maxSize),_data(maxSize,maxSize,anpi::DoNotInitialize) {

    size_t idx=0;
    for (size_t r=0;r<_maxSize;++r) {
      for (size_t c=0;c<_maxSize;++c) {
        _data(r,c)=idx++;
      }
    }
  }

  /// Prepare the evaluation of given size
  void prepare(const size_t size) {
    assert (size<=this->_maxSize);
    this->_a=std::move(anpi::Matrix<T>(size,size,_data.data()));
    this->_b=this->_a;
  }
};


/// Evaluation method for the Crout algorithm 
template<typename T>
class benchLUCrout : public benchLU<T> {
public:
  /// Constructor
  benchLUCrout(const size_t n) : benchLU<T>(n) { }
  
  // Evaluate add on-copy
  inline void eval() {
  	anpi::Matrix<T> LU;
  	std::vector<size_t> permut;
  	anpi::luCrout(this->_a, LU, permut);
  }
};


/// Evaluation method for the Doolittle algorithm
template<typename T>
class benchLUDoolittle : public benchLU<T> {
public:
  /// Constructor
  benchLUDoolittle(const size_t n) : benchLU<T>(n) { }
  
  // Evaluate add on-copy
  inline void eval() {
  	anpi::Matrix<T> LU;
  	std::vector<size_t> permut;
  	anpi::luDoolittle(this->_a, LU, permut);
  }
};

BOOST_AUTO_TEST_SUITE( LU )

/**
 * Instantiate and test the methods of the Matrix class
 */
BOOST_AUTO_TEST_CASE( Decomposition ) {

  std::vector<size_t> sizes = {  24,  32,  48,  64,
                                 96, 128, 192, 256,
                                384, 512, 768,1024,
                               1536,2048,3072,4096};

  const size_t n=sizes.back();
  const size_t repetitions=100;
  std::vector<anpi::benchmark::measurement> times;

  {
    benchLUCrout<double>  baoc(n);

    // Measure on-copy add
    ANPI_BENCHMARK(sizes,repetitions,times,baoc);
    
    //::anpi::benchmark::write("LU_Crout_double.txt",times);
    //::anpi::benchmark::plotRange(times,"Crout (double)","r");
  }

  /*{
    benchLUDoolittle<double>  baoc(n);

    // Measure on-copy add
    ANPI_BENCHMARK(sizes,repetitions,times,baoc);
    
    ::anpi::benchmark::write("LU_Doolittle_double.txt",times);
    ::anpi::benchmark::plotRange(times,"Doolittle (double)","g");
  }*/
  
  ::anpi::benchmark::show();
}
  
BOOST_AUTO_TEST_SUITE_END()
