/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @author Pablo Alvarado
 * @date   29.12.2017
 */


#include <chrono>
#include <iostream>
#include <ostream>
#include <fstream>
#include <limits>

#include <PlotPy.hpp>



namespace anpi {

  namespace plot{
		
		template<typename T>
		struct measurements {
			std::vector<T> eps;
  		std::vector<size_t> calls;
		};

		template<typename T>
		struct benchmarkData {
			inline benchmarkData<T>() : measures(4) {};
			std::string method;
			std::vector<std::string> functions = {"t1", "t2", "t3", "t4"};
  		std::vector<measurements<T>> measures;
		};
  }

  namespace benchmarkPlot {


    /**
     * Save a file with each measurement in a row.
     *
     * The meaning of the columns is as follows:
     * # Size
     * # Average
     * # Standard deviation
     * # Minimum
     * # Maximum  
     */
		template<typename T>
    void write(std::ostream& stream,
               const std::vector<anpi::plot::benchmarkData<T>>& data) {
      for (T i=0; i<data.size(); ++i) {
				anpi::plot::benchmarkData<T> tempData = data[i];
				stream << tempData.method << std::endl;
					
				for (T n=0; n<tempData.functions.size(); ++n) {
					anpi::plot::measurements<T> tempMeasures = tempData.measures[n];
					stream << tempData.functions[n] << std::endl;

					for (T m=0; m<tempMeasures.eps.size(); ++m) {
						stream << tempMeasures.eps[m] << " \t";
						stream << tempMeasures.calls[m] << std::endl;					
					}
				}
				stream << " ------------------------ \n" << std::endl;
      }
    }

    /**
     * Save a file with each measurement in a row
     */
		template<typename T>
    void write(const std::string& filename,
               const std::vector<anpi::plot::benchmarkData<T>>& data) {
      std::ofstream os(filename.c_str());
      write(os,data);
      os.close();
    }

		/*
		template<typename T>
		void plot(const anpi::plot::benchmarkData<T>& benchData,
							const std::string& legend,
							const std::string& color) {
			
			static anpi::Plot2d<double> plotter;
			plotter.initialize(1);
			plotter.plot(benchData.measures, y, legend, color)
		}*/

    

  } // namespace benchmarkPlot
} // namespace anpi

