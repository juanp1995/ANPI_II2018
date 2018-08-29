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
  } // namespace plot

  namespace benchmarkPlot {

		const std::vector<std::string> colors = {"b", "g", "r", "c"};


    /**
     * Save a file with the benchmark data of the root finders
     *
     * 
     * File structure:
     * -> Method name
     *	  ->Function 
     *			->- eps
     *  		-> number of call to operator()
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


		template<typename T>
		std::vector<std::vector<T>> extractData (const anpi::plot::measurements<T> dataStruct){
			
			std::vector<std::vector<T>> data(2);
			std::vector<T> eps(dataStruct.eps.size()), calls(dataStruct.calls.size());
	
			for (size_t i=0; i<eps.size(); ++i){
				eps[i] = dataStruct.eps[i];
				calls[i] = dataStruct.calls[i];
			}

			data[0] = eps;
			data[1] = calls;

			return data; 
		}

		
		template<typename T>
		void plot(const anpi::plot::benchmarkData<T>& data,
							const std::string& subTitle,
							const std::vector<T> xRange,
							const std::vector<T> yRange,
							int id,
							const int subPlot) {

			std::vector<std::string> functions = data.functions;
			static anpi::Plot2d<T> plotter;
			plotter.initialize(id);
			plotter.supTitle(data.method);
			plotter.subPlot(subPlot);
			plotter.setTitle(subTitle);
			plotter.setXLabel("eps");
			plotter.setYLabel("calls");
			plotter.subPlotHeight(0.7);
			//plotter.setXRange(xRange[0], xRange[1]);
			//plotter.setYRange(yRange[0], yRange[1]);
			

			for (size_t i=0; i<functions.size(); ++i){
				std::vector<std::vector<T>> benchData = extractData(data.measures[i]);
				plotter.plot(benchData[0], benchData[1], functions[i],  colors[i]);

			}
		}


		void show() {
			static anpi::Plot2d<double> plotter;
			plotter.show();
		}

    

  } // namespace benchmarkPlot
} // namespace anpi

