
#include <cstdlib>
#include <string>
#include <iostream>
#include <algorithm>

#include <boost/program_options.hpp>
#include <boost/type_traits/is_complex.hpp>


/// Convert the given string to lowercase
std::string tolower(std::string s) {
  std::transform(s.begin(),s.end(),s.begin(),
                 [](unsigned char c) -> unsigned char {
                   return std::tolower(c);
                 });
  return s;
}


namespace bmt=boost::math::tools; // for polynomial

namespace po = boost::program_options;

////////////////////////////////////////////////////////////////////////////
//  Main program
////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {

  try {
    po::options_description desc("Allowed options");
    desc.add_options()
      ("topTemperature,t",po::value<std::string>()->default_value(""),
       "Temperature in the upper border")
       
     	("bottomTemperature, b", po::value<std::string>()->default_value(""),
     	"Temperature in the lower border")
     	
     	("leftTemperature, l", po::value<std::string>()->default_value(""),
     	"Temperature in the left border")
     	
     	("rightTemperature, d", po::value<std::string>()->default_value(""),
     	"Temperature in the right border")
     	
     	("isolate, i", "Isolate borders")
     	
     	("termalProfile, p", po::value<std::string>()->default_value(""),
     	"File with a termal profile")
     	
     	("horizontalPix, h", po::value<std::string>()->default_value(""),
     	"Number of horizontal pixels of the solution")
     	
     	("verticalPix, v", po::value<std::string>()->default_value(""),
     	"Number of vertical pixels of the solution")
     	
     	("deactivateVisual, q", "Deactivate all forms of visualitations")
     	
     	("heatFlow, f", "Activate the visualization of the heat flow")
     	
     	("heatGrid, g", po::value<std::string>()->default_value(""),
     	"Size of the grid for the visualization of the heat flow")
      ;

    po::variables_map vm;

    po::store(parse_command_line(argc, argv, desc), vm);
    
    if (vm.count("help")) {
      std::cout
        << desc << "\n\n";
      std::cout
        << "For the coefficient and root types, use one of the following:\n"
           "  float    real numbers with single precision \n"
           "  double   real numbers with double precision \n"
           "  fcomplex complex numbers with float components\n"
           "  dcomplex complex numbers with double components\n\n"
        << "The polynomial formulae are composed by one or more \n"
           "polynomial terms:\n"
           "  polynomial  := term [ {'+'|'-'} term ]*\n"
           "  term        := coefficient 'x^' exponent\n"
           "  coefficient := double | complex\n"
           "  complex     := { '(' double ',' double ')' } | (double)\n"
           "  exponent    := unsigned integral\n\n"
        << "Examples of valid polynomials:\n"
           "  x^3 + 2x + 1\n"
           "  -5x^4 + 2.5x^2 + x\n"
           "  -3x + 5x^4 + 1 -2x^2\n"
           "  (0,1)x^4 + (5,2)x^2 + 1.5x + (1.5,2)\n"
           "The last example has some complex coefficients"
        << std::endl;
      
      return EXIT_SUCCESS;
    }

    po::notify(vm);

    std::string poly = vm["poly"].as<std::string>();

    
    if (vm.count("coefficients")) {
      coefType = as_type(vm["coefficients"].as<std::string>());
    }
    
    if (vm.count("roots")) {
      rootType = as_type(vm["roots"].as<std::string>());
    }

    if (vm.count("muller")) {
      config.method = Muller;
    }

    if (vm.count("jenkinstraub")) {
      config.method = JenkinsTraub;
    }

    if (vm.count("polish")) {
      config.polish = anpi::PolishRoots;
    }

		if (vm.count("start")) {
			std::string start = vm["start"].as<std::string>();
			config.start = boost::lexical_cast<std::complex<double> >('('+start+')');
			std::cout << "initial point: " << start << std::endl;			
		}
    
    // Dispatch with the proper types to call the real workers
    dispatch(config,poly,coefType,rootType);
    
  } catch(po::error& e) { 
    std::cerr << "Error:\n  " << e.what() << std::endl << std::endl; 
    return EXIT_FAILURE;
  } catch(std::exception& e)  { 
    std::cerr << "Unhandled Exception reached the top of main:\n  " 
              << e.what() << "\nApplication will now exit" << std::endl; 
    return EXIT_FAILURE;
  } 
  
  return EXIT_SUCCESS;
}
