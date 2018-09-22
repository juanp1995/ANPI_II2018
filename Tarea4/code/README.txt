-----------------------------------------------
               COMPILATION
-----------------------------------------------

Create a directory build:

> mkdir build;

Go into that directory

> cd build;

You can choose to build a release version with:

> cmake ../ -DCMAKE_BUILD_TYPE=Release

or a debug version with

> cmake ../ -DCMAKE_BUILD_TYPE=Debug

And build everything with

> make

-----------------------------------------------
            EXECUTION OF TESTS
-----------------------------------------------

Go into the directory "bin" created after
the compilation

> cd build/bin/

To see the available tests, run

> ./tester --list_content

To run all tests

> ./tester

To run a specific suite of tests

> ./tester -t <suite> -r detailed

  Example: > ./tester -t LU -r detailed
	
To run a specific test case

> ./tester --run_test=<suite>/<testCase>

  Example: > ./tester --run_test=LU/Crout
	
-----------------------------------------------
         EXECUTION OF BENCHMARKS
-----------------------------------------------

Go into the directory "bin" created after
the compilation

> cd build/bin

To see available benchmarks

> ./benchmark --list_content

To run all two benchmarks

> ./benchmark -t Matrix -r detailed

  Note: The graph with the results of the 
        benchmark for the addition of matrices
        will show up first, after close it, the
        execution will continue and the graph
        with the results of the benchmark
        for the decomposition methods will
        show up.

To run the benchmark for the addition of
matrices

> ./benchmark --run_test=Matrix/Add

To run the benchmark for the LU
decomposition methods

> ./benchmark --run_test=Matrix/Decomposition
