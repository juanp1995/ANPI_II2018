/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Juan Pablo Brenes Coto
 * @Date  : 07.10.2018
 */

#ifndef ANPI_RESISTORGRID_HPP
#define ANPI_RESISTORGRID_HPP

#include <iostream>
#include "Exception.hpp"
#include <stdlib.h>
#include <vector>

#include "Matrix.hpp"

namespace anpi {

//Pack a pair of indices of the nodes of a resistor
struct IndexPair{
	std::size_t row1;
	std::size_t col1;
	std::size_t row2;
	std::size_t col2;
  };

  class ResistorGrid {
  private:
    //Matrix of the current equation system
    Matrix<float> A_;
    //Vector of the current equation system
    std::vector<float> b_;
    //Vector of solutions of the current equation system
    std::vector<float> x_;

    //Raw map data
    Matrix<float> rawMap_;

    /**
     * Load the image of a map and converts it to anpi Matrix
     */
    void loadMap(std::string filename);

  public:
    //Constructors
    ResistorGrid();

    /**
     * Maps the index of the terminal nodes of a resistor to a linear index
     */
    std::size_t nodesToIndex(const std::size_t row1, const std::size_t col1, const std::size_t row2,
        const std::size_t col2);

    /**
     * Convert and index to the pair of node coordinates
     */
    IndexPair indexToNodes(const std::size_t idx);

    /**
     * Construct the grid from the given file
     * @return true if successful or false otherwise
     */
    bool build(const std::string filename);

    /**
     * Compute the internal data to navigate between the given nodes
     */
    bool navigate(const IndexPair& nodes);
  };
}

#include "ResistorGrid.cpp"

#endif
