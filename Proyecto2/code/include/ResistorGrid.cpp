/*
 * Copyright (C) 2017
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Juan Pablo Brenes Coto
 * @Date:   7.10.2018
 */

#include <opencv2/core.hpp>    // For cv::Mat
#include <opencv2/highgui.hpp> // For cv::imread/imshow
#include <string>
#include <cstdlib>

#include "ResistorGrid.hpp"
#include "Exception.hpp"


namespace anpi {

  ResistorGrid::ResistorGrid() {
  }


  /**
   * Loads the map image and converts it into an anpi Matrix
   *
   * @param[in] filename Name of the map image
   */
  void ResistorGrid::loadMap(std::string filename) {
    std::string mapPath = std::string( ANPI_DATA_PATH) + "/" + filename;

    cv::Mat_<float> map;

    cv::imread(mapPath.c_str(), CV_LOAD_IMAGE_GRAYSCALE).convertTo(map, CV_32FC1);

    if (!map.data) throw anpi::Exception("Error loading the image");
    map /= 255.0f;

    anpi::Matrix<float, std::allocator<float> > amapTmp(map.rows, map.cols, map.ptr<float>());
    // And transform it to a SIMD-enabled matrix
    anpi::Matrix<float> tmp(amapTmp);
    this->rawMap_ = tmp;
  }


  /**
   * Maps two nodes of the resistor grid to the corresponding resistor index
   *
   * @param[in] row1 Row of the first node
   * @param[in] col1 Column of the first node
   * @param[in] row2 Row of the second node
   * @param[in] col2 Column of the second node
   *
   * @return index of the resistor between the two nodes
   */
  std::size_t ResistorGrid::nodesToIndex(const std::size_t row1,
      const std::size_t col1, const std::size_t row2, const std::size_t col2) {

    std::size_t resistorIndx = 0;
    if ((row1 == row2) && (col1 == col2))
      throw anpi::Exception("Invalid node index: Same node");

    //Horizontal resistor
    if ((row1 == row2)
        && ((std::max(col1, col2) - std::min(col1, col2)) == 1)) {
      resistorIndx = row1 * (this->rawMap_.cols() - 1) + std::min(col1, col2);
    }

    //Vertical resistor
    else if (((std::max(row1, row2) - std::min(row1, row2)) == 1)
        && (col1 == col2)) {
      resistorIndx = this->rawMap_.rows() * (this->rawMap_.cols() - 1)
          + std::min(row1, row2) * this->rawMap_.cols() + col2;
    }

    //Invalid mapping of nodes
    else {
      throw anpi::Exception("Invalid mapping of nodes");
    }
    return resistorIndx;
  } // nodesToIndex


  /**
   * Maps the index of a resistor to the indexes of the two nodes where it's
   * connected in the grid
   *
   * @param[in] idx Index of the resistor
   *
   * @return Indexpair struct with the indexes grid nodes
   */
  IndexPair ResistorGrid::indexToNodes(const std::size_t idx) {

    if (idx > this->rawMap_.rows())
      throw anpi::Exception("Invalid resistor index: outside grid");

    IndexPair resistor = { };

    //Index corresponds to an horizontal resistor
    if (idx < ((this->rawMap_.cols() - 1) * this->rawMap_.rows())) {
      for (size_t i = 0; i < this->rawMap_.rows(); ++i) {
        if (i * (this->rawMap_.cols() - 1) <= idx
            && idx < (this->rawMap_.cols() - 1) * (i + 1)) {
          resistor.row1 = i;
          break;
        }
        resistor.row2 = resistor.row1;
        resistor.col1 = idx - resistor.row1 * (this->rawMap_.cols() - 1);
        resistor.col2 = resistor.col1 + 1;
      }
    }

    //Index corresponds to a vertical resistor
    else {
      size_t verticalR = this->rawMap_.rows() * (this->rawMap_.cols() - 1);
      for (std::size_t i = 0; i < this->rawMap_.rows(); ++i) {
        if (idx >= (verticalR + i * this->rawMap_.cols())
            && idx < (verticalR + (i + 1) * this->rawMap_.cols())) {
          resistor.row1 = 1;
          break;
        }
      }
      resistor.row2 = resistor.row1 + 1;
      resistor.col1 = idx - verticalR - resistor.row1 * this->rawMap_.cols();
      resistor.col2 = resistor.col1;
    }

    return resistor;
  } // indexToNodes



  /**
   * Construct the grid from the given file
   *
   * @return true if successful or false otherwise
   */
  bool ResistorGrid::build(const std::string filename) {
    throw anpi::Exception("To be implemented");
  }

  /**
   * Compute the internal data to navigate between the given nodes
   */
  bool ResistorGrid::navigate(const IndexPair& nodes) {
    throw anpi::Exception("To be implemented");
  }



}



