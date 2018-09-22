/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Juan Pablo Brenes
 * @Date  : 03.03.2018
 */
 
 #include "LUDoolittle.hpp"
 
#ifndef ANPI_LU_HPP
#define ANPI_LU_HPP
 
 namespace anpi{
 
 	
 	/*
 	 *	Calls the fastest LU decomposition method
 	 *
 	 */
	template<typename T>
	inline void lu(const anpi::Matrix<T>& A,
								anpi::Matrix<T>& LU,
								std::vector<size_t>& p){
								
		luDoolittle(A, LU, p);
	}
 
 }//anpi
 
 #endif