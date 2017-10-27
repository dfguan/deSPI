/*
 * =====================================================================================
 *
 *       Filename:  special_kmers.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/27/2017 03:16:49 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include "special_kmers.hpp"

bool compare_value(kmersSpchar r, kmersSpchar q)
{
	return r.value < q.value;
}	
