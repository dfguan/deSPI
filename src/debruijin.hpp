// =====================================================================================
//
//       Filename:  debruijin.hpp
//
//    Description:  header file of debruijin cpp
//
//        Version:  1.0
//        Created:  10/29/2015 02:37:30 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Dengfeng Guan, <dfguan@hit.edu.cn>
//   Organization:  Harbin Institute of Technology
//
// =====================================================================================
#ifndef _DEBRUIJIN_HPP
#define _DEBRUIJIN_HPP

#include <cstdio>
#include <iostream>
#include <vector>
//#include <cstdint>
using namespace std;

#include <stdint.h>
#include <zlib.h>
#include "error.hpp"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);



typedef struct _kmersSpchar{
	uint64_t value;
	//uint64_t lastKmer;
	//uint8_t  spPos;//indicate how many nucleotides before '#'
	uint8_t  infor;//fore 5bp indicate sp pos and last 3 kmer indicate lastchar
	bool operator<(const struct _kmersSpchar& r ) const{
		uint8_t m,rm;
		m = infor >> 3;
		rm = r.infor >> 3;
			
		if (m != rm ) {
			if (m < rm ) {
				uint8_t move = (rm - m)<< 1;
				return 	value > (r.value>>move);
			} 
			else {
				uint8_t move = (m - rm) << 1;
				return (value >> move) > r.value;
			}

		} else 
			return value > r.value;
	}
}kmersSpchar;

//int build_deb(char *path, char *refPath, vector<uint64_t> &kmerValue, vector<uint16_t> &kmerInfo);//the path indicate a file containing all sorted k+1mers  

int build_deb(char *path, char *refPath, vector<uint64_t> &kmerValue, vector<uint16_t> &kmerInfo, vector<uint64_t>& p_heads, vector<uint64_t>&  p_tails); //the path indicate a file containing all sorted k+1mers  

uint8_t setInEdge(uint8_t info, uint8_t edge);

uint8_t setOutEdge(uint8_t info, uint8_t edge);

uint64_t binSearch(vector<uint64_t> &kmersValue, uint64_t lowerBound, uint64_t upperBound, uint64_t key);


uint64_t transIntoBits(char *str_kmer, uint8_t len);

string transIntoChars(uint64_t v); 

extern uint8_t _kmer;

extern const uint8_t Bit[];

extern uint64_t *counter;

#define PREINDEXLEN 13

#endif

