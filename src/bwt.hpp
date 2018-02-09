// =====================================================================================
//
//       Filename:  bwt.hpp
//
//    Description:  header file of bwt
//
//        Version:  1.0
//        Created:  11/03/2015 03:42:23 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Dengfeng Guan, <dfguan@hit.edu.cn>
//   Organization:  Harbin Institute of Technology
//
// =====================================================================================
#ifndef _BWT_HPP
#define _BWT_HPP
#include <string>
#include <iostream>
#include <vector>
using namespace std;

#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "error.hpp"

#define PREINDEXLEN 13

class bwt{
	
	uint8_t* 		bwt_occ;
	
	uint64_t 		len_bwt_occ;

	uint8_t* 		AGCTCounter;

	uint64_t		rank[5];
	

	int 			threshold;

	
	//uint64_t* 		occCheck;

	//uint64_t 		occCount;


	uint64_t LFC(uint64_t r, uint8_t c);
	uint64_t occ(uint64_t r, uint8_t c);
	uint32_t locate(uint64_t sp);
	uint32_t locate(uint64_t sp, uint8_t *bytes, char *qual, int loc, int &extend);
public:	
	char*			bwt_str;
		
	uint64_t 		len_bwt_str;
	
	uint64_t 		bwt_str_ptr;

	uint32_t 		*taxonIDTab;
	
	uint64_t 		tidSize;	
	
	uint64_t 		tid_ptr;
	
	uint64_t 		*hash_index;		
	bwt(){
		bwt_str = NULL;
		taxonIDTab = NULL;
		hash_index = NULL;
		bwt_occ = NULL;
		AGCTCounter = NULL;
		uint64_t hash_index_size = (uint64_t)1 <<((PREINDEXLEN<<1) + 1);
		hash_index = (uint64_t *)calloc(hash_index_size, sizeof(uint64_t));
		for (int i=0; i<5; rank[i++] = 0);
	}
	bwt(char *bt, uint64_t len, uint64_t *p_hash_index): bwt_str(bt),  len_bwt_str(len) {
		hash_index = p_hash_index;
		for (int i=0; i<5; rank[i++] = 0);
	}
	bwt(int thres){
		threshold = thres;
		for(int i=0; i<5; rank[i++]=0);
	}
	~bwt() {
		if (bwt_str) delete[] bwt_str;
		if (taxonIDTab) delete[] taxonIDTab;
		if (hash_index) delete[] hash_index;
		if (bwt_occ) delete[] bwt_occ;
		if (AGCTCounter) delete[] AGCTCounter;
	}	
	uint64_t transIntoBits(uint8_t *bytes, uint8_t len);

	
	//int exactMatch(char *str, int len, uint64_t& sp, uint64_t& ep, int& match_len);
	int exactMatch(uint8_t *bytes, int len, int& match_len, uint32_t* assignedTID); 
	
	int exactMatch(uint8_t  *bytes, char *qual, int len, int& match_len, uint32_t* assignedTID);
	int dump_index(const char *dirPath, vector<uint32_t>& p_nkmerTID);
	int dump_index(const char *dirPath);
	
	int load_index(const char *dirPath);
	
	int rmdup(uint32_t *temp, int counter, uint32_t* assignedTID);
};


#endif
