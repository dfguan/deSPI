// =====================================================================================
//
//       Filename:  preprocess.hpp
//
//    Description:  header file of Preprocess.cpp
//
//        Version:  1.0
//        Created:  10/21/2015 09:38:34 AM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Dengfeng Guan, <dfguan@hit.edu.cn>
//   Organization:  Harbin Institute of Technology
//
// =====================================================================================
#ifndef _PROCESS_HPP
#define _PROCESS_HPP
//#include "debruijin.hpp"
#include <iostream>
//#include <string>
#include <algorithm>
//#include <string.h>
#include <cstring>
#include <string>
#include <cstdlib>
#include <map>
#include <set>
#include <fstream>
#include <vector>
using namespace std;

//#include <string.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#include <stdint.h>
#include <zlib.h>
#include "error.hpp"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);
//#include <iostream>
//#include <cstdint>

#define PREINDEXLEN 13



typedef struct _kmersSpchar{
	uint64_t value;
	//uint64_t lastKmer;
	//uint8_t  spPos;//indicate how many nucleotides before '#'
	uint8_t  infor;//fore 5bp indicate sp pos and last 3 kmer indicate lastchar
	uint32_t taxID; //store tid
	bool operator<(const struct _kmersSpchar& r ) const{
		uint8_t m,rm;
		m = infor >> 3;
		rm = r.infor >> 3;
			
		if (m != rm ) {
			if (m < rm ) {
				uint8_t move = (rm - m)<< 1;
				return 	value <= (r.value>>move);
			} else {
				uint8_t move = (m - rm) << 1;
				return (value >> move) < r.value;
			}

		} else 
			return value < r.value;
	}
	struct _kmersSpchar & operator=(const struct _kmersSpchar& r)  {
		value = r.value;
		infor = r.infor;
		taxID = r.taxID;
		return *this;	
	}
}kmersSpchar;

typedef struct {
	uint32_t gID;
	uint32_t taxonID;
}gid_taxid;



//uint8_t Bit[];

//extern const uint8_t Bit[];

extern map<uint32_t, uint32_t> taxonomyTree;
//int build_deb(char *path, char *refPath, vector<uint64_t> &kmerValue, vector<uint16_t> &kmerInfo);//the path indicate a file containing all sorted k+1mers  

//int build_deb(char *path, char *refPath, vector<uint64_t> &kmerValue, vector<uint16_t> &kmerInfo, vector<uint64_t>& p_heads, vector<uint64_t>&  p_tails); //the path indicate a file containing all sorted k+1mers  
//int build_deb(char *path, char *refPath, vector<uint64_t> &kmerValue, vector<uint16_t> &kmerInfo, vector<uint32_t>& kmerTID, vector<uint64_t>& p_heads, vector<uint64_t>&  p_tails, vector<gid_taxid>& p_taxonIDTab); //the path indicate a file containing all sorted k+1mers  

int build_deb(const char *refPath, uint64_t  *kmerValue, uint16_t *kmerInfo, uint32_t *kmerTID, uint64_t kmerNum, vector<uint64_t>& p_heads, vector<uint64_t>&  p_tails); //the path indicate a file containing all sorted k+1mers  

uint8_t setInEdge(uint8_t info, uint8_t edge);

uint8_t setOutEdge(uint8_t info, uint8_t edge);

//uint64_t binSearch(vector<uint64_t> &kmersValue, uint64_t lowerBound, uint64_t upperBound, uint64_t key);

uint64_t binSearch(uint64_t *kmersValue,  uint64_t lowerBound, uint64_t upperBound,uint64_t key);
uint64_t transIntoBits(char *str_kmer, uint8_t len);

uint64_t transRCIntoBits(char *str_kmer, uint8_t len);
string transIntoChars(uint64_t v, uint8_t len); 


int cutOffMulEdges(uint64_t *p_kmerValue, uint16_t *p_kmerInfo, uint64_t kmerNum);
//int cutOffMulEdges(vector<uint64_t>& p_kmerValue, vector<uint16_t>& p_kmerInfo);

//int handleFrstLastKmer(vector<uint64_t> & p_kmerValue, vector<uint16_t> &p_kmerInfo, vector<uint64_t> &h, vector<uint64_t> &t);

int handleFrstLastKmer(uint64_t  *p_kmerValue, uint16_t *p_kmerInfo, vector<uint64_t> &h, vector<uint64_t> &t);

int output(uint64_t* p_kmerValue, uint64_t kmerNum, uint16_t* p_kmerInfo, vector<kmersSpchar>& p_2kmers, vector<uint8_t>& p_2kmers_0p, uint32_t* p_taxonIDTab, vector<uint32_t>& nkmerTID );//initate last char and produce unsorted 2kmers
//int output(vector<uint64_t>& p_kmerValue, vector<uint16_t>& p_kmerInfo, vector<kmersSpchar>& p_2kmers, vector<kmersSpchar>& p_2kmers_0p, vector<uint32_t>& p_taxonIDTab, vector<uint32_t>& nkmerTID );//initate last char and produce unsorted 2kmers

int genSpKmers(uint64_t kvalue, uint32_t p_tid, vector<kmersSpchar> & p_2k, vector<uint8_t>& p_2k_0p);

bool isEnd(uint16_t info);

bool isStart(uint16_t info);

//int setLabel(vector<uint64_t> &kmerValue, vector<uint16_t> &kmerInfo, vector<uint64_t>& p_heads, vector<uint64_t>& p_tails);//start end start&end transit

int setLabel(uint64_t *kmerValue, uint16_t *kmerInfo, uint64_t kmerNum, vector<uint64_t>& p_heads, vector<uint64_t>& p_tails);//start end start&end transit

uint8_t outEdgesNum(uint8_t info);

uint8_t inEdgesNum(uint8_t info);

int setStart(uint16_t & marker) ;

int setEnd(uint16_t & marker);





int uid2taxonID(char *path, vector<uint64_t> p_fkmer, vector<uint32_t>& correspond, vector<gid_taxid>& p_taxonID);

int taxonTree(const char *taxonomyNodesPath);

int taxonID(const char *giTaxidPath, vector<gid_taxid> & taxonIDTab);

int preprocess(const char* refPath, const char *kmerPath, const char *taxonomyNodesPath,const char *giTaxidPath, string& bwt_s, vector<uint32_t>& nkmerTID, uint64_t *hash_index);

//int mergeSort(vector<uint64_t>& p_kmersValue, vector<uint16_t>& p_kmersInfo, vector<uint32_t>& taxonIDTab, vector<uint32_t>& nkmerTID, vector<kmersSpchar>& p_2kmers, string& p_bwt, uint64_t *hash_index, uint64_t index_start_point);

int mergeSort(uint64_t* p_kmersValue, uint64_t kmerNum, uint16_t* p_kmersInfo, uint32_t* taxonIDTab, vector<uint32_t>& nkmerTID, vector<kmersSpchar>& p_2kmers, string& p_bwt, uint64_t *hash_index, uint64_t index_start_point);
//int mergeSort(uint64_t* p_kmersValue, uint64_t kmerNum, uint16_t* p_kmersInfo, vector<uint32_t>& taxonIDTab, vector<uint32_t>& nkmerTID, vector<kmersSpchar>& p_2kmers, string& p_bwt, uint64_t *hash_index, uint64_t index_start_point);

uint64_t binSearch_pre(vector<uint64_t>& values, uint64_t lowerBound, uint64_t upperBound, uint64_t key );

int sortKmers(vector<kmersSpchar>& p_kmers);

uint64_t findInsertPos(uint64_t* p_kmersValue, uint64_t p_s, uint64_t p_e, uint64_t key);
//uint64_t findInsertPos(vector<uint64_t>& p_kmersValue, uint64_t p_s, uint64_t p_e, uint64_t key);

uint32_t LCA(uint32_t tid1, uint32_t tid2);

uint32_t getTID(char *seq_name) ;

uint32_t binSearch(vector<gid_taxid>& p_taxonIDTab,uint64_t lowerBound, uint64_t upperBound, uint32_t key ); 

int uid2taxonID(char *path, vector<uint64_t> p_fkmer, vector<uint32_t>& correspond, vector<gid_taxid>& p_taxonID);

int _getline(FILE *in, char *ln);

void *sort2kmers(void *v_2kmers);

void *assignTID(void * ar);
#endif
