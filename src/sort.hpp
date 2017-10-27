// =====================================================================================
//
//       Filename:  sort.hpp
//
//    Description:  header file of genKmers
//
//        Version:  1.0
//        Created:  10/22/2015 01:58:21 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Dengfeng Guan, <dfguan@hit.edu.cn>
//   Organization:  Harbin Institute of Technology
//
// =====================================================================================
#ifndef _SORT_HPP
#define _SORT_HPP

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <time.h>
#include <fstream>
#include "special_kmers.hpp"
#include "se_node.hpp"

using namespace std;

//#define PREINDEXLEN 13
int genSortedKmers(char *refListPath, const uint8_t kmerSize, char *disorderedKmerpath, char *evo_tree_path);

int mergeSort(char *fpkmers, const uint8_t _kmer, se_node *_se_,  kmersSpchar* p_2kmers, uint64_t kmersSpcharN, string& p_bwt, vector<uint32_t>& nkmerTID, uint64_t *hash_index, uint64_t index_start_point);
int mergeSort(char *fpkmers, const uint8_t _kmer, se_node *_se_,  kmersSpchar* p_2kmers, uint64_t kmersSpcharN, bwt *c_bwt);
int mergeSort(char *fpkmers, const uint8_t _kmer, se_node *_se_,  uint64_t kmersSpcharN, bwt *c_bwt);
int mergeSort(char *fpkmers, const uint8_t _kmer, se_node *_se_, uint64_t kmersSpcharN, string& p_bwt, vector<uint32_t>& nkmerTID, uint64_t *hash_index, uint64_t index_start_point);
int processKmers(const char *refListPath, const uint8_t kmerSize,const  char *disorderedKmerpath,const char *evo_tree_path, string &bwt_s, uint64_t *hash_index, vector<uint32_t> &nkmerTID);
int processKmers(const char *refListPath, const uint8_t kmerSize,const  char *disorderedKmerpath,const char *evo_tree_path, string &bwt_s, uint64_t *hash_index, vector<uint32_t> &nkmerTID, int diff);
int processKmers(const char *refListPath, const uint8_t kmerSize,const  char *disorderedKmerpath,const char *evo_tree_path, bwt *c_bwt);
//int processKmers(char *refListPath, const uint8_t kmerSize, char *disorderedKmerpath, char *evo_tree_path, string &bwt_s, uint64_t *hash_index, vector<uint32_t> &nkmerTID);
//int mergeKmers(char *sortedKmersPath, const uint8_t _kmer_len, kmersSpchar *_kspchar, uint64_t kmersSpcharN, map<uint64_t, struct end_tid> &pairs, string &bwt_s, uint64_t *hash_index, vector<uint32_t> nkmerTID, uint64_t index_start_point);
//int mergeKmers(char *sortedKmersPath, const uint8_t _kmer_len, kmersSpchar *_kspchar, uint64_t kmersSpcharN, kbtree_t(uint64_t) *bt, string &bwt_s, uint64_t *hash_index, vector<uint32_t> nkmerTID, uint64_t index_start_point);
int mergeKmers(char *sortedKmersPath, const uint8_t _kmer_len, kmersSpchar *_kspchar, uint64_t kmersSpcharN, se_node *_se_, string &bwt_s, uint64_t *hash_index, vector<uint32_t> &nkmerTID, uint64_t index_start_point);
int mergeKmers(char *sortedKmersPath, const uint8_t _kmer_len, kmersSpchar *_kspchar, uint64_t kmersSpcharN, se_node *_se_, bwt *c_bwt);
int processKmers(const char *refListPath, const uint8_t kmerSize,const  char *disorderedKmerpath,const char *evo_tree_path, bwt *c_bwt);
//int mergeSort(char *fpkmers, const uint8_t _kmer, kbtree_t(uint64_t) *bt,  kmersSpchar* p_2kmers, uint64_t kmersSpcharN, string& p_bwt, vector<uint32_t>& nkmerTID, uint64_t *hash_index, uint64_t index_start_point);
//int mergeSort(char *fpkmers, const uint8_t _kmer, kbtree_t(uint64_t) *bt, uint64_t kmersSpcharN, string& p_bwt, vector<uint32_t>& nkmerTID, uint64_t *hash_index, uint64_t index_start_point);
void *call_special_kmers(void *arg);
void* sort_thread(void *k);
void *multiThreadSort(void *arg);
int sort_thread_new(uint64_t *writeBuf,  uint64_t bufferSize,int file_ind, uint8_t kmerSize);	
void *keyway_sort_thread(void *);
void quickSort(uint64_t *base, uint64_t num, int comp(const void *, const void *));

void quickSort(uint64_t *base, uint64_t num, int comp(const void *, const void *), int type);

void swap(const void *a, const void *b);

void swap(const void *a, const void *b, int type);
int cmp(const void *a, const void *b);

uint64_t rcKmer(uint64_t p_kmerValue, uint8_t p_kmerSize); 
uint64_t rcOcc(uint64_t );
uint64_t binarySearch(uint64_t mk, uint64_t *target, int64_t up);//return the upper bound
#endif
