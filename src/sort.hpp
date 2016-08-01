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
using namespace std;
int genSortedKmers(char *disorderedKmerpath, const uint8_t kmerSize);

void *multiThreadSort(void *arg);

void quickSort(uint64_t *base, uint64_t num, int comp(const void *, const void *));

void swap(const void *a, const void *b);

int cmp(const void *a, const void *b);

uint64_t rcKmer(uint64_t p_kmerValue, uint8_t p_kmerSize); 

uint64_t binarySearch(uint64_t mk, uint64_t *target, int64_t up);//return the upper bound
#endif
