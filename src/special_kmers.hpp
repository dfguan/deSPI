/*
 * =====================================================================================
 *
 *       Filename:  special_kmers.hpp
 *
 *    Description:  generate special kmers
 *
 *        Version:  1.0
 *        Created:  04/10/2017 02:01:19 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#ifndef _SPECIAL_KMERS_H
#define _SPECIAL_KMERS_H

#include <stdint.h>
#include <stdio.h>


#include <string>
//#include <sort>
#include <iostream>
#include <algorithm>

//#include "evo_tree.hpp"
#include "error.hpp"
#include "bwt.hpp"
#include <queue>


using namespace std;

#define F_ARY_BUFFER_SIZE 0X100000 //1M

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

typedef struct _f_spChar{
	kmersSpchar kchar;
	int ind;
	
	bool operator<(const struct _f_spChar& r ) const{
		uint8_t m,rm;
		m = kchar.infor >> 3;
		rm = r.kchar.infor >> 3;
			
		if (m < rm ) {
			uint8_t move = (rm - m)<< 1;
			return 	kchar.value > (r.kchar.value>>move);
		} else {
			uint8_t move = (m - rm) << 1;
			return (kchar.value >> move) >= r.kchar.value;
		}

	}
}f_spChar;

bool compare_value(kmersSpchar r, kmersSpchar q);


class sp_kmers{
	private:	 	
	uint8_t kmer_len;	
	uint64_t ind;	
	bool useFile;
	FILE **f_ary;
	kmersSpchar *f_ary_buf;
	bool use_f_ary_buf;
	uint64_t buf_ind;
	public:

	kmersSpchar *_kspchar;
	uint64_t _kspchar_counter;
	
	sp_kmers( uint8_t _kmer_len){
		kmer_len = _kmer_len;
		//_evo_tree = _evo_tree_; 
		useFile = false;
		use_f_ary_buf = false;
		buf_ind = 0;
		ind = 0;
		_kspchar = NULL;
		f_ary = NULL;
		f_ary_buf = NULL;
	}
	int initiate_space(size_t kspCounter)
	{
		//be careful of error may happen
		_kspchar_counter = kspCounter * (kmer_len - 1);	
		_kspchar = (kmersSpchar *)malloc(sizeof(kmersSpchar) * _kspchar_counter);
		if (!_kspchar) {
			fprintf(stderr,"Fail to allocate memory space, Required:%lu\n, use disk as alternative", _kspchar_counter);
			f_ary = (FILE **)malloc(sizeof(FILE *) * (kmer_len - 1));
			if (!f_ary) {
				fprintf(stderr,"Fail to allocate disk space, now quit\n");
				exit(1);		
			
			} else {
				for (int i=0; i < kmer_len - 1; ++i) {
					char fileName[32];
					sprintf(fileName, ".ksp.%d", i);
					f_ary[i] = fopen(fileName, "wb");	
				}
			}
			useFile = true;
			if ((f_ary_buf = new kmersSpchar[sizeof(kmersSpchar) *(kmer_len - 1) * F_ARY_BUFFER_SIZE])) {
				use_f_ary_buf = true;	
			}	
		}

		return 0;
	}	
	string transIntoChars(uint64_t v, uint8_t len) 
	{
		string s(len, 0);
		
		char Chars[] = {'A','C','G','T'};
		
		
		for (uint8_t i=len-1;i!=0;--i) {
			s[i] = Chars[v&0x3];
			v = v>>2;
		}
		s[0] = Chars[v&0x3];

		return s;
	}
	int gen_sp_kmers(uint64_t kvalue, uint32_t tid, bwt *c_bwt) {
		char dna[] = {'A','C','G','T'};			
		uint8_t	move = (kmer_len - 1) << 1;

		uint64_t mask = (uint64_t)-1 >> (64 -( ((kmer_len-1))<<1));
		
		//uint64_t mask = (0x11)<<((_kmer-1)<<1);
		//uint8_t lastChar =  last>>((_kmer-1) << 1);
		
		kmersSpchar assignVar;
		
		if (!useFile) {
			for (uint8_t i = move; i>0;) {
			
				uint8_t lastChar = (kvalue >> i) & 0x3;//there if 
				
				assignVar.infor  = lastChar;
				
				assignVar.value = kvalue & mask;
				
				assignVar.infor |= (i >> 1) << 3;
				
				assignVar.taxID = tid;

				mask >>= 2;
				i -= 2;
				_kspchar[ind++] = assignVar;

			} 
		} else {
			if (!use_f_ary_buf) {
			
				for (uint8_t i = move; i>0;) {
				
					uint8_t lastChar = (kvalue >> i) & 0x3;//there if 
					
					assignVar.infor  = lastChar;
					
					assignVar.value = kvalue & mask;
					
					assignVar.infor |= (i >> 1) << 3;
					
					assignVar.taxID = tid;

					mask >>= 2;
					i -= 2;
					fwrite(&assignVar, sizeof(_kmersSpchar), 1, f_ary[i >> 1]);
				} 
			} else {
				for (uint8_t i = move; i>0;) {
				
					uint8_t lastChar = (kvalue >> i) & 0x3;//there if 
					
					assignVar.infor  = lastChar;
					
					assignVar.value = kvalue & mask;
					
					assignVar.infor |= (i >> 1) << 3;
					
					assignVar.taxID = tid;

					mask >>= 2;
					i -= 2;
					f_ary_buf[F_ARY_BUFFER_SIZE * ((uint64_t)(i >> 1)) + buf_ind] = assignVar;
				} 
				buf_ind++;
				if (buf_ind >= F_ARY_BUFFER_SIZE) {
					for (uint64_t i = 0; i < kmer_len - 1; ++i) {
						fwrite(f_ary_buf + i * F_ARY_BUFFER_SIZE, sizeof(kmersSpchar), F_ARY_BUFFER_SIZE, f_ary[i]);
					}
					buf_ind = 0;	
				}
			}
		}
		uint8_t lastChar = kvalue & 0x3;
	
		//cout<<s_ind<<endl;
		//bwt_s += dna[lastChar];
		//tids.push_back(tid);
	
		c_bwt->bwt_str[c_bwt->bwt_str_ptr++] = dna[lastChar];
		c_bwt->taxonIDTab[c_bwt->tid_ptr++] = tid;		
		return NORMAL_EXIT;
		//2k.push_back();
	
	}
	int gen_sp_kmers(uint64_t kvalue, uint32_t tid, string &bwt_s, vector<uint32_t> &tids) {
		char dna[] = {'A','C','G','T'};			
		uint8_t	move = (kmer_len - 1) << 1;

		uint64_t mask = (uint64_t)-1 >> (64 -( ((kmer_len-1))<<1));
		
		//uint64_t mask = (0x11)<<((_kmer-1)<<1);
		//uint8_t lastChar =  last>>((_kmer-1) << 1);
		
		kmersSpchar assignVar;
		
		if (!useFile) {
			for (uint8_t i = move; i>0;) {
			
				uint8_t lastChar = (kvalue >> i) & 0x3;//there if 
				
				assignVar.infor  = lastChar;
				
				assignVar.value = kvalue & mask;
				
				assignVar.infor |= (i >> 1) << 3;
				
				assignVar.taxID = tid;

				mask >>= 2;
				i -= 2;
				_kspchar[ind++] = assignVar;

			} 
		} else {
			if (!use_f_ary_buf) {
			
				for (uint8_t i = move; i>0;) {
				
					uint8_t lastChar = (kvalue >> i) & 0x3;//there if 
					
					assignVar.infor  = lastChar;
					
					assignVar.value = kvalue & mask;
					
					assignVar.infor |= (i >> 1) << 3;
					
					assignVar.taxID = tid;

					mask >>= 2;
					i -= 2;
					fwrite(&assignVar, sizeof(_kmersSpchar), 1, f_ary[i >> 1]);
				} 
			} else {
				for (uint8_t i = move; i>0;) {
				
					uint8_t lastChar = (kvalue >> i) & 0x3;//there if 
					
					assignVar.infor  = lastChar;
					
					assignVar.value = kvalue & mask;
					
					assignVar.infor |= (i >> 1) << 3;
					
					assignVar.taxID = tid;

					mask >>= 2;
					i -= 2;
					f_ary_buf[F_ARY_BUFFER_SIZE * ((uint64_t)(i >> 1)) + buf_ind] = assignVar;
				} 
				buf_ind++;
				if (buf_ind >= F_ARY_BUFFER_SIZE) {
					for (uint64_t i = 0; i < kmer_len - 1; ++i) {
						fwrite(f_ary_buf + i * F_ARY_BUFFER_SIZE, sizeof(kmersSpchar), F_ARY_BUFFER_SIZE, f_ary[i]);
					}
					buf_ind = 0;	
				}
			}
		}
		uint8_t lastChar = kvalue & 0x3;
	
		//cout<<s_ind<<endl;
		bwt_s += dna[lastChar];
		tids.push_back(tid);
		
		return NORMAL_EXIT;
		//2k.push_back();
	
	}
	
	int flush_sp_kmers()
	{
		if (buf_ind) {
			for (uint64_t i = 0; i < kmer_len - 1; ++i) {
				fwrite(f_ary_buf + i * F_ARY_BUFFER_SIZE, sizeof(kmersSpchar), buf_ind, f_ary[i]);
			}
		}	
	
		return NORMAL_EXIT;	
	}

	//	close opened files and release  f_ary_buf	
	int closeFiles()
	{
		
		if (useFile) {
			for (int i = 0; i < kmer_len-1; ++i) fclose(f_ary[i]);
			if (f_ary_buf) 
				delete[] f_ary_buf;
		}
		return NORMAL_EXIT;	
	}
	//sort files and merge to one file with name "ksp.srt"
	int externalSort()
	{
		for (int i=0; i < kmer_len - 1; ++i) {
			char fileName[32];
			sprintf(fileName, ".ksp.%d", i);
			f_ary[i] = fopen(fileName, "rb+");	
		}
		
		uint64_t buf_len = _kspchar_counter / (kmer_len - 1);  // in case it is 1 
		f_ary_buf = new kmersSpchar[buf_len + 1];
		uint64_t left_ele[kmer_len-1] ; 
		
		for (int i = 0; i < kmer_len - 1; ++i) {
			left_ele[i] = buf_len;
		}	
			 
		if (!f_ary_buf) {
			fprintf(stderr, "Fail to allocate space for external sort, now exit\n");
			exit(1);
		}	
		for (int i = 0; i < kmer_len - 1; ++i) {
			fread(f_ary_buf, sizeof(kmersSpchar), buf_len, f_ary[i]);	
			stable_sort(f_ary_buf, f_ary_buf + buf_len, compare_value);
			fseek(f_ary[i], 0, SEEK_SET);
			fwrite(f_ary_buf, sizeof(kmersSpchar), buf_len, f_ary[i]);
		}	
			
		for (int i=0; i < kmer_len - 1; ++i) {
			fseek(f_ary[i], 0, SEEK_SET);
		}
		
		FILE *fp = fopen("ksp.srt", "wb");
		std::priority_queue<f_spChar> prq;
		
		f_spChar f;
		for(int i = 0; i < kmer_len - 1; ++i) {
			if (left_ele[i]) {
				fread(f_ary_buf, sizeof(kmersSpchar), 1, f_ary[i]);
				//cout<<"INSERT"<<transIntoChars(f_ary_buf[0].value, f_ary_buf[0].infor >> 3)<<endl;
				f.kchar = f_ary_buf[0];
				f.ind = i;
				prq.push(f);
				--left_ele[i];
			}
		}		
		
		uint64_t i_buf_ind = 1;
		while (!prq.empty()) {
			f_ary_buf[i_buf_ind++] = prq.top().kchar;
			//cout<<transIntoChars(f_ary_buf[i_buf_ind-1].value, f_ary_buf[i_buf_ind-1].infor >> 3)<<endl;
			int ind = prq.top().ind;	
			if (i_buf_ind > buf_len) {
				fwrite(f_ary_buf+1, sizeof(kmersSpchar), buf_len, fp);
				i_buf_ind = 1;	
			}
			prq.pop();
			
			if (left_ele[ind]) {
				fread(f_ary_buf, sizeof(kmersSpchar), 1, f_ary[ind]);
				f.kchar = f_ary_buf[0];
				f.ind = ind;
				prq.push(f);
				--left_ele[ind];
			}
		}	


		if (i_buf_ind > 1) {
				fwrite(f_ary_buf+1, sizeof(kmersSpchar), i_buf_ind - 1, fp);
		}
		fclose(fp);
					
		for (int i=0; i < kmer_len - 1; ++i) {
			fclose(f_ary[i]);
		}
			
		for (int i=0; i < kmer_len - 1; ++i) {
			char fileName[32];
			sprintf(fileName, ".ksp.%d", i);
			remove(fileName);
		}
		
		
		return NORMAL_EXIT;			
	
	}
	int sortKmers()
	{	
		if (!useFile) {
			kmersSpchar* start = _kspchar;
			kmersSpchar* end = _kspchar + _kspchar_counter;
			stable_sort(start, end);
		} else {
			externalSort();
		
		}
		return NORMAL_EXIT;
	}
	
	//void init(se_node *_se_node_, string & bwt_s)
	//{
		
	//}

	
	~sp_kmers() {
		if (_kspchar) 
			free(_kspchar);
		if (f_ary) 
			delete[] f_ary;
		if (f_ary_buf) 
			delete[] f_ary_buf;

		//if (_kspchar_counter) 
			//free(_kspchar_counter);
	}
};

#endif


