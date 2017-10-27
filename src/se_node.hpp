/*
 * =====================================================================================
 *
 *       Filename:  se_node.hpp
 *
 *    Description:  header of start and end node
 *
 *        Version:  1.0
 *        Created:  04/10/2017 01:39:23 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#ifndef _SE_NODE_H
#define _SE_NODE_H

#include <stdio.h>
#include <zlib.h>
#include <stdint.h>


#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <string>
using namespace std;

#include "kseq.h"
#include "error.hpp"
#include "evo_tree.hpp"
#include "khash.h"
#include "kbtree.h"

#define SE_THREAD_NUM 8
#define _SE_PAIR_CMP(a, b) (((a).s > (b).s) - ((a).s < (b).s))

KSEQ_INIT(gzFile, gzread);
KHASH_SET_INIT_INT64(64)




const uint8_t Bit[] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

typedef struct _se_pair{
	uint64_t s,e;
	uint32_t tid;	

}se_pair;



class se_node{
	private:	
	int file_id; 	
	vector<const char *> fileNameList;
	int 		fileNum;

	string          fileNames;
	
	vector<uint32_t> taxons;	
	//set<uint64_t> start_nodes;
	//set<uint64_t> end_nodes;
	pthread_rwlock_t rwlock;
	
	khash_t(64) *end_nodes;

	uint8_t kmer_len;
	uint64_t out_mask;
	uint64_t in_shift;	
	evo_tree  *_evo_tree;
	int countNonZero(int i) {
		int j = -1;
		for (int z = 0; z < 4; ++z, i >>= 1) 
			if (i & 0x1) 
				++j;
		return j;	
	}
	void error(char *er)
	{
		cerr<<er<<endl;
		exit(1);
	}

	uint64_t transIntoBits(char *str_kmer, uint8_t len)
	{
		uint64_t value = 0;
		
		for(uint8_t i=0; i<len; ++i) value = (value<<2)|Bit[str_kmer[i]];			
		
		return value;

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

	uint64_t transRCIntoBits(char *str_kmer, uint8_t len) 
	{
		uint64_t value = 0;
		
		for (uint8_t i=0; i<len;++i) value = (value<<2)|(Bit[str_kmer[len - i - 1]]^0x3);

		return value;
	}
	static void *static_manager(void *context) {
				//fprintf(stderr,"%s", ((se_node *)context)->fileNames.c_str());
		return (((se_node *)context)->manager());
	}
	void *manager() {
				//fprintf(stderr,"%s", fileNames.c_str());
		while (true) {
			int id;
			pthread_rwlock_wrlock(&rwlock);
			id = file_id++;
			pthread_rwlock_unlock(&rwlock);
			//fprintf(stderr,"%d\t%d", id, fileNum);
			if (id < fileNum) {
				//fprintf(stderr,"%s", fileNames.c_str());
				find_s_e_pair(fileNameList[id], taxons[id]);			
			
			} else {
			
				//fprintf(stderr,"iquit");
				break;	
			} 
				//break;	
		}
		return NULL;
	}

    	inline uint32_t cas(volatile uint32_t *ptr, uint32_t oval, uint32_t nval) {
      		return __sync_val_compare_and_swap(ptr, oval, nval);
    	}
	int find_s_e_pair(const char *_path, uint32_t taxid)
	{
		
		//fprintf(stderr, "find s e pair :%s\t%u\n", _path, taxid);
		gzFile fp_z;
		
		kseq_t *trunk;
		
		fp_z = gzopen(_path,"r");

		if (!fp_z) {
			//fprintf(stderr,"FAIL TO OPEN FILE\n");	
			return FILE_OPEN_ERROR;
		}
		trunk = kseq_init(fp_z);
		

		//uint64_t mask = ~((uint64_t)0x3 <<((_kmer-1)<<1));
//#ifdef CONSIDER_BOTH_ORIENTATION
		//uint64_t mask_rev[] = {(uint64_t)0x3 << ((_kmer - 1)<<1),  (uint64_t)0x2 << ((_kmer - 1)<<1), (uint64_t)0x1 << ((_kmer - 1)<<1), 0};
//#endif
		//uint32_t gid = 0;	
		//uint32_t tid = atoi(taxid);
		
		//uint16_t lastChar;
		bool find_end = false;
		bool find_start = false;
		//set<uint64_t>::iterator it_set;
		khiter_t it_set;	
		//set<uint64_t>::const_iterator it_set_start_end = start_nodes.end();
		//set<uint64_t>::const_iterator it_set_end_end = end_nodes.end();
		//set<uint64_t>::const_iterator it_set_end_end = end_nodes.end();
		khiter_t it_set_end_end = kh_end(end_nodes); 
	
		//kbitr_t bitr;
		se_pair * p, *rc_p;	
		//map<uint64_t, struct end_tid>::iterator it_map, it_map_rc;
		//map<uint64_t, struct end_tid>::const_iterator it_map_end = s_e_pair.end();
		

		uint64_t mask_rev[] = {(uint64_t)0x3 << ((kmer_len - 1)<<1),  (uint64_t)0x2 << ((kmer_len - 1)<<1), (uint64_t)0x1 << ((kmer_len - 1)<<1), 0};
		//struct end_tid n_p, n_p_rc;
		se_pair t = {0,0,0};
		se_pair n_p, n_p_rc;
		while (kseq_read(trunk)>=0) {

			for(uint64_t i=0; i<trunk->seq.l;++i) {
				if (Bit[trunk->seq.s[i]] != 4) {
					uint64_t start = i;
					while(Bit[trunk->seq.s[++i]]!=  4 && i < trunk->seq.l);
					//uint64_t end = i - 1;
					if (i - start > kmer_len) {
						//intiate head 
						uint64_t key = transIntoBits(trunk->seq.s+start, kmer_len);
						uint64_t rcKey = transRCIntoBits(trunk->seq.s+start, kmer_len); 	
						//if (rcKey == 1886573852641024189) {
							//cout<<transIntoChars(rcKey, kmer_len)<<endl;
						//}
						if (!find_start) {
							//it_map = s_e_pair.find(key);
							t.s = key;
							p = kb_getp(uint64_t, node_pair, &t);
						       // key exist
							if (p) {
						      		n_p_rc.e = rcKey;
						       		find_start = true;		
						       
						       }	
							//if (it_map != it_map_end) {
								//n_p_rc.end_node = rcKey;
								//find_start = true;	
							//}
						}
					       	if (!find_end) {
							it_set = kh_get(64, end_nodes, key);
							if (it_set != it_set_end_end) {
								n_p.e = key;
								t.s = rcKey;
								rc_p = kb_getp(uint64_t, node_pair, &t);
								find_end = true;
							}
							//it_set = end_nodes.find(key);
							//if (it_set != it_set_end_end) {
								//n_p.end_node = key;
								//it_map_rc = s_e_pair.find(rcKey);	
								//find_end = true;
							//}
						}
						if (find_start && find_end) {
							//cout<<"find start and end"<<endl;	
							//cout<<"s_e_node1:"<<key<<"\t"<<rcKey<<endl;	
								//record taxid for recom use
							//if (it != s_e_pair.end()) {
								//it->second = taxid;
							//} else 
								//it->second = LCA(it->second, taxid);
							//it_map->second.end_node = n_p.end_node;
							p->e = n_p.e;
							//it_map_rc->second.end_node = n_p_rc.end_node;
							rc_p->e = n_p_rc.e;
							uint32_t *old_tid = &p->tid;
							uint32_t *old_rc_tid = &rc_p->tid;
							//uint32_t *old_tid = &it_map->second.taxon_id;
							//uint32_t *old_rc_tid = &it_map_rc->second.taxon_id;	
							uint32_t _old_tid = *old_tid, otid;
							uint32_t _old_rc_tid = *old_rc_tid;	
							do {
								otid = _old_tid;
								uint32_t ntid = _evo_tree->LCA(taxid, otid);
								_old_tid = cas(old_tid,otid, ntid); 
							} while (_old_tid != otid);		

							do {
								otid = _old_rc_tid;
							uint32_t ntid = _evo_tree->LCA(taxid, otid);
								_old_rc_tid = cas(old_rc_tid,otid, ntid); 
							} while (_old_rc_tid != otid);	
							find_start = false;
							find_end = false;
						}	
						
						for (uint64_t j = start + 1; j < i - kmer_len + 1; ++j) {
							//cout<<trunk->seq.s[j+kmer_len - 1]<<endl;
							key = ((key & out_mask) << 2)| Bit[trunk->seq.s[j+kmer_len-1]];
							rcKey = (rcKey >> 2) | mask_rev[Bit[trunk->seq.s[j + kmer_len - 1]]];
							
							//cout<<key<<"\t"<<rcKey<<endl;	
							//cout<<111<<key<<"\t"<<rcKey<<endl;	
							//if (rcKey == 1886573852641024189) {
								//cout<<transIntoChars(rcKey, kmer_len)<<endl;
							//}
							if (!find_start) {
								//it_map = s_e_pair.find(key);
								//if (it_map != it_map_end) {
									//n_p_rc.end_node = rcKey;
									//find_start = true;	
								
								//}
								t.s = key;
								p = kb_getp(uint64_t, node_pair, &t);
							       // key exist
								if (p) {
									n_p_rc.e = rcKey;
									find_start = true;		
				       
					       			}	
							}
							if (!find_end) {
								it_set = kh_get(64, end_nodes, key);
								if (it_set != it_set_end_end) {
									n_p.e = key;
									t.s = rcKey;
									rc_p = kb_getp(uint64_t, node_pair, &t);
									find_end = true;
								}
								//it_set = end_nodes.find(key);
								//if (it_set != it_set_end_end) {
									//n_p.end_node = key;
									//it_map_rc = s_e_pair.find(rcKey);	
									//find_end = true;
								//}
							}
							if (find_start && find_end) {
								//cout<<"find start and end"<<endl;	
								//cout<<"s_e_node2:"<<key<<"\t"<<rcKey<<endl;	
								//it_map->second.end_node = n_p.end_node;
								//it_map_rc->second.end_node = n_p_rc.end_node;
								//uint32_t *old_tid = &it_map->second.taxon_id;
								//uint32_t *old_rc_tid = &it_map_rc->second.taxon_id;	
								p->e = n_p.e;
								//it_map_rc->second.end_node = n_p_rc.end_node;
								rc_p->e = n_p_rc.e;
								uint32_t *old_tid = &p->tid;
								uint32_t *old_rc_tid = &rc_p->tid;
								uint32_t _old_tid = *old_tid, otid;
								uint32_t _old_rc_tid = *old_rc_tid;	
								do {
									otid = _old_tid;
									uint32_t ntid = _evo_tree->LCA(taxid, otid);
									_old_tid = cas(old_tid,otid, ntid); 
								} while (_old_tid != otid);	

								do {
									otid = _old_rc_tid;
									uint32_t ntid = _evo_tree->LCA(taxid, otid);
									_old_rc_tid = cas(old_rc_tid,otid, ntid); 
								} while (_old_rc_tid != otid);		
								find_start = false;
								find_end = false;
							}	
						}
							//key = ((key & out_mask) << 2)| Bit[trunk->seq.s[i-1]];
							
					} else {
						
						if (i - start == kmer_len) {
							uint64_t key = transIntoBits(trunk->seq.s + start, kmer_len);
							uint64_t rcKey = transRCIntoBits(trunk->seq.s+start, kmer_len); 	
								
							//cout<<111<<key<<"\t"<<rcKey<<endl;	
							//cout<<"s_e_node3:"<<key<<"\t"<<rcKey<<endl;	
							//cout<<key<<"\t"<<rcKey<<endl;	
							//it_map = s_e_pair.find(key);
							t.s = key;
							p = kb_getp(uint64_t, node_pair, &t);
							t.s = rcKey;
							rc_p = kb_getp(uint64_t, node_pair, &t);
							//it_map_rc = s_e_pair.find(rcKey);	
							//it_map->second.end_node = key;
							p->e = key;
							//it_map_rc->second.end_node = rcKey;
							rc_p->e = rcKey;
							//uint32_t *old_tid = &it_map->second.taxon_id;
							//uint32_t *old_rc_tid = &it_map_rc->second.taxon_id;	
							uint32_t *old_tid = &p->tid;
							uint32_t *old_rc_tid = &rc_p->tid;
							uint32_t _old_tid = *old_tid, otid;
							uint32_t _old_rc_tid = *old_rc_tid;	
							do {
								otid = _old_tid;
								uint32_t ntid = _evo_tree->LCA(taxid, otid);
								_old_tid = cas(old_tid,otid, ntid); 
							} while (_old_tid != otid);		

							do {
								otid = _old_rc_tid;
								uint32_t ntid = _evo_tree->LCA(taxid, otid);
								_old_rc_tid = cas(old_rc_tid,otid, ntid); 
							} while (_old_rc_tid != otid);		

							find_start = false;
							find_end = false;
						}	
					
					}

				}
			}
		}
	
		kseq_destroy(trunk);
		gzclose(fp_z);	
		//fprintf(stderr,"find -s -e endingggggg");
		return NORMAL_EXIT;
	}	

			

	public:
	//map<uint64_t, struct end_tid> s_e_pair;
	KBTREE_INIT(uint64_t, se_pair, _SE_PAIR_CMP)
	kbtree_t(uint64_t) *node_pair;	
	bool isNotSingle[32];
	se_node(const char *refList, const uint8_t _kmer_len, evo_tree *_evo_tree_) {
	    
	    fileNum = 0;
	    FILE *fp = fopen(refList, "r");
	    if (fp) {
		char s_taxid[1024];
		char path[1025] = {0};	
		while (fscanf(fp, "%s %s", path,s_taxid) != EOF) {
		
			//file_arg.push_back(path_start + path_len);
			fileNames += path;
			fileNames += '\0';
			//path_len += strlen(path);	
			
			taxons.push_back((uint32_t)strtoull(s_taxid, NULL, 10));
		}
		const char *path_start = fileNames.c_str();	
		const char *path_end = path_start + fileNames.size() - 1;
		
		fileNameList.push_back(path_start);
		fileNum++;	
		for (; path_start < path_end; ++path_start) {
			if ( 0 == *path_start ) {
				fileNameList.push_back(path_start + 1);	
				fileNum++;
			}
		}
			
		fclose(fp);	
	    } else
		   error("failed to open file list");

			//fprintf(stderr,"%s", fileNames.c_str());
	    for (int i = 0; i < 32; ++i) {
	    	if ((i & 0x10) || countNonZero(i)) 
			isNotSingle[i] = true;
	       else 
		       isNotSingle[i] = false;	       
	    }

	    kmer_len = _kmer_len;
		//uint64_t mask = ~((uint64_t)0x3 <<((kmer_len-1)<<1));
		out_mask = ~((uint64_t)0x3 <<((kmer_len-1)<<1));
	    //out_mask = (((uint64_t) 1) << (2*kmer_len- 1)) - 1;
	    in_shift = (kmer_len - 1) << 1; 
	    _evo_tree = _evo_tree_;
		
	    node_pair = kb_init(uint64_t, KB_DEFAULT_SIZE);
	    if (!node_pair) {
	    	fprintf(stderr, "fail to init nodes");
		exit(1);
	    } 
	
	    end_nodes = kh_init(64);
	    if (!end_nodes) {
	    	fprintf(stderr, "fail to init nodes");
		exit(1);
	    } 
	}
	
	int insert_start(uint64_t k)
	{
		se_pair n = {k, 0, 0};
		if(!kb_getp(uint64_t, node_pair, &n))
		       kb_putp(uint64_t, node_pair, &n); 	
		return NORMAL_EXIT;
	}

	int insert_end(uint64_t k)
	{
		int ret;
		kh_put(64, end_nodes, k, &ret);
		return ret;
	}
	int pro_se_node(uint64_t kmer, uint64_t degree, uint64_t rcKmer, uint64_t rc_degree) 
	{
		uint8_t in_degree = degree >> 5;
		uint8_t out_degree = degree & 0x1F;
		uint8_t rc_in_degree = rc_degree >> 5;
		uint8_t rc_out_degree = rc_degree & 0x1F;
		if (isNotSingle[in_degree]) {
			//fprintf(stderr,"insert kmer as start node\n");
			insert_start(kmer);
			insert_end(rcKmer);
			for (uint64_t i = 0; i < 4; ++i) {
				if (in_degree & 0x1) {
					insert_end((kmer >> 2) | (i << in_shift));	
				}
				in_degree >>= 1;	
			}
			for (uint64_t i = 0; i < 4; ++ i) {
				if (rc_out_degree & 0x1) {
					insert_start(((rcKmer & out_mask ) << 2)| i );	
				}
				rc_out_degree >>= 1;
			}
		} 
		if (isNotSingle[out_degree]) {
			//fprintf(stderr,"insert kmer as end node\n");
			insert_start(rcKmer);
			insert_end(kmer);
			for (uint64_t i = 0; i < 4; ++i) {
				if (out_degree & 0x1) {
					insert_start(((kmer & out_mask ) << 2)| i );	
				
				}	
				out_degree >>= 1; 
					
			}
			for (uint64_t i = 0; i < 4; ++ i) {
				if (rc_in_degree & 0x1) {
					insert_end((rcKmer >> 2) | (i << in_shift));	
				}
				rc_in_degree >>= 1;
			}

		}	
		return NORMAL_EXIT;	
	}
	
	void run(void *)
	{
		//struct end_tid et = {0, 0};
		//set<uint64_t>::iterator it;
		
		//fprintf(stderr,"%lu\t%lu\n", start_nodes.size(), end_nodes.size());		
		//for ( it = start_nodes.begin(); it != start_nodes.end(); ++it) {
			//fprintf(stderr,"%lu", *it);
			//s_e_pair[*it] = et;
		//}	
		//start_nodes.clear();
		
		cerr<<kb_size(node_pair)<<endl;
		cerr<<kh_size(end_nodes)<<endl;
		//kbitr_t tmp;
		       //se_pair *p;	
		//kb_itr_first(uint64_t, node_pair, &tmp);
		//while (kb_itr_valid(&tmp)) {
			//p = &kb_itr_key(se_pair, &tmp);
			//printf("%lu\t%u\n", p->e, p->tid);
			//kb_itr_next(uint64_t, node_pair, &tmp);	
		//}

		//for ( it = end_nodes.begin(); it != end_nodes.end(); ++it) {
			//fprintf(stderr,"endNode:%lu\n", *it);
			//s_e_pair[*it] = et;
		//}	
		//start_nodes.clear();
		pthread_t thread_t[SE_THREAD_NUM];
		pthread_rwlock_init(&rwlock, NULL);
		for (int i = 0; i < SE_THREAD_NUM; ++i) {
			pthread_create(&thread_t[i], NULL, static_manager, this);	
		}
		//static_manager(this);		
		for (int i = 0; i < SE_THREAD_NUM; ++i) {
			//fprintf(stderr,"hahahh");
			pthread_join(thread_t[i], NULL);	
		}

		if (end_nodes) kh_destroy(64, end_nodes);
		if (fileNames.size()) fileNames.clear();
		if (fileNameList.size()) fileNameList.clear();
		if (taxons.size()) taxons.clear();
	}

};

#endif
