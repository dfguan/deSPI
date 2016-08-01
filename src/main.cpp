// =====================================================================================
//
//       Filename:  main.cpp
//
//    Description:  index reference and classify reads
//
//        Version:  1.0
//        Created:  11/23/2015 08:51:10 AM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Dengfeng Guan, <dfguan@hit.edu.cn>
//   Organization:  Harbin Institute of Technology
//
// =====================================================================================

#include <sys/time.h>
#include <pthread.h>
#include "preprocess.hpp"

#include "bwt.hpp"

#include "ui.hpp"

#define LEN_LIMIT 2000

#define N_NEEDED 50000

map<uint32_t, uint32_t> taxonomyTree;

extern uint8_t _kmer;

int read_seq;
pthread_rwlock_t rwlock;
const uint8_t Bit3[] = {
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

//const uint8_t rev[]={
	//78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,
	//78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,
	//78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,
	//78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,
	//78,84,78,71,78,78,78,67,78,78,78,78,78,78,67,78,
	//78,78,78,78,65,78,78,78,78,78,78,78,78,78,78,78,
	//78,84,78,71,78,78,78,67,78,78,78,78,78,78,67,78,
	//78,78,78,78,65,78,78,78,78,78,78,78,78,78,78,78
//};

uint64_t total_sequences = 0;
uint64_t total_bps = 0;

//int revComRead(char* str, int len, char *rstr)
//{
	//for (int i=0; i<len; ++i) rstr[i] = rev[str[len- 1 - i]];
	//return NORMAL_EXIT;	

//}

typedef struct {
	uint32_t 	tid;
	char 		*read_Name;
	char 		isClassified;
	char 		algn_r[2048];	
} cly_r;

typedef struct {
	uint8_t threadid;
	uint8_t *_byteformat;
	//int trunkNum;	
	int*    _interv;
	uint32_t* _assignedTID;
	int 	n_seqs;
	kseq_t * seqs_poi;
	bwt 	*p_bwt;
	cly_r  *res;
} thread_aux;
uint32_t resolve_tree(map<uint32_t, uint32_t> &hit_times)
 {
	 //if (hit_times.size() == 1) return hit_times.begin()->first;
	 set<uint32_t> max_taxa;
	 uint32_t max_taxon = 0, max_score = 0;
	 map<uint32_t, uint32_t>::iterator it = hit_times.begin();
	 
	 // Sum each taxon's LTR path
	 while (it != hit_times.end()) {
		 uint32_t taxon = it->first;
		 uint32_t node = taxon;
		 uint32_t score = 0;
		 while (node > 0) {
			 //cout<<node<<endl;
			 score += hit_times[node];
			 node = taxonomyTree[node];
		 }

		 if (score > max_score) {
			 max_taxa.clear();
			 max_score = score;
			 max_taxon = taxon;
		 }
		 else if (score == max_score) {
			 if (max_taxa.empty())
				 max_taxa.insert(max_taxon);
			 max_taxa.insert(taxon);
		 }

		 ++it;
	 }

	 // If two LTR paths are tied for max, return LCA of all
	 if (! max_taxa.empty()) {
		 set<uint32_t>::iterator sit = max_taxa.begin();
		 max_taxon = *sit;
		 for (sit++; sit != max_taxa.end(); sit++)
			 max_taxon = LCA(max_taxon, *sit);
	 }

	 return max_taxon;
}
uint32_t max_hit(map<uint32_t, uint32_t> &hit_times) {

	 set<uint32_t> max_taxa;
	 uint32_t max_taxon = 0, max_score = 0;
	 map<uint32_t, uint32_t>::iterator it = hit_times.begin();
	 
	 // Sum each taxon's LTR path
	 while (it != hit_times.end()) {
		 uint32_t taxon = it->first;
		 uint32_t score = hit_times[taxon];

		 if (score > max_score) {
			 max_taxa.clear();
			 max_score = score;
			 max_taxon = taxon;
		 }
		 else if (score == max_score) {
			 if (max_taxa.empty())
				 max_taxa.insert(max_taxon);
			 max_taxa.insert(taxon);
		 }

		 ++it;
	 }
	 // If two LTR paths are tied for max, return LCA of all
	 if (! max_taxa.empty()) {
		 set<uint32_t>::iterator sit = max_taxa.begin();
		 max_taxon = *sit;
		 for (sit++; sit != max_taxa.end(); sit++)
			 max_taxon = LCA(max_taxon, *sit);
	 }

	 return max_taxon;

}
int _interval; 
//int sec_intervals[1024];
int transIntoBytes( char *str, int len, uint8_t *byteFormat) 
{
	for(int i=0; i <len; ++i) {
		byteFormat[i] = Bit3[str[i]];
	}
	byteFormat[len] = 4;
	return 0;
}	
int splitBytes(uint8_t *bytes, int len, int thres, int* intervals)
{
	int index = 0;
	int counter = 0;
	int start = 0;
	for (int i=0; i<len; ++i) {
		if (bytes[i] != 4) {
			start = i;
			while (bytes[++i]!=4 && i < len);
				
			int t_len = i - start;
			if (t_len > thres) {
				intervals[index++] = start;
				intervals[index++] = t_len;	
				++counter;
		}
		}

	}

	return counter;
}

int classify_seq_2(kseq_t *trunk, int trunkNum,  bwt* bt, cly_r *results, int *intervals, uint8_t *byteFormat, uint32_t* assignedTID)
{
		
	int read_len;
	
	char *useRead;

	char *qual;

	int match_len;
	int  matchNum;	
	int counterIntervals;
	_interval =1; 	
	for (int i=0; i<trunkNum; ++i) {
		read_len = trunk[i].seq.l;
		useRead = trunk[i].seq.s;
		qual = trunk[i].qual.s;
		map<uint32_t, uint32_t> hit_times; 
		transIntoBytes(useRead, read_len, byteFormat);
		//cout<<useRead<<endl;i
		counterIntervals = splitBytes(byteFormat, read_len, _kmer, intervals);
		//cout<<counterIntervals<<endl;
		for (int j=0; j<counterIntervals;++j) {
			int start = intervals[(j<<1)];
			int p_len = intervals[(j<<1) + 1];
			//cout<<p_len<<"\t"<<start<<endl;
			do {
				matchNum = bt->exactMatch(byteFormat + start,   p_len, match_len, assignedTID);
				//cout<<match_len<<"\t"<<mismatch<<endl;	
				if (matchNum) {
					while (matchNum) hit_times[assignedTID[--matchNum]] += match_len;
					p_len -= (match_len - _kmer);
				} else 
					p_len -= _interval;
				//cout<<match_len<<"\t"<<p_len<<"\t"<<spacedLen<<endl;
			} while (p_len > _kmer);
		
		}	
			//cout<<endl;
		cly_r *results_bulk = results + i;
		int buff_pos = 0;
		
		results_bulk->read_Name = trunk[i].name.s;


		map<uint32_t, uint32_t>::iterator it = hit_times.begin();
		while (it != hit_times.end()) {
			buff_pos += sprintf(results_bulk->algn_r + buff_pos, "%u:%u ", it->first, it->second );
			++it;
		}
		results_bulk->algn_r[buff_pos] = '\0';
		//cout<<endl;
		uint32_t lab = resolve_tree(hit_times);

		if (lab) {
			results_bulk->tid = lab;
			results_bulk->isClassified = 'C';
		} else {
			results_bulk->tid = 0;
			results_bulk->isClassified = 'U';
		}	
	}
	//fprintf(stderr,"seconds %f\n",(float)tk/CLOCKS_PER_SEC);

	return NORMAL_EXIT;
}

//int spacedLeng[] = {0,5,7,11,17};


//int classify_seq(kseq_t *trunk, int trunkNum,  bwt& bt, cly_r *results, int offset)

int classify_seq(kseq_t *trunk, int trunkNum,  bwt* bt, cly_r *results, int *intervals, uint8_t *byteFormat, uint32_t* assignedTID)
{
		
	int read_len;
	
	char *useRead;

	char *qual;

	int match_len;
	int  matchNum;	
	//int spacedLen;
	int counterIntervals;
	int hit_counter;
	//fprintf(stderr,"dddddd");	
	for (int i=0; i<trunkNum; ++i) {
		read_len = trunk[i].seq.l;
		useRead = trunk[i].seq.s;
		//fprintf(stderr,"%d\n",read_len);
		qual = trunk[i].qual.s;
		map<uint32_t, uint32_t> hit_times; 
		transIntoBytes(useRead, read_len, byteFormat);
		//cout<<useRead<<endl;i
		counterIntervals = splitBytes(byteFormat, read_len, _kmer, intervals);
		//cout<<counterIntervals<<endl;
		int spacedLen = 20;
		hit_counter = 0;
		//int ind = 4;
		while (spacedLen) {
			//spacedLen = spacedLeng[ind];
			for (int j=0; j<counterIntervals;++j) {
				int start = intervals[(j<<1)];
				int p_len = intervals[(j<<1) + 1];
				//cout<<p_len<<"\t"<<start<<endl;
				do {                                    
					matchNum = bt->exactMatch(byteFormat + start,   p_len, match_len, assignedTID);
					//cout<<match_len<<"\t"<<mismatch<<endl;	
						
					//cout<<match_len<<"\t"<<p_len<<"\t"<<spacedLen<<endl;
					//if (matchNum) {
						//++hit_counter;
						
						while (matchNum) {
							++hit_counter;
							hit_times[assignedTID[--matchNum]] += match_len;
							//cout<<assignedTID[matchNum]<<"\t";
						}
					
					//}
					//cout<<endl;
					if (match_len < spacedLen) p_len -= spacedLen;
					else p_len -= (match_len/spacedLen)*spacedLen;
				} while (p_len > _kmer);
			}	
			if (hit_counter)
				break;
			else 
				//--ind;
				spacedLen -= 5;
		
		}	
			//cout<<endl;
		cly_r *results_bulk = results + i;
		int buff_pos = 0;
		
		results_bulk->read_Name = trunk[i].name.s;

		//hit_times[1286170] = 29;
		map<uint32_t, uint32_t>::iterator it = hit_times.begin();
		while (buff_pos < 2026 && it != hit_times.end() ) {
			buff_pos += sprintf(results_bulk->algn_r + buff_pos, "%u:%u ", it->first, it->second );
			++it;
		}
		results_bulk->algn_r[buff_pos] = '\0';
		//cout<<endl;
		uint32_t lab = resolve_tree(hit_times);

		if (lab) {
			results_bulk->tid = lab;
			results_bulk->isClassified = 'C';
		} else {
			results_bulk->tid = 0;
			results_bulk->isClassified = 'U';
		}	
	}
	//fprintf(stderr,"seconds %f\n",(float)tk/CLOCKS_PER_SEC);

	return NORMAL_EXIT;
}

void report_stats(struct timeval time1, struct timeval time2) {
	  time2.tv_usec -= time1.tv_usec;
	  time2.tv_sec -= time1.tv_sec;
	  if (time2.tv_usec < 0) {
	    time2.tv_sec--;
	    time2.tv_usec += 1000000;
	  }
	  double seconds = time2.tv_usec;
	  seconds /= 1e6;
	  seconds += time2.tv_sec;
	  //cerr << "\r";
	  fprintf(stderr, 
		  "%llu sequences (%.2f Mbp) processed in %.3fs (%.1f Kseq/m, %.2f Mbp/m).\n",
		  (unsigned long long) total_sequences, 0 / 1.0e6, seconds,
		  total_sequences / 1.0e3 /(seconds / 60),0);
          //total_bases / 1.0e6 / (seconds / 60) );
 // fprintf(stderr, "  %llu sequences classified (%.2f%%)\n",
    //      (unsigned long long) total_classified, total_classified * 100.0 / total_sequences);
  //fprintf(stderr, "  %llu sequences unclassified (%.2f%%)\n",
      //    (unsigned long long) (total_sequences - total_classified),
       //   (total_sequences - total_classified) * 100.0 / total_sequences);
}

int read_reads(kstream_t *_fp, kseq_t *_seqs, int n_needed)
{
	kseq_t *temp = _seqs;

	int i = 0;
	//temp[i].f = _fp;
	//temp[i].last_char = temp[n_needed - 1].last_char;
	
	while(i <n_needed && (temp[i].f = _fp) && kseq_read(temp+i)>=0 ) 
		//int z = 0;
		//fprintf(stderr, "%s\t%d\n",temp[i].name.s, i);
		++i;
		//fprintf(stderr,"%d\n",temp[i].last_char);
	
	return i;

}

int output_results(cly_r *results, int n_results)
{
	int ind;

	for(int i=0; i<n_results; ++i) {
		ind = i;
		cout<<results[ind].isClassified<<"\t"<<results[ind].read_Name<<"\t"<<results[ind].tid<<"\t"<<results[ind].algn_r<<endl;
	}
	return NORMAL_EXIT;


}
int build_index(opts* p_opt)
{

	const char *kmerPath = p_opt->sortedKmer.c_str();
	const char *refPath = p_opt->ref.c_str();

	const char *taxonomyNodesPath = p_opt->tids.c_str();
	const char *giTaxidPath = p_opt->gids.c_str();
	const char *dirPath = p_opt->output.c_str();
	

	//vector<uint64_t> fKmer;
	_kmer = p_opt->kmer;

	string  bwt_s;
	vector<uint32_t> nKmerTaxonID;	
	
	fprintf(stderr,"start preprocessing......\n");
	

	uint64_t hash_index_size = (uint64_t)1 <<((PREINDEXLEN<<1) + 1);
	
	uint64_t *hash_index = new uint64_t[hash_index_size]();

	preprocess(refPath, kmerPath, taxonomyNodesPath, giTaxidPath,  bwt_s, nKmerTaxonID, hash_index);

	fprintf(stderr,"writing index....\n");
	//char *dirPath = ".";
	//
	bwt bt(bwt_s.c_str(), bwt_s.length(),hash_index);

	bt.bwt_init();

	bt.write_info(dirPath, nKmerTaxonID);
	
	return NORMAL_EXIT;


}

int thread_initiate(thread_aux *p_aux,int num_threads, bwt *_bt, kseq_t *p_seqs, cly_r *_res, uint8_t *p_byteformat, int *p_interv, uint32_t *p_tids)
{
	for(int i=0; i<num_threads; ++i) {
		p_aux[i].threadid = i;
		p_aux[i]._byteformat = p_byteformat;
		p_aux[i]._interv = p_interv;
		p_aux[i]._assignedTID = p_tids;		
		p_aux[i].seqs_poi = p_seqs;
		p_aux[i].p_bwt = _bt;
		p_aux[i].res = _res;
	}
	return NORMAL_EXIT;
}

static void *thread_worker(void *data)
{
	thread_aux *aux = (thread_aux*)data;
	int _read_seq;
	while (1) {
		pthread_rwlock_wrlock(&rwlock);
		_read_seq = read_seq++;
		pthread_rwlock_unlock(&rwlock);
		if (_read_seq < aux->n_seqs) {
			uint16_t off = (uint16_t)(aux->threadid) <<10;
			classify_seq(aux->seqs_poi+ _read_seq, 1,aux->p_bwt, aux->res+_read_seq, aux->_interv + off, aux->_byteformat+off,aux->_assignedTID + off );
		} else break; 
			
	}
	//uint16_t off = (uint16_t)(aux->threadid) <<10;
	//classify_seq(aux->seqs_poi,aux->trunkNum ,aux->p_bwt, aux->res, aux->_interv + off, aux->_byteformat+off,aux->_assignedTID + off );

}

int classify(opts *p_opt)
{
	
	string taxonomyTreePath = p_opt->lib + "/map";

	_kmer = --p_opt->seed;
	
	bwt *bt =new bwt(_kmer);
	
	fprintf(stderr,"loading index\n");

	bt->load_info(p_opt->lib.c_str());
	
	taxonTree(taxonomyTreePath.c_str());
	//map<uint32_t, uint32_t>::iterator it = taxonomyTree.begin();
	//while (it!=taxonomyTree.end()) {
	//	cout<<it->first<<"\t"<<it->second<<endl;

	//	++it;
	//}
	//cout<<taxonomyTree.size()<<endl;
	fprintf(stderr,"classifying...\n");
	
	gzFile fp;

	kstream_t *_fp;

	//uint32_t *nKmerTaxonID = bt.taxonIDTab;	
	

	kseq_t *seqs = (kseq_t *) calloc(N_NEEDED, sizeof(kseq_t));

	if (!seqs) return MEM_ALLOCATE_ERROR;
	

	cly_r *results = (cly_r *)calloc(N_NEEDED, sizeof(cly_r));

	if (!results) return MEM_ALLOCATE_ERROR;

	struct timeval tv1, tv2;
	
	int n_seqs;
	total_sequences = 0;
	gettimeofday(&tv1,NULL);
	


	int num_reads = p_opt->reads.size();	
	
	uint32_t assignedTID[p_opt->num_threads * 1024];
	uint8_t byteFormat[p_opt->num_threads * 1024]; 
	int intervals[p_opt->num_threads * 1024];
	
	if (p_opt->num_threads >1) {
		thread_aux *aux = new thread_aux[p_opt->num_threads];
		thread_initiate(aux,p_opt->num_threads,bt, seqs, results, byteFormat, intervals, assignedTID);
		for (int i=0; i<num_reads; ++i) {
			
			fp = gzopen(p_opt->reads[i].c_str(), "r");

			if (!fp) return FILE_OPEN_ERROR;

			_fp = ks_init(fp);
			
			while ((n_seqs = read_reads(_fp, seqs, N_NEEDED) )> 0) 	{
				read_seq = 0;
				//int average_read_num = n_seqs/p_opt->num_threads;
				//int pointer = 0;
				pthread_t *tid = (pthread_t *)calloc(p_opt->num_threads, sizeof(pthread_t));
				for (int j=0; j < p_opt->num_threads;++j) {
					//aux[j].trunkNum = average_read_num;
					//aux[j].seqs_poi = seqs + pointer;
					//aux[j].res = results + pointer;
					//pointer += average_read_num;
					aux[j].n_seqs = n_seqs;
					pthread_create(&tid[j], NULL, thread_worker, aux + j);
				}
					//aux[j].trunkNum = n_seqs - j*average_read_num;
					//aux[j].seqs_poi = seqs + pointer;
					//aux[j].res = results + pointer;
					//pointer += average_read_num;
					//pthread_create(&tid[j], NULL, thread_worker, aux + j);
				
				//classify_seq_2(seqs, n_seqs , bt, results );
				for (int j = 0; j < p_opt->num_threads; ++j) pthread_join(tid[j], 0);
				free(tid);
				output_results(results, n_seqs);
				total_sequences += n_seqs;

			}	
			
			gzclose(fp);	
				
		}
	} else {
		for (int i=0; i<num_reads; ++i) {
			
			fp = gzopen(p_opt->reads[i].c_str(), "r");

			if (!fp) return FILE_OPEN_ERROR;

			_fp = ks_init(fp);
			
			while ((n_seqs = read_reads(_fp, seqs, N_NEEDED) )> 0) 	{
			
				classify_seq(seqs, n_seqs, bt, results, intervals, byteFormat , assignedTID);
				//classify_seq_2(seqs, n_seqs , bt, results );
				output_results(results, n_seqs);
				total_sequences += n_seqs;

			}	
			
			gzclose(fp);	
				
		}
	
	}	

	
	gettimeofday(&tv2,NULL);
	report_stats(tv1,tv2);

	//fprintf(stderr,"%f seconds\n",((float)t)/CLOCKS_PER_SEC);
	if (results) free(results);
	if (seqs) free(seqs);
	
	return NORMAL_EXIT;

}

int main(int argc, char *argv[])
{
	opts *_opt = new opts;

	UI ui(_opt);	
	
	if (ui.opt_parse(argc,argv,_opt) == ERROR_PARSE_PARAMS)
		exit(1);

	if (_opt->isClassify) classify(_opt);
	else build_index(_opt);
		
	//char *st;
	
	//uint32_t *uts;

	//uint64_t z;
	//load_index(dirPath, st, uts, &z );		
	//cerr<<<<endl;
	//fprintf(stderr,"%s\n",st);
	
	//uint64_t sp, ep;

	//bwt b(st, uts, z);
	
	//b.bwt_init();
	//char *read = "GGCT";
	//cout<<b.exactMatch(read,4,sp,ep )<<endl;
	//cout<<sp<<"\t"<<ep<<endl;	
	//read reads file
	//output(kmerValue, kmerInfo, _2kmers);
	

	return NORMAL_EXIT;

}
