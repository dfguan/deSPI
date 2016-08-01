// =====================================================================================
//
//       Filename:  index.cpp
//
//    Description:  index reference
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

#include "preprocess.hpp"

#include "bwt.hpp"

#define LEN_LIMIT 2000

#define N_NEEDED 50000

map<uint32_t, uint32_t> taxonomyTree;

extern uint8_t _kmer;

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

const uint8_t rev[]={
	78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,
	78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,
	78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,
	78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,
	78,84,78,71,78,78,78,67,78,78,78,78,78,78,67,78,
	78,78,78,78,65,78,78,78,78,78,78,78,78,78,78,78,
	78,84,78,71,78,78,78,67,78,78,78,78,78,78,67,78,
	78,78,78,78,65,78,78,78,78,78,78,78,78,78,78,78
};

uint64_t total_sequences = 0;
uint64_t total_bps = 0;

int revComRead(char* str, int len, char *rstr)
{
	for (int i=0; i<len; ++i) rstr[i] = rev[str[len- 1 - i]];
	return NORMAL_EXIT;	

}

typedef struct {
	uint32_t 	tid;
	char 		*read_Name;
	char 		isClassified;
	char 		algn_r[2048];	
} cly_r;


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
uint32_t max_hit_tid(map<uint32_t, uint32_t> &hit_times) {

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
uint32_t assignedTID[1024];
uint8_t byteFormat[1024]; 
int _interval; 
int intervals[1024];
int sec_intervals[1024];
int transIntoBytes( char *str, int len) 
{
	for(int i=0; i <len; ++i) {
		byteFormat[i] = Bit3[str[i]];
	}
	byteFormat[len] = 4;
	return 0;
}	
int splitBytes(uint8_t *bytes, int len, int thres)
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

int classify_seq_2(kseq_t *trunk, int trunkNum,  bwt& bt, cly_r *results)
{
		
	int read_len;
	
	char *useRead;

	char *qual;

	int match_len;
	int  matchNum;	
	int counterIntervals;
	for (int i=0; i<trunkNum; ++i) {
		read_len = trunk[i].seq.l;
		useRead = trunk[i].seq.s;
		qual = trunk[i].qual.s;
		map<uint32_t, uint32_t> hit_times; 
		transIntoBytes(useRead, read_len);
		//cout<<useRead<<endl;i
		counterIntervals = splitBytes(byteFormat, read_len, _kmer);
		//cout<<counterIntervals<<endl;
		for (int j=0; j<counterIntervals;++j) {
			int start = intervals[(j<<1)];
			int p_len = intervals[(j<<1) + 1];
			//cout<<p_len<<"\t"<<start<<endl;
			do {
				matchNum = bt.exactMatch(byteFormat + start,   p_len, match_len, assignedTID);
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



int classify_seq(kseq_t *trunk, int trunkNum,  bwt& bt, cly_r *results, int offset)
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
		qual = trunk[i].qual.s;
		map<uint32_t, uint32_t> hit_times; 
		transIntoBytes(useRead, read_len);
		//cout<<useRead<<endl;i
		counterIntervals = splitBytes(byteFormat, read_len, _kmer);
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
					matchNum = bt.exactMatch(byteFormat + start,   p_len, match_len, assignedTID);
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
		cly_r *results_bulk = results + (i << 1) + offset;
		int buff_pos = 0;
		
		results_bulk->read_Name = trunk[i].name.s;

		//hit_times[1286170] = 29;
		map<uint32_t, uint32_t>::iterator it = hit_times.begin();
		while (buff_pos < 2026 && it != hit_times.end() ) {
			//buff_pos += sprintf(results_bulk->algn_r + buff_pos, "%u:%u ", it->first, it->second );
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

int output_results(cly_r *results, int n_results, bool isPairEnd)
{
	if (isPairEnd) {
	
		for(int i=0; i<n_results; ++i) {
			int ind = i << 1;
			cout<<results[ind].isClassified<<"\t"<<results[ind].read_Name<<"\t"<<results[ind].tid<<"\t"<<results[ind].algn_r<<endl;
			++ind;
			cout<<results[ind].isClassified<<"\t"<<results[ind].read_Name<<"\t"<<results[ind].tid<<"\t"<<results[ind].algn_r<<endl;
		}
	
	} else {
		for(int i=0; i<n_results; ++i) {
			int ind = i << 1;
			cout<<results[ind].isClassified<<"\t"<<results[ind].read_Name<<"\t"<<results[ind].tid<<"\t"<<results[ind].algn_r<<endl;
		}
	
	}
	
	return NORMAL_EXIT;


}
int main(int argc, char *argv[])
{
	int choice = atoi(argv[1]);
	bool isPairEnd = false;

	if (choice) {
	
		char *kmerPath = argv[2];
		char *refPath = argv[3];

		char *taxonomyNodesPath = argv[4];
		char *giTaxidPath = argv[5];
		char *dirPath = argv[6];
		

		//vector<uint64_t> fKmer;
		_kmer = (uint8_t) atoi(argv[7]); 
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
		//uint64_t sp,ep;	
		//bt.exactMatch("GCTTCGCTGTTATTGGCACCAATTGGATCAC",31, sp,ep);
	} else {
	
		char *readPath;
	
		
		char *taxonomyNodesPath;
		
		char * dirPath;
		
		char *readPath_s;	
		
		
		fprintf(stderr,"%d", argc);	
		if (argc > 7 ) {
			isPairEnd = true;
		
			readPath = argv[2];
		
			readPath_s = argv[3];	
			
			taxonomyNodesPath = argv[4];
			
			dirPath = argv[5];
			_kmer = (uint8_t) atoi(argv[6]);
			
			_interval = atoi(argv[7]);
	
		} else {
		
			readPath = argv[2];
			
			taxonomyNodesPath = argv[3];
			
			
			dirPath = argv[4];
			
			_kmer = (uint8_t) atoi(argv[5]);
			
			_interval = atoi(argv[6]);
		
		}


		--_kmer;
		
		bwt bt(_kmer);
		
		fprintf(stderr,"loading index\n");

		bt.load_info(dirPath);
		
		taxonTree(taxonomyNodesPath);
		//map<uint32_t, uint32_t>::iterator it = taxonomyTree.begin();
		//while (it!=taxonomyTree.end()) {
		//	cout<<it->first<<"\t"<<it->second<<endl;

		//	++it;
		//}
		//cout<<taxonomyTree.size()<<endl;
		fprintf(stderr,"classifying...\n");
		
		gzFile fp;
		
		//uint32_t *nKmerTaxonID = bt.taxonIDTab;	
		
		fp = gzopen(readPath, "r");

		if (!fp) return FILE_OPEN_ERROR;

		kstream_t *_fp = ks_init(fp);

		kseq_t *seqs = (kseq_t *) calloc(N_NEEDED, sizeof(kseq_t));

		if (!seqs) return MEM_ALLOCATE_ERROR;
		//
		//parameters for pair-end reads
		gzFile fp_s;
		kstream_t *_fp_s;
		kseq_t *seqs_s;

		if (isPairEnd) {
		
			fp_s = gzopen(readPath_s, "r");

			if (!fp_s) return FILE_OPEN_ERROR;

			_fp_s  = ks_init(fp_s);

			seqs_s = (kseq_t *) calloc(N_NEEDED, sizeof(kseq_t));
			
			if (!seqs_s) return MEM_ALLOCATE_ERROR;
		
		
		}

		if (!seqs) return MEM_ALLOCATE_ERROR;
		
		cly_r *results = (cly_r *)calloc(2 * N_NEEDED, sizeof(cly_r));




  		struct timeval tv1, tv2;
		
		int n_seqs;
		total_sequences = 0;
		gettimeofday(&tv1,NULL);
		if (isPairEnd) {
		
			while ((n_seqs = read_reads(_fp, seqs, N_NEEDED) )> 0 && read_reads(_fp_s, seqs_s, N_NEEDED)) {
				
				classify_seq(seqs, n_seqs , bt, results, 0);
				classify_seq(seqs_s, n_seqs, bt, results, 1);
				output_results(results, n_seqs, isPairEnd);
				total_sequences += n_seqs;
			}	
		
		
		} else {
		
			while ((n_seqs = read_reads(_fp, seqs, N_NEEDED) )> 0) {
			
				classify_seq(seqs, n_seqs , bt, results, 0);
				output_results(results, n_seqs, isPairEnd);
				total_sequences += n_seqs;

			}	
		
		}	
		gettimeofday(&tv2,NULL);
		report_stats(tv1,tv2);

		//fprintf(stderr,"%f seconds\n",((float)t)/CLOCKS_PER_SEC);
		if (results) free(results);
		if (seqs) free(seqs);

	}
	
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
