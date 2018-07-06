// =====================================================================================
//
//       Filename:  sort.cpp
//
//    Description: generate sorted Kmers
//
//        Version:  1.0
//        Created:  10/22/2015 03:23:13 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Dengfeng Guan, <dfguan@hit.edu.cn>
//   Organization:  Harbin Institute of Technology
//
// =====================================================================================


#include "index.hpp"
#include "jreader.hpp"

//#include "se_node.hpp"
#include "evo_tree.hpp"

#include <queue>


#define BUFFERSIZE (1<<20)
//#define KMER_LENGTH_PlusOne 32
#define BUCKET_LENGTH 12 //need less than KMER_LENGTH_PlusOne
#define BUCKET_CAPACITY (1<<(BUCKET_LENGTH<<1))
#define THREAD_NUM 32 //should not surpass 2^8-1
#define STKSIZ 40
#define FILE_NUM_LIMIT 256 
#define DG_BUFFERSIZE (1<<27) //128M

//#define KMER_LENGTH 17

uint64_t *countKmer = NULL;
uint64_t *hashKmer = NULL;
uint64_t trans[256];
uint64_t totalKmerNum;
uint64_t totalCounterNum;
uint64_t segCount[THREAD_NUM];
uint64_t multiFlag=1;
bool     *isStartNode;
int 	useFileNum;

uint64_t set_bit [] = {8, 4, 2 , 1, 16}; 
struct thread_aux{
	se_node *_s_e;
	sp_kmers *_s_p;	
	string  &bwt_s;
	vector<uint32_t> &tids;

};
struct thread_aux_2{
	se_node *_s_e;
	sp_kmers *_s_p;	
	bwt *c_bwt;
	uint8_t kmersize;
	uint64_t kmerN;
};
typedef struct element{
	uint64_t key;
	uint64_t val;	
	int      ind;
	bool operator<(const element &e) const{
		return key > e.key;	
	}
}ele;

//KBTREE_INIT(uint64_t, se_pair, _SE_PAIR_CMP)

//int main(int argc, char *argv[]) {
		
	//char *refListPath = argv[1];
	//char *s_kmer_len = argv[2];
	//char *kmersPath = argv[3];
	//char *evo_tree_path = argv[4];
	 //return genSortedKmers(refListPath, uint8_t(atoi(s_kmer_len)), kmersPath, evo_tree_path);
//}
//processKmers(refListpath, kmer, kmerpath, bwt_s, hash_index);

int printDegrees(uint64_t degrees)
{
	
	char bitChars[32] = {'A','C','G','T'};
	//bitChars[1] = 'A';
	//bitChars[2] = 'C';
	//bitChars[4] = 'G';
	//bitChars[8] = 'T';V	

	uint8_t out_degree = degrees & 0x1F;
	uint8_t in_degree = degrees >> 5;
	for (int i=0; i < 4; ++i) {
		if (in_degree & 0x1) {
			cout<<bitChars[i]<<",";
		}
			in_degree >>= 1;
	}
	cout<<"\t";

	for (int i=0; i < 4; ++i) {
		if (out_degree & 0x1) {
			cout<<bitChars[i]<<",";
		}
			out_degree >>= 1;
	}
	cout<<endl;
	return NORMAL_EXIT;	
}

int processKmers(const char *refListPath, const uint8_t kmerSize,const  char *disorderedKmerpath,const char *evo_tree_path, string &bwt_s, uint64_t *hash_index, vector<uint32_t> &nkmerTID)
{
    	//////////////////////////////////////kmercounting///////////////////////////////////////////
	countKmer=(uint64_t *)calloc(BUCKET_CAPACITY,sizeof(uint64_t));
	uint64_t bufferSize=BUFFERSIZE;//<<1;
	//uint64_t occBuf;
	uint64_t *writeBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t bufPoint=0;
	FILE *fpKmer=NULL;
	char *kmerPath = "database.srt";
	fpKmer=fopen(kmerPath ,"wb");
	uint64_t tempMove=(kmerSize - BUCKET_LENGTH)<<1; //32 defines kmer_plus_one
	fprintf(stderr,"start kmercounting\n");
	//jellyfish db reader from kraken src
	JReader jr(disorderedKmerpath);
	uint64_t val_len = jr.get_val_len();
	uint64_t key_ct = jr.get_key_ct();
	uint64_t key_len = jr.get_key_len();
	uint64_t pair_size = key_len + val_len;
	char pair[pair_size];


	ifstream input_file(jr.get_db_name().c_str(), std::ifstream::binary);
	if (!input_file) return FILE_OPEN_ERROR;
	
	input_file.seekg(jr.header_size(), ios_base::beg);
	//cout<<key_len<<"\t"<<val_len<<"\t"<<jr.header_size()<<"\t"<<key_ct<<endl;
	// Create a copy of the offsets array for use as insertion positions
	//FILE * input_file = fopen(jr.get_db_name().c_str(),"rb");
	
	evo_tree *_evo_tree = new evo_tree(evo_tree_path);
	se_node *_se_node = new se_node(refListPath, kmerSize, _evo_tree);
	//cout<<pair_size<<endl;
	for (uint64_t i = 0; i < key_ct; i++) {
		input_file.read(pair, pair_size);
		//for (int i = 0; i < pair_size; ++i)
			//cout<<std::hex<<pair[i];
		//cout<<endl;
		uint64_t kmerValue = 0;
		uint64_t rcKmerValue = 0;
		uint64_t occBuf = 0;
		uint64_t rcoccBuf = 0;
		memcpy(&kmerValue, pair, key_len);
		//cout<<kmerValue<<endl;
		//cout<<std::hex<<kmerValue<<endl;
				
		memcpy(&occBuf, pair + key_len, val_len);//copy occ buff	
		//if (kmerValue == 0)
			//fprintf(stderr,"%lx",occBuf);
		//cout<<std::hex<<kmerValue<<"\t"<<std::hex<<occBuf<<endl;
		//cout<<kmerValue<<"\t"<<occBuf<<endl;
		rcKmerValue = rcKmer(kmerValue,kmerSize);
		
		rcoccBuf = rcOcc(occBuf);	
		//if (kmerValue == 45587144286002457 )
			//cout<<transIntoChars(kmerValue, kmerSize)<<endl;
		//fprintf(stderr,"got it %lx", occBuf >> 32);		
		//if (rcKmerValue == 45587144286002457)	
	//V		fprintf(stderr,"got it rc: %lx", rcoccBuf >> 32);		
		if (kmerValue > rcKmerValue) {
			fprintf(stderr,"false format, use jellyfish regenerate kmers");
			exit(1); 
		}
		//if (rcKmerValue == 0)
			//fprintf(stderr,"%lx",occBuf);
		
		++countKmer[kmerValue >> tempMove];
		//cout<<kmerSize<<endl;	
		++countKmer[rcKmerValue >> tempMove];
		
		//fprintf(stderr, "%lu\t%lx\t%lx\t%lx\t%lx\n", i, kmerValue, rcKmerValue,occBuf, rcoccBuf );
		//cout<<kmerValue<<"\t"<<transIntoChars(kmerValue, kmerSize)<<endl;
		//cout<<rcKmerValue<<"\t"<<transIntoChars(rcKmerValue,kmerSize)<<endl;
		_se_node->pro_se_node( kmerValue, occBuf >> 32, rcKmerValue, rcoccBuf >> 32);	
		//cout<<kmerValue<<"\t"<<transIntoChars(kmerValue, kmerSize)<<endl;
		//cout<<rcKmerValue<<"\t"<<transIntoChars(rcKmerValue,kmerSize)<<endl;
				
		writeBuf[bufPoint++] = kmerValue;
		writeBuf[bufPoint++] = occBuf;
		writeBuf[bufPoint++] = rcKmerValue;
		writeBuf[bufPoint++] = rcoccBuf; 
		if(bufPoint>=bufferSize) {
			fwrite(writeBuf,sizeof(uint64_t),bufferSize,fpKmer);
			bufPoint=0;
		}
		//cout<<transIntoChars(kmer,21)<<endl;

	}

	if(bufPoint) 	fwrite(writeBuf,sizeof(uint64_t),bufPoint,fpKmer);

	fprintf(stderr,"done with cycling\n");
	/*
	while(fscanf(fp1,"%s%lu",seqBuf,&occBuf)!=EOF) {
		uint64_t seq=0;
		for (uint8_t i=0; i<BUCKET_LENGTH;++i)  seq = (seq << 2)| trans[seqBuf[i]];
		countKmer[seq]++;
		for (uint8_t i=BUCKET_LENGTH;i<kmerSize;++i) seq = (seq<<2) | trans[seqBuf[i]];

	}
	*/
	input_file.close();
	
	if (writeBuf)
		free(writeBuf);
	
	fclose(fpKmer);
	//fclose(fp1);
	//start sort and find start end point
	
	sp_kmers *_sp_kmers = new sp_kmers(kmerSize);

	pthread_t thread_1, thread_2;

		
	pthread_create(&thread_1, NULL, sort_thread, (void *)&kmerSize);	
	//sort_thread((void *)&kmerSize);	
	struct thread_aux ta = {_se_node, _sp_kmers, bwt_s, nkmerTID};
	pthread_create(&thread_2, NULL, call_special_kmers,(void *)&ta);	
	//call_special_kmers(&ta);	
	pthread_join(thread_1, NULL);
	pthread_join(thread_2, NULL);
	
	isStartNode = &_se_node->isNotSingle[0];
	mergeKmers("database.srt", kmerSize, _sp_kmers->_kspchar,_sp_kmers->_kspchar_counter, _se_node, bwt_s, hash_index, nkmerTID, bwt_s.size()-1); 
	
	if (_sp_kmers)
		delete _sp_kmers;
	if (_se_node)
		delete _se_node;
		
	return 0;
}

int index(const char *refListPath, const uint8_t kmerSize,const  char *disorderedKmerpath,const char *evo_tree_path, bwt *c_bwt)
{
    	//////////////////////////////////////kmercounting///////////////////////////////////////////
	countKmer=(uint64_t *)calloc(BUCKET_CAPACITY,sizeof(uint64_t));
	

	uint64_t tempMove=(kmerSize - BUCKET_LENGTH)<<1; //32 defines kmer_plus_one
	fprintf(stderr,"start kmercounting\n");




	//jellyfish db reader from kraken src
	JReader jr(disorderedKmerpath);
	uint64_t val_len = jr.get_val_len();
	uint64_t key_ct = jr.get_key_ct();
	uint64_t key_len = jr.get_key_len();
	uint64_t pair_size = key_len + val_len;
	char pair[pair_size];


	ifstream input_file(jr.get_db_name().c_str(), std::ifstream::binary);
	if (!input_file) return FILE_OPEN_ERROR;
	input_file.seekg(jr.header_size(), ios_base::beg);
	
	
	//cout<<key_len<<"\t"<<val_len<<"\t"<<jr.header_size()<<"\t"<<key_ct<<endl;
	// Create a copy of the offsets array for use as insertion positions
	//FILE * input_file = fopen(jr.get_db_name().c_str(),"rb");
	uint64_t bufferSize= ((key_ct  + FILE_NUM_LIMIT - 1) / (FILE_NUM_LIMIT)) << 3;//<<1;
	uint64_t *writeBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	if (!writeBuf) {
		fprintf(stderr, "Require at least %lu bytes, Memory isn't enough\n", bufferSize);
	
	}
	uint64_t bufPoint=0;
	uint64_t bufLimit = bufferSize >> 1;
	hashKmer = writeBuf + bufLimit; 


	evo_tree *_evo_tree = new evo_tree(evo_tree_path);
	se_node *_se_node = new se_node(refListPath, kmerSize, _evo_tree);
	//cout<<pair_size<<endl;
	int file_ind = 0;

	for (uint64_t i = 0; i < key_ct; i++) {
		input_file.read(pair, pair_size);
		//for (int i = 0; i < pair_size; ++i)
			//cout<<std::hex<<pair[i];
		//cout<<endl;
		uint64_t kmerValue = 0;
		uint64_t rcKmerValue = 0;
		uint64_t occBuf = 0;
		uint64_t rcoccBuf = 0;
		memcpy(&kmerValue, pair, key_len);
		memcpy(&occBuf, pair + key_len, val_len);//copy occ buff	
		
		rcKmerValue = rcKmer(kmerValue,kmerSize);
		rcoccBuf = rcOcc(occBuf);	
		//if (kmerValue == 45587144286002457 )
			//cout<<transIntoChars(kmerValue, kmerSize)<<endl;
		if (kmerValue > rcKmerValue) {
			fprintf(stderr,"false format, use jellyfish regenerate kmers");
			exit(1); 
		}
		//if (rcKmerValue == 0)
			//fprintf(stderr,"%lx",occBuf);
		
		++countKmer[kmerValue >> tempMove];
		++countKmer[rcKmerValue >> tempMove];
		
		//fprintf(stderr, "%lu\t%lx\t%lx\t%lx\t%lx\n", i, kmerValue, rcKmerValue,occBuf, rcoccBuf );
		//cout<<kmerValue<<"\t"<<transIntoChars(kmerValue, kmerSize)<<endl;
		//cout<<rcKmerValue<<"\t"<<transIntoChars(rcKmerValue,kmerSize)<<endl;
		_se_node->pro_se_node( kmerValue, occBuf >> 32, rcKmerValue, rcoccBuf >> 32);	
		//cout<<kmerValue<<"\t"<<transIntoChars(kmerValue, kmerSize)<<endl;
		//cout<<rcKmerValue<<"\t"<<transIntoChars(rcKmerValue,kmerSize)<<endl;
				
		writeBuf[bufPoint++] = kmerValue;
		writeBuf[bufPoint++] = occBuf;
		writeBuf[bufPoint++] = rcKmerValue;
		writeBuf[bufPoint++] = rcoccBuf; 
		if(bufPoint>=bufLimit) {
			sort_thread_new(writeBuf, bufLimit, file_ind, kmerSize);	
			//fwrite(writeBuf,sizeof(uint64_t),bufferSize,fpKmer);
			memset(countKmer, 0, BUCKET_CAPACITY * 8);
			++file_ind;
			bufPoint=0;
		}
		//cout<<transIntoChars(kmer,21)<<endl;

	}

	if(bufPoint)  {
			sort_thread_new(writeBuf, bufPoint,file_ind++, kmerSize);	
	}	
		
	useFileNum = file_ind;	
		//fwrite(writeBuf,sizeof(uint64_t),bufPoint,fpKmer);

	fprintf(stderr,"done with cycling\n");
	/*
	while(fscanf(fp1,"%s%lu",seqBuf,&occBuf)!=EOF) {
		uint64_t seq=0;
		for (uint8_t i=0; i<BUCKET_LENGTH;++i)  seq = (seq << 2)| trans[seqBuf[i]];
		countKmer[seq]++;
		for (uint8_t i=BUCKET_LENGTH;i<kmerSize;++i) seq = (seq<<2) | trans[seqBuf[i]];

	}
	*/
	input_file.close();
	
	if (writeBuf)
		free(writeBuf);
	
	//fclose(fpKmer);
	//fclose(fp1);
	//start sort and find start end point
	
	sp_kmers *_sp_kmers = new sp_kmers(kmerSize);

	pthread_t thread_1, thread_2;

		
	pthread_create(&thread_1, NULL, keyway_sort_thread, NULL);	
	//sort_thread((void *)&kmerSize);	
	struct thread_aux_2 ta = {_se_node, _sp_kmers, c_bwt, kmerSize, (key_ct << 1)};
	pthread_create(&thread_2, NULL, call_special_kmers,(void *)&ta);	
	//call_special_kmers(&ta);	
	pthread_join(thread_1, NULL);
	pthread_join(thread_2, NULL);
	
	isStartNode = &_se_node->isNotSingle[0];
	fprintf(stderr, "start merging...");
	mergeKmers("database.srt", kmerSize, _sp_kmers->_kspchar,_sp_kmers->_kspchar_counter, _se_node, c_bwt); 
	fprintf(stderr, "end merging...");
	if (_sp_kmers)
		delete _sp_kmers;
	if (_se_node)
		delete _se_node;
		
	return 0;
}
int processKmers(const char *refListPath, const uint8_t kmerSize,const  char *disorderedKmerpath,const char *evo_tree_path, string &bwt_s, uint64_t *hash_index, vector<uint32_t> &nkmerTID, int differ)
{
    	//////////////////////////////////////kmercounting///////////////////////////////////////////
	countKmer=(uint64_t *)calloc(BUCKET_CAPACITY,sizeof(uint64_t));
	

	uint64_t tempMove=(kmerSize - BUCKET_LENGTH)<<1; //32 defines kmer_plus_one
	fprintf(stderr,"start kmercounting\n");




	//jellyfish db reader from kraken src
	JReader jr(disorderedKmerpath);
	uint64_t val_len = jr.get_val_len();
	uint64_t key_ct = jr.get_key_ct();
	uint64_t key_len = jr.get_key_len();
	uint64_t pair_size = key_len + val_len;
	char pair[pair_size];


	ifstream input_file(jr.get_db_name().c_str(), std::ifstream::binary);
	if (!input_file) return FILE_OPEN_ERROR;
	input_file.seekg(jr.header_size(), ios_base::beg);
	
	
	//cout<<key_len<<"\t"<<val_len<<"\t"<<jr.header_size()<<"\t"<<key_ct<<endl;
	// Create a copy of the offsets array for use as insertion positions
	//FILE * input_file = fopen(jr.get_db_name().c_str(),"rb");
	uint64_t bufferSize= ((key_ct  + FILE_NUM_LIMIT - 1) / (FILE_NUM_LIMIT)) << 3;//<<1;
	uint64_t *writeBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	if (!writeBuf) {
		fprintf(stderr, "Require at least %lu bytes, Memory isn't enough\n", bufferSize);
	
	}
	uint64_t bufPoint=0;
	uint64_t bufLimit = bufferSize >> 1;
	hashKmer = writeBuf + bufLimit; 


	evo_tree *_evo_tree = new evo_tree(evo_tree_path);
	se_node *_se_node = new se_node(refListPath, kmerSize, _evo_tree);
	//cout<<pair_size<<endl;
	int file_ind = 0;

	for (uint64_t i = 0; i < key_ct; i++) {
		input_file.read(pair, pair_size);
		//for (int i = 0; i < pair_size; ++i)
			//cout<<std::hex<<pair[i];
		//cout<<endl;
		uint64_t kmerValue = 0;
		uint64_t rcKmerValue = 0;
		uint64_t occBuf = 0;
		uint64_t rcoccBuf = 0;
		memcpy(&kmerValue, pair, key_len);
		memcpy(&occBuf, pair + key_len, val_len);//copy occ buff	
		
		rcKmerValue = rcKmer(kmerValue,kmerSize);
		rcoccBuf = rcOcc(occBuf);	
		//if (kmerValue == 45587144286002457 )
			//cout<<transIntoChars(kmerValue, kmerSize)<<endl;
		if (kmerValue > rcKmerValue) {
			fprintf(stderr,"false format, use jellyfish regenerate kmers");
			exit(1); 
		}
		//if (rcKmerValue == 0)
			//fprintf(stderr,"%lx",occBuf);
		
		++countKmer[kmerValue >> tempMove];
		++countKmer[rcKmerValue >> tempMove];
		
		//fprintf(stderr, "%lu\t%lx\t%lx\t%lx\t%lx\n", i, kmerValue, rcKmerValue,occBuf, rcoccBuf );
		//cout<<kmerValue<<"\t"<<transIntoChars(kmerValue, kmerSize)<<endl;
		//cout<<rcKmerValue<<"\t"<<transIntoChars(rcKmerValue,kmerSize)<<endl;
		_se_node->pro_se_node( kmerValue, occBuf >> 32, rcKmerValue, rcoccBuf >> 32);	
		//cout<<kmerValue<<"\t"<<transIntoChars(kmerValue, kmerSize)<<endl;
		//cout<<rcKmerValue<<"\t"<<transIntoChars(rcKmerValue,kmerSize)<<endl;
				
		writeBuf[bufPoint++] = kmerValue;
		writeBuf[bufPoint++] = occBuf;
		writeBuf[bufPoint++] = rcKmerValue;
		writeBuf[bufPoint++] = rcoccBuf; 
		if(bufPoint>=bufLimit) {
			sort_thread_new(writeBuf, bufLimit, file_ind, kmerSize);	
			//fwrite(writeBuf,sizeof(uint64_t),bufferSize,fpKmer);
			memset(countKmer, 0, BUCKET_CAPACITY * 8);
			++file_ind;
			bufPoint=0;
		}
		//cout<<transIntoChars(kmer,21)<<endl;

	}

	if(bufPoint)  {
			sort_thread_new(writeBuf, bufPoint,file_ind++, kmerSize);	
	}	
		
	useFileNum = file_ind;	
		//fwrite(writeBuf,sizeof(uint64_t),bufPoint,fpKmer);

	fprintf(stderr,"done with cycling\n");
	/*
	while(fscanf(fp1,"%s%lu",seqBuf,&occBuf)!=EOF) {
		uint64_t seq=0;
		for (uint8_t i=0; i<BUCKET_LENGTH;++i)  seq = (seq << 2)| trans[seqBuf[i]];
		countKmer[seq]++;
		for (uint8_t i=BUCKET_LENGTH;i<kmerSize;++i) seq = (seq<<2) | trans[seqBuf[i]];

	}
	*/
	input_file.close();
	
	if (writeBuf)
		free(writeBuf);
	
	//fclose(fpKmer);
	//fclose(fp1);
	//start sort and find start end point
	
	sp_kmers *_sp_kmers = new sp_kmers(kmerSize);

	pthread_t thread_1, thread_2;

		
	pthread_create(&thread_1, NULL, keyway_sort_thread, NULL);	
	//sort_thread((void *)&kmerSize);	
	struct thread_aux ta = {_se_node, _sp_kmers, bwt_s, nkmerTID};
	pthread_create(&thread_2, NULL, call_special_kmers,(void *)&ta);	
	//call_special_kmers(&ta);	
	pthread_join(thread_1, NULL);
	pthread_join(thread_2, NULL);
	
	isStartNode = &_se_node->isNotSingle[0];
	mergeKmers("database.srt", kmerSize, _sp_kmers->_kspchar,_sp_kmers->_kspchar_counter, _se_node, bwt_s, hash_index, nkmerTID, bwt_s.size()-1); 
	
	if (_sp_kmers)
		delete _sp_kmers;
	if (_se_node)
		delete _se_node;
		
	return 0;
}
uint64_t findInsertPos(uint64_t* p_kmersValue, uint64_t p_s, uint64_t p_e, uint64_t key)
{
	for(uint64_t i=p_s; i < p_e; i+=2) {
		if (key <= p_kmersValue[i])
			return i;
	}
	return p_e;

}

//*********************mergesort with file special kmers and bwt ***************//

int mergeSort(char *fpkmers, const uint8_t _kmer, se_node *_se_,  uint64_t kmersSpcharN, bwt *c_bwt)
{


	fprintf(stderr, "merging using disk space...\n");
	uint64_t tidPtr = c_bwt->tid_ptr;
	uint64_t ind = c_bwt->bwt_str_ptr;
	
	char *bwt_str = c_bwt->bwt_str;
	uint32_t *nkmerTID = c_bwt->taxonIDTab;	

	uint64_t *hash_index = c_bwt->hash_index;
			
	uint64_t kmerNum;
	uint64_t allocatedNum;

	
	FILE *fp = fopen(fpkmers, "rb");
	
	fread(&kmerNum, sizeof(uint64_t), 1 , fp);


	FILE *fp_spchar = fopen("ksp.srt","rb");	
	if (!fp_spchar) {
		fprintf(stderr,"Fail to open sorted special kmers' file, now exit\n");
		exit(1);
	}
	//fprintf(stderr,"%lu\n",kmerNum);	
	//allocatedNum = kmerNum << 1;
	uint64_t *kmersValue = NULL;
	allocatedNum = DG_BUFFERSIZE;
	
	kmersValue = new uint64_t[allocatedNum];
	if (!kmersValue) {
		fprintf(stderr, "Failed to allocate buffer, requires:%lu\n", allocatedNum * sizeof(uint64_t));
		exit(1);	
	}
	uint64_t ksp_allocateNum = allocatedNum >> 1;
	kmersSpchar *_ksp = new kmersSpchar[ksp_allocateNum];
	
	if (!_ksp) {
		fprintf(stderr, "Failed to allocate buffer, requires:%lu\n", ksp_allocateNum * sizeof(kmersSpchar));
		exit(1);	
	}
	//do {
		//kmersValue = new uint64_t[allocatedNum];
		//allocatedNum >>= 1;
	//}while (!kmersValue);
	
	//allocatedNum <<= 1;
		

	uint64_t start;
	uint64_t end;
	start = 0;
	end = start + fread(kmersValue, sizeof(uint64_t), allocatedNum, fp);
	char Chars[] = {'A','C','G','T','#','$'};
	char bitChars[32] = {0};
	bitChars[1] = 'A';
	bitChars[2] = 'C';
	bitChars[4] = 'G';
	bitChars[8] = 'T';	
	
	//map<uint64_t, struct end_tid>::iterator it_map;
	//map<uint64_t, struct end_tid>::const_iterator it_end = se_pair.end();	
	//string bwt = "";
	se_pair *s,t;	
	//fprintf(stderr,"%lu\t%lu\n", end, p_2kmers.size());	
	int mov = (_kmer - PREINDEXLEN) << 1;
	uint64_t pre = -1;
	uint64_t suffix_v = 0;
	bool isEnd = false;
	//vector<kmersSpchar>::iterator it = p_2kmers.begin();
	uint64_t input_num;
	for (uint64_t z=0; z<kmersSpcharN; z += input_num) {
		input_num = fread(_ksp, sizeof(kmersSpchar), ksp_allocateNum, fp_spchar);
		for (uint64_t q = 0 ; q < input_num; ++q) {
			kmersSpchar* it = _ksp + q; 
			uint8_t off = it->infor >> 3;
			//cout<<it->value<<endl;
			uint64_t key = it->value << ((_kmer - off)<< 1); 
			uint64_t loc;	
			do {
				loc = findInsertPos(kmersValue, start,end, key);
				//cout<<"start:"<<start<<"\t"<<"loc:"<<loc<<"\t"<<"end:"<<end<<endl;
				for (uint64_t i = start; i < loc;i += 2){
					//cout<<kmersValue[i]<<endl;
					t.s = kmersValue[i];
					if (isStartNode[kmersValue[i+1]>>37] || (_se_->kb_getp(uint64_t, _se_->node_pair, &t))) 
						bwt_str[ind] = '#';
						//p_bwt += '#';
					else 
						
						bwt_str[ind] = bitChars[kmersValue[i+1] >> 37];
						//p_bwt += bitChars[kmersValue[i+1] >> 37];
					//p_bwt += Chars[p_kmersInfo[i]>>10];
					//++ind;
					//cout<<transIntoChars(kmersValue[i], _kmer)<<"\t"<<ind<<endl;
					//if (!(ind & 0x1F)) nkmerTID.push_back(kmersValue[i+1]&0xFFFFFFFF);
					if (!(ind & 0x1F)) nkmerTID[tidPtr++] = (kmersValue[i+1]&0xFFFFFFFF);
					if ((suffix_v = (kmersValue[i]>>mov)) != pre) {
						//fprintf(stderr,"%lu\t",suffix_v);
						hash_index[suffix_v<<1] = ind;
						hash_index[(suffix_v<<1)+1] = ind + 1;
						pre = suffix_v;
					} else 	
						++hash_index[(suffix_v<<1)+1];
					++ind;
				}
			
				if (loc == end && !isEnd) {
					start = 0;
					end = start + fread(kmersValue, sizeof(uint64_t), allocatedNum, fp);
					if (end == 0) {
						isEnd = true;	
						break;
					}

				} else 
					break;	
			} while (true);
			
			bwt_str[ind] = Chars[it->infor & 0x7]; 
		
			if (!(ind & 0x1F)) nkmerTID[tidPtr++] = (it->taxID);
			
			//cout<<transIntoChars(it->value, off)<<"\t"<<ind<<endl;
			if (off >= PREINDEXLEN) {
			
				if ((suffix_v = (it->value >>((off-PREINDEXLEN)<<1))) != pre) {
					//fprintf(stderr,"%lu\t",suffix_v);
					hash_index[suffix_v<<1] = ind;
					hash_index[(suffix_v<<1)+1] = ind + 1;
					pre = suffix_v; 
				} else 	
					++hash_index[(suffix_v<<1)+1];	
			}
			start = loc;
			++ind;
			//cout<<"stringLength"<<p_bwt.size()<<endl;	
		}
		
	}
	//fprintf(stderr, "hahah%lu\t%lu\n",start, end);
	//cout<<"wwwwww"<<endl;
	do {
			//cout<<"start:"<<start<<"\t"<<"end:"<<end<<endl;
		for (uint64_t i = start; i < end; i += 2)	{
			t.s = kmersValue[i];
			if (isStartNode[kmersValue[i+1]>>37] || (_se_->kb_getp(uint64_t, _se_->node_pair, &t))) 
				bwt_str[ind] = '#';
				//p_bwt += '#';
			else 
				bwt_str[ind] = bitChars[kmersValue[i+1] >> 37];
				//p_bwt += bitChars[kmersValue[i+1] >> 37];
			//p_bwt += Chars[p_kmersInfo[i]>>10];
			//cout<<transIntoChars(kmersValue[i], _kmer)<<"\t"<<ind<<endl;
			if (!(ind & 0x1F)) nkmerTID[tidPtr++] = (kmersValue[i+1]&0XFFFFFFFF);
			
			if ((suffix_v = (kmersValue[i]>>mov)) != pre) {
				//fprintf(stderr,"%lu\t",suffix_v);
				hash_index[suffix_v<<1] = ind;
				hash_index[(suffix_v<<1)+1] = ind + 1;
				pre = suffix_v;
			} else 	
				++hash_index[(suffix_v<<1)+1];
			++ind;
		}
		
		if (!isEnd) {
			start = 0;
			end = start + fread(kmersValue, sizeof(uint64_t), allocatedNum, fp);
			if (start == end) {
				isEnd = true;	
				break;
			}
		} else 
			break;
	
	} while (true);
	if (kmersValue) delete[] kmersValue;
	fclose(fp_spchar);
	fclose(fp);
	fprintf(stderr, "end merging with disk space\n");
	return NORMAL_EXIT;

}

int mergeSort(char *fpkmers, const uint8_t _kmer, se_node *_se_,  kmersSpchar* p_2kmers, uint64_t kmersSpcharN, bwt *c_bwt)
{
	fprintf(stderr, "start merging using memory\n");		
	uint64_t tidPtr = c_bwt->tid_ptr;
	uint64_t ind = c_bwt->bwt_str_ptr;
	
	char *bwt_str = c_bwt->bwt_str;
	uint32_t *nkmerTID = c_bwt->taxonIDTab;	

	uint64_t *hash_index = c_bwt->hash_index;
	
	
	uint64_t kmerNum;
	uint64_t allocatedNum = 2;

	uint64_t start;
	uint64_t end;
	
	FILE *fp = fopen(fpkmers, "rb");
	
	fread(&kmerNum, sizeof(uint64_t), 1 , fp);
	
	//fprintf(stderr,"%lu\n",kmerNum);	
	allocatedNum <<= 30; // 
	uint64_t *kmersValue = NULL;
	fprintf(stderr, "break1\n");		
	while (!kmersValue) {
		kmersValue = new(std::nothrow) uint64_t[allocatedNum];
		allocatedNum >>= 1;
	}
	fprintf(stderr, "break2\n");		
	allocatedNum <<= 1;

	start = 0;
	end = start + fread(kmersValue, sizeof(uint64_t), allocatedNum, fp);
	char Chars[] = {'A','C','G','T','#','$'};
	char bitChars[32] = {0};
	bitChars[1] = 'A';
	bitChars[2] = 'C';
	bitChars[4] = 'G';
	bitChars[8] = 'T';	

	//se_pair *s;
	se_pair t;	
	//map<uint64_t, struct end_tid>::iterator it_map;
	//map<uint64_t, struct end_tid>::const_iterator it_end = se_pair.end();	
	//string bwt = "";
	
	//fprintf(stderr,"%lu\t%lu\n", end, p_2kmers.size());	
	int mov = (_kmer - PREINDEXLEN) << 1;
	uint64_t pre = -1;
	uint64_t suffix_v = 0;
	//uint64_t ind = index_start_point;
	bool isEnd = false;
	//vector<kmersSpchar>::iterator it = p_2kmers.begin();
	fprintf(stderr, "break3\n");		
	for (uint64_t z=0; z<kmersSpcharN; ++z) {
		kmersSpchar* it = p_2kmers + z; 
		uint8_t off = it->infor >> 3;
		//cout<<it->value<<endl;
		uint64_t key = it->value << ((_kmer - off)<< 1); 
		uint64_t loc;	
		do {
			loc = findInsertPos(kmersValue, start,end, key);
			//cout<<"start:"<<start<<"\t"<<"loc:"<<loc<<"\t"<<"end:"<<end<<endl;
			for (uint64_t i = start; i < loc;i += 2){
				//cout<<kmersValue[i]<<endl;
				t.s = kmersValue[i];
				if (isStartNode[kmersValue[i+1]>>37] || (_se_->kb_getp(uint64_t, _se_->node_pair, &t))) 
					bwt_str[ind] = '#';
					//p_bwt += '#';
				else 
					
					bwt_str[ind] = bitChars[kmersValue[i+1] >> 37];
					//p_bwt += bitChars[kmersValue[i+1] >> 37];
				//p_bwt += Chars[p_kmersInfo[i]>>10];
				//++ind;
				//cout<<transIntoChars(kmersValue[i], _kmer)<<"\t"<<ind<<endl;
				//if (!(ind & 0x1F)) nkmerTID.push_back(kmersValue[i+1]&0xFFFFFFFF);
				if (!(ind & 0x1F)) nkmerTID[tidPtr++] = (kmersValue[i+1]&0xFFFFFFFF);
				if ((suffix_v = (kmersValue[i]>>mov)) != pre) {
					//fprintf(stderr,"%lu\t",suffix_v);
					hash_index[suffix_v<<1] = ind;
					hash_index[(suffix_v<<1)+1] = ind + 1;
					pre = suffix_v;
				} else 	
					++hash_index[(suffix_v<<1)+1];
				++ind;
			}
		
			if (loc == end && !isEnd) {
				start = 0;
				end = start + fread(kmersValue, sizeof(uint64_t), allocatedNum, fp);
				if (end == 0) {
					isEnd = true;	
					break;
				}

			} else 
				break;	
		} while (true);
		
		bwt_str[ind] = Chars[it->infor & 0x7]; 
	
		if (!(ind & 0x1F)) nkmerTID[tidPtr++] = (it->taxID);
		
		//cout<<transIntoChars(it->value, off)<<"\t"<<ind<<endl;
		if (off >= PREINDEXLEN) {
		
			if ((suffix_v = (it->value >>((off-PREINDEXLEN)<<1))) != pre) {
				//fprintf(stderr,"%lu\t",suffix_v);
				hash_index[suffix_v<<1] = ind;
				hash_index[(suffix_v<<1)+1] = ind + 1;
				pre = suffix_v; 
			} else 	
				++hash_index[(suffix_v<<1)+1];	
		}
		start = loc;
		++ind;
		//cout<<"stringLength"<<p_bwt.size()<<endl;	
	}
	//fprintf(stderr, "hahah%lu\t%lu\n",start, end);
	//cout<<"wwwwww"<<endl;
	do {
			//cout<<"start:"<<start<<"\t"<<"end:"<<end<<endl;
		for (uint64_t i = start; i < end; i += 2)	{
			t.s = kmersValue[i];
			if (isStartNode[kmersValue[i+1]>>37] || (_se_->kb_getp(uint64_t, _se_->node_pair, &t))) 
				bwt_str[ind] = '#';
				//p_bwt += '#';
			else 
				bwt_str[ind] = bitChars[kmersValue[i+1] >> 37];
				//p_bwt += bitChars[kmersValue[i+1] >> 37];
			//p_bwt += Chars[p_kmersInfo[i]>>10];
			//cout<<transIntoChars(kmersValue[i], _kmer)<<"\t"<<ind<<endl;
			if (!(ind & 0x1F)) nkmerTID[tidPtr++] = (kmersValue[i+1]&0XFFFFFFFF);
			
			if ((suffix_v = (kmersValue[i]>>mov)) != pre) {
				//fprintf(stderr,"%lu\t",suffix_v);
				hash_index[suffix_v<<1] = ind;
				hash_index[(suffix_v<<1)+1] = ind + 1;
				pre = suffix_v;
			} else 	
				++hash_index[(suffix_v<<1)+1];
			++ind;
		}
		
		if (!isEnd) {
			start = 0;
			end = start + fread(kmersValue, sizeof(uint64_t), allocatedNum, fp);
			if (start == end) {
				isEnd = true;	
				break;
			}
		} else 
			break;
	
	} while (true);
	
	fprintf(stderr, "break4\n");		
	if (kmersValue) delete   []kmersValue;
	fprintf(stderr, "end merging using memory\n");		
	return NORMAL_EXIT;




}

//**********************without using bwt ********************/

int mergeSort(char *fpkmers, const uint8_t _kmer, se_node *_se_,  uint64_t kmersSpcharN, string& p_bwt, vector<uint32_t>& nkmerTID, uint64_t *hash_index, uint64_t index_start_point)
{
	uint64_t kmerNum;
	uint64_t allocatedNum;

	
	FILE *fp = fopen(fpkmers, "rb");
	
	fread(&kmerNum, sizeof(uint64_t), 1 , fp);


	FILE *fp_spchar = fopen("ksp.srt","rb");	
	if (!fp_spchar) {
		fprintf(stderr,"Fail to open sorted special kmers' file, now exit\n");
		exit(1);
	}
	//fprintf(stderr,"%lu\n",kmerNum);	
	//allocatedNum = kmerNum << 1;
	uint64_t *kmersValue = NULL;
	allocatedNum = DG_BUFFERSIZE;
	
	kmersValue = new uint64_t[allocatedNum];
	if (!kmersValue) {
		fprintf(stderr, "Failed to allocate buffer, requires:%lu\n", allocatedNum * sizeof(uint64_t));
		exit(1);	
	}
	uint64_t ksp_allocateNum = allocatedNum >> 1;
	kmersSpchar *_ksp = new kmersSpchar[ksp_allocateNum];
	
	if (!_ksp) {
		fprintf(stderr, "Failed to allocate buffer, requires:%lu\n", ksp_allocateNum * sizeof(kmersSpchar));
		exit(1);	
	}
	//do {
		//kmersValue = new uint64_t[allocatedNum];
		//allocatedNum >>= 1;
	//}while (!kmersValue);
	
	//allocatedNum <<= 1;
		

	uint64_t start;
	uint64_t end;
	start = 0;
	end = start + fread(kmersValue, sizeof(uint64_t), allocatedNum, fp);
	char Chars[] = {'A','C','G','T','#','$'};
	char bitChars[32] = {0};
	bitChars[1] = 'A';
	bitChars[2] = 'C';
	bitChars[4] = 'G';
	bitChars[8] = 'T';	
	
	//map<uint64_t, struct end_tid>::iterator it_map;
	//map<uint64_t, struct end_tid>::const_iterator it_end = se_pair.end();	
	//string bwt = "";
	se_pair *s,t;	
	//fprintf(stderr,"%lu\t%lu\n", end, p_2kmers.size());	
	int mov = (_kmer - PREINDEXLEN) << 1;
	uint64_t pre = -1;
	uint64_t suffix_v = 0;
	uint64_t ind = index_start_point;
	bool isEnd = false;
	//vector<kmersSpchar>::iterator it = p_2kmers.begin();
	uint64_t input_num;
	for (uint64_t z=0; z<kmersSpcharN; z += input_num) {
		input_num = fread(_ksp, sizeof(kmersSpchar), ksp_allocateNum, fp_spchar);
		for (uint64_t q = 0 ; q < input_num; ++q) {
			kmersSpchar* it = _ksp + q; 
			uint8_t off = it->infor >> 3;
			//cout<<it->value<<endl;
			uint64_t key = it->value << ((_kmer - off)<< 1); 
			uint64_t loc;	
			do {
				loc = findInsertPos(kmersValue, start,end, key);
				//cout<<"start:"<<start<<"\t"<<"loc:"<<loc<<"\t"<<"end:"<<end<<endl;
				for (uint64_t i = start; i < loc;i += 2){
					//cout<<kmersValue[i]<<endl;
					t.s = kmersValue[i];
					if (isStartNode[kmersValue[i+1]>>37] || (_se_->kb_getp(uint64_t, _se_->node_pair, &t))) 
						p_bwt += '#';
					else 
						p_bwt += bitChars[kmersValue[i+1] >> 37];
					//p_bwt += Chars[p_kmersInfo[i]>>10];
					++ind;
					//cout<<transIntoChars(kmersValue[i], _kmer)<<"\t"<<ind<<endl;
					if (!(ind & 0x1F)) nkmerTID.push_back(kmersValue[i+1]&0xFFFFFFFF);
					if ((suffix_v = (kmersValue[i]>>mov)) != pre) {
						//fprintf(stderr,"%lu\t",suffix_v);
						hash_index[suffix_v<<1] = ind;
						hash_index[(suffix_v<<1)+1] = ind + 1;
						pre = suffix_v;
					} else 	
						++hash_index[(suffix_v<<1)+1];
				}
			
				if (loc == end && !isEnd) {
					start = 0;
					end = start + fread(kmersValue, sizeof(uint64_t), allocatedNum, fp);
					if (end == 0) {
						isEnd = true;	
						break;
					}

				} else 
					break;	
			} while (true);
			
			p_bwt += Chars[it->infor & 0x7]; 
			++ind;
		
			if (!(ind & 0x1F)) nkmerTID.push_back(it->taxID);
			
			//cout<<transIntoChars(it->value, off)<<"\t"<<ind<<endl;
			if (off >= PREINDEXLEN) {
			
				if ((suffix_v = (it->value >>((off-PREINDEXLEN)<<1))) != pre) {
					//fprintf(stderr,"%lu\t",suffix_v);
					hash_index[suffix_v<<1] = ind;
					hash_index[(suffix_v<<1)+1] = ind + 1;
					pre = suffix_v; 
				} else 	
					++hash_index[(suffix_v<<1)+1];	
			}
			start = loc;
			//cout<<"stringLength"<<p_bwt.size()<<endl;	
		}
		
	}
	//fprintf(stderr, "hahah%lu\t%lu\n",start, end);
	//cout<<"wwwwww"<<endl;
	do {
			//cout<<"start:"<<start<<"\t"<<"end:"<<end<<endl;
		for (uint64_t i = start; i < end; i += 2)	{
			t.s = kmersValue[i];
			if (isStartNode[kmersValue[i+1]>>37] || (_se_->kb_getp(uint64_t, _se_->node_pair, &t))) 
				p_bwt += '#';
			else 
				p_bwt += bitChars[kmersValue[i+1] >> 37];
			//p_bwt += Chars[p_kmersInfo[i]>>10];
			++ind;
			//cout<<transIntoChars(kmersValue[i], _kmer)<<"\t"<<ind<<endl;
			if (!(ind & 0x1F)) nkmerTID.push_back(kmersValue[i+1]&0XFFFFFFFF);
			
			if ((suffix_v = (kmersValue[i]>>mov)) != pre) {
				//fprintf(stderr,"%lu\t",suffix_v);
				hash_index[suffix_v<<1] = ind;
				hash_index[(suffix_v<<1)+1] = ind + 1;
				pre = suffix_v;
			} else 	
				++hash_index[(suffix_v<<1)+1];
		}
		
		if (!isEnd) {
			start = 0;
			end = start + fread(kmersValue, sizeof(uint64_t), allocatedNum, fp);
			if (start == end) {
				isEnd = true;	
				break;
			}
		} else 
			break;
	
	} while (true);
	
	if (kmersValue) delete[] kmersValue;
	fclose(fp_spchar);
	fclose(fp);
	return NORMAL_EXIT;

}
/* 
int mergeSort(char *fpkmers, const uint8_t _kmer, se_node *_se_,  uint64_t kmersSpcharN, string& p_bwt, vector<uint32_t>& nkmerTID, uint64_t *hash_index, uint64_t index_start_point)
{
	uint64_t kmerNum;
	uint64_t allocatedNum;

	
	FILE *fp = fopen(fpkmers, "rb");
	
	fread(&kmerNum, sizeof(uint64_t), 1 , fp);


	FILE *fp_spchar = fopen("ksp.srt","rb");	
	if (!fp_spchar) {
		fprintf(stderr,"Fail to open sorted special kmers' file, now exit\n");
		exit(1);
	}
	//fprintf(stderr,"%lu\n",kmerNum);	
	//allocatedNum = kmerNum << 1;
	uint64_t *kmersValue = NULL;
	allocatedNum = DG_BUFFERSIZE;
	
	kmersValue = new uint64_t[allocatedNum];
	if (!kmersValue) {
		fprintf(stderr, "Failed to allocate buffer, requires:%lu\n", allocatedNum * sizeof(uint64_t));
		exit(1);	
	}
	uint64_t ksp_allocateNum = allocatedNum >> 1;
	kmersSpchar *_ksp = new kmersSpchar[ksp_allocateNum];
	
	if (!_ksp) {
		fprintf(stderr, "Failed to allocate buffer, requires:%lu\n", ksp_allocateNum * sizeof(kmersSpchar));
		exit(1);	
	}
	//do {
		//kmersValue = new uint64_t[allocatedNum];
		//allocatedNum >>= 1;
	//}while (!kmersValue);
	
	//allocatedNum <<= 1;
		

	uint64_t start;
	uint64_t end;
	start = 0;
	end = start + fread(kmersValue, sizeof(uint64_t), allocatedNum, fp);
	char Chars[] = {'A','C','G','T','#','$'};
	char bitChars[32] = {0};
	bitChars[1] = 'A';
	bitChars[2] = 'C';
	bitChars[4] = 'G';
	bitChars[8] = 'T';	
	
	//map<uint64_t, struct end_tid>::iterator it_map;
	//map<uint64_t, struct end_tid>::const_iterator it_end = se_pair.end();	
	//string bwt = "";
	se_pair *s,t;	
	//fprintf(stderr,"%lu\t%lu\n", end, p_2kmers.size());	
	int mov = (_kmer - PREINDEXLEN) << 1;
	uint64_t pre = -1;
	uint64_t suffix_v = 0;
	uint64_t ind = index_start_point;
	bool isEnd = false;
	//vector<kmersSpchar>::iterator it = p_2kmers.begin();
	uint64_t input_num;
	for (uint64_t z=0; z<kmersSpcharN; z += input_num) {
		input_num = fread(_ksp, sizeof(kmersSpchar), ksp_allocateNum, fp_spchar);
		for (uint64_t q = 0 ; q < input_num; ++q) {
			kmersSpchar* it = _ksp + q; 
			uint8_t off = it->infor >> 3;
			//cout<<it->value<<endl;
			uint64_t key = it->value << ((_kmer - off)<< 1); 
			uint64_t loc;	
			do {
				loc = findInsertPos(kmersValue, start,end, key);
				//cout<<"start:"<<start<<"\t"<<"loc:"<<loc<<"\t"<<"end:"<<end<<endl;
				for (uint64_t i = start; i < loc;i += 2){
					//cout<<kmersValue[i]<<endl;
					t.s = kmersValue[i];
					if (isStartNode[kmersValue[i+1]>>37] || (_se_->kb_getp(uint64_t, _se_->node_pair, &t))) 
						p_bwt += '#';
					else 
						p_bwt += bitChars[kmersValue[i+1] >> 37];
					//p_bwt += Chars[p_kmersInfo[i]>>10];
					++ind;
					cout<<transIntoChars(kmersValue[i], _kmer)<<"\t"<<ind<<endl;
					if (!(ind & 0x1F)) nkmerTID.push_back(kmersValue[i+1]&0xFFFFFFFF);
					if ((suffix_v = (kmersValue[i]>>mov)) != pre) {
						//fprintf(stderr,"%lu\t",suffix_v);
						hash_index[suffix_v<<1] = ind;
						hash_index[(suffix_v<<1)+1] = ind + 1;
						pre = suffix_v;
					} else 	
						++hash_index[(suffix_v<<1)+1];
				}
			
				if (loc == end && !isEnd) {
					start = 0;
					end = start + fread(kmersValue, sizeof(uint64_t), allocatedNum, fp);
					if (end == 0) {
						isEnd = true;	
						break;
					}

				} else 
					break;	
			} while (true);
			
			p_bwt += Chars[it->infor & 0x7]; 
			++ind;
		
			if (!(ind & 0x1F)) nkmerTID.push_back(it->taxID);
			
			cout<<transIntoChars(it->value, off)<<"\t"<<ind<<endl;
			if (off >= PREINDEXLEN) {
			
				if ((suffix_v = (it->value >>((off-PREINDEXLEN)<<1))) != pre) {
					//fprintf(stderr,"%lu\t",suffix_v);
					hash_index[suffix_v<<1] = ind;
					hash_index[(suffix_v<<1)+1] = ind + 1;
					pre = suffix_v; 
				} else 	
					++hash_index[(suffix_v<<1)+1];	
			}
			start = loc;
			//cout<<"stringLength"<<p_bwt.size()<<endl;	
		}
		
	}
	//fprintf(stderr, "hahah%lu\t%lu\n",start, end);
	//cout<<"wwwwww"<<endl;
	do {
			//cout<<"start:"<<start<<"\t"<<"end:"<<end<<endl;
		for (uint64_t i = start; i < end; i += 2)	{
			t.s = kmersValue[i];
			if (isStartNode[kmersValue[i+1]>>37] || (s = kb_getp(uint64_t, bt, &t))) 
				p_bwt += '#';
			else 
				p_bwt += bitChars[kmersValue[i+1] >> 37];
			//p_bwt += Chars[p_kmersInfo[i]>>10];
			++ind;
			cout<<transIntoChars(kmersValue[i], _kmer)<<"\t"<<ind<<endl;
			if (!(ind & 0x1F)) nkmerTID.push_back(kmersValue[i+1]&0XFFFFFFFF);
			
			if ((suffix_v = (kmersValue[i]>>mov)) != pre) {
				//fprintf(stderr,"%lu\t",suffix_v);
				hash_index[suffix_v<<1] = ind;
				hash_index[(suffix_v<<1)+1] = ind + 1;
				pre = suffix_v;
			} else 	
				++hash_index[(suffix_v<<1)+1];
		}
		
		if (!isEnd) {
			start = 0;
			end = start + fread(kmersValue, sizeof(uint64_t), allocatedNum, fp);
			if (start == end) {
				isEnd = true;	
				break;
			}
		} else 
			break;
	
	} while (true);
	
	if (kmersValue) delete[] kmersValue;
	fclose(fp_spchar);
	fclose(fp);
	return NORMAL_EXIT;

}
*/
/*  
int mergeSort(char *fpkmers, const uint8_t _kmer, kbtree_t(uint64_t) *bt,  uint64_t kmersSpcharN, string& p_bwt, vector<uint32_t>& nkmerTID, uint64_t *hash_index, uint64_t index_start_point)
{
	uint64_t kmerNum;
	uint64_t allocatedNum;

	
	FILE *fp = fopen(fpkmers, "rb");
	
	fread(&kmerNum, sizeof(uint64_t), 1 , fp);


	FILE *fp_spchar = fopen("ksp.srt","rb");	
	if (!fp_spchar) {
		fprintf(stderr,"Fail to open sorted special kmers' file, now exit\n");
		exit(1);
	}
	//fprintf(stderr,"%lu\n",kmerNum);	
	//allocatedNum = kmerNum << 1;
	uint64_t *kmersValue = NULL;
	allocatedNum = DG_BUFFERSIZE;
	
	kmersValue = new uint64_t[allocatedNum];
	if (!kmersValue) {
		fprintf(stderr, "Failed to allocate buffer, requires:%lu\n", allocatedNum * sizeof(uint64_t));
		exit(1);	
	}
	uint64_t ksp_allocateNum = allocatedNum >> 1;
	kmersSpchar *_ksp = new kmersSpchar[ksp_allocateNum];
	
	if (!_ksp) {
		fprintf(stderr, "Failed to allocate buffer, requires:%lu\n", ksp_allocateNum * sizeof(kmersSpchar));
		exit(1);	
	}
	//do {
		//kmersValue = new uint64_t[allocatedNum];
		//allocatedNum >>= 1;
	//}while (!kmersValue);
	
	//allocatedNum <<= 1;
		

	uint64_t start;
	uint64_t end;
	start = 0;
	end = start + fread(kmersValue, sizeof(uint64_t), allocatedNum, fp);
	char Chars[] = {'A','C','G','T','#','$'};
	char bitChars[32] = {0};
	bitChars[1] = 'A';
	bitChars[2] = 'C';
	bitChars[4] = 'G';
	bitChars[8] = 'T';	
	
	//map<uint64_t, struct end_tid>::iterator it_map;
	//map<uint64_t, struct end_tid>::const_iterator it_end = se_pair.end();	
	//string bwt = "";
	se_pair *s,t;	
	//fprintf(stderr,"%lu\t%lu\n", end, p_2kmers.size());	
	int mov = (_kmer - PREINDEXLEN) << 1;
	uint64_t pre = -1;
	uint64_t suffix_v = 0;
	uint64_t ind = index_start_point;
	bool isEnd = false;
	//vector<kmersSpchar>::iterator it = p_2kmers.begin();
	uint64_t input_num;
	for (uint64_t z=0; z<kmersSpcharN; z += input_num) {
		input_num = fread(_ksp, sizeof(kmersSpchar), ksp_allocateNum, fp_spchar);
		for (uint64_t q = 0 ; q < input_num; ++q) {
			kmersSpchar* it = _ksp + q; 
			uint8_t off = it->infor >> 3;
			//cout<<it->value<<endl;
			uint64_t key = it->value << ((_kmer - off)<< 1); 
			uint64_t loc;	
			do {
				loc = findInsertPos(kmersValue, start,end, key);
				//cout<<"start:"<<start<<"\t"<<"loc:"<<loc<<"\t"<<"end:"<<end<<endl;
				for (uint64_t i = start; i < loc;i += 2){
					//cout<<kmersValue[i]<<endl;
					t.s = kmersValue[i];
					if (isStartNode[kmersValue[i+1]>>37] || (s = kb_getp(uint64_t, bt, &t))) 
						p_bwt += '#';
					else 
						p_bwt += bitChars[kmersValue[i+1] >> 37];
					//p_bwt += Chars[p_kmersInfo[i]>>10];
					++ind;
					cout<<transIntoChars(kmersValue[i], _kmer)<<"\t"<<ind<<endl;
					if (!(ind & 0x1F)) nkmerTID.push_back(kmersValue[i+1]&0xFFFFFFFF);
					if ((suffix_v = (kmersValue[i]>>mov)) != pre) {
						//fprintf(stderr,"%lu\t",suffix_v);
						hash_index[suffix_v<<1] = ind;
						hash_index[(suffix_v<<1)+1] = ind + 1;
						pre = suffix_v;
					} else 	
						++hash_index[(suffix_v<<1)+1];
				}
			
				if (loc == end && !isEnd) {
					start = 0;
					end = start + fread(kmersValue, sizeof(uint64_t), allocatedNum, fp);
					if (end == 0) {
						isEnd = true;	
						break;
					}

				} else 
					break;	
			} while (true);
			
			p_bwt += Chars[it->infor & 0x7]; 
			++ind;
		
			if (!(ind & 0x1F)) nkmerTID.push_back(it->taxID);
			
			cout<<transIntoChars(it->value, off)<<"\t"<<ind<<endl;
			if (off >= PREINDEXLEN) {
			
				if ((suffix_v = (it->value >>((off-PREINDEXLEN)<<1))) != pre) {
					//fprintf(stderr,"%lu\t",suffix_v);
					hash_index[suffix_v<<1] = ind;
					hash_index[(suffix_v<<1)+1] = ind + 1;
					pre = suffix_v; 
				} else 	
					++hash_index[(suffix_v<<1)+1];	
			}
			start = loc;
			//cout<<"stringLength"<<p_bwt.size()<<endl;	
		}
		
	}
	//fprintf(stderr, "hahah%lu\t%lu\n",start, end);
	//cout<<"wwwwww"<<endl;
	do {
			//cout<<"start:"<<start<<"\t"<<"end:"<<end<<endl;
		for (uint64_t i = start; i < end; i += 2)	{
			t.s = kmersValue[i];
			if (isStartNode[kmersValue[i+1]>>37] || (s = kb_getp(uint64_t, bt, &t))) 
				p_bwt += '#';
			else 
				p_bwt += bitChars[kmersValue[i+1] >> 37];
			//p_bwt += Chars[p_kmersInfo[i]>>10];
			++ind;
			cout<<transIntoChars(kmersValue[i], _kmer)<<"\t"<<ind<<endl;
			if (!(ind & 0x1F)) nkmerTID.push_back(kmersValue[i+1]&0XFFFFFFFF);
			
			if ((suffix_v = (kmersValue[i]>>mov)) != pre) {
				//fprintf(stderr,"%lu\t",suffix_v);
				hash_index[suffix_v<<1] = ind;
				hash_index[(suffix_v<<1)+1] = ind + 1;
				pre = suffix_v;
			} else 	
				++hash_index[(suffix_v<<1)+1];
		}
		
		if (!isEnd) {
			start = 0;
			end = start + fread(kmersValue, sizeof(uint64_t), allocatedNum, fp);
			if (start == end) {
				isEnd = true;	
				break;
			}
		} else 
			break;
	
	} while (true);
	
	if (kmersValue) delete[] kmersValue;
	fclose(fp_spchar);
	fclose(fp);
	return NORMAL_EXIT;

}
*/
int mergeSort(char *fpkmers, const uint8_t _kmer, se_node *_se_,  kmersSpchar* p_2kmers, uint64_t kmersSpcharN, string& p_bwt, vector<uint32_t>& nkmerTID, uint64_t *hash_index, uint64_t index_start_point)
{
	
	uint64_t kmerNum;
	uint64_t allocatedNum;

	uint64_t start;
	uint64_t end;
	
	FILE *fp = fopen(fpkmers, "rb");
	
	fread(&kmerNum, sizeof(uint64_t), 1 , fp);
	
	//fprintf(stderr,"%lu\n",kmerNum);	
	allocatedNum = kmerNum << 1;
	uint64_t *kmersValue = NULL;
	do {
		kmersValue = new uint64_t[allocatedNum];
		allocatedNum >>= 1;
	}while (!kmersValue);
	
	allocatedNum <<= 1;

	start = 0;
	end = start + fread(kmersValue, sizeof(uint64_t), allocatedNum, fp);
	char Chars[] = {'A','C','G','T','#','$'};
	char bitChars[32] = {0};
	bitChars[1] = 'A';
	bitChars[2] = 'C';
	bitChars[4] = 'G';
	bitChars[8] = 'T';	

	//se_pair *s;
	se_pair t;	
	//map<uint64_t, struct end_tid>::iterator it_map;
	//map<uint64_t, struct end_tid>::const_iterator it_end = se_pair.end();	
	//string bwt = "";
	
	//fprintf(stderr,"%lu\t%lu\n", end, p_2kmers.size());	
	int mov = (_kmer - PREINDEXLEN) << 1;
	uint64_t pre = -1;
	uint64_t suffix_v = 0;
	uint64_t ind = index_start_point;
	bool isEnd = false;
	//vector<kmersSpchar>::iterator it = p_2kmers.begin();
	for (uint64_t z=0; z<kmersSpcharN; ++z) {
		kmersSpchar* it = p_2kmers + z; 
		uint8_t off = it->infor >> 3;
		//cout<<it->value<<endl;
		uint64_t key = it->value << ((_kmer - off)<< 1); 
		uint64_t loc;	
		do {
			loc = findInsertPos(kmersValue, start,end, key);
			//cout<<"start:"<<start<<"\t"<<"loc:"<<loc<<"\t"<<"end:"<<end<<endl;
			for (uint64_t i = start; i < loc;i += 2){
				//cout<<kmersValue[i]<<endl;
				t.s = kmersValue[i];
				if (isStartNode[kmersValue[i+1]>>37] ||(_se_->kb_getp(uint64_t, _se_->node_pair, &t))) 
					p_bwt += '#';
				else 
					p_bwt += bitChars[kmersValue[i+1] >> 37];
				//p_bwt += Chars[p_kmersInfo[i]>>10]
				++ind;
				//cout<<transIntoChars(kmersValue[i], _kmer)<<"\t"<<ind<<endl;
				if (!(ind & 0x1F)) nkmerTID.push_back(kmersValue[i+1]&0xFFFFFFFF);
				if ((suffix_v = (kmersValue[i]>>mov)) != pre) {
					//fprintf(stderr,"%lu\t",suffix_v);
					hash_index[suffix_v<<1] = ind;
					hash_index[(suffix_v<<1)+1] = ind + 1;
					pre = suffix_v;
				} else 	
					++hash_index[(suffix_v<<1)+1];
			}
		
			if (loc == end && !isEnd) {
				start = 0;
				end = start + fread(kmersValue, sizeof(uint64_t), allocatedNum, fp);
				if (end == 0) {
					isEnd = true;	
					break;
				}

			} else 
				break;	
		} while (true);
		
		p_bwt += Chars[it->infor & 0x7]; 
		++ind;
	
		if (!(ind & 0x1F)) nkmerTID.push_back(it->taxID);
		
		//cout<<transIntoChars(it->value, off)<<"\t"<<ind<<endl;
		if (off >= PREINDEXLEN) {
		
			if ((suffix_v = (it->value >>((off-PREINDEXLEN)<<1))) != pre) {
				//fprintf(stderr,"%lu\t",suffix_v);
				hash_index[suffix_v<<1] = ind;
				hash_index[(suffix_v<<1)+1] = ind + 1;
				pre = suffix_v; 
			} else 	
				++hash_index[(suffix_v<<1)+1];	
		}
		start = loc;
		//cout<<"stringLength"<<p_bwt.size()<<endl;	
	}
	//fprintf(stderr, "hahah%lu\t%lu\n",start, end);
	//cout<<"wwwwww"<<endl;
	do {
			//cout<<"start:"<<start<<"\t"<<"end:"<<end<<endl;
		for (uint64_t i = start; i < end; i += 2)	{
			t.s = kmersValue[i];
			if (isStartNode[kmersValue[i+1]>>37] || (_se_->kb_getp(uint64_t, _se_->node_pair, &t))) 
				p_bwt += '#';
			else 
				p_bwt += bitChars[kmersValue[i+1] >> 37];
			//p_bwt += Chars[p_kmersInfo[i]>>10];
			++ind;
			//cout<<transIntoChars(kmersValue[i], _kmer)<<"\t"<<ind<<endl;
			if (!(ind & 0x1F)) nkmerTID.push_back(kmersValue[i+1]&0XFFFFFFFF);
			
			if ((suffix_v = (kmersValue[i]>>mov)) != pre) {
				//fprintf(stderr,"%lu\t",suffix_v);
				hash_index[suffix_v<<1] = ind;
				hash_index[(suffix_v<<1)+1] = ind + 1;
				pre = suffix_v;
			} else 	
				++hash_index[(suffix_v<<1)+1];
		}
		
		if (!isEnd) {
			start = 0;
			end = start + fread(kmersValue, sizeof(uint64_t), allocatedNum, fp);
			if (start == end) {
				isEnd = true;	
				break;
			}
		} else 
			break;
	
	} while (true);
	
	if (kmersValue) delete[] kmersValue;
	return NORMAL_EXIT;




}

int mergeKmers(char *sortedKmersPath, const uint8_t _kmer_len, kmersSpchar *_kspchar, uint64_t kmersSpcharN, se_node *_se_, bwt *c_bwt)
{
	if (_kspchar) 	
       		mergeSort(sortedKmersPath,_kmer_len , _se_,  _kspchar, kmersSpcharN, c_bwt);
	else 
			mergeSort(sortedKmersPath, _kmer_len, _se_, kmersSpcharN, c_bwt);
	return NORMAL_EXIT;
}

int mergeKmers(char *sortedKmersPath, const uint8_t _kmer_len, kmersSpchar *_kspchar, uint64_t kmersSpcharN, se_node *_se_, string &bwt_s, uint64_t *hash_index, vector<uint32_t> &nkmerTID, uint64_t index_start_point)
{
	if (_kspchar) 	
       		mergeSort(sortedKmersPath,_kmer_len , _se_,  _kspchar, kmersSpcharN,  bwt_s, nkmerTID, hash_index, index_start_point);
	else 
		mergeSort(sortedKmersPath, _kmer_len, _se_, kmersSpcharN, bwt_s, nkmerTID, hash_index, index_start_point);
	return NORMAL_EXIT;
}
/*  
int mergeKmers(char *sortedKmersPath, const uint8_t _kmer_len, kmersSpchar *_kspchar, uint64_t kmersSpcharN, kbtree_t(uint64_t) *bt, string &bwt_s, uint64_t *hash_index, vector<uint32_t> nkmerTID, uint64_t index_start_point)
{
		
	if (_kspchar) 	
       		mergeSort(sortedKmersPath,_kmer_len , bt,  _kspchar, kmersSpcharN,  bwt_s, nkmerTID, hash_index, index_start_point);
	else 
		mergeSort(sortedKmersPath, _kmer_len,  bt, kmersSpcharN, bwt_s, nkmerTID, hash_index, index_start_point);
	
	
		
	
	//cout<<bwt_s<<endl;
	
	FILE *bwt_string_fp = fopen("bwt_string","w");
	
	fwrite(bwt_s.c_str(), sizeof(char), bwt_s.size(), bwt_string_fp);	
	
	fclose(bwt_string_fp);

	return NORMAL_EXIT;




}
*/
/*
int mergeKmers(char *sortedKmersPath, const uint8_t _kmer_len, kmersSpchar *_kspchar, uint64_t kmersSpcharN, map<uint64_t, struct end_tid> &pairs, string &bwt_s, uint64_t *hash_index, vector<uint32_t> nkmerTID, uint64_t index_start_point)
{
		
	if (_kspchar) 	
       		mergeSort(sortedKmersPath,_kmer_len , pairs,  _kspchar, kmersSpcharN,  bwt_s, nkmerTID, hash_index, index_start_point);
	else 
		mergeSort(sortedKmersPath, _kmer_len,  pairs, kmersSpcharN, bwt_s, nkmerTID, hash_index, index_start_point);
	
	
		
	
	//cout<<bwt_s<<endl;
	
	FILE *bwt_string_fp = fopen("bwt_string","w");
	
	fwrite(bwt_s.c_str(), sizeof(char), bwt_s.size(), bwt_string_fp);	
	
	fclose(bwt_string_fp);

	return NORMAL_EXIT;




}
*/

void* call_special_kmers(void *arg) 
{
	se_node *_se_ = ((struct thread_aux_2 *)arg)->_s_e;
	sp_kmers *_sp_kmers = ((struct thread_aux_2 *)arg)->_s_p;
	bwt * c_bwt = ((struct thread_aux_2 *)arg)->c_bwt;

	uint64_t kmersize = ((struct thread_aux_2 *)arg)->kmersize;
	uint64_t key_ct = ((struct thread_aux_2 *)arg)->kmerN;
	
	fprintf(stderr,"start run se_node\n");
	_se_->run(NULL);
	fprintf(stderr,"start initiate specail kmers' space\n");
	uint64_t node_size = kb_size(_se_->node_pair);
	_sp_kmers->initiate_space(node_size);//be careful of space
	//map<uint64_t, struct end_tid> &T = _se_->s_e_pair; 
	//map<uint64_t, struct end_tid>::iterator it = T.begin();
	
	//map<uint64_t, struct end_tid>::const_iterator it_end = T.end();
	//start to init bwt
	c_bwt->len_bwt_str = node_size * kmersize + key_ct;
	
	
	c_bwt->tidSize = node_size - ((node_size - 1) >> 5) + ((c_bwt->len_bwt_str - 1) >> 5);	
	

	if (!(c_bwt->bwt_str = (char *)malloc(sizeof(char)*c_bwt->len_bwt_str))) {
		fprintf(stderr, "Fail to allocate space for bwt string\n");	
		exit(1);	
	}
	if (!(c_bwt->taxonIDTab = (uint32_t *)malloc(sizeof(uint32_t)*c_bwt->tidSize))) {
		fprintf(stderr, "Fail to allocate space for bwt string\n");	
		exit(1);	
	}
	
	c_bwt->tid_ptr = 0;
	c_bwt->bwt_str_ptr = 0;

	
	
	fprintf(stderr,"start generate specail kmers space\n");
	
	kbitr_t itr;
	se_pair *t;
	se_node::kbtree_t(uint64_t) *h = _se_->node_pair;
	//_se_->kb_itr_first(uint64_t, h, &itr); 
	se_node::kb_itr_first(uint64_t, h, &itr);
	while (kb_itr_valid(&itr)) {
		//cout<<it->first<<"\t"<<it->second.end_node<<endl;
		t = &kb_itr_key(se_pair, &itr);
		
		_sp_kmers->gen_sp_kmers(t->e, t->tid, c_bwt);
		se_node::kb_itr_next(uint64_t, h, &itr);	
	}
	
	fprintf(stderr,"flush out left kmers\n");
	_sp_kmers->flush_sp_kmers();	
	fprintf(stderr,"close files and remove useless variables\n");
	_sp_kmers->closeFiles();
	fprintf(stderr,"sort kmers\n");
	_sp_kmers->sortKmers();

	fprintf(stderr,"end sort kmers\n");
	return NULL;
}
/*  
 *
void* call_special_kmers(void *arg) 
{
	se_node *_se_ = ((struct thread_aux *)arg)->_s_e;
	sp_kmers *_sp_kmers = ((struct thread_aux *)arg)->_s_p;
	string &bwt_s = ((struct thread_aux *)arg)->bwt_s;
	vector<uint32_t> &tids = ((struct thread_aux *)arg)->tids;
	fprintf(stderr,"start run se_node\n");
	_se_->run(NULL);
	fprintf(stderr,"start initiate specail kmers' space\n");
	_sp_kmers->initiate_space(kb_size(_se_->node_pair));//be careful of space
	//map<uint64_t, struct end_tid> &T = _se_->s_e_pair; 
	//map<uint64_t, struct end_tid>::iterator it = T.begin();
	
	//map<uint64_t, struct end_tid>::const_iterator it_end = T.end();
	fprintf(stderr,"start generate specail kmers space\n");
	
	kbitr_t itr;
	se_pair *t;
	se_node::kbtree_t(uint64_t) *h = _se_->node_pair;
	//_se_->kb_itr_first(uint64_t, h, &itr); 
	se_node::kb_itr_first(uint64_t, h, &itr);
	while (kb_itr_valid(&itr)) {
		//cout<<it->first<<"\t"<<it->second.end_node<<endl;
		t = &kb_itr_key(se_pair, &itr);
		
		_sp_kmers->gen_sp_kmers(t->e, t->tid, bwt_s, tids);
		se_node::kb_itr_next(uint64_t, h, &itr);	
	}
	
	fprintf(stderr,"flush out left kmers\n");
	_sp_kmers->flush_sp_kmers();	
	fprintf(stderr,"close files and remove useless variables\n");
	_sp_kmers->closeFiles();
	fprintf(stderr,"sort kmers\n");
	_sp_kmers->sortKmers();

	fprintf(stderr,"end sort kmers\n");
	return NULL;
}
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *  */
int sort_thread_new(uint64_t *readBuf, uint64_t bufferSize, int file_ind, uint8_t kmerSize)
{
	//const uint8_t kmerSize = *(uint8_t *) k;
	//char *kmerPath = "database.srt";
	//uint64_t bufferSize=BUFFERSIZE;//<<1;
	fprintf(stderr,"sort %d file\n", file_ind);
	uint64_t tempMove=(kmerSize - BUCKET_LENGTH)<<1; //32 defines kmer_plus_one
	//FILE *fpKmer=NULL;
	//fprintf(stderr,"countKmer[0]=%lu\n",countKmer[0]);

	for(uint64_t i=1;i<BUCKET_CAPACITY;i++)	countKmer[i]=countKmer[i-1]+countKmer[i];
		
	totalKmerNum=countKmer[BUCKET_CAPACITY-1];
	//fprintf(stderr, "%lu element\n", totalKmerNum);	
	totalCounterNum = totalKmerNum;
	segCount[0]=0;
	for(uint64_t i=1;i<THREAD_NUM;i++) {
		uint64_t locateNum=i*(totalKmerNum/THREAD_NUM);
		uint64_t locate=binarySearch(locateNum,countKmer,BUCKET_CAPACITY-1);
		segCount[i]=locate;
	}

	//uint64_t *readBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	//uint64_t readNum;
	
	//hashKmer=(uint64_t *)calloc(totalKmerNum + totalCounterNum,sizeof(uint64_t));

	//fpKmer=fopen(kmerPath ,"rb");
	fprintf(stderr,"start kmer distributing\n");
	//while((readNum=fread(readBuf,sizeof(uint64_t),bufferSize,fpKmer))>0) { 
		//fprintf(stderr,"readNum=%lu\n",readNum );
	for(uint64_t i=0;i<bufferSize;i+=2) {
		uint64_t tempSeq=readBuf[i]>>(tempMove);
		//uint64_t tempOcc=readBuf[i+1];
		//countKmer[tempSeq] stands for the position
		countKmer[tempSeq]--;//because it is upper bound(cannot reach)
	//fprintf(stderr,"%lu\n", 2*countKmer[tempSeq]);	
		hashKmer[2*countKmer[tempSeq]] = readBuf[i];
		hashKmer[2*countKmer[tempSeq]+1] = readBuf[i+1];
		//	
	//	hashKmer[countKmer[tempSeq]][1]=readBuf[i+1];
	}
	//}
	
	//for (uint64_t i = 2*countKmer[5749]; i < 2*countKmer[5750];i += 2)
	//{
		//cout<<transIntoChars(hashKmer[i], 31)<<"\t"<<hashKmer[i]<<endl;
	
	//}	
	//fclose(fpKmer);
	fprintf(stderr,"start sorting\n");
	time_t start=time(0);
	clock_t mul_start=clock();
	//quickSort(hashKmer+countKmer[20] ,countKmer[21]-countKmer[20],cmp);
	 
	//////////////////////////////////////multiThread//////////////////////////////////////////////////
	pthread_t myThread[THREAD_NUM];
	//pthread_attr_t attr[THREAD_NUM];
	for(uint64_t i=0;i<THREAD_NUM;i++) {
		//fprintf(stderr,"i=%lu\n",i );
		//pthread_attr_init(&attr[i]);
		//pthread_attr_setscope(&attr[i], PTHREAD_SCOPE_SYSTEM);
		int check=pthread_create( &myThread[i], NULL, multiThreadSort, (void*) i);
		//multiThreadSort((void *)i);
		if(check) {
		fprintf(stderr,"threadNum:%lu, Error - pthread_create() return code: %d\n",i,check);
		exit(EXIT_FAILURE);
		    }
    	//pthread_join( myThread[i], NULL);	
	}
	
	for(uint64_t i=0;i<THREAD_NUM;i++) {
		multiFlag=THREAD_NUM;
		pthread_join( myThread[i], NULL);
	}
	
	time_t finish=time(0);
	clock_t mul_end=clock();
	long duration=(finish-start);
	double mul_dur=(double)(mul_end-mul_start)/CLOCKS_PER_SEC;
	fprintf(stderr,"sort time is %lf (%ld)\n",mul_dur,duration );
	/////////////////////////////////////end///////////////////////////////////////////////////////////
	
	/*
	for (i = 0; i < totalKmerNum-1; ++i)
	{
		if(hashKmer[i][0]>hashKmer[i+1][0]) 
		{
			fprintf(stderr,"hashKmer[%lu]=%lu, hashKmer[%lu]=%lu\n",i,hashKmer[i][0],i+1,hashKmer[i+1][0]);
		}
	}
	*/
	char file_name[20] ={0};
	sprintf(file_name, ".%d", file_ind);

	FILE *fpKmer = fopen(file_name, "wb");	
	 
	//fpKmer=fopen(,"wb");
	//for (uint64_t i=0;i<2*totalKmerNum;i+=2) {
		//fprintf(stderr,"%lu\n",hashKmer[i]);
		//cout<<transIntoChars(hashKmer[i], kmerSize)<<"\t";
		//printDegrees(hashKmer[i+1]>>32);
		//getchar();
	//}
	
	fwrite(&totalKmerNum, sizeof(uint64_t), 1, fpKmer);

	fwrite(hashKmer,sizeof(uint64_t),totalKmerNum + totalCounterNum,fpKmer);
	
	fclose(fpKmer);
	//free(kmerPath);
	//if (readBuf)
		//free(readBuf);
	//if (hashKmer)
		//free(hashKmer);
	//if (countKmer)
		//free(countKmer);
	return 0;
}

void* sort_thread(void *k)
{
	const uint8_t kmerSize = *(uint8_t *) k;
	char *kmerPath = "database.srt";
	uint64_t bufferSize=BUFFERSIZE;//<<1;
	uint64_t tempMove=(kmerSize - BUCKET_LENGTH)<<1; //32 defines kmer_plus_one
	FILE *fpKmer=NULL;
	fprintf(stderr,"countKmer[0]=%lu\n",countKmer[0]);

	for(uint64_t i=1;i<BUCKET_CAPACITY;i++)	countKmer[i]=countKmer[i-1]+countKmer[i];
		
	totalKmerNum=countKmer[BUCKET_CAPACITY-1];
	totalCounterNum = totalKmerNum;
	segCount[0]=0;
	for(uint64_t i=1;i<THREAD_NUM;i++) {
		uint64_t locateNum=i*(totalKmerNum/THREAD_NUM);
		uint64_t locate=binarySearch(locateNum,countKmer,BUCKET_CAPACITY-1);
		segCount[i]=locate;
	}

	uint64_t *readBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t readNum;
	
	hashKmer=(uint64_t *)calloc(totalKmerNum + totalCounterNum,sizeof(uint64_t));

	fpKmer=fopen(kmerPath ,"rb");
	fprintf(stderr,"start kmer distributing\n");
	while((readNum=fread(readBuf,sizeof(uint64_t),bufferSize,fpKmer))>0) { 
		//fprintf(stderr,"readNum=%lu\n",readNum );
		for(uint64_t i=0;i<readNum;i+=2) {
			uint64_t tempSeq=readBuf[i]>>(tempMove);
			//uint64_t tempOcc=readBuf[i+1];
			//countKmer[tempSeq] stands for the position
			countKmer[tempSeq]--;//because it is upper bound(cannot reach)
			
			hashKmer[2*countKmer[tempSeq]] = readBuf[i];
			hashKmer[2*countKmer[tempSeq]+1] = readBuf[i+1];
			//	
		//	hashKmer[countKmer[tempSeq]][1]=readBuf[i+1];
		}
	}
	
	//for (uint64_t i = 2*countKmer[5749]; i < 2*countKmer[5750];i += 2)
	//{
		//cout<<transIntoChars(hashKmer[i], 31)<<"\t"<<hashKmer[i]<<endl;
	
	//}	
	fclose(fpKmer);
	fprintf(stderr,"start sorting\n");
	time_t start=time(0);
	clock_t mul_start=clock();
	//quickSort(hashKmer+countKmer[20] ,countKmer[21]-countKmer[20],cmp);
	 
	//////////////////////////////////////multiThread//////////////////////////////////////////////////
	pthread_t myThread[THREAD_NUM];
	//pthread_attr_t attr[THREAD_NUM];
	for(uint64_t i=0;i<THREAD_NUM;i++) {
		//fprintf(stderr,"i=%lu\n",i );
		//pthread_attr_init(&attr[i]);
		//pthread_attr_setscope(&attr[i], PTHREAD_SCOPE_SYSTEM);
		int check=pthread_create( &myThread[i], NULL, multiThreadSort, (void*) i);
		//multiThreadSort((void *)i);
		if(check) {
		fprintf(stderr,"threadNum:%lu, Error - pthread_create() return code: %d\n",i,check);
		exit(EXIT_FAILURE);
		    }
    	//pthread_join( myThread[i], NULL);	
	}
	
	for(uint64_t i=0;i<THREAD_NUM;i++) {
		multiFlag=THREAD_NUM;
		pthread_join( myThread[i], NULL);
	}
	
	time_t finish=time(0);
	clock_t mul_end=clock();
	long duration=(finish-start);
	double mul_dur=(double)(mul_end-mul_start)/CLOCKS_PER_SEC;
	fprintf(stderr,"sort time is %lf (%ld)\n",mul_dur,duration );
	/////////////////////////////////////end///////////////////////////////////////////////////////////
	
	/*
	for (i = 0; i < totalKmerNum-1; ++i)
	{
		if(hashKmer[i][0]>hashKmer[i+1][0]) 
		{
			fprintf(stderr,"hashKmer[%lu]=%lu, hashKmer[%lu]=%lu\n",i,hashKmer[i][0],i+1,hashKmer[i+1][0]);
		}
	}
	*/
	 
	fpKmer=fopen(kmerPath,"wb");
	//for (uint64_t i=0;i<2*totalKmerNum;i+=2) {
		//fprintf(stderr,"%lu\n",hashKmer[i]);
		//cout<<transIntoChars(hashKmer[i], kmerSize)<<"\t";
		//printDegrees(hashKmer[i+1]>>32);
		//getchar();
	//}
	
	fwrite(&totalKmerNum, sizeof(uint64_t), 1, fpKmer);

	fwrite(hashKmer,sizeof(uint64_t),totalKmerNum + totalCounterNum,fpKmer);
	
	fclose(fpKmer);
	//free(kmerPath);
	if (readBuf)
		free(readBuf);
	if (hashKmer)
		free(hashKmer);
	if (countKmer)
		free(countKmer);
	return 0;
}
uint64_t rcKmer(uint64_t p_kmerValue, uint8_t p_kmerSize) {
	uint64_t rckv = -1;
	uint64_t mask = 0x3;

	for (uint8_t i=0; i<p_kmerSize; ++i) {
		rckv = (rckv<<2)|(p_kmerValue & mask );   
		p_kmerValue >>= 2;	
	}
	return (~rckv);

}

uint64_t rcOcc(uint64_t occ)
{
	uint64_t in_out_degree = occ >> 32;
	
	uint64_t taxid = occ & 0XFFFFFFFF;

	uint64_t check_bit = 0x1; 
	uint64_t rv = 0;	
	for (int i = 0; i < 5; ++i) {
		if (in_out_degree & check_bit) {
			rv |= set_bit[i];			
		}
		in_out_degree >>= 1;
	}
	rv <<= 5;
	for (int i = 0; i < 5; ++i) {
		if (in_out_degree & check_bit) {
			rv |= set_bit[i];			
		}
		in_out_degree >>= 1;
	}
	return (rv << 32) | taxid;
}


void *multiThreadSort(void *arg)
{
	clock_t start=clock();
	uint64_t num=(uint64_t) arg;
	uint64_t i;
	uint64_t low,up;
	low=segCount[num];
	if(num<THREAD_NUM-1)	up=segCount[num+1];

	else 			up=BUCKET_CAPACITY-1;
	
	//fprintf(stderr,"thread:%lu  [%lu,%lu]\n",num,low,up);
	for(i=low;i<up;i++) {
		//fprintf(stderr,"thread:%lu  %lu\n",num,countKmer[i+1] - countKmer[i]);
		quickSort(hashKmer+2 *countKmer[i] ,countKmer[i+1]-countKmer[i],cmp,1);
		//quickSort((void *)&hashKmer[countKmer[i]],countKmer[i+1]-countKmer[i],2*sizeof(uint64_t),cmp);
		//qsort(hashKmer[countKmer[i]],(countKmer[i+1]-countKmer[i]),2*sizeof(uint64_t),cmp);
	}
	if(i==BUCKET_CAPACITY-1) {
		quickSort(hashKmer+ 2 * countKmer[i] ,totalKmerNum-countKmer[i],cmp,1);
		//quickSort((void *)&hashKmer[countKmer[i]],totalKmerNum-countKmer[i],2*sizeof(uint64_t),cmp);
		//qsort(hashKmer[countKmer[i]],(totalKmerNum-countKmer[i]),2*sizeof(uint64_t),cmp);
	}
	
	clock_t finish=clock();
	double duration=(double)(finish-start)/CLOCKS_PER_SEC;
	duration=duration/multiFlag;
	//fprintf(stderr,"thread:%lu  multiFlag=%lu cost:%lf\n",num,multiFlag,duration );
	return NULL;
}
void quickSort(uint64_t *base, uint64_t num, int comp(const void *, const void *), int type)
{
	if(num<2) return;
	uint64_t lo=0, hi=num-1, mid;
	uint64_t d_lo, d_hi, d_mid;
	uint64_t loguy, higuy;
	uint64_t d_loguy, d_higuy;

	uint64_t lostk[STKSIZ], histk[STKSIZ];
	int stkptr=0;
	size_t size;
	while(1) {
		size=hi-lo+1;
		mid=lo+(size>>1);
		d_lo = lo << 1;
		d_mid = mid << 1;
		d_hi = hi << 1;
		if(comp(base+d_lo, base+d_mid) > 0) swap(base+d_lo, base+d_mid, 2);
		if(comp(base+d_lo, base+d_hi) > 0) swap(base+d_lo, base+d_hi, 2);
		if(comp(base+d_mid, base+d_hi) > 0) swap(base+d_mid, base+d_hi, 2);	
		loguy=lo;
		higuy=hi;
		while(1) {
			do {
				loguy++;
				d_loguy = loguy << 1;
			}while(loguy<=hi && comp(base+d_loguy, base+d_mid)<=0);
			do {
				higuy--;
				d_higuy = higuy << 1;
			}while(higuy>mid && comp(base+d_higuy, base+d_mid)>0);
			if(higuy<loguy) break;

			swap(base+d_loguy, base+d_higuy,2);
			if(mid==higuy) {
				mid=loguy;
				d_mid = mid << 1;
			}
		}
		//swap((void *)base[mid],(void *)base[higuy]);
		while(higuy>lo&&comp(base+d_higuy, base+d_mid)==0)	{higuy--; d_higuy = higuy << 1;}

		if ( higuy - lo >= hi - loguy )	{ //deal with the smaller block firstx
			if (lo < higuy){
				lostk[stkptr] = lo;  
				histk[stkptr] = higuy;  
				++stkptr;  
			}                           // save big recursion for later   
			if (loguy < hi) {  
				lo = loguy;  
				continue;           // do small recursion  
			}
		}  else {
			if (loguy < hi) {  
				lostk[stkptr] = loguy;  
				histk[stkptr] = hi;  
				++stkptr;               // save big recursion for later   
			} 
			if (lo < higuy) {  
				hi = higuy;  
				continue;           // do small recursion  
			} 
		}
		--stkptr;  
		if (stkptr >= 0) {  
			lo = lostk[stkptr];  
			hi = histk[stkptr];  
			continue;           // pop subarray from stack  
		}  else  return; 
	}
	
}

void quickSort(uint64_t *base, uint64_t num, int comp(const void *, const void *))
{
	if(num<2) return;
	uint64_t lo=0, hi=num-1, mid;
	int64_t loguy, higuy;
	uint64_t lostk[STKSIZ], histk[STKSIZ];
	int stkptr=0;
	size_t size;
	while(1) {
		size=hi-lo+1;
		mid=lo+(size>>1);
		if(comp(base+lo, base+mid) > 0) swap(base+lo, base+mid);
		if(comp(base+lo, base+hi) > 0) swap(base+lo, base+hi);
		if(comp(base+mid, base+hi) > 0) swap(base+mid, base+hi);	
		loguy=lo;
		higuy=hi;
		while(1) {
			do {
				loguy++;
			}while(loguy<=hi && comp( base+loguy, base+mid)<=0);
			do {
				higuy--;
			}while(higuy>mid && comp(base+higuy, base+mid)>0);
			if(higuy<loguy) break;

			swap(base+loguy, base+higuy);
			if(mid==higuy) mid=loguy;
		}
		//swap((void *)base[mid],(void *)base[higuy]);
		while(higuy>lo&&comp(base+higuy, base+mid)==0)	higuy--;

		if ( higuy - lo >= hi - loguy )	{ //deal with the smaller block firstx
			if (lo < higuy){
				lostk[stkptr] = lo;  
				histk[stkptr] = higuy;  
				++stkptr;  
			}                           // save big recursion for later   
			if (loguy < hi) {  
				lo = loguy;  
				continue;           // do small recursion  
			}
		}  else {
			if (loguy < hi) {  
				lostk[stkptr] = loguy;  
				histk[stkptr] = hi;  
				++stkptr;               // save big recursion for later   
			} 
			if (lo < higuy) {  
				hi = higuy;  
				continue;           // do small recursion  
			} 
		}
		--stkptr;  
		if (stkptr >= 0) {  
			lo = lostk[stkptr];  
			hi = histk[stkptr];  
			continue;           // pop subarray from stack  
		}  else  return; 
	}
	
}
void swap(const void *a, const void *b, int type)
{
	uint64_t *temp_a = (uint64_t *) a;
	uint64_t *temp_b = (uint64_t *) b;	
	//fprintf(stderr, "run %lu\t",temp_a - hashKmer);
	uint64_t v = *temp_a ; 	
	*temp_a = *temp_b;
	*temp_b = v;

	++temp_a;
	++temp_b;

	v = *temp_a ; 	
	*temp_a = *temp_b;
	*temp_b = v;
	//uint64_t (*va=((uint64_t (*)[2])a),(*vb)[2]=((uint64_t (*)[2])b);
	
	
	//temp=va[0][0];
	//va[0][0]=vb[0][0];
	//vb[0][0]=temp;
	//temp=va[0][1];
	//va[0][1]=vb[0][1];
	//vb[0][1]=temp;
}
void swap(const void *a, const void *b)
{
	uint64_t temp = *(uint64_t *) a;
	*(uint64_t *) a = *(uint64_t *) b;
	*(uint64_t *) b = temp;
	//uint64_t (*va=((uint64_t (*)[2])a),(*vb)[2]=((uint64_t (*)[2])b);
	
	
	//temp=va[0][0];
	//va[0][0]=vb[0][0];
	//vb[0][0]=temp;
	//temp=va[0][1];
	//va[0][1]=vb[0][1];
	//vb[0][1]=temp;
}


int cmp(const void *a, const void *b)
{
	uint64_t va = *(uint64_t *) a,vb= *(uint64_t *) b;

	//fprintf(stderr,"va=%lu; vb=%lu\n",va,vb );
	if(va<vb) return -1;
	else if(va==vb) return 0;
	else return 1;
}


uint64_t binarySearch(uint64_t mk, uint64_t *target, int64_t up)//return the upper bound
{
 	 int64_t low=0,mid=0;
 	 while(low<=up) {
     		mid=(low+up)>>1;
     		if(mk<target[mid])   		up=mid-1;
     		else if(mk>target[mid])  	low=mid+1;
     		else 				return mid;
  	}
 	return low;
}


void *keyway_sort_thread(void *)
{
	FILE *fp[useFileNum];
	uint64_t left_kmers[useFileNum]; //counter of kmers left in a file
	uint64_t *buffer = new uint64_t[BUFFERSIZE];
	bool useBuffer = false;
	if (buffer) useBuffer = true; 	
	for (int i = 0; i < useFileNum; ++i) {
		char file_name[10];
		sprintf(file_name, ".%d", i);
		fp[i] = fopen(file_name, "rb");
		if (!fp[i]) {
			fprintf(stderr,"Failed to open .%d file", i);
			exit(1);
		} else {
			fread(left_kmers+i, 8, 1, fp[i]);
		}	
	}
		
	//heapsort
		
	std::priority_queue<ele> prq;
	
	FILE *fpKmer=NULL;
	char *kmerPath = "database.srt";
	fpKmer=fopen(kmerPath ,"wb");
	uint64_t zong = 0;	
	fwrite(&zong, 8, 1, fpKmer);

	ele e;
	for (int i = 0 ; i < useFileNum; ++i) {
		if (left_kmers[i]) {
			fread(&e.key, 8, 1, fp[i]);
			fread(&e.val, 8, 1, fp[i]);
			e.ind = i;
			--left_kmers[i];
			prq.push(e);
		}
	}	
	
	if (!useBuffer) {
		while (!prq.empty()) {
			uint64_t key = prq.top().key;
			uint64_t val = prq.top().val;
			int ind = prq.top().ind;
		
			prq.pop();
			++zong;
			if (left_kmers[ind]) {
				fread(&e.key, 8, 1, fp[ind]);
				fread(&e.val, 8, 1, fp[ind]);
				e.ind = ind;
				--left_kmers[ind];
				prq.push(e);
			}
			fwrite(&key, 8, 1, fpKmer);	
			fwrite(&val, 8, 1, fpKmer);	
		}	
	} else {
		uint64_t bufPoint = 0;
		while (!prq.empty()) {
			uint64_t key = prq.top().key;
			uint64_t val = prq.top().val;
			int ind = prq.top().ind;
		
			prq.pop();
			++zong;
			if (left_kmers[ind]) {
				fread(&e.key, 8, 1, fp[ind]);
				fread(&e.val, 8, 1, fp[ind]);
				e.ind = ind;
				--left_kmers[ind];
				prq.push(e);
			}
			buffer[bufPoint++] = key;
			buffer[bufPoint++] = val;
			if (bufPoint >= BUFFERSIZE) {
				fwrite(buffer, sizeof(uint64_t), BUFFERSIZE, fpKmer);	
				bufPoint = 0;	
			}
		}	
		if (bufPoint) {
				fwrite(buffer, sizeof(uint64_t), bufPoint, fpKmer);	
		}
		delete[] buffer;	
	}
	
	fseek(fpKmer, 0, SEEK_SET);
	fwrite(&zong, 8 , 1, fpKmer);	
	fclose(fpKmer);
	
	for (int i = 0; i < FILE_NUM_LIMIT; ++i) {
		char file_name[10];
		sprintf(file_name, ".%d", i);
		remove(file_name);
	}
	
	return NULL;
}


