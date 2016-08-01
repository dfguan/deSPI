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

#include "sort.hpp"
#include "jreader.hpp"
#define BUFFERSIZE (1<<20)
//#define KMER_LENGTH_PlusOne 32
#define BUCKET_LENGTH 12 //need less than KMER_LENGTH_PlusOne
#define BUCKET_CAPACITY (1<<(BUCKET_LENGTH<<1))
#define THREAD_NUM 32 //should not surpass 2^8-1
#define STKSIZ 40
//#define KMER_LENGTH 17

uint64_t *countKmer = NULL;
uint64_t *hashKmer = NULL;
uint64_t trans[256];
uint64_t totalKmerNum;
uint64_t segCount[THREAD_NUM];
uint64_t multiFlag=1;


int main(int argc, char *argv[]) {
 return genSortedKmers(argv[1], uint8_t(atoi(argv[2])));
}
int genSortedKmers(char *disorderedKmerpath, const uint8_t kmerSize)
{
    	//////////////////////////////////////kmercounting///////////////////////////////////////////
	countKmer=(uint64_t *)calloc(BUCKET_CAPACITY,sizeof(uint64_t));
	uint64_t bufferSize=BUFFERSIZE;//<<1;
	uint64_t occBuf;
	uint64_t *writeBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t bufPoint=0;
	FILE *fpKmer=NULL;
	char *kmerPath = "sortedKmers";
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
	input_file.seekg(jr.header_size(), ios_base::beg);

	// Create a copy of the offsets array for use as insertion positions
	//FILE * input_file = fopen(jr.get_db_name().c_str(),"rb");
	
	//if (!input_file) return FILE_OPEN_ERROR;

	for (uint64_t i = 0; i < key_ct; i++) {
		input_file.read(pair, pair_size);
		uint64_t kmerValue,rcKmerValue = 0;
		memcpy(&kmerValue, pair, key_len);
		
		++countKmer[kmerValue >> tempMove];
		//cout<<kmerSize<<endl;	
		rcKmerValue = rcKmer(kmerValue,kmerSize);
		++countKmer[rcKmerValue >> tempMove];
		//cout<<kmerValue<<"\t"<<transIntoChars(kmerValue, kmerSize)<<endl;
		//cout<<rcKmerValue<<"\t"<<transIntoChars(rcKmerValue,kmerSize)<<endl;
				
		writeBuf[bufPoint++] = kmerValue;
		writeBuf[bufPoint++] = rcKmerValue;
		//writeBuf[bufPoint++]=occBuf;
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
	free(writeBuf);
	fclose(fpKmer);
	//fclose(fp1);
	/////////////////////////////////////////kmer distributing///////////////////////////////////////////////
	fprintf(stderr,"countKmer[0]=%lu\n",countKmer[0]);

	for(uint64_t i=1;i<BUCKET_CAPACITY;i++)	countKmer[i]=countKmer[i-1]+countKmer[i];
		
	totalKmerNum=countKmer[BUCKET_CAPACITY-1];
	segCount[0]=0;
	for(uint64_t i=1;i<THREAD_NUM;i++) {
		uint64_t locateNum=i*(totalKmerNum/THREAD_NUM);
		uint64_t locate=binarySearch(locateNum,countKmer,BUCKET_CAPACITY-1);
		segCount[i]=locate;
	}

	uint64_t *readBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t readNum;
	hashKmer=(uint64_t *)calloc(totalKmerNum,sizeof(uint64_t));
	fpKmer=fopen(kmerPath ,"rb");
	fprintf(stderr,"start kmer distributing\n");
	while((readNum=fread(readBuf,sizeof(uint64_t),bufferSize,fpKmer))>0) { 
		//fprintf(stderr,"readNum=%lu\n",readNum );
		for(uint64_t i=0;i<readNum;++i) {
			uint64_t tempSeq=readBuf[i]>>(tempMove);
			//uint64_t tempOcc=readBuf[i+1];
			//countKmer[tempSeq] stands for the position
			countKmer[tempSeq]--;//because it is upper bound(cannot reach)
			hashKmer[countKmer[tempSeq]] = readBuf[i];
		//	hashKmer[countKmer[tempSeq]][1]=readBuf[i+1];
		}
	}
	fclose(fpKmer);
	fprintf(stderr,"start sorting\n");
	time_t start=time(0);
	clock_t mul_start=clock();
	//quickSort(hashKmer+countKmer[20] ,countKmer[21]-countKmer[20],cmp);
	 
	//////////////////////////////////////multiThread//////////////////////////////////////////////////
	pthread_t myThread[THREAD_NUM];
	//pthread_attr_t attr[THREAD_NUM];
	for(uint64_t i=0;i<THREAD_NUM;i++) {
		fprintf(stderr,"i=%lu\n",i );
		//pthread_attr_init(&attr[i]);
		//pthread_attr_setscope(&attr[i], PTHREAD_SCOPE_SYSTEM);
		int check=pthread_create( &myThread[i], NULL, multiThreadSort, (void*) i);
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
//	for (uint64_t i=0;i<totalKmerNum;++i) {
//		fprintf(stderr,"%lu\n",hashKmer[i]);
	//	getchar();
//	}
	
	fwrite(&totalKmerNum, sizeof(uint64_t), 1, fpKmer);

	fwrite(hashKmer,sizeof(uint64_t),totalKmerNum,fpKmer);
	
	fclose(fpKmer);
	//free(kmerPath);
	free(readBuf);
	free(hashKmer);
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
void *multiThreadSort(void *arg)
{
	clock_t start=clock();
	uint64_t num=(uint64_t) arg;
	uint64_t i;
	uint64_t low,up;
	low=segCount[num];
	if(num<THREAD_NUM-1)	up=segCount[num+1];

	else 			up=BUCKET_CAPACITY-1;
	
	fprintf(stderr,"thread:%lu  [%lu,%lu]\n",num,low,up);
	for(i=low;i<up;i++) {
		quickSort(hashKmer+countKmer[i] ,countKmer[i+1]-countKmer[i],cmp);
		//quickSort((void *)&hashKmer[countKmer[i]],countKmer[i+1]-countKmer[i],2*sizeof(uint64_t),cmp);
		//qsort(hashKmer[countKmer[i]],(countKmer[i+1]-countKmer[i]),2*sizeof(uint64_t),cmp);
	}
	if(i==BUCKET_CAPACITY-1) {
		quickSort(hashKmer+ countKmer[i] ,totalKmerNum-countKmer[i],cmp);
		//quickSort((void *)&hashKmer[countKmer[i]],totalKmerNum-countKmer[i],2*sizeof(uint64_t),cmp);
		//qsort(hashKmer[countKmer[i]],(totalKmerNum-countKmer[i]),2*sizeof(uint64_t),cmp);
	}
	
	clock_t finish=clock();
	double duration=(double)(finish-start)/CLOCKS_PER_SEC;
	duration=duration/multiFlag;
	fprintf(stderr,"thread:%lu  multiFlag=%lu cost:%lf\n",num,multiFlag,duration );
	return NULL;
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


