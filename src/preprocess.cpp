// =====================================================================================
//
//       Filename:  preprocess.cpp
//
//    Description:  Generate and Sort Kmers contain "$#" 
//
//        Version:  1.0
//        Created:  10/20/2015 09:25:48 AM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Dengfeng Guan, <dfguan@hit.edu.cn>
//   Organization:  Harbin Institute of Technology
//
// =====================================================================================

#include "preprocess.hpp"


//quick sort kmers carrying specail characters
//uint8_t _kmer;

//uint64_t *counter;

///#include "debruijin.hpp"
uint8_t _kmer;

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

uint64_t *counter;

uint32_t *gidCounter;

int outputDeb(uint64_t *p_kmerValue, uint16_t *p_kmerInfo, uint64_t p_kmerNum) 
{
	uint64_t p_size = p_kmerNum - 1;
	char Chars[] = {'A','C','G','T'};
	cout<<"debruijn:" <<endl;
	for(uint64_t i=0; i<p_size;++i) {
		cout<<transIntoChars(p_kmerValue[i], _kmer)<<"\t";
		uint8_t inEdges =((p_kmerInfo[i]>>4)&0xF); 
		uint8_t outEdges = ((p_kmerInfo[i])&0xF);
		for(uint8_t j=0;j<4;++j,inEdges >>= 1) if(inEdges&1) 	cout<<Chars[j]<<",";
		
		cout<<"\t";

		for(uint8_t j=0;j<4;++j,outEdges >>= 1) if(outEdges&1) 	cout<<Chars[j]<<","; 
		cout<<endl;	
	}
	cout<<"end"<<endl;
	return NORMAL_EXIT;


}
//fore part of kmer will be use to limit the range 
int print(char *str) 
{
	for(uint8_t i=0; i <_kmer; ++i) {
		fprintf(stderr,"%c",str[i]);

	}

	fprintf(stderr,"\n");
	return NORMAL_EXIT;

}

int build_deb(const char *refPath, uint64_t  * kmerValue, uint16_t *kmerInfo, uint32_t *kmerTID, uint64_t kmerNum, vector<uint64_t>& p_heads, vector<uint64_t>&  p_tails) //the path indicate a file containing all sorted k+1mers  
{
	//uint64_t prevValue = -1;

	//uint64_t *values = new uint64_t[0x2710];//store 10000 value;
	
	//uint64_t useless;
	uint64_t counterSize = (1<<(PREINDEXLEN<<1)) + 1;

	counter = new uint64_t[counterSize]();
		
	

	uint16_t move = (_kmer - PREINDEXLEN) << 1;


	//uint64_t nodeCounter = 0;//index

	// initiate out edges

/*	
	while(!feof(fp)) {
		uint32_t b_size = fread(values,sizeof(uint64_t),0x2710,fp);
		//fprintf(stderr,"%lu\t%lu\n",value[0],useless);
		for(uint32_t i=0; i<b_size;++i) {
			++counter[values[i]>>move]; 
			kmerValue.push_back(values[i]);
			kmerInfo.push_back(0);
			kmerTID.push_back(0);
		
		}
	} 
*/
	//fprintf(stderr,"gogo\n"); 
		/*	
		uint64_t kValue = value[0] >> 2;
		uint8_t outEdge = value[0] & 3;

		if (prevValue != kValue ) {
			kmerValue.push_back(kValue);
			kmerInfo.push_back(0);
			//p_taxon_ID.push_back(0);
			
			prevValue = kValue;
			++nodeCounter;

		} 
		kmerInfo[nodeCounter - 1] = setOutEdge(uint8_t(kmerInfo[nodeCounter-1]),outEdge);
		*/
	


	//fclose(fp);
	//initiate ranges 
	uint64_t sum = 0;
	uint64_t temp = 0;

	for (uint64_t i=0; i<kmerNum; ++i) ++counter[kmerValue[i]>>move];

	for(uint64_t i=0; i<counterSize;++i) {
		temp = counter[i];
		counter[i] = sum;
		sum = sum + temp;
	} 
	// initate preparations
	// start to initiate inout degrees this
	gzFile fp_z;
	
	kseq_t *trunk;
	
	fp_z = gzopen(refPath,"r");

	if (!fp_z)
		return FILE_OPEN_ERROR;
	trunk = kseq_init(fp_z);
	

	uint64_t mask = ~((uint64_t)0x3 <<((_kmer-1)<<1));
#ifdef CONSIDER_BOTH_ORIENTATION
	uint64_t mask_rev[] = {(uint64_t)0x3 << ((_kmer - 1)<<1),  (uint64_t)0x2 << ((_kmer - 1)<<1), (uint64_t)0x1 << ((_kmer - 1)<<1), 0};
#endif
	//uint32_t gid = 0;	
	uint32_t tid = 0;
	
	uint16_t lastChar;
	while (kseq_read(trunk)>=0) {
		//fprintf(stderr,"haha\n");
		tid = getTID(trunk->name.s);
		//fprintf(stderr,"kago\n");
		//uint32_t gid_loc = binSearch(p_taxonIDTab, gidCounter[gid>>10], gidCounter[(gid>>10)+1] - 1, gid);
		//if (~gid_loc)  tid = p_taxonIDTab[gid_loc].taxonID;   
		//else 
		//tid = 0;
		for(uint64_t i=0; i<trunk->seq.l;++i) {
			if (Bit[trunk->seq.s[i]] != 4) {
				uint64_t start = i;
				while(Bit[trunk->seq.s[++i]]!=  4 && i < trunk->seq.l);
				//uint64_t end = i - 1;
				if (i - start > _kmer) {
					//intiate head 
					uint64_t key = transIntoBits(trunk->seq.s+start, _kmer);
					
					uint64_t loc = binSearch(kmerValue, counter[key>>move],counter[(key>>move)+1]-1, key);
					
					kmerInfo[loc] = setOutEdge((uint8_t) kmerInfo[loc], Bit[trunk->seq.s[start + _kmer]]);
					
					kmerTID[loc] = LCA(kmerTID[loc], tid); 
					
					p_heads.push_back(loc);
					
#ifdef CONSIDER_BOTH_ORIENTATION
					uint64_t keyR = transRCIntoBits(trunk->seq.s+start, _kmer);

					loc = binSearch(kmerValue, counter[keyR>>move],counter[(keyR>>move)+1]-1, keyR);
					
					//lastChar = Bit[trunk->seq.s[start + _kmer]]^0x3;					
					//kmerInfo[loc] = setInEdge((uint8_t) kmerInfo[loc],(uint8_t)lastChar);
					
					kmerInfo[loc] = setInEdge((uint8_t) kmerInfo[loc],Bit[trunk->seq.s[start + _kmer]]^0x3);
					//kmerInfo[loc] |= (lastChar << 10);
					
					kmerTID[loc] = LCA(kmerTID[loc], tid); 
					
					p_tails.push_back(loc);
#endif
					for (uint64_t j = start + 1; j < i - _kmer; ++j) {
						key = ((key & mask) << 2)| Bit[trunk->seq.s[j+_kmer-1]];
						//print(trunk->seq.s+j);		
						loc = binSearch(kmerValue, counter[key>>move],counter[(key>>move)+1]-1, key);
					
						//lastChar = Bit[trunk->seq.s[j-1]];
						kmerInfo[loc] = setInEdge((uint8_t) kmerInfo[loc], Bit[trunk->seq.s[j-1]]);
						//kmerInfo[loc] = setInEdge((uint8_t) kmerInfo[loc],(uint8_t)lastChar);
						
						//kmerInfo[loc] |= (lastChar << 10);
						
						kmerInfo[loc] = setOutEdge((uint8_t) kmerInfo[loc], Bit[trunk->seq.s[j+_kmer]]);
						//fprintf(stderr,"cycling ....\n");	
						kmerTID[loc] = LCA(kmerTID[loc], tid); 
#ifdef CONSIDER_BOTH_ORIENTATION
						
						keyR = (keyR >> 2) | mask_rev[Bit[trunk->seq.s[j + _kmer - 1]]];

						loc = binSearch(kmerValue, counter[keyR>>move],counter[(keyR>>move)+1]-1, keyR);
						
						kmerInfo[loc] = setOutEdge((uint8_t) kmerInfo[loc], (Bit[trunk->seq.s[j-1]]^0x3));
						
						//lastChar = Bit[trunk->seq.s[j + _kmer]]^0x3;					
							
						kmerInfo[loc] = setInEdge((uint8_t) kmerInfo[loc], Bit[trunk->seq.s[j + _kmer]]^0x3);
						//kmerInfo[loc] = setInEdge((uint8_t) kmerInfo[loc],(uint8_t)lastChar);
						//kmerInfo[loc] |= (lastChar << 10);
						kmerTID[loc] = LCA(kmerTID[loc], tid); 
#endif						
					} 
					
					key = ((key & mask) << 2)| Bit[trunk->seq.s[i-1]];
				
					loc = binSearch(kmerValue, counter[key>>move],counter[(key>>move)+1]-1, key);
					//lastChar = Bit[trunk->seq.s[i - _kmer -1]];
					kmerInfo[loc] = setInEdge((uint8_t) kmerInfo[loc],Bit[trunk->seq.s[i - _kmer -1]]);
					//kmerInfo[loc] = setInEdge((uint8_t) kmerInfo[loc],(uint8_t)lastChar);
					
					//kmerInfo[loc] |= (lastChar<<10);	
					kmerTID[loc] = LCA(kmerTID[loc], tid); 

					p_tails.push_back(loc);	
					//fprintf(stderr,"done\n");	
#ifdef CONSIDER_BOTH_ORIENTATION
					keyR = (keyR >> 2) | mask_rev[Bit[trunk->seq.s[i - 1]]];

					loc = binSearch(kmerValue, counter[keyR>>move],counter[(keyR>>move)+1]-1, keyR);
					
					kmerInfo[loc] = setOutEdge((uint8_t) kmerInfo[loc], (Bit[trunk->seq.s[i - _kmer -1]]^0x3));
					
					kmerTID[loc] = LCA(kmerTID[loc], tid); 
					p_heads.push_back(loc);	
					
#endif
				} else {
					if (i - start == _kmer) {
						uint64_t key = transIntoBits(trunk->seq.s + start, _kmer);
						
						uint64_t loc = binSearch(kmerValue, counter[key>>move],counter[(key>>move)+1]-1, key);
						
						p_heads.push_back(loc);
			
						kmerTID[loc] = LCA(kmerTID[loc], tid); 
					
						p_tails.push_back(loc);	
				
#ifdef CONSIDER_BOTH_ORIENTATION
						uint64_t keyR = transRCIntoBits(trunk->seq.s+start, _kmer);
				
						loc = binSearch(kmerValue, counter[keyR>>move],counter[(keyR>>move)+1]-1, keyR);

						p_heads.push_back(loc);
						kmerTID[loc] = LCA(kmerTID[loc], tid); 
						p_tails.push_back(loc);
#endif
					}
				
				}

			}
		}
	}
	kseq_destroy(trunk);
	gzclose(fp_z);	

	


	return NORMAL_EXIT;
 
	//kmerValue.push_back(-1);
	//kmerInfo.push_back(0);
	//p_taxon_ID.push_back(0);
	//initiate In edges
	/*
	fprintf(stderr,"initiate in\n");
	uint64_t mask = ~((uint64_t)0x3<<((_kmer-1)<<1));
	vector<uint64_t>::iterator itKmerValue = kmerValue.begin();
       	vector<uint16_t>::iterator itKmerInfo = kmerInfo.begin();
	//vector<uint32_t>::iterator itTaxonID = p_taxon_ID.begin();	
	

	for(uint64_t i=0; i<nodeCounter; ++i) {
		
		uint8_t outEdges = (uint8_t)kmerInfo[i];
		uint8_t inEdge = (uint8_t)(kmerValue[i] >>((_kmer-1) << 1)&0x3);
		uint64_t p_value = (kmerValue[i] & mask) << 2; 
		for(uint8_t j=0;j<4;++j){
			if (outEdges & 1) {
				uint64_t kValue = p_value | j;
				uint64_t pos = binSearch( kmerValue, 0, nodeCounter - 1, kValue );
				if (kmerValue[pos] != kValue) {
					fprintf(stderr,"not equal\n");
					kmerValue.insert(itKmerValue+pos,kValue); //repetive initiation wastes times
					kmerInfo.insert(itKmerInfo+pos,0);
					//p_taxon_ID.insert(itTaxonID + pos, 0);
					++nodeCounter;
				}
				kmerInfo[pos] = setInEdge((uint8_t)kmerInfo[pos],inEdge);		
			}
			outEdges >>= 1;
		} 
	}
	*/
	//done
	//
	//
	//outputDeb(kmerValue, kmerInfo);

//	return NORMAL_EXIT;

}

uint8_t setInEdge(uint8_t info, uint8_t edge)
{
	return (info|(1<<(edge + 4)));	

}

uint8_t setOutEdge(uint8_t info, uint8_t edge)
{
	return (info | (1<<edge)); 

}

uint64_t binSearch(uint64_t *kmersValue,  uint64_t lowerBound, uint64_t upperBound,uint64_t key)
{
	//if (key < kmersValue[lowerBound])		return lowerBound;
	//if (key > kmersValue[upperBound]) 		return upperBound+1;

	uint64_t mid; 
	while (lowerBound <= upperBound) {
		mid = lowerBound + ((upperBound - lowerBound) >> 1);
		
		if (kmersValue[mid] < key) {
		lowerBound = mid + 1;
		} else if (kmersValue[mid] > key) {
			upperBound = mid - 1;
		
		} else
			return mid;
	}

	return -1;
}
uint64_t transIntoBits(char *str_kmer, uint8_t len)
{
	uint64_t value = 0;
	
	for(uint8_t i=0; i<len; ++i) value = (value<<2)|Bit[str_kmer[i]];			
	
	return value;

}

uint64_t transRCIntoBits(char *str_kmer, uint8_t len) 
{
	uint64_t value = 0;
	
	for (uint8_t i=0; i<len;++i) value = (value<<2)|(Bit[str_kmer[len - i - 1]]^0x3);

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

int sortKmers(kmersSpchar* start, kmersSpchar* end)
{
	stable_sort(start, end);
	return NORMAL_EXIT;
}

uint64_t findInsertPos(uint64_t* p_kmersValue, uint64_t p_s, uint64_t p_e, uint64_t key)
{
	for(uint64_t i=p_s; i < p_e; ++i) {
		if (key <= p_kmersValue[i])
			return i;
	}
	return p_e;

}

int mergeSort(uint64_t* p_kmersValue, uint64_t kmerNum, uint16_t* p_kmersInfo, uint32_t* taxonIDTab, vector<uint32_t>& nkmerTID, kmersSpchar* p_2kmers, uint64_t kmersSpcharN, string& p_bwt, uint64_t *hash_index, uint64_t index_start_point)
{
	char Chars[] = {'A','C','G','T','#','$'};
	
	//string bwt = "";
	
	uint64_t start = 0;
	uint64_t end = kmerNum;
	//fprintf(stderr,"%lu\t%lu\n", end, p_2kmers.size());	
	int mov = (_kmer - PREINDEXLEN) << 1;
	uint64_t pre = -1;
	uint64_t suffix_v = 0;
	uint64_t ind = index_start_point;
	//vector<kmersSpchar>::iterator it = p_2kmers.begin();
	for (uint64_t z=0; z<kmersSpcharN; ++z) {
		kmersSpchar* it = p_2kmers + z; 
		uint8_t off = it->infor >> 3;
		//cout<<it->value<<endl;
		uint64_t key = it->value << ((_kmer - off)<< 1); 
		uint64_t loc = findInsertPos(p_kmersValue, start,end, key);
		for (uint64_t i = start; i < loc;++i)	{
			p_bwt += Chars[p_kmersInfo[i]>>10];
			++ind;
			//cout<<transIntoChars(p_kmersValue[i], _kmer)<<"\t"<<ind<<endl;
			if (!(ind & 0x1F)) nkmerTID.push_back(taxonIDTab[i]);
			if ((suffix_v = (p_kmersValue[i]>>mov)) != pre) {
				//fprintf(stderr,"%lu\t",suffix_v);
				hash_index[suffix_v<<1] = ind;
				hash_index[(suffix_v<<1)+1] = ind + 1;
				pre = suffix_v;
			} else 	
				++hash_index[(suffix_v<<1)+1];
			
	
			
		}
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
	}
	fprintf(stderr, "hahah%lu\t%lu\n",start, end);
	//cout<<"wwwwww"<<endl;
	for (uint64_t i = start; i < end; ++i)	{
		p_bwt += Chars[p_kmersInfo[i]>>10];
		++ind;
		//cout<<transIntoChars(p_kmersValue[i], _kmer)<<"\t"<<ind<<endl;
		if (!(ind & 0x1F)) nkmerTID.push_back(taxonIDTab[i]);
		
		if ((suffix_v = (p_kmersValue[i]>>mov)) != pre) {
			//fprintf(stderr,"%lu\t",suffix_v);
			hash_index[suffix_v<<1] = ind;
			hash_index[(suffix_v<<1)+1] = ind + 1;
			pre = suffix_v;
		} else 	
			++hash_index[(suffix_v<<1)+1];
	}

	return NORMAL_EXIT;

}

uint32_t LCA(uint32_t tid1, uint32_t tid2)
{
	if (tid1 == 0 || tid2 == 0)
      	return tid1 ? tid1 : tid2;

	set<uint32_t> pool;
	while (tid1 > 0) {
		pool.insert(tid1);
		tid1 = taxonomyTree[tid1];
	}
	while (tid2 > 0) {
		if (pool.count(tid2) > 0)
			return tid2;
		tid2 = taxonomyTree[tid2];
	}
	return 1;
	


}



uint32_t getTID(char *seq_name) 
{
	char *tokens = strtok(seq_name,"|");

	tokens = strtok(NULL,"|");

	return strtoul(tokens,NULL, 10);


}
// gid_taxonID can't supass 2^32-1
uint32_t binSearch(vector<gid_taxid>& p_taxonIDTab,uint64_t lowerBound, uint64_t upperBound, uint32_t key ) 
{
	if (p_taxonIDTab[lowerBound].gID > key || p_taxonIDTab[upperBound].gID < key) return -1;
	
	uint32_t mid;
	
	while (lowerBound <= upperBound) {
		mid = lowerBound + ((upperBound - lowerBound) >> 1);

		if (p_taxonIDTab[mid].gID < key) lowerBound = mid + 1;
		else if (p_taxonIDTab[mid].gID>key) upperBound = mid - 1;
		else return mid;
	}	
	return -1;
}

//path: reference 
//p_fkmer : first kmer of a path
//p_taxonID: corresponding between gid and uid
//correspond: taxonID
/*  
int uid2taxonID(char *path, vector<uint64_t> p_fkmer, vector<uint32_t>& correspond, vector<gid_taxid>& p_taxonID)
{
	gzFile fp;
	kseq_t *trunk;
	
	fp = gzopen(path, "r");
	
	if (!fp) return FILE_OPEN_ERROR;
	
	trunk = kseq_init(fp);

	//when read into several chromosome may be different
	uint32_t pos;
	//uint64_t ind;
	uint32_t end;
	
	uint64_t p_fkmer_size = p_fkmer.size();
	
	uint64_t p_taxID_size = p_taxonID.size(); 


	uint32_t tID;
	
	uint64_t mask =~((uint64_t)0x3 <<((_kmer-1)<<1)); 

	while(kseq_read(trunk)>=0) {
	//find a position kmers were not N;
		//get tid through gid;
		uint32_t gID = getGID(trunk->name.s);
		
		uint32_t gID_loc = binSearch(p_taxonID, 0, p_taxID_size-2, gID);

		if (~gID_loc) { 
		
			tID = p_taxonID[gID_loc].taxonID;

			end = trunk->seq.l - _kmer + 1;

			pos = findKmerPos(trunk->seq.s, 0, end);

			while(pos < end) {

				uint64_t value = transIntoBits(trunk->seq.s+pos,_kmer);
				uint64_t loca = binSearch_pre(p_fkmer, 0, p_fkmer_size-2,value);
				if (~loca) correspond[loca] = LCA(correspond[loca],tID);
				for (uint32_t i=pos+_kmer;i<end;++i){
					if (trunk->seq.s[i] == 'N') {
						pos = i + 1;
						break;
					} else {
						value = (value & mask)<<2 | Bit[trunk->seq.s[i]];
						uint64_t loca = binSearch_pre(p_fkmer, 0, p_fkmer_size-2,value);
						if (~loca) correspond[loca] = LCA(correspond[loca],tID);
					}
				} 


			}		
		}
	
	
	
	
	}	
	return NORMAL_EXIT;


}
*/
//the species number can surpass 2^32 - 2

// taxoi
uint64_t binSearch_pre(vector<uint64_t>& values, uint64_t lowerBound, uint64_t upperBound, uint64_t key )
{
 	if (values[lowerBound] > key || values[upperBound] < key) return -1;
	
	uint32_t mid;
	
	while (lowerBound <= upperBound) {
		mid = lowerBound + ((upperBound - lowerBound) >> 1);

		if (values[mid] < key) lowerBound = mid + 1;
		else if (values[mid] > key) upperBound = mid - 1;
		else return mid;
	}	
	return -1;

}
/*
int taxonID(const char *giTaxidPath, vector<gid_taxid> & taxonIDTab)
{	
	FILE *fp = fopen(giTaxidPath, "r");
	
	if (!fp)
		return FILE_OPEN_ERROR;

	gid_taxid ins;
	
	uint32_t gID;
	uint32_t taxID;
	
	gidCounter = new uint32_t[1000000](); //should be changed if gid is larger
	while(EOF != fscanf(fp, "%u %u", &gID, &taxID)) {
		ins.gID = gID;
		ins.taxonID = taxID;
		taxonIDTab.push_back(ins);
		//fprintf(stderr,"%u\t%u\n", gID, taxID);
		++gidCounter[(gID >> 10)];
	}
	uint32_t sum = 0;
	uint32_t temp = 0;
	for(uint32_t i=0; i<1000000;++i) {
		temp = gidCounter[i];
		gidCounter[i] = sum;
		sum = sum + temp;
	}


	fclose(fp);
	return NORMAL_EXIT;

}
*/

int _getline(FILE *in, char *ln)
{	
	char c;
	int count = 0;

	while ((c = getc(in))!= EOF) {
		if (c != '\n')
			ln[count++] = c;
		else 
			break;
	}

	ln[count] = '\0';
	
	return count;
}

int taxonTree(const char *taxonomyNodesPath)
{
	FILE *fp = fopen(taxonomyNodesPath,"r");

	char line[1024];
	//istring token;
	//stringstream iss;
	char *token;
	
	uint32_t tid;
	uint32_t p_tid;

	while(_getline(fp,line)) {
		token = strtok(line,"\t|");
		tid = strtoul(token,NULL,10);
		token = strtok(NULL,"\t|");
		p_tid = strtoul(token,NULL,10);
		taxonomyTree[tid] = p_tid;
	}
	//end mark
	taxonomyTree[1] = 0;
	return NORMAL_EXIT;

}


/* 
void *assignTID(void * ar)
{
	arr_p * art = (arr_p *)ar;
	uid2taxonID(art->path, art->p_fk, art->corres, art->gID_taxID);
	return NULL;	
}
*/

uint8_t outEdgesNum(uint8_t info)
{	
	uint8_t outCount = 0;
	for (uint8_t i=0; i<4; ++i) {
		if ((info>>i)&1) ++outCount;
	}	
	return outCount;
}

uint8_t inEdgesNum(uint8_t info)
{      
	info >>= 4;
	uint8_t inCount = 0;
	for (uint8_t i=0; i<4; ++i) {
		if ((info>>i)&1) ++inCount;
	}	
	return inCount;
}

int setStart(uint16_t& marker) 
{
	marker |= (1 << 9);

	return NORMAL_EXIT;

}

int setEnd(uint16_t& marker)
{	
	marker |= (1 << 8);

	return NORMAL_EXIT;
}



int cutOffMulEdges(uint64_t *p_kmerValue, uint16_t *p_kmerInfo, uint64_t kmerNum)
{
	uint64_t size = kmerNum;
	//uint64_t mask = 
	
	uint16_t move = (_kmer - PREINDEXLEN) << 1;
	uint64_t mask = ~((uint64_t)0x3 <<((_kmer-1)<<1));
	for (uint64_t i=0; i<size;++i) {
		//bool isMulIn = mulitpleIn(p_kmerInfo[i]);
		//bool isMulOut = multipleOut(p_kmerInfo[i]);
		uint8_t inCount = inEdgesNum(p_kmerInfo[i]);
		uint8_t outCount = outEdgesNum(p_kmerInfo[i]);
		//fprintf(stderr,"%lu",i);	
		if ( 1 != inCount) {
			setStart( p_kmerInfo[i]);
			if (inCount) {
				uint8_t inEdges = p_kmerInfo[i] >> 4;
				for (uint8_t j=0; j< 4; ++j,inEdges >>= 1) {
					if (inEdges&1) {
						uint64_t key = (p_kmerValue[i] >> 2)|((uint64_t)j<<((_kmer-1)<<1));
						uint64_t loc = binSearch(p_kmerValue, counter[key>>move], counter[(key>>move)+1]-1,key );
						setEnd( p_kmerInfo[loc]);
					} 
				}	
			}
		}

		if ( 1 != outCount) {	
			setEnd( p_kmerInfo[i]);
			if (outCount) {
				uint8_t outEdges = p_kmerInfo[i];
				for (uint8_t j=0; j< 4; ++j,outEdges >>= 1) {
					if (outEdges&1) {
						uint64_t key =((p_kmerValue[i]&mask)<<2)|j;//j is uint8_t,A mistakes 
						uint64_t loc = binSearch(p_kmerValue, counter[key>>move], counter[(key>>move)+1] - 1, key);
						setStart(p_kmerInfo[loc]);
					} 
				}	
			}
				
		}
			
	}

	return NORMAL_EXIT;

}
 
int handleFrstLastKmer(uint64_t  *p_kmerValue, uint16_t *p_kmerInfo, vector<uint64_t> &h, vector<uint64_t> &t)
{
	uint64_t hSize = h.size();
	//uint64_t p_kmerValue_size = p_kmerValue.size() - 1;
	
	uint16_t move = (_kmer - PREINDEXLEN) << 1;
	
	for (uint64_t i=0;i<hSize;++i) {
	 	//uint64_t loc = binSearch(p_kmerValue,0,p_kmerValue_size-1,h[i]);
		setStart(p_kmerInfo[h[i]]);
		uint8_t inCount = inEdgesNum(p_kmerInfo[h[i]]);
		if (inCount) {
			uint8_t inEdges = p_kmerInfo[h[i]] >> 4;
			for (uint8_t j=0; j< 4; ++j, inEdges >>= 1) {
				if (inEdges&1) {
					uint64_t key = (p_kmerValue[h[i]] >> 2)|((uint64_t)j<<((_kmer-1)<<1));
					uint64_t loca = binSearch(p_kmerValue,counter[key>>move], counter[(key>>move)+1] - 1, key );
					setEnd(p_kmerInfo[loca]);
				} 
				
			}
		}
	}
	
	uint64_t tSize = t.size();
	uint64_t mask = ~((uint64_t)0x3 <<((_kmer-1)<<1));

	for (uint64_t i=0;i<tSize;++i) {
	 	//uint64_t loc = binSearch(p_kmerValue, counter[],p_kmerValue_size-1,t[i]);
		setEnd(p_kmerInfo[t[i]]);
		uint8_t outCount = outEdgesNum(p_kmerInfo[t[i]]);
		if (outCount) {
			uint8_t outEdges = p_kmerInfo[t[i]];
			for (uint8_t j=0; j< 4; ++j,outEdges >>= 1) {
				if (outEdges&1) {
					uint64_t key = (p_kmerValue[t[i]] & mask) << 2 | j;
					uint64_t loca = binSearch(p_kmerValue,counter[key>>move], counter[(key>>move) + 1] - 1,key );
					setStart(p_kmerInfo[loca]);
				} 
				
			}
		}
	}

	return 0;
}

uint64_t countEndNum(uint16_t* p_kmerInfo,uint64_t p_kmerNum)
{
	uint64_t _count = 0;
	for (uint64_t i=0; i<p_kmerNum; ++i) {
		if(isEnd(p_kmerInfo[i])) ++_count; 
	//	fprintf(stderr,"%x\n",p_kmerInfo[i]);
	}
	return _count;
}	

int setLabel(uint64_t *kmerValue, uint16_t *kmerInfo, uint64_t kmerNum, vector<uint64_t>& p_heads, vector<uint64_t>& p_tails, uint64_t& countEndN)//start end start&end transit
{
	cutOffMulEdges(kmerValue, kmerInfo, kmerNum);			
	fprintf(stderr,"stop cutting off\n");
	handleFrstLastKmer(kmerValue, kmerInfo, p_heads, p_tails);	
	countEndN = countEndNum(kmerInfo, kmerNum);	
	return NORMAL_EXIT;

}	



bool isStart(uint16_t info) 
{
	return (( info >> 9 ) & 1 );

}

bool isEnd(uint16_t info) 
{
	return ((info >> 8) & 1);

}

int genSpKmers(uint64_t kvalue, uint32_t p_tid, kmersSpchar* p_2k, uint8_t* p_2k_0p, uint64_t& s_ind)
{	
	uint8_t	move = (_kmer - 1) << 1;

	uint64_t mask = (uint64_t)-1 >> (64 -( ((_kmer-1))<<1));
	
	//uint64_t mask = (0x11)<<((_kmer-1)<<1);
	//uint8_t lastChar =  last>>((_kmer-1) << 1);
	
	kmersSpchar assignVar;

	for (uint8_t i = move; i>0;) {
	
		uint8_t lastChar = (kvalue >> i) & 0x3;//there if 
		
		assignVar.infor  = lastChar;
		
		assignVar.value = kvalue & mask;
		
		assignVar.infor |= (i >> 1) << 3;
		
		assignVar.taxID = p_tid;

		mask >>= 2;
		i -= 2;
		p_2k[s_ind++] = assignVar;

	} 
	
	
	uint8_t lastChar = kvalue & 0x3;
	
	//cout<<s_ind<<endl;
	p_2k_0p[s_ind/(_kmer-1) - 1] = lastChar;


	return NORMAL_EXIT;
	//2k.push_back();

}


int output(uint64_t* p_kmerValue, uint64_t kmerNum, uint16_t* p_kmerInfo, kmersSpchar* p_2kmers, uint8_t* p_2kmers_0p, uint32_t* p_taxonIDTab, vector<uint32_t>& nkmerTID )//initate last char and produce unsorted 2kmers
{
	
	uint64_t p_kmerValue_size = kmerNum;

	uint16_t lastChar = 0x5 << 10;
	
	uint16_t move = (_kmer - PREINDEXLEN) << 1;

	uint64_t mask = ~((uint64_t)0x3 <<((_kmer-1)<<1));
	
	//bool isPair = false;
	
	//uint64_t lastKmerValue;

	uint8_t bitMove = ((_kmer-1)<<1);
	
	uint32_t uid = 0;	
#ifdef OUTPUT_STRING_PATH

	char Chars[] = {'A','C','G','T'};
	
	uint64_t totalLen = 0;
	string path;
#endif
	//nkmerTID.push_back(0);
	uint32_t ttid;
	bool isprint;
	uint64_t _ind = 0;

	for (uint64_t i=0; i < p_kmerValue_size; ++i) {
		if (isStart(p_kmerInfo[i])) {
			//path = "";
			//preValue = p_kmerValue[i];
			isprint = false;

			p_kmerInfo[i] |= lastChar;
			uint64_t loc = i;
			
			ttid = p_taxonIDTab[loc];

			nkmerTID.push_back(p_taxonIDTab[loc]);
#ifdef OUTPUT_STRING_PATH
			path = transIntoChars(p_kmerValue[i], _kmer);
			totalLen += _kmer;
			//cout<<loc<<"->";
#endif
			while(!isEnd(p_kmerInfo[loc])) {
				
				lastChar = (p_kmerValue[loc] >>	bitMove) << 10;
				//fprintf(stderr,"%x\t",p_kmerInfo[loc]&0xF);	
				uint8_t j = 0; while (!((p_kmerInfo[loc]>>j)&1)) ++j;
				
				//uint64_t key = ((p_kmerValue[loc]&mask)<<2)|(p_kmerInfo[loc]>>10);
				
				uint64_t key = ((p_kmerValue[loc]&mask)<<2)|j;
				loc = binSearch(p_kmerValue, counter[key>>move], counter[(key>>move)+1]-1, key);
#ifdef OUTPUT_STRING_PATH

				if (p_taxonIDTab[loc]!=ttid) {
					isprint = true;
					path += "1";
				}
#endif
				p_kmerInfo[loc] |= lastChar;

				//p_taxonID[loc] = uid;
#ifdef OUTPUT_STRING_PATH
				path += Chars[j];
				++totalLen;
			//	cout<<loc<<"->";
#endif
			}
			
			
#ifdef OUTPUT_STRING_PATH
			//cout<<endl;
			//if (isprint)
				cout<<uid<<"\t"<<path<<endl;
#endif
			//lastKmerValue = p_kmerValue[loc];
			genSpKmers(p_kmerValue[loc], p_taxonIDTab[loc],  p_2kmers, p_2kmers_0p, _ind);	//#1 < #2 .... < $; they are all less than all AGCT
			
			lastChar = 0x4<<10;

			++uid;
			
		}
	
	}
#ifdef OUTPUT_STRING_PATH
	cout<<uid<<"\t"<<totalLen<<endl;
#endif
		
	return NORMAL_EXIT;
}

int preprocess(const char* refPath, const char *kmerPath, const char *taxonomyNodesPath, const char *giTaxidPath, string& bwt_s, vector<uint32_t>& nkmerTID, uint64_t *hash_index)
{
	


	
	uint64_t kmerNum;


	FILE *fp = fopen(kmerPath, "rb");
	
	fread(&kmerNum, sizeof(uint64_t), 1 , fp);
	
	fprintf(stderr,"%lu\n",kmerNum);	
	uint64_t *kmersValue = new uint64_t[kmerNum];
	
	uint16_t *kmersInfo = new uint16_t[kmerNum]();

	uint32_t *kmersTID = new uint32_t[kmerNum]();

	if (!kmersValue| !kmersInfo | !kmersTID) return MEM_ALLOCATE_ERROR;

	fread(kmersValue, sizeof(uint64_t), kmerNum, fp);
	
	
	
	
	vector<uint64_t> heads;
	vector<uint64_t> tails;

	vector<gid_taxid> taxonIDTab;
		
	fprintf(stderr,"taxonomy tree built....\n");
	taxonTree(taxonomyNodesPath);

	//fprintf(stderr, "scan gid and tid\n");
	//taxonID(giTaxidPath, taxonIDTab);	

	fprintf(stderr,"building debruijin\n");

	build_deb(refPath, kmersValue, kmersInfo, kmersTID, kmerNum,  heads, tails);
	
	//taxonIDTab.clear();	
	
	fprintf(stderr,"set labels...\n");

	//outputDeb(kmersValue, kmersInfo, kmerNum);
	uint64_t endNodeNum;
	setLabel(kmersValue, kmersInfo, kmerNum, heads, tails, endNodeNum);
	//sort2kmers((void *)p_2kmers);
	heads.clear();
	tails.clear();
	//getchar();
	fprintf(stderr, "end node is %lu\n",endNodeNum);

	kmersSpchar* _2kmers = new kmersSpchar[endNodeNum*(_kmer - 1)];

	uint8_t* _2kmers_0p = new uint8_t[endNodeNum];

	fprintf(stderr,"outputing... \n");
	output(kmersValue, kmerNum, kmersInfo, _2kmers, _2kmers_0p, kmersTID, nkmerTID);
	
	fprintf(stderr,"sorting.... \n");
	sortKmers(_2kmers, _2kmers + endNodeNum*(_kmer - 1));
	
	//for (uint64_t i=0; i <endNodeNum*(_kmer-1);++i) cout<<transIntoChars(_2kmers[i].value, _2kmers[i].infor>>3)<<endl;	
	char Chars[] = {'A','C','G','T','#','$'};
	fprintf(stderr,"after soriting ...\n");
	bwt_s = "";

	fprintf(stderr,"initiate 0 positions\n");
	
	for (uint64_t i=0; i<endNodeNum; ++i) 	bwt_s += Chars[_2kmers_0p[i] & 0x7];
	fprintf(stderr,"merging... \n");
	
	mergeSort(kmersValue, kmerNum,kmersInfo, kmersTID, nkmerTID,  _2kmers,endNodeNum*(_kmer-1),  bwt_s, hash_index, endNodeNum-1);
	
	if (kmersInfo) delete kmersInfo;
      	if (kmersValue) delete kmersValue;
	if (kmersTID) delete kmersTID;
	//cout<<bwt_s<<endl;
/*	
	FILE *bwt_string_fp = fopen("bwt_string","w");
	
	fwrite(bwt_s.c_str(), sizeof(char), bwt_s.size(), bwt_string_fp);	
	
	fclose(bwt_string_fp);
*/
	return NORMAL_EXIT;

}


	
