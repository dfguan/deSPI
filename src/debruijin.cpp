// =====================================================================================
//
//       Filename:  debruijin.cpp
//
//    Description:  construct debruijin graph and output sorted kmers without specail character and kmers that involved specail char not sorted
//
//        Version:  1.0
//        Created:  10/29/2015 08:38:41 AM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Dengfeng Guan, <dfguan@hit.edu.cn>
//   Organization:  Harbin Institute of Technology
//
// =====================================================================================

#include "debruijin.hpp"

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



int outputDeb(vector<uint64_t> &p_kmerValue, vector<uint16_t> &p_kmerInfo) 
{
	uint64_t p_size = p_kmerValue.size() - 1;
	char Chars[] = {'A','C','G','T'};
	cout<<"debruijn:" <<endl;
	for(uint64_t i=0; i<p_size;++i) {
		cout<<transIntoChars(p_kmerValue[i])<<"\t";
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

int build_deb(char *path, char *refPath, vector<uint64_t> &kmerValue, vector<uint16_t> &kmerInfo, vector<uint64_t>& p_heads, vector<uint64_t>&  p_tails) //the path indicate a file containing all sorted k+1mers  
{
	//uint64_t prevValue = -1;

	FILE *fp = fopen(path, "rb");
	
	uint64_t *values = new uint64_t[0x2710];//store 10000 value;
	
	//uint64_t useless;
	uint64_t counterSize = (1<<(PREINDEXLEN<<1)) + 1;

	counter = new uint64_t[counterSize]();
	


	uint16_t move = (_kmer - PREINDEXLEN) << 1;


	//uint64_t nodeCounter = 0;//index

	// initiate out edges
	while(!feof(fp)) {
		uint32_t b_size = fread(values,sizeof(uint64_t),0x2710,fp);
		//fprintf(stderr,"%lu\t%lu\n",value[0],useless);
		for(uint32_t i=0; i<b_size;++i) {
			++counter[values[i]>>move]; 
			kmerValue.push_back(values[i]);
			kmerInfo.push_back(0);
		
		}
	} 
 
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
	


	fclose(fp);
	//initiate ranges 
	uint64_t sum = 0;
	uint64_t temp = 0;
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
	
	while (kseq_read(trunk)>=0) {
		for(uint64_t i=0; i<trunk->seq.l;++i) {
			if (Bit[trunk->seq.s[i]] != 4) {
				uint64_t start = i;
				while(Bit[trunk->seq.s[++i]]!=  4 && i < trunk->seq.l);
				//uint64_t end = i - 1;
				if (i - start > _kmer) {
					//intiate head 
					uint64_t key = transIntoBits(trunk->seq.s+start, _kmer);
					//print(trunk->seq.s+start);
					//fprintf(stderr, "\t%lu\t%lu\n",key,counter[key>>move]);
					
					uint64_t loc = binSearch(kmerValue, counter[key>>move],counter[(key>>move)+1]-1, key);
					
					kmerInfo[loc] = setOutEdge((uint8_t) kmerInfo[loc], Bit[trunk->seq.s[start + _kmer]]);
					
					p_heads.push_back(loc);
					
					for (uint64_t j = start + 1; j < i - _kmer; ++j) {
						key = ((key & mask) << 2)| Bit[trunk->seq.s[j+_kmer-1]];
						//print(trunk->seq.s+j);		
						loc = binSearch(kmerValue, counter[key>>move],counter[(key>>move)+1]-1, key);
						kmerInfo[loc] = setInEdge((uint8_t) kmerInfo[loc], Bit[trunk->seq.s[j-1]]);
						kmerInfo[loc] = setOutEdge((uint8_t) kmerInfo[loc], Bit[trunk->seq.s[j+_kmer]]);
					} 
					
					key = ((key & mask) << 2)| Bit[trunk->seq.s[i-1]];
				
					loc = binSearch(kmerValue, counter[key>>move],counter[(key>>move)+1]-1, key);
					
					kmerInfo[loc] = setInEdge((uint8_t) kmerInfo[loc], Bit[trunk->seq.s[i - _kmer -1]]);
					
					p_tails.push_back(loc);		
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

uint64_t binSearch(vector<uint64_t> &kmersValue,  uint64_t lowerBound, uint64_t upperBound,uint64_t key)
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

string transIntoChars(uint64_t v) 
{
	string s(_kmer, 0);
	
	char Chars[] = {'A','C','G','T'};
	
	
	for (uint8_t i=_kmer-1;i!=0;--i) {
		s[i] = Chars[v&0x3];
		v = v>>2;
	}
	s[0] = Chars[v&0x3];

	return s;
}

/*
int main(int argc, char *argv[])
{
	
 vector<uint64_t> kmerValue;
 
 vector<uint16_t> kmerInfo;//lower 8 bit indicates in edges and outedges, upper 8 bits only use
 
 vector<uint32_t> taxonID;//first will be used to store UID
 


 vector<kmersSpchar> _2kmers;

 vector<uint64_t> heads;
 vector<uint64_t> tails;

 vector<uint64_t> fKmer;
 _kmer = (uint8_t) atoi(argv[3]); 
 
 fprintf(stderr,"build debruijin.....\n");
 
 build_deb(argv[1], argv[2], kmerValue, kmerInfo, heads, tails);

 
 fprintf(stderr,"set labels...\n");
 
 setLabel(kmerValue, kmerInfo, heads, tails);
 
 fprintf(stderr,"output \n");

 output(kmerValue, kmerInfo, _2kmers, taxonID, fKmer);
 
 return NORMAL_EXIT;
}

*/

