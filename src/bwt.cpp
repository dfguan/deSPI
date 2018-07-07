// =====================================================================================
//
//       Filename:  bwt.cpp
//
//    Description:  conduct BWT to debruijin
//
//        Version:  1.0
//        Created:  10/20/2015 09:22:52 AM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Dengfeng Guan, <dfguan@hit.edu.cn>
//   Organization:  Harbin Institute of Technology
//
// =====================================================================================
#include "bwt.hpp"
//#define PREINDEXLEN 13

//char *bwt_str;
const uint8_t Bit2[] = {
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
//uint8_t _count[] = {0,0,0,0,0};
//const uint8_t Bit[256];
int print_str(uint8_t *bytes, int len )
{	
	char alphabet[] = {'A','C','G','T'};

	for (int i=0; i<len; ++i) cout<<alphabet[bytes[i]];
	cout<<endl;
	return 0;
}
uint64_t bwt::occ(uint64_t r, uint8_t c)
{	
	uint16_t mask[] = {0xFFFF, 0xFFF0, 0xFF00, 0xF000};//constricted by operating system only for linux
	uint64_t pos = r >> 8;
	
	uint64_t base = *(uint64_t *)(bwt_occ + pos * 168 + (c<<3));
	

	uint64_t start = pos * 168 + 40;
	
	uint64_t extra = r & 0x3;
	
	uint64_t end = (pos + 1)*40 + ((r - extra) >> 1);
	
	uint64_t count = 0;
	//fprintf(stderr,"dddddddd------%x\n", c);
	for (uint64_t i = start; i < end; i += 2) {
		//uint32_t ind = (uint16_t *)
		//fprintf(stderr,"%x\n",*(uint16_t *)(bwt_occ+i));
		count += AGCTCounter[(*(uint16_t *)(bwt_occ+i))*5 + c];//here won't surpass the range
	}
	
	//fprintf(stderr,"count is:%lu\n",count);	
	count += AGCTCounter[(*(uint16_t *)(bwt_occ + end) | mask[extra])*5+c];
	//fprintf(stderr,"count is:%lu\n",count);	
	return base + count;

}

uint64_t bwt::LFC(uint64_t r, uint8_t c) 
{
	//fprintf(stderr,"%hu\n",c);
	//uint64_t ra = rank[c];
	//uint64_t oc = occ(r,c);
	//fprintf(stderr,"%lu\t%lu\n",ra,oc);
	//return ra+oc;
	return rank[c] + occ(r,c);
	

}

uint64_t bwt::transIntoBits(uint8_t *bytes_kmer, uint8_t len)
{
	uint64_t value = 0;
	
	for(uint8_t i=0; i<len; ++i) value = (value<<2)|bytes_kmer[i];			
	
	return value;

}




//give a sp get the unitig until a # is met
int bwt::get_utg(uint64_t sp, string &utg)
{
	uint64_t transInd; // = ((sp>>8)+1)*40 + (sp >> 1);
	char nucl_table[] = {'A','C','G','T'};
	uint8_t c; //= (bwt_occ[transInd] >> ((sp & 0x1) << 2))& 0xF;
	int len_utg = 0;
	while (true) {
		//fprintf(stderr,"%lu\t%x\n",sp,c);
		transInd = ((sp>>8)+1)*40 + (sp >> 1);
		c = (bwt_occ[transInd] >> ((sp & 0x1) << 2))& 0xF;
		if (c == 4) break;//break when c == 4 is met
		utg += nucl_table[c];
		++len_utg;//utg_len could be better ? 
		sp = LFC(sp, c);
	}
	//flip the unitig
	for ( int i = 0; i < (len_utg >> 1); ++i) {
		char c = utg[i];
		utg[i] = utg[len_utg - 1 - i];
		utg[len_utg - 1 - i] = c;
	}
	return len_utg;//return length ?  
}



uint32_t bwt::locate(uint64_t sp)
{
	uint64_t transInd = ((sp>>8)+1)*40 + (sp >> 1);
	//fprintf(stderr,"%lu\n",transInd);
	uint8_t c = (bwt_occ[transInd] >> ((sp & 0x1) << 2))& 0xF;
	//fprintf(stderr,"%x\n",bwt_occ[transInd]);
	//char ts[] = {'A','C','G','T','#'};
	//string bush = "";
	//bush += ts[c];
	while (c != 4 && (sp & 0x1F)) {
		sp = LFC(sp, c);
		//fprintf(stderr,"%lu\t%x\n",sp,c);
		transInd = ((sp>>8)+1)*40 + (sp >> 1);
		c = (bwt_occ[transInd] >> ((sp & 0x1) << 2))& 0xF;
		//bush += ts[c];
	}
	
	if (c == 4) {
		//cout<<"11\t"<<bush<<endl;
		return taxonIDTab[occ(sp, c)];
	} else {
		//cout<<"12\t"<<bush<<"\t"<<sp<<endl;
		return taxonIDTab[(sp>>5) - ((rank[0]+31)>>5) + rank[0]];
	}


}

uint32_t bwt::locate(uint64_t sp, uint8_t *bytes, char *qual, int loc, int &extend)
{
	uint64_t transInd = ((sp>>8)+1)*40 + (sp >> 1);
	uint8_t c = (bwt_occ[transInd] >> ((sp & 0x1) << 2))& 0xF;
//	char ts[] = {'A','C','G','T','#'};
//	string bush = "";
//	bush += ts[c];
	int mismatch = 0;
	while (c != 4 && (sp & 0x1F)) {
		sp = LFC(sp, c);
		transInd = ((sp>>8)+1)*40 + (sp >> 1);
		c = (bwt_occ[transInd] >> ((sp & 0x1) << 2))& 0xF;
		--loc;
		if (loc > -1) {
			if (c!= bytes[loc] && qual[loc] > 20 ) ++mismatch;
		       	if (mismatch < 2) ++extend;	
		}			
//		bush += ts[c];
	}
	
	if (c == 4) {
//		cout<<"11\t"<<bush<<endl;
		return taxonIDTab[occ(sp, c)];
	} else {
	//	cout<<"12\t"<<bush<<"\t"<<sp<<endl;
		return taxonIDTab[(sp>>5) - ((rank[0]+31)>>5) + rank[0]];
	}


}




int bwt::exactMatch(uint8_t *bytes, char *qual, int len, int& match_len, uint32_t* assignedTID) 
{
	uint64_t sp,ep;	
	//int match_len;
	//mismatch = 0;
	//fwrite(str,sizeof(char),len,stdout);
	//cout<<endl;
	uint64_t prefixValue = transIntoBits(bytes + len - PREINDEXLEN, PREINDEXLEN ) << 1;
	sp = hash_index[prefixValue];
	ep = hash_index[prefixValue + 1];
	int i = len - PREINDEXLEN - 1;
	uint8_t c ;//= str[len-1];
	//sp = rank[Bit2[c]];
	//ep = rank[Bit2[c]+1];
	//int i = len - 2;
	//cout<<len<<"---------------------"<<endl;
	//cout<<sp<<"\t"<<ep<<endl;
	uint64_t sp_pre = sp;
	uint64_t ep_pre = ep;
	//uint64_t ep_pre = ep;
	//char ts[] = {'A','C','G','T','#'};
	//string li = "";	
	//fprintf(stderr,"%lu\t%lu\n", sp,ep);	
		//uint64_t transInd = ((sp>>8)+1)*40 + (sp >> 1);
		//uint8_t ck = (bwt_occ[transInd] >> ((sp & 0x1) << 2))& 0xF;
		//li += ts[c];
	while (sp < ep && i >= 0) {
		c = bytes[i];
		
		sp_pre = sp;
		ep_pre = ep;
		//ep_pre = ep;

		sp = LFC(sp, c);
		ep = LFC(ep, c);
	
		//transInd = ((sp>>8)+1)*40 + (sp >> 1);
		//ck = (bwt_occ[transInd] >> ((sp & 0x1) << 2))& 0xF;
		//li += ts[ck];	
		--i;

		//fprintf(stderr,"%lu\t%lu\n", sp,ep);	
		//cout<<sp<<"\t"<<ep<<endl;
	}
	//cout<<endl;
	if (sp < ep) {
		match_len = len;
		i = i + 1;
		//cout<<sp<<"\t"<<ep<<endl;
	
	} else {
		
		match_len = len - i - 2;
	       	sp = sp_pre;
		ep = ep_pre;
		//cout<<match_len<<"\t"<<sp<<"\t2"<<ep<<endl;
		i = i + 2;
	}

	int counterTID = 0;
	int extend = 0;
	int max_extend = -1;
	uint32_t temp_tid;
	if (match_len > threshold) {
		
		while (sp < ep) {
			temp_tid = locate(sp, bytes, qual, i, extend);
			if (extend > max_extend) {
				counterTID = 0;
				assignedTID[counterTID++] = temp_tid;
				max_extend = extend;
			} else 
				if (extend == max_extend) 
					assignedTID[counterTID++] = temp_tid;
			
			++sp;
		}
	//cout<<len - match_len - 1<<"\t"<<match_len<<"\t";
	//cout<<li<<"\t";
	//fwrite(str+ i,sizeof(char), match_len, stdout);	
	//cout<<"\t"<<assignedTID[0]<<endl;

	
	}
	match_len += max_extend;

	return counterTID;
		//return false;
	/*
	int extra_step = 0;
	clock_t t = clock();	
	if (sp < ep) {
		match_len = len;
		c = bwt_str[sp];
		while (Bit2[c]!=4) {
			++extra_step;
			sp = LFC(sp, c);
			c = bwt_str[sp];
		}
		sp = occ(sp, c);
	time_consume = clock() - t;
	cout<<"extrac:\t"<<extra_step<<endl;
		return true;
		
	} else {
		if (sp_pre == ep_pre - 1) {
			i += 2;
			c = bwt_str[sp_pre];
			while (i >= 0 && Bit2[c] != 4) {
				
				if (c != str[i]) ++mismatch;
				//if (mismatch > 3) break;
				cout<<c<<"\t"<<str[i]<<endl;
				sp_pre = LFC(sp_pre, c);
				c = bwt_str[sp_pre];
				-- i;
			}
			if (i>=0) {
				match_len = len - i - 1;
			} else match_len = len; 
			c = bwt_str[sp_pre];
			while (Bit2[c]!=4) {
				++extra_step; 
				sp_pre = LFC(sp_pre, c);
				c = bwt_str[sp_pre];
				//cout<<sp_pre<<endl;
			} 
			
			//cout<<sp_pre<<"\t"<<c<<endl;	
			sp = occ(sp_pre, c);
	time_consume = clock() - t;
	cout<<"extrac:\t"<<extra_step<<endl;
			return true;
		} else  {

	time_consume = clock() - t;
			return false;	
		}
	
	} 
	*/
	//cout<<time_consume<<endl;

}
int bwt::rmdup(uint32_t *tempTids, int counter, uint32_t* assignedTID)
{
	int counter_assignedTID = 1;
	
	for(int i=1; i<counter; ++i) {
		int j;
		for(j=0; j<counter_assignedTID; ++j) {
			if (tempTids[i] == assignedTID[j]) break;
		}
		if (j == counter_assignedTID) assignedTID[counter_assignedTID++] =  tempTids[i];
	}
	return counter_assignedTID;

}


int bwt::exactMatch(uint8_t *bytes, int len, int& match_len, uint32_t* assignedTID) 
{
	uint64_t sp,ep;	
	//int match_len;
	//mismatch = 0;
	//fwrite(str,sizeof(char),len,stdout);
	//cout<<endl;
	uint32_t tempTids[1024];
	uint64_t prefixValue = transIntoBits(bytes + len - PREINDEXLEN, PREINDEXLEN ) << 1;
	sp = hash_index[prefixValue];
	ep = hash_index[prefixValue + 1] ;
	int i = len - PREINDEXLEN - 1;
	uint8_t c ;//= str[len-1];
	//sp = rank[Bit2[c]];
	//ep = rank[Bit2[c]+1];
	//int i = len - 2;
	//cout<<len<<"---------------------"<<endl;
	//cout<<sp<<"\t"<<ep<<endl;
	uint64_t sp_pre = sp;
	uint64_t ep_pre = ep;
	//uint64_t ep_pre = ep;
	//char ts[] = {'A','C','G','T','#'};
	//string li = "";	
	//fprintf(stderr,"%lu\t%lu\n", sp,ep);	
		//uint64_t transInd = ((sp>>8)+1)*40 + (sp >> 1);
		//uint8_t ck = (bwt_occ[transInd] >> ((sp & 0x1) << 2))& 0xF;
		//li += ts[c];
	if (sp >= ep) {
		i = len - 3;
	}
	while (sp < ep && i >= 0) {
		c = bytes[i];
		
		sp_pre = sp;
		ep_pre = ep;
		//ep_pre = ep;

		sp = LFC(sp, c);
		ep = LFC(ep, c);
	
		//transInd = ((sp>>8)+1)*40 + (sp >> 1);
		//ck = (bwt_occ[transInd] >> ((sp & 0x1) << 2))& 0xF;
		//li += ts[ck];	
		--i;

		//fprintf(stderr,"%lu\t%lu\n", sp,ep);	
		//cout<<sp<<"\t"<<ep<<endl;
	}
	//cout<<endl;
	if (sp < ep) {
		match_len = len;
		i = i + 1;
		//cout<<sp<<"\t"<<ep<<endl;
	
	} else {
		
		match_len = len - i - 2;
	       	sp = sp_pre;
		ep = ep_pre;
		//cout<<match_len<<"\t"<<sp<<"\t2"<<ep<<endl;
		i = i + 2;
	}
	int counterTID = 0;

	//cout<<i<<"\t"<<match_len<<"\t";
	//print_str(bytes + i,match_len);
	int count_times = ep - sp;
	if (match_len > threshold  && count_times < 1024) {
		//fprintf(stderr,"sdddddd");	
		//print_str(bytes + i, match_len);
		//fprintf(stderr,"%lu,%lu\n",sp,ep);
		while (sp < ep) {
			tempTids[counterTID++] = locate(sp);
			++sp;
			//cout<<"\t"<<assignedTID[counterTID - 1];
		}
		//cout<<endl;
		
		//cout<<len - match_len - 1<<"\t"<<match_len<<"\t";
		//cout<<li<<"\t";
		//fwrite(str+ i,sizeof(char), match_len, stdout);	
		//cout<<"\t"<<assignedTID[0]<<endl;
		
		assignedTID[0] = tempTids[0];
	
	}
	if (counterTID > 1) 
		return rmdup(tempTids, counterTID, assignedTID);
	else 
	       return counterTID;

		//return false;
	/*
	int extra_step = 0;
	clock_t t = clock();	
	if (sp < ep) {
		match_len = len;
		c = bwt_str[sp];
		while (Bit2[c]!=4) {
			++extra_step;
			sp = LFC(sp, c);
			c = bwt_str[sp];
		}
		sp = occ(sp, c);
	time_consume = clock() - t;
	cout<<"extrac:\t"<<extra_step<<endl;
		return true;
		
	} else {
		if (sp_pre == ep_pre - 1) {
			i += 2;
			c = bwt_str[sp_pre];
			while (i >= 0 && Bit2[c] != 4) {
				
				if (c != str[i]) ++mismatch;
				//if (mismatch > 3) break;
				cout<<c<<"\t"<<str[i]<<endl;
				sp_pre = LFC(sp_pre, c);
				c = bwt_str[sp_pre];
				-- i;
			}
			if (i>=0) {
				match_len = len - i - 1;
			} else match_len = len; 
			c = bwt_str[sp_pre];
			while (Bit2[c]!=4) {
				++extra_step; 
				sp_pre = LFC(sp_pre, c);
				c = bwt_str[sp_pre];
				//cout<<sp_pre<<endl;
			} 
			
			//cout<<sp_pre<<"\t"<<c<<endl;	
			sp = occ(sp_pre, c);
	time_consume = clock() - t;
	cout<<"extrac:\t"<<extra_step<<endl;
			return true;
		} else  {

	time_consume = clock() - t;
			return false;	
		}
	
	} 
	*/
	//cout<<time_consume<<endl;

}
//load bwt and taxonID
int bwt::load_index(const char *dirPath)
{
	//fprintf(stderr,"%s",dirPath);	
	string bwtFilePath(dirPath);
	string tidFIlePath(dirPath);
	if (bwtFilePath[bwtFilePath.length()-1] == '/') {
		bwtFilePath += "despi.bwt";
		tidFIlePath += "despi.tid";
	
	} else {
		bwtFilePath += "/despi.bwt";
		tidFIlePath += "/despi.tid";

	}
	
	FILE *bwt_fp = fopen(bwtFilePath.c_str(), "rb");
	FILE *tid_fp = fopen(tidFIlePath.c_str(), "rb");

	if (NULL == bwt_fp || NULL == tid_fp)
		return FILE_OPEN_ERROR;
	
	
	//uint64_t bwt_len; 
      	
	fread(&len_bwt_occ, 8, 1, bwt_fp);

	uint8_t* temp_bwt_occ  = new uint8_t [len_bwt_occ];
	
	fread(temp_bwt_occ, sizeof(uint8_t), len_bwt_occ, bwt_fp);
	
	bwt_occ = temp_bwt_occ;	
/*
	for (uint64_t i=0; i < len_bwt_occ; i += 168) {
		for (uint8_t z=0; z<5;++z) 
		fprintf(stderr, "%lu\t%lu\n",i + z*8, *(uint64_t *)(bwt_occ+i + z*8));
	
	
	}
*/
	/*	
	fread(&occCount, 8, 1, bwt_fp);

	occCheck = new uint64_t[occCount];
	

	fread(occCheck, sizeof(uint64_t), occCount, bwt_fp);
	*/
	fread(rank, sizeof(uint64_t), 5, bwt_fp);


	uint64_t hash_index_size = (uint64_t)1 <<((PREINDEXLEN<<1) + 1);
	
	hash_index = new uint64_t[hash_index_size]();

	fread(hash_index,sizeof(uint64_t),hash_index_size,bwt_fp);
	
	uint64_t AGCTCounterSize;

	fread(&AGCTCounterSize, sizeof(uint64_t), 1, bwt_fp);

	AGCTCounter = new uint8_t[AGCTCounterSize];

	fread(AGCTCounter, sizeof(uint8_t), AGCTCounterSize, bwt_fp);
/*	
	for (uint64_t i=0; i<hash_index_size;++i) {
		hash_index[i] += rank[0];
	
	
	}
*/
	fread(&tidSize, sizeof(uint64_t), 1, tid_fp);

	taxonIDTab = new uint32_t[tidSize];	
	
	fread(taxonIDTab, sizeof(uint32_t), tidSize, tid_fp);

	fclose(bwt_fp);
	fclose(tid_fp);
	return NORMAL_EXIT;
}

uint8_t countZero(uint64_t v)
{
	uint8_t counter = 0;
	for (uint8_t i=0; i<4;++i) {
		if ((v&0xf) == 0) ++counter;
		v >>= 4;
	
	}
	return counter;
}
int bwt::dump_index(const char *dirPath)
{
	
	string bwtFilePath(dirPath);
	string tidFIlePath(dirPath);
	if (bwtFilePath[bwtFilePath.length()-1] == '/') {
		bwtFilePath += "despi.bwt";
		tidFIlePath += "despi.tid";
	
	} else {
		bwtFilePath += "/despi.bwt";
		tidFIlePath += "/despi.tid";

	}
	
	FILE *bwt_fp = fopen(bwtFilePath.c_str(), "wb");
	FILE *tid_fp = fopen(tidFIlePath.c_str(), "wb");

	if (NULL == bwt_fp || NULL == tid_fp)
		return FILE_OPEN_ERROR;
	
	
	//here combine occ and bwt_str
	uint32_t bufferSize = 168 << 8;
      	
	uint8_t *buffer = new uint8_t[bufferSize + 5]; // ensure that buffer is long enough acutally it won't surpass buffersize + 2

	if (buffer == NULL) {
		fprintf(stderr, "Fail to allocate space, exit \n");
		exit(1);
	}	
	//uint64_t theEnd = (len_bwt_str >> 8 );
	fprintf(stderr,"start to dump bwt idx to disk\n");	
	uint64_t occCount = ((len_bwt_str + 255)>>8)*5;
	uint64_t byteLen = ((len_bwt_str + 1) >> 1) + 2 + (occCount << 3);
	
	fwrite(&byteLen, 8, 1, bwt_fp);

	uint64_t leftOneChar = len_bwt_str & 0x1;
	
	uint32_t bufferPoint = 0;
	
	//uint8_t fillOcc = 0;
	
	uint64_t occCheck[] = {0,0,0,0,0};	
	for(uint64_t i=0; i< len_bwt_str - leftOneChar; i += 2) {

		if (0 == (i&0xFF)) {
			//fprintf(stderr,"%u\n",bufferPoint);
			for(uint8_t j=0; j< 5; ++j) {
				memcpy((uint64_t *)(buffer + bufferPoint),occCheck + j, sizeof(uint64_t));
				bufferPoint += 8;	
			}
		}
		//each character takes 4 bits
		++occCheck[Bit2[bwt_str[i]]];
		++occCheck[Bit2[bwt_str[i+1]]];
	
		buffer[bufferPoint++] = (Bit2[bwt_str[i+1]]<<4)|(Bit2[bwt_str[i]]);

		if (bufferPoint >= bufferSize) {
			fwrite(buffer, sizeof(uint8_t), bufferSize, bwt_fp);
			bufferPoint = 0;
		}  
	}
	
	if (leftOneChar) {
		buffer[bufferPoint++] = 0xF0|Bit2[bwt_str[len_bwt_str-1]];
		++occCheck[Bit2[bwt_str[len_bwt_str-1]]];
	} 
	
	buffer[bufferPoint++] = 0xFF;
	buffer[bufferPoint++] = 0xFF;

	fwrite(buffer,sizeof(uint8_t), bufferPoint, bwt_fp);
	/*
      	fwrite(&len_bwt_str, 8, 1, bwt_fp);

	//p_bwt_s  = new char [bwt_len];
	fprintf(stderr,"1\n");

	fwrite(bwt_str, sizeof(char), len_bwt_str, bwt_fp);
	
	fwrite(&occCount, 8, 1, bwt_fp);

	fwrite(occCheck, sizeof(uint64_t),occCount,bwt_fp);
	*/
	rank[0] = occCheck[4];
	for (uint8_t i=1; i<5; ++i) rank[i] = rank[i-1] + occCheck[i-1];
	fwrite(rank, sizeof(uint64_t),5,bwt_fp);	
	//*table_len = bwt_len;
	
	uint64_t hash_index_size = (uint64_t)1 <<((PREINDEXLEN<<1) + 1);
	
	fwrite(hash_index,sizeof(uint64_t),hash_index_size,bwt_fp);

	//p_nkmerTID = new uint32_t[bwt_len];
	uint64_t AGCTCounterSize = (1 << 16) * 5;
	fwrite(&AGCTCounterSize, 8, 1, bwt_fp);
	//fprintf(stderr,"%lutexst1\n",AGCTCounterSize);
	//if (AGCTCounter != NULL) fprintf(stderr,"haha\n");
	//fprintf(stderr,"texst33\n");
	uint64_t mask[] = {0, 0x1111, 0x2222, 0x3333, 0x4444};
	AGCTCounter = new uint8_t[AGCTCounterSize];
	if (AGCTCounter == NULL) {
		uint8_t ACGT_Array[] = {0,0,0,0,0};
		for (uint64_t i=0; i <= 0xFFFF; ++i) {
			for (uint8_t j=0; j < 5; ++j) {
				ACGT_Array[j] = countZero(i^mask[j]);
			}	
			fwrite(ACGT_Array, sizeof(uint8_t), 5, bwt_fp);
		}
	} else {
		for (uint64_t i=0; i <= 0xFFFF; ++i) {
			for (uint8_t j=0; j < 5; ++j) AGCTCounter[i*5 + j] = countZero(i^mask[j]);	
		}
		fwrite(AGCTCounter, sizeof(uint8_t), AGCTCounterSize, bwt_fp);
	} 
	
	fprintf(stderr,"finish dumping index to disk\n");
	fprintf(stderr,"start dumping SA to disk\n");
	fwrite(&tidSize, sizeof(uint64_t), 1, tid_fp );

	fwrite(taxonIDTab, sizeof(uint32_t), tidSize, tid_fp);
	
	fprintf(stderr,"finish dumping SA to disk\n");

	fclose(bwt_fp);
	fclose(tid_fp);
	
	return NORMAL_EXIT;
}


//int bwt::dump_index(const char *dirPath, vector<uint32_t>& p_nkmerTID)
//{
	
	//string bwtFilePath(dirPath);
	//string tidFIlePath(dirPath);
	//if (bwtFilePath[bwtFilePath.length()-1] == '/') {
		//bwtFilePath += "despi.bwt";
		//tidFIlePath += "despi.tid";
	
	//} else {
		//bwtFilePath += "/despi.bwt";
		//tidFIlePath += "/despi.tid";

	//}
	
	//FILE *bwt_fp = fopen(bwtFilePath.c_str(), "wb");
	//FILE *tid_fp = fopen(tidFIlePath.c_str(), "wb");

	//if (NULL == bwt_fp || NULL == tid_fp)
		//return FILE_OPEN_ERROR;
	
	
	//here combine occ and bwt_str
	//uint32_t bufferSize = 168 << 8;
          
	//uint8_t *buffer = new uint8_t[bufferSize + 5]; // just ensure that buffer is long enough acutally it won't surpass buffersize + 2

	
	//uint64_t theEnd = (len_bwt_str >> 8 );
	//fprintf(stderr,"eere\n");	
	//uint64_t byteLen = ((len_bwt_str + 1) >> 1) + 2 + (occCount << 3);
	
	//fwrite(&byteLen, 8, 1, bwt_fp);

	//uint64_t leftOneChar = len_bwt_str & 0x1;
	
	//fprintf(stderr,"eere1\n");	
	//uint32_t bufferPoint = 0;
	
	//uint8_t fillOcc = 0;
	
	//uint64_t occPoint = 0;

	//for(uint64_t i=0; i< len_bwt_str - leftOneChar; i += 2) {
		//if (0 == (i&0xFF)) {
			//fprintf(stderr,"%u\n",bufferPoint);
			//for(uint8_t j=0; j< 5; ++j) {
				//memcpy((uint64_t *)(buffer + bufferPoint),occCheck+occPoint*5 + j, sizeof(uint64_t));
				//bufferPoint += 8;	
			//}
			//++occPoint;
		//}
		//buffer[bufferPoint++] = (Bit2[bwt_str[i+1]]<<4)|(Bit2[bwt_str[i]]);

		//if (bufferPoint >= bufferSize) {
			//fprintf(stderr,"%lu\n",i);
			//fwrite(buffer, sizeof(uint8_t), bufferSize, bwt_fp);
			//bufferPoint = 0;
		//}  
	//}
	
	//if (leftOneChar) {
		//buffer[bufferPoint++] = 0xF0|Bit2[bwt_str[len_bwt_str-1]];
	//} 
	
	//buffer[bufferPoint++] = 0xFF;
	//buffer[bufferPoint++] = 0xFF;

	//fwrite(buffer,sizeof(uint8_t), bufferPoint, bwt_fp);
	//fprintf(stderr,"end dumping bwt...\n");
	
          //fwrite(&len_bwt_str, 8, 1, bwt_fp);

	//p_bwt_s  = new char [bwt_len];
	//fprintf(stderr,"1\n");

	//fwrite(bwt_str, sizeof(char), len_bwt_str, bwt_fp);
	
	//fwrite(&occCount, 8, 1, bwt_fp);

	//fwrite(occCheck, sizeof(uint64_t),occCount,bwt_fp);
	//*/
	//fprintf(stderr, "start to dump SA...\n");	
	//fwrite(rank, sizeof(uint64_t),5,bwt_fp);	
	//*table_len = bwt_len;
	//uint64_t hash_index_size = (uint64_t)1 <<((PREINDEXLEN<<1) + 1);
	
	//fwrite(hash_index,sizeof(uint64_t),hash_index_size,bwt_fp);

	//p_nkmerTID = new uint32_t[bwt_len];
	//uint64_t AGCTCounterSize = (1 << 16) * 5;
	//fprintf(stderr,"%lutexst1\n",AGCTCounterSize);
	//if (AGCTCounter != NULL) fprintf(stderr,"haha\n");
	//AGCTCounter = new uint8_t[AGCTCounterSize];
	//fprintf(stderr,"texst33\n");
	//if (AGCTCounter == NULL) fprintf(stderr, "fail to allocate space\n");

	//fwrite(&AGCTCounterSize, 8, 1, bwt_fp);
	//uint64_t mask[] = {0, 0x1111, 0x2222, 0x3333, 0x4444};

	//for (uint64_t i=0; i <= 0xFFFF; ++i) {
		//for (uint8_t j=0; j < 5; ++j) AGCTCounter[i*5 + j] = countZero(i^mask[j]);	
	//}
	
	//fwrite(AGCTCounter, sizeof(uint8_t), AGCTCounterSize, bwt_fp);
	
	//fprintf(stderr,"texst2\n");
	//tidSize = p_nkmerTID.size();

	//fwrite(&tidSize, sizeof(uint64_t), 1, tid_fp );

	//fwrite(&p_nkmerTID[0], sizeof(uint32_t), tidSize, tid_fp);
	
	//fprintf(stderr,"\n");

	//fclose(bwt_fp);
	//fclose(tid_fp);
	
	//return NORMAL_EXIT;

//}
/*
//int bwt::bwt_init()
{
	//256 is interval
	occCount = ((len_bwt_str + 255)>>8)*5;
	
	occCheck = new uint64_t[occCount];
	
	uint64_t ind = 0;
	
	uint64_t tempOcc[] = {0,0,0,0,0};

	uint8_t kago = 0;
	
	//uint64_t endSpot = (((len_bwt_str+255) >> 8) << 8) + 1;

	for(uint64_t i=0; i < len_bwt_str; ++i) {
		if (kago  == 0) {
			occCheck[(ind*5)] = tempOcc[0];
			occCheck[(ind*5)+1] = tempOcc[1];
			
			occCheck[(ind*5)+2] = tempOcc[2];
			occCheck[(ind*5)+3] = tempOcc[3];	
			occCheck[(ind*5)+4] = tempOcc[4];
			
			++ind;
			//for (int z=0; z < 5; ++z) {
			//	fprintf(stderr, "%lu\n", tempOcc[z]);
			//}	
				
		}
		++tempOcc[Bit2[bwt_str[i]]];
		++kago;	
	}

	//for(uint64_t i=endSpot; i<len_bwt_str; ++i) ++tempOcc[Bit2[bwt_str[i]]]; 

	rank[0] = tempOcc[4];
	for (uint8_t i=1; i<5; ++i) rank[i] = rank[i-1] + tempOcc[i-1];
	//for (uint8_t i=0;i<5;++i) cout<<rank[i]<<endl;
	//rank[5] = rank[4];			
	return NORMAL_EXIT;
}
*/
