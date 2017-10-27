
// =====================================================================================
//
//       Filename:  jreader.cpp
//
//    Description:  jellyfish 1.1.11 reader 
//
//        Version:  1.0
//        Created:  12/08/2015 04:52:21 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Dengfeng Guan, <dfguan@hit.edu.cn>
//   Organization:  Harbin Institute of Technology
//
// =====================================================================================
#include "jreader.hpp"

JReader::JReader(string filename_str, string mode, size_t size) 
{
  open_file(filename_str, mode, size);
}

void JReader::open_file(string filename_str, string mode, size_t size) 
{
	const char *filename = filename_str.c_str();
	
	Input_DB_filename = filename_str;
	
	FILE *fp = fopen(filename, mode.c_str());
	if (!fp) {
		fprintf(stderr, "Fail to open file, Please check file path\n");
		exit(1);	
	}
	char data[1024];

	char *DATABASE_FILE_TYPE = "JFLISTDN";
	
	char *fptr = data;
	
	fread(data, sizeof(char), JR_HEADER_DATA_SIZE, fp);	
	
	if (strncmp(fptr, DATABASE_FILE_TYPE, strlen(DATABASE_FILE_TYPE))) {
		fprintf(stderr," not in proper format\n");
		exit(1);
	}
	
	memcpy(&key_bits, fptr + 8, 8);
	memcpy(&val_len, fptr + 16, 8);
	memcpy(&key_ct, fptr + 48, 8);
	//if (val_len != 4) {
		//fprintf(stderr,"can only handle 4 byte DB values\n");
		//exit(1);
	//}
	k = key_bits / 2;
	key_len = key_bits / 8 + !! (key_bits % 8);
	fclose(fp);
}
/*  
void JReader::open_file(string filename_str, string mode, size_t size) 
{
	const char *filename = filename_str.c_str();
	
	Input_DB_filename = filename_str;
	
	int o_flags = O_RDONLY;
	
	fd = open(filename, o_flags, 0666);
	// Second try for R/W if failure was due to non-existence
	if (fd < 0) {
		fprintf(stderr,"can't open file %s\n",filename);
		exit(1);
	}


	struct stat sb;
	if (fstat(fd, &sb) < 0) {

		fprintf(stderr,"can't stat file %s\n",filename);
		exit(1);

	}
	size_t filesize = sb.st_size;

	int m_flags = MAP_PRIVATE;


	fptr = (char *)mmap(0, filesize, PROT_READ | PROT_WRITE, m_flags, fd, 0);
	if (fptr == MAP_FAILED) {
		fprintf(stderr,"can't mmap %s\n",filename);
		exit(1);

	}
	//valid = true;
	//fptr = ptr;
	char *DATABASE_FILE_TYPE = "JFLISTDN";

	if (strncmp(fptr, DATABASE_FILE_TYPE, strlen(DATABASE_FILE_TYPE))) {
		fprintf(stderr," not in proper format\n");
		exit(1);
	
	}
	memcpy(&key_bits, fptr + 8, 8);
	memcpy(&val_len, fptr + 16, 8);
	memcpy(&key_ct, fptr + 48, 8);
	//if (val_len != 4) {
		//fprintf(stderr,"can only handle 4 byte DB values\n");
		//exit(1);
	//}
	k = key_bits / 2;
	key_len = key_bits / 8 + !! (key_bits % 8);
}
*/

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
uint64_t JReader::get_key_ct()
{
	return key_ct;
}

uint64_t JReader::get_val_len()
{
	return val_len;


}
string JReader::get_db_name()
{

	return Input_DB_filename;
}
uint64_t JReader::get_key_len()
{
	return key_len;
}
void JReader::reader(char *data) 
{
  uint64_t pair_size = key_len + val_len;
  char pair[pair_size];

  ifstream input_file(Input_DB_filename.c_str(), std::ifstream::binary);
  input_file.seekg(header_size(), ios_base::beg);

  // Create a copy of the offsets array for use as insertion positions
  for (uint64_t i = 0; i < key_ct; i++) {
    input_file.read(pair, pair_size);
    uint64_t kmer = 0;
    memcpy(&kmer, pair, key_len);
    cout<<transIntoChars(kmer,21)<<endl;
  
  }
  input_file.close();

}

size_t JReader::header_size() 
{ 
	return 72 + 2 * (4 + 8 * key_bits); 
}
/*
int main(int argc, char *argv[])
{
	JReader jr(argv[1]);
	
	jr.reader(NULL);




}
*/
