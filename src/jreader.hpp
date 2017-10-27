// =====================================================================================
//
//       Filename:  jreader.hpp
//
//    Description:  header of jreader
//
//        Version:  1.0
//        Created:  12/11/2015 04:42:52 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Dengfeng Guan, <dfguan@hit.edu.cn>
//   Organization:  Harbin Institute of Technology
//
// =====================================================================================

#ifndef _JREADER_HPP
#define _JREADER_HPP
#include <string>
#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;

#include <stdint.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>

#define JR_HEADER_DATA_SIZE 56

class JReader{
	private:
    		char 		*fptr;
		int 		fd;
		uint8_t 	k;
		uint64_t 	key_bits;
		uint64_t 	key_len;
		uint64_t 	val_len;
		uint64_t 	key_ct;
		string 		Input_DB_filename;
	public:

		void reader(char *data) ;
		JReader(string filename_str, string mode="r", size_t size=0); 
		void open_file(string filename_str, string mode, size_t size); 
		size_t header_size();
		uint64_t get_key_len();
		uint64_t get_val_len();
		uint64_t get_key_ct();
		string get_db_name();
};

string transIntoChars(uint64_t v, uint8_t len);
#endif
