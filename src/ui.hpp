// =====================================================================================
//
//       Filename:  classify.cpp
//
//    Description:  classify reads
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

//#include "desc.h"
#include <string>
#include <vector>

using namespace std;

#include <stdint.h>
#include <string.h>
#define PACKAGE_NAME "gsw"

typedef struct options {
	bool isClassify;
	uint8_t kmer;
	string sortedKmer;
	string gids;
	string tids;
	string ref;
	string output;

	string lib;
	uint8_t seed;
	int inv;
	int num_threads;
	vector<string> reads;
	int argc;
	char **argv;
}opts;


class UI {
public:
	UI(opts *opt);
	int usage();
	int ind_usage();
	int classify_usage();
	int opt_parse(int argc, char *argv[], opts* opt);
};
