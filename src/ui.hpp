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
#include "desc.h"

#define INDEX 1
#define CLASSIFY 2
#define VIEW 3
//should separate different paramaters for different opitons 
//
typedef struct options {
	uint8_t kmer;//both 
	string sortedKmer; //index
	//string gids;
	string ref;//used?
	string output; //index output dir
	
	int		option; // 1 build index 2 classify 3 view;	
	bool isPaired;//classify
	string index_dir;//still used ? maybe not
	uint8_t seed;//classify 
	int intv;//clasisfy 
	int iter;//classify not used 
	int n_thread; //classify 
	vector<string> read_fns; //classify
	string evo_tree_path;
	
	//view
	bool out_all;
	string taxids;
	
	int argc; //all
	char **argv;//all
}opts;


class UI {
public:
	UI(opts *opt);
	int usage();
	int ind_usage();
	int classify_usage();
	int view_usage();
	int opt_parse(int argc, char *argv[], opts* opt);
};
