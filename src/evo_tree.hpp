// =====================================================================================
//
//       Filename:  evo_tree.hpp
//
//    Description:  evolutionary tree 
//
//        Version:  1.0
//        Created:  04/03/2017 08:41:15 AM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Dengfeng Guan, <dfguan@hit.edu.cn>
//   Organization:  Harbin Institute of Technology
//
// =====================================================================================

#ifndef _EVO_TREE_H
#define _EVO_TREE_H

#include <stdint.h>
#include <string.h>

#include <iostream>
#include <map>
#include <string>
#include <set>

//using namespace std; /dont use this it will lead to troubles

class  evo_tree{
private:
	std::map<uint32_t, uint32_t> tree;
	
public:
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
	void error(const char *msg) {
		::std::cerr << "Error: " << msg << ::std::endl;
		exit(1);
	}
	
	evo_tree(const char *path) {
		FILE *fp = fopen(path,"r");
		
		if (fp) {
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
				tree[tid] = p_tid;
			}
			//end mark
			tree[1] = 0;
		} else {
			error("Fail to open evolutionary tree");	
		}
	
	}
	uint32_t LCA(uint32_t tid1, uint32_t tid2) {
	
		if (tid1 == 0 || tid2 == 0)
		return tid1 ? tid1 : tid2;

		std::set<uint32_t> pool;
		while (tid1 > 0) {
			pool.insert(tid1);
			tid1 = tree[tid1];
		}
		while (tid2 > 0) {
			if (pool.count(tid2) > 0)
				return tid2;
			tid2 = tree[tid2];
		}
		return 1;
	}
};


#endif


