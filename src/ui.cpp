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
#include "ui.hpp"
#include "error.hpp"

#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>



char *const short_options = "m:k:g:a:r:o:s:i:t:h";
struct option long_options[] = {
    { "kmer_len",     1,   NULL,    'm'   },
    { "kmers",     1,   NULL,    'k' },
    { "gids",     1,   NULL,    'g'   },
    { "tids",     1,   NULL,    'a'   },
    { "ref",	1,NULL,	'r'},
    { "output", 1,  NULL, 'o'},
    { "seed_len", 1, NULL, 's'},
    { "inv", 1, NULL, 'i'},
    { "threads",1, NULL,'t'},
    //{"gapextended", 1,  NULL,'e'},
    //{"match",   1,  NULL,'c'},
    { "help",    0,  NULL,'h'},
    { 0,     0,   0,    0   }
};


UI::UI(opts *opt)
{
	opt->seed = 31;
	opt->kmer = 31;
	opt->inv = 0;
	opt->num_threads = 1;

}

int UI::ind_usage()
{
        fprintf(stderr, "Usage:     %s [Options] <SortedKmer> <gidTID> <Taxonomy> <Reference> <IndexDir>\n\n", PACKAGE_NAME); 
        fprintf(stderr, "<SortedKmer>           Sorted kmer in binary file\n");
        fprintf(stderr, "<gidTid>               correspondence between gid and tid\n");
        fprintf(stderr, "<Taxonomy>             taxonomy tree\n");
        fprintf(stderr, "<Reference>             References file, in FASTQ/FASTA format\n");
        fprintf(stderr, "<IndexDir>             The directory storing RHT index\n");
        
        fprintf(stderr, "Options:   -k, --kmer       <uint8_t>          kmer length to build index[31]\n"); 
        fprintf(stderr, "           -h, --help                          help\n");
        fprintf(stderr, "\n"); 
	return ERROR_PARSE_PARAMS;

}

int UI::classify_usage()
{
        fprintf(stderr, "Usage:     %s [Options] <IndexDir> <ReadFiles>\n\n", PACKAGE_NAME); 
        fprintf(stderr, "<IndexDir>             The directory storing RHT index\n");
        fprintf(stderr, "<ReadFiles>             Reads files, in FASTQ/FASTA format\n");
        
        fprintf(stderr, "Options:   -s, --seed_len      <uint8_t>          lower bound length to find a seed [31]\n"); 
        fprintf(stderr, "           -t, --threads       <int>              number of threads [1]\n");
        fprintf(stderr, "           -h, --help                             help\n");
        fprintf(stderr, "\n"); 
	return ERROR_PARSE_PARAMS;

}

int UI::usage()
{
        fprintf(stderr, "Some description here\n\n");
        fprintf(stderr, "     Usage:\n");
	fprintf(stderr, "           %s -h/--help\n", PACKAGE_NAME);
	fprintf(stderr, "           %s command [options] <arguments>\n", PACKAGE_NAME);
        fprintf(stderr, "     Examples:\n");
	fprintf(stderr, "           %s -h/--help\n", PACKAGE_NAME);
	fprintf(stderr, "           %s command [options] <arguments>\n", PACKAGE_NAME);
        fprintf(stderr, "     Further help:\n");
	fprintf(stderr, "           %s help commands\n", PACKAGE_NAME);
	fprintf(stderr, "           %s help <command>\n", PACKAGE_NAME);
	return ERROR_PARSE_PARAMS;
}
int cmd_usage()
{
	fprintf(stderr,"%s commands are:\n\n", PACKAGE_NAME);
	fprintf(stderr,"                index      build index of reference file\n");
	fprintf(stderr,"                classify   classify reads with given index\n");
	fprintf(stderr,"                help       list further help information\n");
	fprintf(stderr,"\n");
	return ERROR_PARSE_PARAMS;

}

int UI::opt_parse(int argc, char *argv[], opts* opt)
{
	int c; 
    	int option_index=0;

	//fprintf(stderr,"%denter\t%s\n", argc, argv[0]);
	if (argc < 2) return usage();
	if (strcmp(argv[1],"help")==0) {
		if (argc < 3) return usage();

		if (strcmp(argv[2], "index") == 0) {
			return ind_usage();
		} else if (strcmp(argv[2], "classify")== 0) 
				return classify_usage();
			else 
				if (strcmp(argv[2],"commands") == 0) return cmd_usage();
				else return usage();

	} else {
		if (strcmp(argv[1],"index") == 0) {
			opt->isClassify = false;
		} else if (strcmp(argv[1], "classify") == 0) opt->isClassify = true;	
			else return usage();
	}
	while((c = getopt_long(argc, argv, short_options, long_options, &option_index))>=0){
		switch(c){
		    case 's':
			opt->seed = (uint8_t)atoi(optarg);
			break;
		    case 'h':
			return usage();
			break;
		    case 'w':
			//opt->hit_limit = atoi(optarg);
			break;
		    case 'k':
			opt->kmer = (uint8_t)atoi(optarg);
			break;
		    case 't':
			opt->num_threads = atoi(optarg);
			break;
		    default:
			fprintf(stderr,"inappropriate parameters\n");
			return usage();
		      
		}
	}
	if (opt->isClassify) {
		if (optind + 3 > argc) {
			fprintf(stderr, "[opt_parse]: arguments can't be omited!\n"); 
			return classify_usage(); 
		} else {
			++optind;	
			opt->lib = argv[optind++];
			//fprintf(stderr,"%d\t%s\n",optind, argv[optind++]);
			while (optind < argc) opt->reads.push_back(argv[optind++]);
		}
			
	} else {
	
		if (optind + 6 != argc) {
			fprintf(stderr, "[opt_parse]: arguments can't be omited!\n"); 
			return ind_usage(); 
		
		} else {
			++optind;
			opt->sortedKmer = argv[optind++];
			opt->gids = argv[optind++];
			opt->tids = argv[optind++];
			opt->ref = argv[optind++];
			opt->output = argv[optind++];
		
		} 
	}
	opt->argv = argv;
	opt->argc = argc;
	return 0;  

}
