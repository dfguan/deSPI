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



char *const short_options = "k:x:r:s:i:t:pah";
struct option long_options[] = {
    //{ "kmer_len",     1,   NULL,    'm'   },
    { "kmers",     1,   NULL,    'k' },
    //{ "gids",     1,   NULL,    'g'   },
    //{ "tids",     1,   NULL,    'a'   },
    //{ "ref",	1,NULL,	'r'},
    //{ "output", 1,  NULL, 'o'},
    {"max_it", 1, NULL, 'r'},
    { "seed_len", 1, NULL, 's'},
    { "interv", 1, NULL, 'i'},
    { "threads",1, NULL,'t'},
    {"paired", 0, NULL, 'p'},
    //{"gapextended", 1,  NULL,'e'},
    //{"match",   1,  NULL,'c'},
	{"all", 0, NULL, 'a'},
	{"taxid", 1, NULL, 'x'},
	{ "help",    0,  NULL,'h'},
    { 0,     0,   0,    0   }
};


UI::UI(opts *opt)
{
	opt->seed = 24;
	opt->kmer = 31;
	//opt->inv = 5;
	opt->n_thread = 1;
	//opt->iteration = 4;
	opt->intv= 4;
	opt->iter = 5;
	opt->isPaired = false;
	opt->out_all = true;
	opt->taxids = "";
}

int UI::view_usage()
{
	fprintf(stderr, "\n"); 
	fprintf(stderr, "Program:   %s\n", PACKAGE_NAME); 
	fprintf(stderr, "Version:   %s\n", PACKAGE_VERSION); 
	fprintf(stderr, "Contact:   %s\n\n", CONTACT); 
	fprintf(stderr, "Usage:     %s  view  [Options] <IndexDir>\n\n", PACKAGE_NAME); 
	fprintf(stderr, "<IndexDir>              the directory contains deSPI index\n");
	
	fprintf(stderr, "           -x, --taxid      <uint32_t>+          output unitigs corresponding to the taxonomy id, use comma to join multiple ids, output all unitigs if not set\n"); 
	fprintf(stderr, "           -h, --help                           help\n");
	fprintf(stderr, "\n"); 
	return ERROR_PARSE_PARAMS;

}




int UI::ind_usage()
{
	fprintf(stderr, "\n"); 
	fprintf(stderr, "Program:   %s\n", PACKAGE_NAME); 
	fprintf(stderr, "Version:   %s\n", PACKAGE_VERSION); 
	fprintf(stderr, "Contact:   %s\n\n", CONTACT); 
        fprintf(stderr, "Usage:     %s  index  [Options] <SortedKmer> <Taxonomy> <Reference> <IndexDir>\n\n", PACKAGE_NAME); 
        fprintf(stderr, "<SortedKmer>            sorted kmers\n");
        //fprintf(stderr, "<gidTid>               correspondence between gid and tid\n");
        fprintf(stderr, "<Taxonomy>              taxonomy tree\n");
        fprintf(stderr, "<Reference>             references file, in FASTQ/FASTA format\n");
        fprintf(stderr, "<IndexDir>              the directory used to store deSPI index\n\n");
        
        fprintf(stderr, "Options:   -k, --kmer       <uint8_t>          kmer size of de Bruijin graph [31]\n"); 
        fprintf(stderr, "           -h, --help                          help\n");
        fprintf(stderr, "\n"); 
	return ERROR_PARSE_PARAMS;

}

int UI::classify_usage()
{
	fprintf(stderr, "\n"); 
	fprintf(stderr, "Program:   %s\n", PACKAGE_NAME); 
	fprintf(stderr, "Version:   %s\n", PACKAGE_VERSION); 
	fprintf(stderr, "Contact:   %s\n\n", CONTACT); 
        fprintf(stderr, "Usage:     %s  classify  [Options] <IndexDir> <ReadFiles>\n\n", PACKAGE_NAME); 
        fprintf(stderr, "<IndexDir>              the directory contains deSPI index\n");
        fprintf(stderr, "<ReadFiles>             reads files, in FASTQ/FASTA format (separated by space)\n\n");
        
        fprintf(stderr, "Options:   -s, --seed_len      <uint8_t>          lower bound of seed length from 24 [24]\n"); 
        fprintf(stderr, "           -t, --threads       <int>              number of threads [1]\n");
	//fprintf(stderr, "           -r, --max_it        <uint8_t>          maximal iteration times [4]\n");
		fprintf(stderr, "           -i, --seed_step     <int8_t>           seed step size [4]\n");
		fprintf(stderr, "           -p, --paired        <int8_t>           set for paired end reads [False]\n");
		fprintf(stderr, "           -h, --help                             help\n");
        fprintf(stderr, "\n"); 
	return ERROR_PARSE_PARAMS;

}

int UI::usage()
{
	fprintf(stderr, "\n"); 
	fprintf(stderr, "     Program:   %s\n", PACKAGE_NAME); 
	fprintf(stderr, "     Version:   %s\n", PACKAGE_VERSION); 
	fprintf(stderr, "     Contact:   %s\n\n", CONTACT); 
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
	fprintf(stderr,"                index      build index for a reference library\n");
	fprintf(stderr,"                classify   classify metagenomics reads\n");
	fprintf(stderr,"                view       extract sequences in database\n");
	fprintf(stderr,"                help       list further help information\n");
	fprintf(stderr,"\n");
	return ERROR_PARSE_PARAMS;

}

int UI::opt_parse(int argc, char *argv[], opts* opt)
{
	int c; 
    	int option_index=0;
	
	if (argc < 2) return usage();
	if (strcmp(argv[1],"help")==0) {
		if (argc < 3) return usage();
		if (strcmp(argv[2], "index") == 0) 
			return ind_usage();
		else if (strcmp(argv[2], "classify")== 0) 
			return classify_usage();
		else if (strcmp(argv[2], "view") == 0) 
			return view_usage();
		else if (strcmp(argv[2],"commands") == 0) 
			return cmd_usage();
		else 
			return usage();
	} else {
		if (strcmp(argv[1],"index") == 0) 
			opt->option = INDEX;
		else if (strcmp(argv[1], "classify") == 0) 
			opt->option = CLASSIFY;
		else if (strcmp(argv[1], "view") == 0)
			opt->option = VIEW;	
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
		    case 'r':
				opt->iter = atoi(optarg);
				break;
		    case 'i':
				opt->intv = atoi(optarg);
				break;
		    case 'k':
				opt->kmer = (uint8_t)atoi(optarg);
				break;
		    case 't':
				opt->n_thread = atoi(optarg);
				break;
		    case 'p':
				opt->isPaired = true;
				break;
			case 'a':
				opt->out_all = true;
				break;
			case 'x':
				opt->taxids = optarg;
				break;
		    default:
				fprintf(stderr,"inappropriate parameters\n");
				return usage();
		      
		}
	}
	if (opt->option == CLASSIFY) {
	//for (int i=optind; i < argc; ++i)
			//fprintf(stderr,"%d\t%s\t",i,argv[i]);
		if (optind + 3 > argc) {
			fprintf(stderr, "[opt_parse]: arguments can't be omited!\n"); 
			return classify_usage(); 
		} else {
			++optind;	
			opt->index_dir = argv[optind++];
			//fprintf(stderr,"%d\t%s\n",optind, argv[optind++]);
			while (optind < argc) opt->read_fns.push_back(argv[optind++]);
		}
			
	} else if (opt->option == INDEX) {
	
		if (optind + 5 > argc) {
			fprintf(stderr, "[opt_parse]: arguments can't be omited!\n"); 
			return ind_usage(); 
		
		} else {
			++optind;
			opt->sortedKmer = argv[optind++];
			//opt->gids = argv[optind++];
			opt->evo_tree_path = argv[optind++];
			opt->ref = argv[optind++];
			opt->output = argv[optind++];
		
		} 
	} else {
		if (optind + 1 > argc) {
			fprintf(stderr, "[opt_parse]: arguments can't be omited!\n"); 
			return view_usage(); 
		}
		opt->index_dir = argv[optind];

	}
	opt->argv = argv;
	opt->argc = argc;
	return 0;  

}
