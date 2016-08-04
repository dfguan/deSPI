### Quick Start 

```
git clone https://github.com/derekguan/deSPI.git

cd deSPI/src

make

cd ..

bin/deSPI-download -o taxonomy taxonomy  

bin/deSPI-index test/refs/ taxonomy/nodes.dmp bin index

bin/deSPI classify index test/reads/reads.fa >Labels
```

---

### Introduction

deSPI is a novel read classification tool which classifies reads by recognizing and analyzing the matches between reads and reference with de Bruijin graph-based lightweight reference indexing.
 
deSPI is mainly designed by Bo Liu and developed by Dengfeng Guan in Center for Bioinformatics, Harbin Institute of Technology, China.

### Parameters
#### deSPI-download
```
deSPI-download [<options>] <database>

ARGUMENT
 <database>        One of refseq, genbank, or taxonomy:
                     - use refseq or genbank for genomic sequences,
                     - taxonomy for taxonomy mappings.

COMMON OPTIONS
 -o <directory>         Folder to which the files are downloaded. Default: '.'.
 -P <# of threads>      Number of processes when downloading (uses xargs). Default: '1'

WHEN USING database refseq OR genbank:
 -d <domain>            What domain to download. One or more of bacteria, viral, archaea (comma separated).
 -a <assembly level>    Only download genomes with the specified assembly level. Default: 'Complete Genome'.
 -c <refseq category>   Only download genomes in the specified refseq category. Default: any.
 -g                     Download GI map.
```
#### deSPI-index
```
deSPI-index <REF_DIR> <NODE_PATH> <DESPI_BIN_DIR> <INDEX_DIR>

ARGUMENT
    <REF_DIR>                directory of reference library.
    <NODE_PATH>              location of nodes.dmp file.
    <DESPI_BIN_DIR>          directory of deSPI execuative files.
    <INDEX_DIR>              directory to store deSPI's index.
```
#### deSPI 
```
Usage:     deSPI classify [Options] <IndexDir> <ReadFiles>

<IndexDir>             The directory storing deSPI index
<ReadFiles>            Reads files, in FASTQ/FASTA format (separated by space)
Options:   -s, --seed_len      <uint8_t>          lower bound of seed length [30]
           -t, --threads       <int>              # of threads [1]
           -h, --help                             help
```

### Memory Requirements

With 1742 refseq complete genomes, deSPI requires 11 gigabyte to complete classification task, however it consumes much more space to build deSPI index. To facilitate 








### User's Guide

#### Database Download
Users could take the following command to build their own reference library. All reference genomes will be download from NCBI genbank or refseq.

```
cd deSPI

bin/deSPI-download -o LOCAL_REF_DIR -P THREAD_NUM -d DOMAIN DATABASE 

```
DOMAIN includes viral, bacteria, archaea

DATABSE includes refseq, genbank, taxonomy



#### Index Construction
After creating a reference library, deSPI index could be constructed by the following command, since deSPI needs nodes.dmp as input, nodes.dmp requires be downloaded first. (Utilize the command above with DATABASE is taxonomy)
```
cd deSPI

bin/deSPI-index LOCAL_REF_DIR taxonmy/nodes.dmp bin INDEX_DIR

```
---

#### Classify Reads
After deSPI is done with index construction. User could run the following command to classify reads files.

```
cd deSPI

bin/deSPI classify <INDEX_DIR> <READS_FILES>

```
---

###Contact
For advising, bug reporting and requiring help, please contact ydwang@hit.edu.cn
