```
 .acgatagcgcgctatata.ACATATATATATATATATATGC|AGAGAGAGAGA.
G                   T                      T           |`-.,_
CTTTTTTTTTTTTTTTTTTTTTCGAGAGAGATATATATATACA|AGAGAGAGAGA|,-'`
 `aaaaaaaaaaaaaaaaaa'                       :    ___   l
                                            /   C   T   \
                                           /     `-'     \
                                         ,'               `.
                                      .aTg.               .aTg.
                                      "AGGGGc..       ,.cGGGGA"
                                        "CAAAAAAAAAAAAAAAAAG"
                                           ""AATTTTTTTGG""
```

# alenhador
Parallel aligner for nucleotide sequences (.fasta) with similarity heuristic filter

# build it

Requires 'g++' with C++14 available and 'make'

```sh
    $ git clone https://github.com/pentalpha/alenhador.git
    $ cd alenhador/
    $ make
```

# --help
```sh
    $ ./bin/alenhador -h
    HELP: 
    ./alenhador [fasta query file] [fasta DB file] n
        [fasta query file] = A .fasta file with sequences of 
                           nucleotides to be searched.
        [fasta DB file]	   = A .fasta file with a lot of sequences,
                           where you want to search the query.
        n                  = The maximum number of results for each 
                           query. Default=6.
```
# Example usage:
```sh
        ./alenhador example-data/query.fa example-data/database.fa 3
        Query sequence:
	       >BC036785.1 Homo sapiens TP53 target 5, mRNA (cDNA clone MGC:46104 IMAGE:5744881), complete cds
	       	Length: 1040bp
        [Profiling threads = 2]
        [Starting sequence profiling thread]
        [Starting sequence profiling thread]
        [Reading DB]
        [Ready to filter sequences by profile]
        [Finished reading DB]
        [Finished profiling sequences]
        [Finished profiling sequences]
        [Finished filtering sequences]

        querySq = '>BC036785.1 Homo sapiens TP53 target 5, mRNA (cDNA clone MGC:46104 IMAGE:5744881), complete cds'
        subject = '>Query-p53'
        Alignment Score = 5200.000000
        Alignment Relative Score = 1.000000

        querySq	0	GCTGGCTGAACTGAGAGGAACAGGGTTGGT	30
        subject	0	GCTGGCTGAACTGAGAGGAACAGGGTTGGT	30

        querySq	30	GCCTGGCACTGGTGTTGCTCCATTCATCTC	60
        subject	30	GCCTGGCACTGGTGTTGCTCCATTCATCTC	60
        
        [Some hundreds of alignment lines later]
        
        querySq	15300	------------------------------	15330
        subject	15300	CCCGGGCAGAGGGATTCCGAACCCGAGAAA	15330

        querySq	15330	----------------------A	15353
        subject	15330	TAAAAGTCTGTTCCACCCCCTGG	15353


        Search for:
	              >BC036785.1 Homo sapiens TP53 target 5, mRNA (cDNA clone MGC:46104 IMAGE:5744881), complete cds
        Overall Results:
        >X15750.1 Rabbit skeletal muscle mRNA for ryanodine receptor
	              Profile Similarity: 1.92308%
	              Alignment Score: 3103
	              Alignment Relative Score: 4.04221%

        >Query-p53
	              Profile Similarity: 94.1346%
	              Alignment Score: 5200
	              Alignment Relative Score: 100%

        >U37667.1 Human BRCA1 gene, partial cds
	              Profile Similarity: 2.42459%
	              Alignment Score: 3050
	              Alignment Relative Score: 16.0695%
```
