```
 .acgatagcgcgctatata.ACATATATATATATATATATGC|AGAGAGAGAGA.
/                   |                      |           |`-.,_
\TTTTTTTTTTTTTTTTTTT|TCGAGAGAGATATATATATACA|AGAGAGAGAGA|,-'`
 `------------------'                       :    ___   l
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
    $ make
```

# --help
```sh
    $ ./bin/alenhador -h
    HELP: 
    ./alenhador [fasta query file] [fasta DB file] n
        [fasta query file]	= A .fasta file with sequences of 
                            nucleotides to be searched.
        [fasta DB file]		= A .fasta file with a lot of sequences,
                            where you want to search the query.
        n			        = The maximum number of results for each 
                            query. Default=6.
    Example usage:
        ./alenhador example-data/query.fa example-data/database.fa 7
```
