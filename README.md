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
		`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"
		   `=`,'=/     `=`,'=/     `=`,'=/     `=`,'=/
		     y==/        y==/        y==/        y==/
		   ,=,-<=`.    ,=,-<=`.    ,=,-<=`.    ,=,-<=`.
		,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_
```
# alenhador
Parallel aligner (Needleman–Wunsch) for nucleotide sequences (.fasta) with similarity heuristic filter.

# what it does?
Searches for query .fasta nucleotide sequences inside a .fasta file with many different sequences, using dynamic and parallel computing.

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

# How it works?
In order to search for query sequences in a big .fasta "database" file, "alenhador" first reads the database and at the same time filters sequences who are not similar enough to the query. Sequences filtered are not stored in memory anymore, and that helps reducing the memory usage. After reducing the possible sequences to align, it does the Needleman–Wunsch algorithm on each sequence concurrently.

## The Heuristic Filter
The project filters database sequences who are not similar to the query, so that alignments are only made with sequences who have a chance of being a match. First, it creates a profile of every sequence in the query and the database.

For a sequence S, the profile contais:
- The number of occurences of a word W in F;
- The sequence length;
- The method "float compare(otherProfile)";
### float compare()
The method "profile.compare(otherProfile)" returns a floating point number 0 <= x <= 1.0, which indicates the chance of "profile" containing "otherProfile". The result is not perfect and precision may vary depending on the words length.
```cpp
float SeqProfile::compare(SeqProfile* other){
    vector<float> comparations;
    for(pair<string, int> entry : other->wordCounts)
    {
        if(wordCounts.find(entry.first) != wordCounts.end()){
            float value = (float)wordCounts[entry.first]/entry.second;            
            comparations.push_back(value);
        }
    }

    if(comparations.size() > 0){
        float result = 0.0;
        for(float comp : comparations){
            result += comp;
        }
        return result/other->length;
    }else{
        return 0.0;
    }
}

```
## The Alignment algorithm
The project uses Needleman–Wunsch, an O(n^2) algorithm for sequence alignment that uses dynamic programming. The result is an optimal alignment. The alignment of each sequence is made concurrently. An absolute score for the alignment is returned along with the alignment itself and the relative score (score/len(subjectSequence)) is displayed too.

## Overall process
```py
for each querySeq in queryFasta:
	filteredSeqs = dbFasta.filterSeqsAccordingTo(querySeq)
	alignments = list()
	for each filteredSeq in filteredSeqs:
		alignment = NeedlemanWunsch(querySeq, filteredSeq)
		alignment.start()
		alignments.add(alignment)
	waitForAllAreFinished(alignments)
	printResults(alignments)
```

### filterSeqsAccordingTo(querySeq)
To maximize the efficiency of this part of the process, it is divided in three parts (each separated and concurrent with the others): reading fasta database, profiling sequences, sorting the most similar.
```py
numberOfProfilingThreads = max(std::thread::hardware_concurrency() - 2, 2)
detach [numberOfProfilingThreads] profiling threads
detach sorting thread
startReadingDB
while(sorting or profiling):
	wait
return filteredSeqs
```
The fasta reading thread (the main thread) reads the file and puts the readen sequences in a profiling queue. The profiling threads pop a sequence from the profiling queue, create of profile of it and the put it in a sorting queue. The sorting thread pops a profiled sequence from the sorting queue and tries to put it in the filteredSeqs list, which has a limited size, and sequences who cannot fit are discarted from the memory.
