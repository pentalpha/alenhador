#ifndef _Fasta_Heuristic_Filter_
#define _Fasta_Heuristic_Filter_

#include <map>
#include <list>
#include <vector>
#include <string>
#include <deque>
#include <mutex>
#include <thread>
#include <math.h>
#include "SeqProfile.h"
#include "NuclSeq.h"

using namespace std;

struct SeqNode {
    NuclSeq* seq;
    float similarity;
};

class FastaHeuristicFilter{
public:
    FastaHeuristicFilter(NuclSeq *querySeq, string dbPath, int maxFiltered);
    list< pair<NuclSeq*, float> > justDoIt();
private:
    int filtered = 0;
    int seqsReadFromDB = 0;
    SeqProfile queryProfile;
    NuclSeq *querySeq;
    string databasePath;
    int maxFiltered;
    static int defaultWordSize;

    bool readingDB = false;
    mutex profileQueueLock;
    deque< pair<string*, string*> > toProfile;

    bool profiling = false;
    mutex sortQueueLock;
    deque<SeqNode> toSort;
    
    bool sorting = false;
    list<SeqNode> filteredNodes;

    void readDB();

    void addToProfilingQueue(string* name, string* content);
    pair<string*, string*> getRawSeq();
    int seqsToProfile();
    void compareSeqs();

    void addToSortingQueue(SeqNode newNode);
    SeqNode getUnsortedNode();
    int seqsToSort();
    void sortSeqs();

};

#endif