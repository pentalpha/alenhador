#ifndef _MSEQ_PROFILE_
#define _MSEQ_PROFILE_

#include <string>
#include <iostream>
#include <map>
#include <iostream>

using namespace std;

class MultiSeqProfile{
public:
    map<string, map<int, int>* > wordCounts;
    map<int, int> seqLength;
    int wordSize;
    MultiSeqProfile(int wordSize);
    void includeSeq(int seqLine, string seq);
    int wordCount(string word, int seqLine);
    void increaceWordCountInSeq(string word, int seqLine);
    int getSeqLength(int index);
    //MultiSeqProfile(string fileName);
    float compare(MultiSeqProfile* other, int seqIndex);
    //void toFile(string fileName, string seqName);
};

#endif