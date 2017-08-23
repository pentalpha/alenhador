#ifndef _SEQ_PROFILE_
#define _SEQ_PROFILE_

#include <string>
#include <iostream>
#include <map>
#include <iostream>

using namespace std;

class SeqProfile{
public:
    map<string, int> wordCounts;
    int length;
    SeqProfile(string *seq, int wordSize);
    int wordCount(string &word);
    void increaceWordCount(string &word);
    float compare(SeqProfile* other);
};

#endif