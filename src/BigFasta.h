#ifndef _BIG_FASTA_
#define _BIG_FASTA_

#include <map>
#include <string>
#include <algorithm>
#include <vector>
#include "NuclSeq.h"

using namespace std;

class BigFasta{
public:
    map<int, NuclSeq*> sequences;
    map<string, int> seqLine;
    BigFasta(string filePath);
    ~BigFasta();

    inline NuclSeq* getNuclSeq(int line){
        return sequences[line];
    }

    inline NuclSeq* getNuclSeq(string name){
        int index = seqLine[name];
        return sequences[index];
    }

    /*inline SeqProfile* profileOf(int startLine){
        return sequences[startLine]->getProfile();
    }*/

    //int wordFrequenceIn(string word, int )

    //vector<pair<int, float> >* similaritiesWith(NuclSeq* querySeq);
    
    void deleteAll(){
        for(pair<int, NuclSeq*> entry : sequences){
            if(entry.second != NULL){
                entry.second->free();
            }
        }
    }

    //void compareWithFasta(BigFasta *queryFasta);
};

#endif