#include <string>
#include <iostream>
#include "Alignment.h"
#include "NuclSeq.h"
#include "BigFasta.h"

using namespace std;



int main(int argc, char *argv[]){
    if(argc >= 3){
        string queryFile(argv[1]);
        string databaseFile(argv[2]);

        makeAlignment(queryFile, databaseFile);
    }else{
        cout << "Please inform the query and database fasta files\n\t$ ./aligner [queryFile] [databaseFile]\n";
    }
}

void makeAlignment(string queryFile, string databaseFile){
    /*string seqA = "GAATTCAGTGATTGAA-GAACGCAGTAACCT";
    string seqB = "---GAACGCAGTAACCT";
    string seqC = "AGAGC";
    string seqD = "CAGC";

    Alignment test(seqA, seqB);
    
    cout << "Initialized alignmenter" << endl;
    auto result = test.align();
    cout << result.first << "\n";
    cout << result.second << "\n";

    NuclSeq seq("../query.fa", 1);
    NuclSeq seq2("../query.fa", 85);
    NuclSeq seq3("../query.fa", 71);
    
    cout << ">" << seq.name << endl << seq.contains(seq3) << endl;
    cout << ">" << seq2.name << endl << seq2.contains(seq3) << endl;
    cout << ">" << seq3.name << endl << seq3.contains(seq3) << endl;

    Alignment test(seq3, seq2);
    auto result = test.align();
    cout << result.first << "\n";
    cout << result.second << "\n";*/

    BigFasta NuclSeq(queryFile, 71);
    BigFasta littleOne(databaseFile);
    cout << "Read entire fasta file" << endl;
    /*for(pair<int, NuclSeq*> entry : littleOne.sequences){
        cout << "Line " << entry.first << endl;
        string* content = entry.second->getContent();
        cout << *content << endl;
    }*/

    vector<pair<int, float> > sims = *littleOne.similaritiesWith(seq2);
    cout << "Query seq: " << seq2->name << endl;
    cout << *seq2->getContent() << endl;
    NuclSeq* bestMatch = NULL;
    float sim = 0.0;
    for(pair<int, float> entry : sims){
        //cout << "Line " << entry.first << endl;
        //cout << "Similarity with\n\t" << littleOne.sequences[entry.first]->name << endl;
        //cout << entry.second << endl;
        bestMatch = littleOne.sequences[entry.first];
        sim = entry.second;
    }

    cout << "Line " << bestMatch->startAtLineXInFasta << endl;
    cout << "Similarity with\n\t" << bestMatch->name << endl;
    cout << sim << endl;

    Alignment test(*seq2, *bestMatch);
    auto result = test.align();
    cout << result.first << "\n";
    cout << result.second << "\n";
}