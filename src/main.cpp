#include <string>
#include <iostream>
#include <map>
#include <list>
#include "Alignment.h"
#include "NuclSeq.h"
#include "FastaHeuristicFilter.h"
#include "BigFasta.h"

using namespace std;

void searchEachInQuery(string query, string db, int max);
void searchSeqInFastaDB(NuclSeq *querySeq, string dbPath, int maxFiltered);
void printHelp();

int main(int argc, char *argv[]){
    if(argc >= 3){
        string queryFile(argv[1]);
        cout << "Query .fasta file: " << queryFile << endl;
        string databaseFile(argv[2]);
        cout << "Database .fasta file: " << queryFile << endl;
        int maxToAlign = 6;
        if(argc >= 4){
            maxToAlign = stoi(argv[3]);
        }
        searchEachInQuery(queryFile, databaseFile, maxToAlign);
    }else if (argc == 2){
        printHelp();
    }else{
        cout << "Please inform the query and database fasta files" << endl;
        printHelp();
    }

    return 0;
}

void printHelp(){
    cout << "HELP: " << endl;
    cout << "./alenhador [fasta query file] [fasta DB file] n" << endl
        << "\t[fasta query file]\t= A .fasta file with sequences of " << endl
        << "\t\t\t\tnucleotides to be searched." << endl
        << "\t[fasta DB file]\t\t= A .fasta file with a lot of sequences," << endl
        << "\t\t\t\twhere you want to search the query." << endl
        << "\tn\t\t\t= The maximum number of results for each " << endl 
        << "\t\t\t\tquery. Default=6." << endl
        << "Example usage:" << endl
        << "\t./alenhador example-data/query.fa example-data/database.fa 7" << endl;
}

void searchEachInQuery(string query, string db, int max){
    BigFasta querySeqs(query);
    for(pair<int, NuclSeq*> entry : querySeqs.sequences){
        NuclSeq *seq = entry.second;
        if(seq == NULL){
            cout << "Query at line " << entry.first << " of " << query 
            << " could not be read correctly, not searching it." << endl;
        }else{
            searchSeqInFastaDB(seq, db, max);
            //delete seq;
        }
    }
}

void searchSeqInFastaDB(NuclSeq *querySeq, string dbPath, int maxFiltered){
    FastaHeuristicFilter filter(querySeq, dbPath, maxFiltered);
    cout << "\n##########################################################" << endl;
    cout << "##########################################################" << endl;
    cout << "##########################################################" << endl;
    cout << "\nQuery sequence:\n\t" << *querySeq->name << endl
        <<"\tLength: " << querySeq->getContent()->length() << "bp" << endl;
    cout << "Similtar to: " << endl;
    list< pair<NuclSeq*, float> > filteredSequences = filter.justDoIt();
    if(filteredSequences.size() > 0){
        for(pair<NuclSeq*, float> seq : filteredSequences){
            cout << "Similarity with\n\t" << *seq.first->name 
            << "(" << seq.first->getContent()->length() <<"bp)" << endl
            << "\n\t" << seq.second << endl;
        }
    }else{
        cout << "\tNone\n";
    }
}

/*void makeAlignment(string queryFile, string databaseFile){
    string seqA = "GAATTCAGTGATTGAA-GAACGCAGTAACCT";
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
    cout << result.second << "\n";

    BigFasta NuclSeq(queryFile, 71);
    BigFasta littleOne(databaseFile);
    cout << "Read entire fasta file" << endl;
    /*for(pair<int, NuclSeq*> entry : littleOne.sequences){
        cout << "Line " << entry.first << endl;
        string* content = entry.second->getContent();
        cout << *content << endl;
    }

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
}*/