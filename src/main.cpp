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
    int count = 0;
    for(pair<int, NuclSeq*> entry : querySeqs.sequences){
        NuclSeq *seq = entry.second;
        if(seq == NULL){
            cout << "Query at line " << entry.first << " of " << query 
            << " could not be read correctly, not searching it." << endl;
        }else{
            if(count > 0){
                cout << "\n##########################################################" << endl;
                cout << "##########################################################" << endl;
                cout << "##########################################################" << endl;
            }
            searchSeqInFastaDB(seq, db, max);
            count++;
            //delete seq;
        }
    }
}

void searchSeqInFastaDB(NuclSeq* querySeq, string dbPath, int maxFiltered){
    FastaHeuristicFilter filter(querySeq, dbPath, maxFiltered);
    
    cout << "\nQuery sequence:\n\t" << *querySeq->name << endl
        <<"\tLength: " << querySeq->getContent()->length() << "bp" << endl;
    pair<list<NuclSeq*>, map<NuclSeq*, float> > filteringResults = filter.justDoIt();
    map<NuclSeq*, float> filteredSequences = filteringResults.second;
    list<NuclSeq*> sortedFilteredSeqs = filteringResults.first;
    map<NuclSeq*, Alignment*> alignments;
    for(NuclSeq* filteredSeq : sortedFilteredSeqs){
        Alignment *alignment = new Alignment(querySeq, filteredSeq);
        alignment->startAlignment();
        alignments[filteredSeq] = alignment;
    }
    int completed = 0;
    map<NuclSeq*, bool> printed;
    for(NuclSeq* filteredSeq : sortedFilteredSeqs){
        printed[filteredSeq] = false;
    }
    while(completed < alignments.size()){
        for(NuclSeq* filteredSeq : sortedFilteredSeqs){
            if(!(printed[filteredSeq]) && alignments[filteredSeq]->isDone()){
                cout << alignments[filteredSeq]->getAlignment() << endl;
                printed[filteredSeq] = true;
                completed++;
            }
        }
    }

    cout << "Search for:\n\t" << *querySeq->name << endl; 
    cout << "Overall Results:" << endl;
    for(NuclSeq* filteredSeq : sortedFilteredSeqs){
        cout << *(filteredSeq->name) << endl;
        cout << "\t" << "Profile Similarity: " << (filteredSequences[filteredSeq])*100 << "%" << endl;
        cout << "\t" << "Alignment Score: " << alignments[filteredSeq]->getAlignmentValue() << endl;
        cout << "\t" << "Alignment Relative Score: " << (alignments[filteredSeq]->getAlignmentRelativeValue())*100 << "%\n\n";
    }
}