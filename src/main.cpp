#include <string>
#include <iostream>
#include <map>
#include <list>
#include <vector>
#include "Alignment.h"
#include "NuclSeq.h"
#include "FastaHeuristicFilter.h"
#include "BigFasta.h"
#include "FileUtils.h"

using namespace std;

void searchEachInQuery(string query, string db, int max, int specificIndex);
void searchSeqInFastaDB(NuclSeq *querySeq, string dbPath, int maxFiltered);
void printHelp();
void printInvalidArgument(string arg);

int main(int argc, char *argv[]){
    vector<string> args;
    bool helpArgFound = false;
    for(int i = 1; i < argc; i++){
        string s(argv[i]);
        if(s == "-h" || s == "--help"){
            helpArgFound = true;
            break;
        }
        args.push_back(s);
    }

    if(helpArgFound){
        printHelp();
    }else{
        string queryFilePath = "";
        string databaseFilePath = "";
        int maxToAlign = 4;
        int specificQuery = -1;

        for(string arg : args){
            if(fileExists(arg.c_str())){
                if(queryFilePath == ""){
                    queryFilePath = arg;
                }else if (databaseFilePath == ""){
                    databaseFilePath = arg;
                }else{
                    printInvalidArgument(arg);
                }
            }else if(arg.length() > 2){
                if(arg[1] == '='){
                    string rawNumber = arg.substr(2,arg.length()-2);
                    try{
                        int value = stoi(rawNumber);
                        //cout << arg << "=" << value << endl;
                        if(arg[0] == 'n'){
                            maxToAlign = value;
                        }else if(arg[0] == 'i'){
                            specificQuery = value;
                        }
                    }catch(...){
                        printInvalidArgument(rawNumber);
                    }
                }else{
                    printInvalidArgument(arg);
                }
            }else{
                printInvalidArgument(arg);
            }
        }

        if(queryFilePath != "" && databaseFilePath != ""){
            cout << "Query\t\t" << queryFilePath << endl;
            cout << "DB\t\t" << databaseFilePath << endl;
            cout << "Max. Results\t" << maxToAlign << endl;
            if(specificQuery == -1){
                cout << "Querys\tALL" << endl;
            }else{
                cout << "Align only\t" << specificQuery+1 << "th" << endl;
            }
            
            searchEachInQuery(queryFilePath, databaseFilePath, maxToAlign, specificQuery);
        }else{
            cout << "Please inform the query and database '*.fasta' files. " << endl;
            printHelp();
        }
    }

    return 0;
}

void printInvalidArgument(string arg){
    cout << "ERROR! Invalid argument: " << arg << endl;
}

void printHelp(){
    cout << "HELP: " << endl;
    cout << "./alenhador <query file name>.fasta <DB file name>.fasta [n=X] [i=Y]" << endl
        << "\t<query file name>.fasta\t= A .fasta file with sequences of " << endl
        << "\t\t\t\tnucleotides to be searched." << endl
        << "\t<DB file name>.fasta\t= A .fasta file with a lot of sequences," << endl
        << "\t\t\t\twhere you want to search the query." << endl
        << "\tn\t\t\t= The maximum number of results for each " << endl 
        << "\t\t\t\tquery. Default=4." << endl
        << "\ti\t\t\t= If specified, will only search for the  " << endl 
        << "\t\t\t\tthe specific query with index 'i'." << endl
        << "Example usage:" << endl
        << "\t./alenhador example-data/query.fa example-data/database.fa n=7 i=0" << endl;
}

void searchEachInQuery(string query, string db, int max, int specificIndex){
    BigFasta querySeqs(query);
    int count = 0;
    for(pair<int, NuclSeq*> entry : querySeqs.sequences){
        NuclSeq *seq = entry.second;
        if(seq == NULL){
            cout << "Query at line " << entry.first << " of " << query 
            << " could not be read correctly, not searching it." << endl;
        }else{
            if(count == specificIndex || specificIndex == -1){
                if(count > 0){
                    cout << "\n##########################################################" << endl;
                    cout << "##########################################################" << endl;
                    cout << "##########################################################" << endl;
                }
                searchSeqInFastaDB(seq, db, max);
            }
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