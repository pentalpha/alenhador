#include "BigFasta.h"

bool lessOnPair(const std::pair<int,float> &left, const std::pair<int,float> &right) {
    return left.second < right.second;
}

BigFasta::BigFasta(string fileName){
    //system("mkdir tmp");
    ifstream fastaFile;
    fastaFile.open(fileName);
    
    if (fastaFile.is_open())
    {
        int nline = 1;
        int lastStart = -1;
        string name = "";
        string *sequence;
        string line;
        while (getline(fastaFile,line))
        {   
            if(line.length() > 0){
                if(line[0] == '<' || line[0] == '>'){
                    if(lastStart > 0){
                        sequences[lastStart] = new NuclSeq(fileName, lastStart, name, sequence);
                        sequence = NULL;
                    }
                    //cout << "Seq start at " << nline << " name:" << endl;
                    name = line;
                    //cout << name << endl;
                    lastStart = nline;
                    sequence = new string("");
                }else if(sequence != NULL){
                    *sequence += line;
                }
            }
            nline++;
        }
        if(sequence != NULL){
            sequences[lastStart] = new NuclSeq(fileName, lastStart, name, sequence);
        }
        fastaFile.close();
    }
}

vector<pair<int, float> >* BigFasta::similaritiesWith(NuclSeq* querySeq){
    vector<pair<int, float> > *similarities = new vector<pair<int, float> >();
    cout << "Started calculating similarities" << endl;
    for(pair<int, NuclSeq*> entry : sequences){
        //cout << "Similarity for seq at " << entry.first << ": " << endl;
        float s = entry.second->contains(querySeq);
        //cout << s << endl;
        pair<int, float> newSim;
        newSim.first = entry.first;
        newSim.second = s;
        similarities->push_back(newSim);
    }
    cout << "Finished calculating similarities." << endl;
    sort(similarities->begin(), similarities->end(), lessOnPair);
    return similarities;
}

BigFasta::~BigFasta(){
    deleteAll();
}