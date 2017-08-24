#include "SeqProfile.h"

SeqProfile::SeqProfile(string *seq, int wordSize){
    if(seq == NULL){
        cout << "Cannot make the profile of a NULL string" << endl;
    }
    length = seq->length();
    string word;
    for(int i = 0; i < seq->length()-wordSize; i++){
        word = seq->substr(i, wordSize);
        increaceWordCount(word);
    }
}

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

void SeqProfile::increaceWordCount(string &word){
    if(wordCounts.find(word) != wordCounts.end()){
        wordCounts[word]++;
    }else{
        wordCounts[word] = 1;
    }
}

int SeqProfile::wordCount(string &word){
    if(wordCounts.find(word) != wordCounts.end()){
        return wordCounts[word];
    }else{
        return 0;
    }
}