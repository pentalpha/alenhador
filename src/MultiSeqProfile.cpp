/*MultiSeqProfile::MultiSeqProfile(string seq, int wordSize){
    length = seq.length();
    string word;
    for(int i = 0; i < seq.length()-wordSize; i++){
        word = seq.substr(i, wordSize);
        if(wordCounts.find(word) != wordCounts.end()){
            wordCounts[word]++;
        }else{
            wordCounts[word] = 1;
        }
    }
}*/

/*MultiSeqProfile::MultiSeqProfile(string fileName){
    ifstream profileFile;
    profileFile.open(fileName);
    int nline = 1;
    if (profileFile.is_open())
    {
        string* contentTmp = new string("");
        string line;
        string word;
        int count;
        while (getline(profileFile,line))
        {   
            if(nline > 2){
                profileFile >> word >> count;
                wordCounts[word] = count;
            }else if (nline == 2){
                profileFile >> length;
            }
            nline++;
        }
        profileFile.close();
    }
}*/

float MultiSeqProfile::compare(MultiSeqProfile* other){
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

/*void MultiSeqProfile::toFile(string fileName, string seqName){
    ofstream profFile;
    profFile.open(fileName);
    if (profFile.is_open())
    {
        //MultiSeqProfile profile(*content, defaultWordSize);
        profFile << seqName << "\n";
        profFile << length << "\n";
        for(pair<string, int> entry : wordCounts){
            profFile << entry.first << "\t" << entry.second << "\n";
        }
        profFile.close();
    }
}*/

MultiSeqProfile::MultiSeqProfile(int wordLength){
    wordSize = wordLength;
    
}

void MultiSeqProfile::includeSeq(int seqLine, string seq){
    string word;
    for(int i = 0; i < seq.length()-wordSize; i++){
        word = seq.substr(i, wordSize);
        increaceWordCountInSeq(word, seqLine);
    }
}

void MultiSeqProfile::increaceWordCountInSeq(string &word, int seqLine){
    if(wordCounts.find(word) != wordCounts.end()){
        wordCounts[word] = new map<int, int>();
        map<int, int> *countsForWord = wordCounts[word];
        *countsForWord[seqLine] = 1;
    }else{
        map<int, int> *countsForWord = wordCounts[word];
        if(countsForWord->find(word) != countsForWord->end()){
            *countsForWord[seqLine] = 1;
        }else{
            *countsForWord[seqLine] = *countsForWord[seqLine] + 1;
        }
    }
}

int MultiSeqProfile::wordCount(string word, int seqLine);
int MultiSeqProfile::getSeqLength(int index);
float MultiSeqProfile::compare(MultiSeqProfile* other, int seqIndex);