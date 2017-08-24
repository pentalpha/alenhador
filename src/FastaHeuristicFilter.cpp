#include "FastaHeuristicFilter.h"

/*  SeqProfile queryProfile;
    NuclSeq *querySeq;
    string databasePath;
    int maxFiltered;

    mutex profileQueueLock;
    deque< pair<string*, string*> > toProfile;

    mutex sortQueueLock;
    deque<SeqNode> toSort;
    
    list<SeqNode> filteredNodes;*/

int FastaHeuristicFilter::defaultWordSize = 11;

FastaHeuristicFilter::FastaHeuristicFilter(NuclSeq *query, string dbPath, int max)
    : queryProfile(query->getContent(), defaultWordSize)
{
    querySeq = query;
    databasePath = dbPath;
    maxFiltered = max;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

list< pair<NuclSeq*, float> > FastaHeuristicFilter::justDoIt(){
    //thread readDBThread(readDB);
    readingDB = true;
    thread compareSeqsThread(&FastaHeuristicFilter::compareSeqs, this);
    thread sortSeqsThread(&FastaHeuristicFilter::sortSeqs, this);
    profiling = true;
    compareSeqsThread.detach();
    sorting = true;
    sortSeqsThread.detach();
    readDB();
    while(profiling || sorting){
        //Nothing to see here
    }
    //cout << "FINISHED THREADS" << endl;
    list< pair<NuclSeq*, float> > finalSeqs;
    for(SeqNode node : filteredNodes){
        pair<NuclSeq*, float> newEntry;
        newEntry.first = node.seq;
        newEntry.second = node.similarity;
        finalSeqs.push_back(newEntry);
    }

    return finalSeqs;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void FastaHeuristicFilter::readDB(){
    cout << "[Reading DB]" << endl;
    ifstream fastaFile;
    fastaFile.open(databasePath);
    
    if (fastaFile.is_open())
    {
        int nline = 1;
        int lastStart = -1;
        string *name;
        string *sequence;
        string line;
        while (getline(fastaFile,line))
        {   
            if(line.length() > 0){
                if(line[0] == '<' || line[0] == '>'){
                    if(lastStart > 0){
                        addToProfilingQueue(name, sequence);
                        sequence = NULL;
                        name = NULL;
                    }
                    name = new string(line);
                    lastStart = nline;
                    sequence = new string("");
                }else if(sequence != NULL){
                    *sequence += line;
                }
            }
            nline++;
        }
        if(sequence != NULL){
            addToProfilingQueue(name, sequence);
        }
        fastaFile.close();
    }else{
        cout << "IO ERROR: Fasta DB file could not be read!" << endl;
        readingDB = false;
        return;
    }
    
    readingDB = false;
    cout << "[Finished reading DB]" << endl;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void FastaHeuristicFilter::addToProfilingQueue(string* name, string* content){
    pair<string*,string*> newProfiling;
    newProfiling.first = name;
    newProfiling.second = content;
    std::lock_guard<std::mutex> lock(profileQueueLock);
    toProfile.push_back(newProfiling);
    seqsReadFromDB++;
}

pair<string*, string*> FastaHeuristicFilter::getRawSeq(){
    std::lock_guard<std::mutex> lock(profileQueueLock);
    if(toProfile.size() > 0){
        pair<string*, string*> rawSeq = toProfile.front();
        //cout << *rawSeq.second << endl;
        toProfile.pop_front();
        return rawSeq;
    }else{
        cout << "There is no sequence in profiling queue" << endl;
    }
}

int FastaHeuristicFilter::seqsToProfile(){
    std::lock_guard<std::mutex> lock(profileQueueLock);
    return toProfile.size();
}

void FastaHeuristicFilter::compareSeqs(){
    cout << "[Ready to make sequence profiles]" << endl;
    while(readingDB || seqsToProfile() > 0){
        if(seqsToProfile() > 0){
            pair<string*, string*> rawSeq = getRawSeq();
            NuclSeq *seq = new NuclSeq(rawSeq.first, rawSeq.second);
            SeqProfile profile(rawSeq.second, defaultWordSize);
            SeqNode node;
            node.seq = seq;
            node.similarity = profile.compare(&queryProfile);
            addToSortingQueue(node);
        }
    }
    profiling = false;
    cout << "[Finished profiling sequences]" << endl;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void FastaHeuristicFilter::addToSortingQueue(SeqNode newNode){
    std::lock_guard<std::mutex> lock(sortQueueLock);
    toSort.push_back(newNode);
}

SeqNode FastaHeuristicFilter::getUnsortedNode(){
    std::lock_guard<std::mutex> lock(sortQueueLock);
    SeqNode node = toSort.front();
    toSort.pop_front();
    return node;
}

int FastaHeuristicFilter::seqsToSort(){
    std::lock_guard<std::mutex> lock(sortQueueLock);
    return toSort.size();
}

void FastaHeuristicFilter::sortSeqs(){
    cout << "[Ready to filter sequences by profile]" << endl;
    
    while(profiling || seqsToSort() > 0){
        if(seqsToSort() > 0){
            SeqNode node = getUnsortedNode();
            int currentFiltered = filteredNodes.size();
            if(currentFiltered > 0){
                if(node.similarity <= filteredNodes.back().similarity){
                    delete node.seq;
                    node.seq = NULL;
                    filtered++;
                }else{
                    list<SeqNode>::iterator index = filteredNodes.begin();
                    for(SeqNode filteredNode : filteredNodes){
                        if(filteredNode.similarity >= node.similarity){
                            index++;
                        }else{
                            break;
                        }
                    }
                    filteredNodes.insert(index, node);
                    if(filteredNodes.size() > maxFiltered){
                        delete filteredNodes.back().seq;
                        filteredNodes.pop_back();
                        filtered++;
                    }
                }
            }else{
                filteredNodes.push_front(node);
            }
        }
    }
    sorting = false;
    cout << string("[Finished filtering sequences]\n\tUnfiltered: ") + to_string(filteredNodes.size())
        + string(", filtered: ") + to_string(filtered) << endl;
}