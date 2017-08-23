#ifndef _NUCL_SEQ_
#define _NUCL_SEQ_

#include <string>
#include <iostream>
#include <map>
#include <vector>
// basic file operations
#include <iostream>
#include <fstream>

#include "FileUtils.h"

using namespace std;

class NuclSeq{
public:
    //string fastaPath/*, profileFileName*/;
    //int startAtLineXInFasta;
    string *name;
    //int defaultWordSize = 9;
    
    NuclSeq(string *name, string *content);
    ~NuclSeq();
    string* getContent(/*bool keep = false*/);

    //float contains(NuclSeq &other/*, bool keep = false*/){
    //    return profile->compare(other.getProfile(/*keep*/));
    //}

    void free();
private:
    string* content;
    //bool alwaysStore;
};

#endif