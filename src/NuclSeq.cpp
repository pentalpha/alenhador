
#include "NuclSeq.h"

NuclSeq::NuclSeq(string *seqName, string *seqContent){
    name = seqName;
    content = seqContent;
}

NuclSeq::~NuclSeq(){
    free();
}

string* NuclSeq::getContent(/*bool keep*/){
    if(content != NULL){
        return content;
    }else{
        return NULL;
    }
}

void NuclSeq::free(){
    if(content != NULL){
        delete content;
        content = NULL;
    }
    if(name != NULL){
        delete name;
        name = NULL;
    }
}


/*NuclSeq::NuclSeq(string path, int startLine, bool keep, bool store)
    : fastaPath(path), profileFileName(path + string(".seqProf")), startAtLineXInFasta(startLine),
    alwaysStore(store)
{
    content = NULL;
    profile = NULL;
    if(startLine != 1 || alwaysStore){
        createFiles(keep);
    }
}*/

/*NuclSeq::NuclSeq(string path, int startLine)
    : fastaPath(path), startAtLineXInFasta(startLine)
{
    content = NULL;
    profile = NULL;
    readFromFile();
    //if(startLine != 1 || alwaysStore){
    //    createFiles(keep);
    //}
}

void NuclSeq::readFromFile(){
    ifstream fastaFile;
    fastaFile.open(fastaPath);
    int nline = 1;
    if (fastaFile.is_open())
    {
        string* contentTmp = new string("");
        string line;
        while (getline(fastaFile,line))
        {   
            if(nline > startAtLineXInFasta){
                if(line.length() > 0){
                    if(line[0] == '<' || line[0] == '>'){
                        break;
                    }
                }
                *contentTmp += line;
            }else if (nline == startAtLineXInFasta){
                name = line;
            }
            nline++;
        }
        fastaFile.close();
        content = contentTmp;
        profile = new SeqProfile(*content, defaultWordSize);
    }
}

NuclSeq::NuclSeq(string path, int startLine, string seqName, string *sequence)
    : fastaPath(path), startAtLineXInFasta(startLine), name(seqName)
{
    content = sequence;
    profile = new SeqProfile(*content, defaultWordSize);
}

/*void NuclSeq::createFiles(bool keep){
    string fileName = string("tmp/") + nameInPath(fastaPath) + to_string(startAtLineXInFasta);
    if (!fileExists(fileName.c_str())){
        getContent(true);
        ofstream fastaFile;
        fastaFile.open(fileName);
        if (fastaFile.is_open())
        {
            fastaFile << name << endl;
            fastaFile << *content << endl;
            fastaFile.close();
        }
        
    }
    fastaPath = fileName;

    profileFileName = fastaPath + string(".seqProf");
    if (!fileExists(profileFileName.c_str())){
        profile = new SeqProfile(*getContent(keep), defaultWordSize);
        profile->toFile(profileFileName, name);
    }
    if(!keep){
        free();
    }
    startAtLineXInFasta = 1;
}*/

/*string* NuclSeq::getContent(bool keep){
    /*if(content != NULL){
        return content;
    }else{
        return NULL;
    }

    ifstream fastaFile;
    fastaFile.open(fastaPath);
    int nline = 1;
    if (fastaFile.is_open())
    {
        string* contentTmp = new string("");
        string line;
        while (getline(fastaFile,line))
        {   
            if(nline > startAtLineXInFasta){
                if(line.length() > 0){
                    if(line[0] == '<' || line[0] == '>'){
                        break;
                    }
                }
                *contentTmp += line;
            }else if (nline == startAtLineXInFasta){
                name = line;
            }
            nline++;
        }
        fastaFile.close();
        if(keep){
            content = contentTmp;
        }
        return contentTmp;
    }else{
        return NULL;
    }
}

/*SeqProfile* NuclSeq::getProfile(bool keep){
    //if(profile != NULL){
    //    return profile;
    //}else{
    //    SeqProfile* profileTmp;
        //if(fileExists(profileFileName.c_str())){
        //    profileTmp = new SeqProfile(profileFileName);
        //}else{
    //        profileTmp = new SeqProfile(*getContent()(keep), defaultWordSize);
        //}
         
        //if(keep){
    //        profile = profileTmp;
        //}
    //    return profile;
    //}
//}*/