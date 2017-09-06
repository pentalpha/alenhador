#ifndef _ALIGNMENT_
#define _ALIGNMENT_

#include <map>
#include <string>
#include <iostream>
#include <deque>
#include <thread>
#include "NuclSeq.h"

using namespace std;

string charDequeToString(deque<char> chars);

class Point{
public:
    int x, y;
    Point(int a, int b)
        : x(a), y(b)
    {
    }

    Point()
        : x(-1), y(-1)
    {
    }
};

class Alignment{
public:
    Alignment(string* a, string* b, string aName, string bName);
    Alignment(NuclSeq* a, NuclSeq* b);
    
    float getAlignmentValue();
    float getAlignmentRelativeValue();
    string getAlignment();
    void startAlignment();
    bool isDone();
    float score(string seqA, string seqB);
private:
    void initCharIndex();
    void align();
    inline void makeNuclSimMatrix(){
        //nuclSimMatrix
    }
    inline short sim(int a, int b){
        return nuclSimMatrix[charIndex[a]][charIndex[b]];
    }
    inline bool elementCalculated(int x, int y);
    int getValueOfXY(int x, int y);

    void allocateBacktrackMatrix();
    void allocateSimMatrix();
    void initMatrix();

    pair<string, string> traceAlignment();
    pair<string, string> alignment;
    float value = 0.0;

    string seqA, seqB;
    string seqAName, seqBName;
    int lenX, lenY;
    //map<char, int> charIndex;
    int charIndex[256];
    short nuclSimMatrix[5][5]= {{5,3,2,4,-2},
                        {3,5,2,2,-2},
                        {2,2,5,1,-2},
                        {4,2,1,5,-2},
                        {-2,-2,-2,-2,-2}};
    int gapStartPenalty = -3;
    int** simMatrix;
    Point** backtrack;

    bool done = false;
};

#endif