#ifndef _ALIGNMENT_
#define _ALIGNMENT_

#include <map>
#include <string>
#include <iostream>
#include <deque>
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
    string seqA, seqB;
    int lenX, lenY;
    //map<char, int> charIndex;
    int charIndex[256];
    short nuclSimMatrix[5][5]= {{5,3,2,4,-2},
                        {3,5,2,2,-2},
                        {2,2,5,1,-2},
                        {4,2,1,5,-2},
                        {-2,-2,-2,-2,-2}};
    int** simMatrix;
    Point** backtrack;

    Alignment(string a, string b);
    Alignment(NuclSeq a, NuclSeq b);

    inline bool elementCalculated(int x, int y);
    int getValueOfXY(int x, int y);
    pair<string, string> align();
private:
    void initCharIndex();

    inline void makeNuclSimMatrix(){
        //nuclSimMatrix
    }

    inline short sim(int a, int b){
        return nuclSimMatrix[charIndex[a]][charIndex[b]];
    }

    void allocateBacktrackMatrix();
    void allocateSimMatrix();
    void initMatrix();
    pair<string, string> traceAlignment();
};

#endif