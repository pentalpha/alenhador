#include "Alignment.h"

using namespace std;

string charDequeToString(deque<char> chars){
    string str = "";
    for(char c : chars){
        str += c;
    }
    return str;
}

Alignment::Alignment(string* a, string* b, string aName, string bName)
    : seqA(*a), seqB(*b), lenX(seqA.length()+1), lenY(seqB.length()+1)
{
    seqAName = aName;
    seqBName = bName;
    makeNuclSimMatrix();
    initCharIndex();
    allocateBacktrackMatrix();
    allocateSimMatrix();
    initMatrix();
}

Alignment::Alignment(NuclSeq* a, NuclSeq* b)
    : seqA(*(a->getContent())), seqB(*(b->getContent())), lenX(seqA.length()+1), lenY(seqB.length()+1)
{
    seqAName = *(a->name);
    seqBName = *(b->name);
    makeNuclSimMatrix();
    initCharIndex();
    allocateBacktrackMatrix();
    allocateSimMatrix();
    initMatrix();
}

void Alignment::initCharIndex(){
    for(int i = 0; i < 256; i++){
        charIndex[i] = -1;
    }
    charIndex[(int)'A'] = 0;
    charIndex[(int)'a'] = 0;
    charIndex[(int)'T'] = 1;
    charIndex[(int)'t'] = 1;
    charIndex[(int)'G'] = 2;
    charIndex[(int)'g'] = 2;
    charIndex[(int)'C'] = 3;
    charIndex[(int)'c'] = 3;
    charIndex[(int)'-'] = 4;
}

void Alignment::allocateBacktrackMatrix(){
    backtrack = new Point*[lenX];
    for(int i = 0; i < lenX; i++){
        backtrack[i] = new Point[lenY];
        for(int j = 0; j < lenY; j++){
            backtrack[i][j] = Point();
        }
    }
}

bool Alignment::elementCalculated(int x, int y){
        //cout << "Verifying if " << x << "_" << y << " has been calculated" << endl;
        //cout << backtrack[1][1].x << endl;
        //cout << backtrack[x-1][y-1].x << endl;
        return (backtrack[x][y].x != -1 && backtrack[x][y].y != -1);
}

void Alignment::allocateSimMatrix(){
    simMatrix = new int*[lenX];
    for(int i = 0; i < lenX; i++){
        simMatrix[i] = new int[lenY];
    }
}

void Alignment::initMatrix(){
    simMatrix[0][0] = 0;
    backtrack[0][0] = Point(0,0);
    for (int i = 1; i < lenX; i++){
        simMatrix[i][0] = i*sim((int)'-',(int)'-');
        backtrack[i][0] = Point(i-1,0);
    }
    for (int i = 1; i < lenY; i++){
        simMatrix[0][i] = i*sim((int)'-',(int)'-');
        backtrack[0][i] = Point(0, i-1);
    }
}

int Alignment::getValueOfXY(int x, int y){
    //cout << "Calculating value of " << x << "_" << y << endl;
    if(elementCalculated(x, y)){
        return simMatrix[x][y];
    }else{
        int res;
        int a = getValueOfXY(x-1,y-1) + sim((int)seqA[x-1], (int)seqB[y-1]);
        int b = getValueOfXY(x,y-1) + sim((int)seqA[x-1], (int)'-');
        int c = getValueOfXY(x-1,y) + sim((int)'-', (int)seqB[y-1]);
        if ((a > b) && (a > c)){
            res = a;
            backtrack[x][y].x = x-1;
            backtrack[x][y].y = y-1;
        }else if ((b > c) && (b > c)){
            res = b;
            backtrack[x][y].x = x;
            backtrack[x][y].y = y-1;
        }else{
            res = c;
            backtrack[x][y].x = x-1;
            backtrack[x][y].y = y;
        }
        simMatrix[x][y] = res;
        //cout << "Calculating value of " << x << "_" << y << ":" << endl;
        //cout << res << endl;
        return res;
    }
}

pair<string, string> Alignment::traceAlignment(){
    Point lastMovement(-1,-1);
    Point current(lenX-1, lenY-1);
    deque<char> traceA, traceB;

    while(current.x > 0 && current.y > 0){
        //cout << "Current element: " << current.x << "_" << current.y << endl;
        if(lastMovement.x == -1 && lastMovement.y == -1){
            traceA.push_front(seqA[current.x-1]);
            traceB.push_front(seqB[current.y-1]);
        }else if(lastMovement.x == -1 && lastMovement.y == 0){
            traceA.push_front(seqA[current.x-1]);
            traceB.push_front('-');
        }else{
            traceA.push_front('-');
            traceB.push_front(seqB[current.y-1]);
        }

        Point nextElement(backtrack[current.x][current.y].x, backtrack[current.x][current.y].y);
        lastMovement.x = nextElement.x-current.x;
        lastMovement.y = nextElement.y-current.y;
        current.x = nextElement.x;
        current.y = nextElement.y;
    }
    pair<string, string> aligns;
    aligns.first = charDequeToString(traceA);
    aligns.second = charDequeToString(traceB);

    return aligns;
}

void Alignment::align(){
    getValueOfXY(lenX-1, lenY-1);
    //cout << "Calculated simMatrix" << endl;
    alignment = traceAlignment();
    value = score(seqA, seqB);
    done = true;
}

void Alignment::startAlignment(){
    thread alignmentThread(&Alignment::align, this);
    done = false;
    alignmentThread.detach();
}

float Alignment::score(string seqA, string seqB){
    bool coveringGap = false;
    int points = 0;
    for (int i = 0; i < seqA.length() && i < seqB.length(); i++){
        char itemA = seqA[i];
        char itemB = seqB[i];
        int s = sim(itemA, itemB);
        if(!coveringGap && s == nuclSimMatrix[4][4]){
            coveringGap = true;
            points += gapStartPenalty;
        }else if(coveringGap && s != nuclSimMatrix[4][4]){
            coveringGap = false;
        }

        points += s;
    }
    return points;
}

bool Alignment::isDone(){
    return done;
}

float Alignment::getAlignmentValue(){
    return value;
}

float Alignment::getAlignmentRelativeValue(){
    return ((getAlignmentValue() / alignment.first.length()) / nuclSimMatrix[0][0]);
}

string Alignment::getAlignment(){
    int lineLength = 30;
    int nextPart = 0;
    string output = "";
    output += string("querySq = '") + seqAName + string("'\n");
    output += string("subject = '") + seqBName + string("'\n");
    output += string("Alignment Score = ") + to_string(getAlignmentValue())
            + string("\n");
    output += string("Alignment Relative Score = ") 
            + to_string(getAlignmentRelativeValue())
            + string("\n\n");
    string querySq = alignment.first;
    string subject = alignment.second;
    int notPrinted = querySq.length();
    while(notPrinted > 0){
        int start = nextPart;
        string qryPart = querySq.substr(start, lineLength);
        string sbjPart = subject.substr(start, lineLength);
        int end = start + qryPart.length();
        nextPart += qryPart.length();
        notPrinted -= qryPart.length();

        output += string("querySq\t") + to_string(start) + string ("\t");
        output += qryPart + string("\t") + to_string(end) + string("\n");
        output += string("subject\t") + to_string(start) + string ("\t");
        output += sbjPart + string("\t") + to_string(end) + string("\n\n");
    }
    return output;
}