#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char *argv[]){
    //read from a file
    ifstream inFile(argv[1]);

    ofstream outFile;
    outFile.open(argv[2]);

    return 0;
}