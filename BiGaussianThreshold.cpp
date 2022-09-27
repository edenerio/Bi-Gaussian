#include <iostream>
#include <fstream>
#include <math.h>
#include <tgmath.h>

using namespace std;

class BiMean {
    public:
        int numRows, numCols, minVal, maxVal, maxHeight;
        double offSet, dividePt;
        int *histAry;
        int *GaussAry;
        char **histGraph;
        char **GaussGraph;

    BiMean(int rows, int cols, int minv, int maxv){
        numRows = rows;
        numCols = cols;
        minVal = minv;
        maxVal = maxv;
        histAry[maxv+1] = {};
        offSet = (maxv-minv)/10;
        dividePt = offSet;
    }

    int loadHist(ifstream& in){
        int retVal = 0;
        while(!in.eof()){
            int x;
            int y;
            in >> x;
            in >> y;
            histAry[x] = y;
            if(retVal < y){
                retVal = y;
            }
        }
        cout << retVal << endl;
        return retVal;
    }

    void allocateArrs(){
        histGraph[maxVal+1][maxHeight+1];
        GaussGraph[maxVal+1][maxHeight+1];
        GaussAry[maxVal+1];
    }

    void plotGraph(int *ary, char **graph, char symbol){
        //maps 1D array onto 2D array with symbol
        //if ary[i] > 0 then graph[i, ary[i]] <- symbol
        //symbol will be * for histGraph and + for GaussGraph
        /*
            FIX THIS FUNCTION!!
        */
        for(int i = 0; i<numRows; i++){
            if(ary[i] > 0){
                graph[i][ary[i]] = symbol;
            }else{
                graph[i][ary[i]] = ' ';
            }
        }
        for(int i = 0; i<numRows; i++){
            for(int j=0; j<numCols; j++){
                cout << graph[i][j];
            }
            cout << endl;
        }
    }

    void addVertical(int thr){
        //add histGraph[thr, j] with | where j=0 to maxHeight
        for(int i=0; i<=maxHeight; i++){
            histGraph[thr][i] = '|';
        }
    }

    double computeMean(int leftIndex, int rightIndex, int x){
        //computes the histogiven portion of from given leftIndex to rightIndex of the histogram
        //and returns the weighted average of the histogram
        //ALGO IN SPECS
        double sum = 0;
        x = 0;
        int numPixels = 0;

        for(int i=leftIndex; i<rightIndex; i++){
            sum += (histAry[i] * i);
            numPixels += histAry[i];
            if(histAry[i] > x){
                x = histAry[i];
            }
        }

        return sum/numPixels;
    }

    double computeVar(int leftIndex, int rightIndex, int x){
        //weight = weight of variance
        //ALGO IN SPECS
        double sum = 0.0;
        int numPixels = 0;
        double sq;
        
        for(int i=leftIndex; i<rightIndex; i++){
            sq = double(i - x) * double(i - x);
            sum += (double(histAry[i]) * sq);
            numPixels += histAry[i];
        }

        return sum/numPixels;
    }

    double modifiedGauss(int x, int mean, int var, int maxheight){
        //instructions at specs
        //NOT DONE
        double n = (x-mean) * (x-mean);
        double c = var * var;
        double e = (-(n/(2*c)));
        return double(maxheight * exp(e));
    }

    int* setZero(int *arr){
        //set array to zero
        int len = *(&arr + 1) - arr;
        for(int i=0; i<len; i++){
            arr[i] = 0;
        }
        return arr;
    }

    double biMeanGauss(int thr, ofstream& outFile){
        //determines the best threshold selection via fitGauss method
        //ALGO IN SPECS
        double sum1, sum2;
        double total = 0;
        double bestThr = dividePt;
        double minSumDiff = 999999.0;

        while(dividePt  < (maxVal - offSet)){
            GaussAry = setZero(GaussAry);
            sum1 = fitGauss(0, dividePt, GaussAry);
            sum2 = fitGauss(dividePt, maxVal, GaussAry);
            total = sum1+sum2;
            if(total < minSumDiff){
                minSumDiff = total;
                bestThr = dividePt;
            }
            outFile << dividePt << " " << sum1 << " " << sum2 << " " << total << " " << minSumDiff << " " << bestThr << endl;
            dividePt++;
        }
        return bestThr;
    }

    double fitGauss(double leftIndex, double rightIndex, int* arr){
        //ALGO IN SPECS
        double mean = computeMean(leftIndex, rightIndex, maxHeight);
        double var = computeVar(leftIndex, rightIndex, mean);
        double sum = 0.0;
        double Gval;
        double maxGval = 0;
        for(int i=leftIndex; i<=rightIndex; i++){
            Gval = modifiedGauss(i, mean, var, maxHeight);
            sum += abs(Gval-(double)histAry[i]);
            GaussAry[i] = int(Gval);
        }
        return sum;
    }

    void bestFitGauss(double bestThr, int* arr){
        //ALGO IN SPECS
        double sum1, sum2;
        arr = setZero(arr);
        sum1 = fitGauss(0, bestThr, arr);
        sum2 = fitGauss(bestThr, maxVal, arr);
    }

    void plotAll(){
        //overlay histGraph with add vertical line and GaussGraph to File1
    }
};

int main(int argc, char *argv[]){
    //read from a file
    if(argc != 4){
        cout << "Invalid arguments" << endl;
        return 1;
    }

    //ifstream inFile(argv[1]);
    ifstream inFile;
    inFile.open(argv[1]);

    ofstream outFile1;
    outFile1.open(argv[2]);

    ofstream outFile2;
    outFile2.open(argv[3]);

    int row, col, minval, maxval;

    if(inFile.is_open()){
        inFile >> row;
        inFile >> col;
        inFile >> minval;
        inFile >> maxval;
    }
    else{
        cout << "Input File not found!" << endl;
        return 2;
    }
    
    BiMean biMean = BiMean(row,col,minval,maxval);
    biMean.maxHeight = biMean.loadHist(inFile);
    
    //Step 3
    biMean.allocateArrs();

    //Step 4
    biMean.plotGraph(biMean.histAry, biMean.histGraph, '*');
    outFile1 << biMean.histGraph;


    //Step 5
    //Done on class BiMean construction
    //offSet and dividePt value assignment

    //Step 6
    double bestThrVal = biMean.biMeanGauss(biMean.dividePt, outFile2);

    //Step 7
    biMean.bestFitGauss(bestThrVal, biMean.GaussAry);

    //Step 8
    biMean.plotGraph(biMean.GaussAry, biMean.GaussGraph, '+');
    outFile1 << biMean.GaussGraph;

    //Step 9
    outFile1 << "Best Threshol Value: " << bestThrVal;

    //Step 10
    biMean.addVertical(bestThrVal);
    outFile1 << biMean.histGraph;

    //Step 11
    //outFile1 << biMean.plotAll()

    //Step 12
    inFile.close();
    outFile1.close();
    outFile2.close();
    return 0;
}