#include <iostream>
#include <fstream>
#include <math.h>
#include <tgmath.h>
//#include <bits/stdc++.h>

using namespace std;

class BiMean {
    public:
        int numRows, numCols, minVal, maxVal, maxHeight;
        double offSet, dividePt;
        int *histAry;
        int *GaussAry;
        char **histGraph;
        char **GaussGraph;
        int *smoothHist;

    BiMean(int rows, int cols, int minv, int maxv){
        numRows = rows;
        numCols = cols;
        minVal = minv;
        maxVal = maxv;
        histAry = new int[maxv+1]();
        GaussAry = new int[maxv+1]();
        smoothHist = new int[maxv+1]();
        offSet = (maxv-minv)/10;
        dividePt = (maxv-minv)/10;
    }

    ~BiMean(){
        delete[] histAry;
        delete[] GaussAry;
        for(int i=0; i<maxHeight+1; i++){
            delete[] histGraph[i];
            delete[] GaussGraph[i];
        }
        delete[] histGraph;
        delete[] GaussGraph;
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
        return retVal;
    }

    void allocateArrs(){
        histGraph = new char*[maxHeight+1];
        GaussGraph = new char*[maxHeight+1];
        for(int i=0; i<maxHeight; i++){
            histGraph[i] = new char[maxVal+1]();
            GaussGraph[i] = new char[maxVal+1]();
        }
    }

    void plotGraph(int *ary, char **graph, char symbol){
        //maps 1D array onto 2D array with symbol
        //if ary[i] > 0 then graph[i, ary[i]] <- symbol
        //symbol will be * for histGraph and + for GaussGraph
        for(int i = 0; i<=maxVal; i++){
            if(ary[i] > 0){
                int point = maxHeight-ary[i];
                for(int j=point; j<maxHeight; j++){
                    graph[j][i] = symbol;
                }
            }
        }
    }

    void addVertical(int thr){
        //add histGraph[thr, j] with | where j=0 to maxHeight
        for(int i=0; i<maxHeight; i++){
            histGraph[i][thr] = '|';
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
            sq = (double)(i - x) * (double)(i - x);
            sum += ((double)histAry[i] * sq);
            numPixels += histAry[i];
        }

        return sum/numPixels;
    }

    double modifiedGauss(int x, int mean, int var, int maxheight){
        //instructions at specs
        //NOT DONE
        double n = (x-mean) * (x-mean);
        double e = (-(n/(2*var)));
        return double(maxheight * exp(e));
    }

    int* setZero(int *arr){
        //set array to zero
        for(int i=0; i<maxVal; i++){
            arr[i] = 0;
        }
        return arr;
    }

    double biMeanGauss(int thr, ofstream& outFile){
        //determines the best threshold selection via fitGauss method
        //ALGO IN SPECS
        double sum1, sum2, total;
        double bestThr = thr;
        double minSumDiff = 999999.0;

        while(thr  < (maxVal - offSet)){
            GaussAry = setZero(GaussAry);
            sum1 = fitGauss(0, thr, GaussAry);
            sum2 = fitGauss(thr, maxVal, GaussAry);
            total = sum1+sum2;
            if(total < minSumDiff){
                minSumDiff = total;
                bestThr = thr;
            }
            outFile << thr << " " << sum1 << " " << sum2 << " " << total << " " << minSumDiff << " " << bestThr << endl;
            thr++;
        }
        return bestThr;
    }

    double fitGauss(double leftIndex, double rightIndex, int* arr){
        //ALGO IN SPECS
        double mean = computeMean(leftIndex, rightIndex, maxHeight);
        double var = computeVar(leftIndex, rightIndex, mean);
        double sum = 0.0;
        double Gval;
        for(int i=leftIndex; i<=rightIndex; i++){
            Gval = modifiedGauss(i, mean, var, maxHeight);
            sum += abs(Gval-(double)histAry[i]);
            arr[i] = (int)Gval;
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

    void makeBorder(char symbol, ofstream& outfile){
        //border for output file
        for(int i=0; i<maxVal*2; i++){
            outfile << symbol;
        }
        outfile << endl;
    }

    void plotAll(ofstream& outFile){
        //overlay histGraph with add vertical line and GaussGraph to File1
        for(int i = 0; i<maxHeight; i++){
            for(int j=0; j<maxVal; j++){
                if((int)histGraph[i][j]=='|'){
                    outFile << histGraph[i][j] << " ";
                }
                else if(((int)histGraph[i][j]>0 && (int)GaussGraph[i][j]>0) && (int)histGraph[i][j]!='|'){
                    outFile << 'x' << " ";
                }else{
                    if((int)histGraph[i][j]>0){
                        outFile << histGraph[i][j] << " ";
                    }
                    if((int)GaussGraph[i][j]>0){
                        outFile << GaussGraph[i][j] << " ";
                    }
                }
                if((int)histGraph[i][j]==0 && (int)GaussGraph[i][j]==0){
                    outFile << "  ";
                }
            }
            outFile << endl;
        }
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
    outFile1 << "Histogram Curve" << endl;
    for(int i = 0; i<biMean.maxHeight; i++){
        for(int j=0; j<biMean.maxVal; j++){
            if((int)biMean.histGraph[i][j]==0) outFile1 << "  ";
            else outFile1 << biMean.histGraph[i][j] << " ";
        }
        outFile1 << endl;
    }
    
    biMean.makeBorder('-', outFile1);
    
    //Step 5
    //Done on class BiMean construction
    //offSet and dividePt value assignment
    
    //Step 6
    double bestThrVal = biMean.biMeanGauss(biMean.dividePt, outFile2);
    
    //Step 7
    biMean.bestFitGauss(bestThrVal, biMean.GaussAry);

    //Step 8
    biMean.plotGraph(biMean.GaussAry, biMean.GaussGraph, '+');
    outFile1 << "Gaussian Curve" << endl;
    for(int i = 0; i<biMean.maxHeight; i++){
        for(int j=0; j<biMean.maxVal; j++){
            if((int)biMean.GaussGraph[i][j]==0) outFile1 << "  ";
            else outFile1 << biMean.GaussGraph[i][j] << " ";
        }
        outFile1 << endl;
    }

    biMean.makeBorder('-', outFile1);

    //Step 9
    outFile1 << "Best Threshol Value: " << bestThrVal << endl;

    biMean.makeBorder('-',outFile1);

    outFile1 << "Histogram curve with a vertical line at the selected threshold" << endl;

    //Step 10
    biMean.addVertical(bestThrVal);
    for(int i = 0; i<biMean.maxHeight; i++){
        for(int j=0; j<biMean.maxVal; j++){
            if((int)biMean.histGraph[i][j]==0) outFile1 << "  ";
            else outFile1 << biMean.histGraph[i][j] << " ";
        }
        outFile1 << endl;
    }

    biMean.makeBorder('-',outFile1);

    outFile1 << "Histogram Curve and two best-fitted Gaussian curves" << endl;

    //Step 11
    biMean.plotAll(outFile1);


    //Step 12
    inFile.close();
    outFile1.close();
    outFile2.close();
    return 0;
}