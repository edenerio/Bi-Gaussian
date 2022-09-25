class BiMean {
    public:
        int numRows, numCols, minVal, maxVal, maxHeight, maxGVal;
        int offSet, dividePt;
    private:
        int *histAry;
        int *GaussAry;
        char **histGraph;
        char **GaussGraph;

    BiMean(int rows, int cols, int minval, int maxval, int maxheight, int maxgval){
        numRows = rows;
        numCols = cols;
        minVal = minval;
        maxVal = maxval;
        maxHeight = maxheight;
        maxGVal = maxgval;
        histAry[maxVal+1] = {};
        GaussAry[maxVal+1] = {};
        histGraph[maxVal+1][maxHeight+1] = {};
        GaussGraph[maxVal+1][maxHeight+1] = {};
        offSet = (maxVal-minVal)/10;
    }

    int loadHist(int *in){
        histAry = in;
        int len = *(&histAry + 1) - histAry;
        int retVal = histAry[0];
        for(int i=0; i<len; i++){
            if(retVal<histAry[i]){
                retVal = histAry[i];
            }
        }
        return retVal;
    }

    void plotGraph(int *ary, int **graph, char symbol){
        //maps 1D array onto 2D array with symbol
        //if ary[i] > 0 then graph[i, ary[i]] <- symbol
        //symbol will be * for histGraph and + for GaussGraph
    }

    void addVertical(){
        //add histGraph[thr, j] with | where j=0 to maxHeight
    }

    double computeMean(int leftIndex, int rightIndex){
        int weight = 0;

        //computes the histogiven portion of from given leftIndex to rightIndex of the histogram
        //and returns the weighted average of the histogram
        //ALGO IN SPECS

        return weight;
    }

    double computeVar(int leftIndex, int rightIndex){
        int weight = 0;

        //weight = weight of variance
        //ALGO IN SPECS

        return weight;
    }

    double modifiedGauss(int x, int mean, int var, int maxHeight){
        //instructions at specs
    }

    void setZero(int *ary){
        //set ary to zero
    }

    int biMeanGauss(int thr){
        //determines the best threshold selection via fitGauss method
        //ALGO IN SPECS
    }

    void fitGauss(){
        //ALGO IN SPECS
    }

    void bestFitGauss(){
        //ALGO IN SPECS
    }

    void plotAll(){
        //overlay histGraph with add vertical line and GaussGraph to File1
    }
};