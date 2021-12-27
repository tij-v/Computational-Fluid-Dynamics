//Name: Shitij Vishwakarma
//Student ID: 902294581
//AEE 5150: CFD
//Project 1: Question 1

//Libraries used in this program
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//Global Definitions
#define thermalDiffusivity 0.1 //ft^2/hour
#define wallThickness 1 //ft
#define initialTemp 100 //degree F
#define finalTemp 300 //degree F

//function declarations
void InitializeArray(double* tempArray);
void FTSCExplicit(double diffusionNumber, double* currWalltemp, double* initWalltemp);
void DuFortFrankel(double diffusionNumber, double* currWalltemp, double* initWalltemp, int n, int caseNum);
void ThomasAlgorithm(double* above, double* below, double* constant, double* diagonal, double* tempHolder);
void FTSCImplicit(double diffusionNumber, double* currWalltemp);
void CrankNicolson(double diffusionNumber, double* currWalltemp);
void printData(int method, double** outputTemp, double caseNumber);

//main subroutine
int main(void) {
    double timeStepOne = 0.01;
    double timeStepTwo = 0.05;
    double runTime = 0.4;
    double deltaX = 0.05; //ft
    double diffusionNumber;

    double* currWalltemp = (double*)malloc(100 * sizeof(double));//stores the nth iteration values
    double* initWalltemp = (double*)malloc(100 * sizeof(double)); //stores the (n+1)th iteration values
    double** outputTemp = (double**)malloc(100 * sizeof(double*));//2D pointer array to point to the results needed at
                                                                  //intervals of 0.1, 0.2, 0.3, 0.4
    double caseArray[2] = { 0.01, 0.05 };

    int iterationCase = 0;
    int fileCounter = 0;
    //initialization the temperature arrays at n = 1;
    currWalltemp[0] = finalTemp;
    currWalltemp[20] = finalTemp;

    for (int i = 0; i < 100; i++) {
        outputTemp[i] = (double*)malloc(5 * sizeof(double));
    }
    for (int caseNum = 0; caseNum < 2; caseNum++) {
        diffusionNumber = ((thermalDiffusivity * caseArray[caseNum]) / (deltaX * deltaX));
        for (int method = 0; method < 4; method++) {
            fileCounter = fileCounter + 1;
            InitializeArray(currWalltemp);
            iterationCase = runTime / caseArray[caseNum];
            for (int n = 1; n < iterationCase + 1; n++) {
                if (method == 0) {
                    FTSCExplicit(diffusionNumber, currWalltemp, initWalltemp);
                }
                else if (method == 1) {
                    DuFortFrankel(diffusionNumber, currWalltemp, initWalltemp, n, caseNum);
                }
                else if (method == 2) {
                    FTSCImplicit(diffusionNumber, currWalltemp);
                }
                else if (method == 3) {
                    CrankNicolson(diffusionNumber, currWalltemp);
                }
                if (n % 10 == 0 && caseNum == 0) {
                    int j = (n / 10) - 1;
                    for (int i = 0; i < 21; i++) {
                        outputTemp[i][j] = currWalltemp[i];
                    }
                    j = 0;
                }
                if (n % 2 == 0 && caseNum == 1) {
                    int j = (n / 2) - 1;
                    for (int i = 0; i < 21; i++) {
                        outputTemp[i][j] = currWalltemp[i];
                    }
                    j = 0;
                }
            }
            printData(method, outputTemp, caseArray[caseNum]);
        }
    }
    for (int i = 0; i < 100; i++) {
        free(outputTemp[i]);
    }
    if (outputTemp) {
        outputTemp = NULL;
    }
    if (currWalltemp) {
        currWalltemp = NULL;
        free(currWalltemp);
    }
    if (initWalltemp) {
        initWalltemp = NULL;
        free(initWalltemp);
    }
    if (outputTemp) {
        outputTemp = NULL;
        free(outputTemp);
    }
}
//initialize array initializes the current array to the inital temperature conditions
void InitializeArray(double* tempArray) {
    for (int i = 1; i < 20; i++) {
        tempArray[i] = initialTemp;
    }
}
//This function solves for the FTSC scheme 
void FTSCExplicit(double diffusionNumber, double* currWalltemp, double* initWalltemp) {
    for (int i = 0; i < 21; i++) {
        initWalltemp[i] = currWalltemp[i];
    }
    for (int i = 1; i < 20; i++) {
        currWalltemp[i] = initWalltemp[i] + (diffusionNumber * (initWalltemp[i + 1] - (2 * initWalltemp[i]) + initWalltemp[i - 1]));
    }
}
//this function solves for the DFF scheme
void DuFortFrankel(double diffusionNumber, double* currWalltemp, double* initWalltemp, int n, int caseNum) {
    double* tempHolder = (double*)malloc(100 * sizeof(double)); //temporary array to hold the (n+1)th value while the other arrays are updated
    //since the DFF scheme uses two initial starting conditions we need to use FTSC to get those, since FTSC is a one-step scheme. 
    if (n == 1) {
        for (int i = 0; i < 21; i++) {
            initWalltemp[i] = currWalltemp[i];
        }
        if (caseNum == 0) {
            FTSCExplicit(diffusionNumber, currWalltemp, initWalltemp);
        }
        else if (caseNum == 1) {
            FTSCExplicit(diffusionNumber, currWalltemp, initWalltemp);
        }
    }
    //at this level we should have two matrices that have the initial values at (n)th and (n+1)th time levels. 
    //we can now start using the DFF scheme to calcualte any further results.
    else if (n > 1) {
        for (int i = 1; i < 20; i++) {
            tempHolder[i] = (((1 - (2 * diffusionNumber)) * initWalltemp[i]) + ((2 * diffusionNumber) * (currWalltemp[i + 1] + currWalltemp[i - 1]))) / (1 + (2 * diffusionNumber));
        }
        for (int i = 1; i < 20; i++) {
            initWalltemp[i] = currWalltemp[i];
            currWalltemp[i] = tempHolder[i];
        }
    }
    if (tempHolder) {
        tempHolder = NULL;
        free(tempHolder);
    }
}
//this function incorporates the thomas algorithm with other functions that might need it. 
void ThomasAlgorithm(double* above, double* below, double* constant, double* diagonal, double* tempHolder) {
    for (int i = 1; i < 19; i++) {
        diagonal[i] = diagonal[i] - ((above[i - 1] * below[i - 1]) / diagonal[i - 1]);
        constant[i] = constant[i] - ((below[i - 1] * constant[i - 1]) / diagonal[i - 1]);
    }
    tempHolder[18] = constant[18] / diagonal[18];
    for (int i = 17; i > -1; i--) {
        tempHolder[i] = (constant[i] - (above[i] * tempHolder[i + 1])) / diagonal[i];
    }
}
//this function solves for the FTSC scheme implicitly. it uses the thomas algorithm
void FTSCImplicit(double diffusionNumber, double* currWalltemp) {
    double* above = (double*)malloc(100 * sizeof(double));
    double* below = (double*)malloc(100 * sizeof(double));
    double* constant = (double*)malloc(100 * sizeof(double));
    double* diagonal = (double*)malloc(100 * sizeof(double));
    double* tempHolder = (double*)malloc(100 * sizeof(double));
    for (int i = 0; i < 18; i++) {
        above[i] = (-1) * diffusionNumber;
        below[i] = (-1) * diffusionNumber;
    }
    for (int i = 0; i < 19; i++) {
        diagonal[i] = 1 + (2 * diffusionNumber);
        constant[i] = currWalltemp[i + 1];
    }
    constant[0] = constant[0] + (diffusionNumber * currWalltemp[0]);
    constant[18] = constant[18] + (diffusionNumber * currWalltemp[20]);

    ThomasAlgorithm(above, below, constant, diagonal, tempHolder);
    for (int i = 0; i < 19; i++) {
        currWalltemp[i + 1] = tempHolder[i];
    }
}
//this function solves for the Cranck-Nicolson scheme, also utilizing the thomas algorithm.
void CrankNicolson(double diffusionNumber, double* currWalltemp) {
    double* above = (double*)malloc(100 * sizeof(double));
    double* below = (double*)malloc(100 * sizeof(double));
    double* constant = (double*)malloc(100 * sizeof(double));
    double* diagonal = (double*)malloc(100 * sizeof(double));
    double* tempHolder = (double*)malloc(100 * sizeof(double));
    for (int i = 0; i < 18; i++) {
        above[i] = diffusionNumber / (-2);
        below[i] = diffusionNumber / (-2);
    }
    for (int i = 0; i < 19; i++) {
        diagonal[i] = (1 + diffusionNumber);
        constant[i] = ((diffusionNumber / 2) * (currWalltemp[i + 2] + currWalltemp[i])) + ((1 - diffusionNumber) * currWalltemp[i + 1]);
    }
    constant[0] = constant[0] + ((diffusionNumber / 2) * currWalltemp[0]);
    constant[18] = constant[18] + ((diffusionNumber / 2) * currWalltemp[20]);
    ThomasAlgorithm(above, below, constant, diagonal, tempHolder);
    for (int i = 0; i < 19; i++) {
        currWalltemp[i + 1] = tempHolder[i];
    }
}
//final function to print the data onto a text file. 
void printData(int method, double** outputTemp, double caseNumber) {
    FILE* filePtr;
    filePtr = fopen("Question1_Output.txt", "a");
    if (method == 0) {
        if (method == 0 && caseNumber == 0.05) {
            fprintf(filePtr, "FTSC Explicit Scheme @delta_t = %.2f\n@t =\t   0.1\t    0.2\t     0.3      0.4\n", caseNumber);
        }
        else {
            fprintf(filePtr, "FTSC Explicit Scheme @delta_t = %.2f\n@t =\t  0.1\t   0.2\t  0.3\t  0.4\n", caseNumber);
        }
    }
    else if (method == 1) {
        fprintf(filePtr, "DuFort-Frankel Scheme @delta_t = %.2f\n@t =\t 0.1\t  0.2\t  0.3\t  0.4\n", caseNumber);
    }
    else if (method == 2) {
        fprintf(filePtr, "FTSC Implicit Scheme @delta_t = %.2f\n@t =\t 0.1\t  0.2\t  0.3\t  0.4\n", caseNumber);
    }
    else if (method == 3) {
        fprintf(filePtr, "Crank-Nicolson Scheme @delta_t = %.2f\n@t =\t 0.1\t  0.2\t  0.3\t  0.4\n", caseNumber);
    }
    if (method == 0 || method == 1 || method == 2 || method == 3) {
        for (int i = 0; i < 21; i++) {
            fprintf(filePtr, "\t"); 
            for (int j = 0; j < 4; j++) {
                if (method == 0 && caseNumber == 0.05) { 
                    if (abs(outputTemp[i][j]) < 1000) {
                        fprintf(filePtr, "%    .2f  ", outputTemp[i][j]);
                    }
                    else if (abs(outputTemp[i][j]) > 1000){
                        fprintf(filePtr, "%    .1e ", outputTemp[i][j]);
                    }
                }
                else {
                    fprintf(filePtr, "%.2f\t", outputTemp[i][j]);
                }
            }
            fprintf(filePtr, "\n");
        }
        fprintf(filePtr, "\n\n");
    }
}