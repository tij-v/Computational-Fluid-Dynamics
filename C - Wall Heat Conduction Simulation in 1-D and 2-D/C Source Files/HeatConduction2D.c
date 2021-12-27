// CLass: CFD 
//Project 1: Question 2    



#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define thermalDiffusivity 0.645 //alpha = ft^2/hr
#define xThickness 3.5 //ft
#define yThickness 3.5 //ft
#define initialTemp 0 
#define finalTemp 200

void arrayInitialization(double** wallTemp);
void ThomasAlgorithm(double* above, double* below, double* constant, double* diagonal, double* tempHolder);
void ADIScheme(int xIter, int yIter, int tIter, double xDiffusionNumber, double yDiffusionNumber, double** wallTemp);
void printData(double** wallTemp, double n);
int main(void) {
    double totalTime = 0.5;
    double** wallTemp = (double**)malloc(36 * sizeof(double*));
    double xDiffusionNumber, yDiffusionNumber;
    double deltaX, deltaY, deltaT;
    
    int xIter, yIter, tIter;

    deltaX = deltaY = 0.1;
    deltaT = 0.01;

    xIter = xThickness / deltaX;
    yIter = yThickness / deltaY;
    tIter = totalTime / deltaT;

    xDiffusionNumber = thermalDiffusivity * (deltaT / (deltaX * deltaX));
    yDiffusionNumber = thermalDiffusivity * (deltaT / (deltaY * deltaY));

    for (int i = 0; i < 36; i++) {
        wallTemp[i] = (double*)malloc(36 * sizeof(double));
    }

    arrayInitialization(wallTemp);
    ADIScheme(xIter, yIter, tIter, xDiffusionNumber, yDiffusionNumber, wallTemp);

    //printf("%d %d %d", xIter, yIter, tIter);

}

void arrayInitialization(double** wallTemp) {
    //Initializing the array. I will use the reverse array methodology for population, thus, implementing
    //a graphical mentality. 

    //"sweeping" through bottom x-axis. 
    for (int i = 0; i < 36; i++) {
        wallTemp[0][i] = finalTemp;
    }
    //"sweeping" through initial y-axis
    for (int j = 0; j < 36; j++) {
        wallTemp[j][0] = finalTemp;
    }
    //maintaining the order of initializing in the reverse array method
    for (int j = 1; j < 36; j++) {
        for (int i = 1; i < 36; i++) {
            wallTemp[j][i] = initialTemp;
        }
    }
}
void ThomasAlgorithm(double* above, double* below, double* constant, double* diagonal, double* tempHolder) {
    for (int i = 1; i < 34; i++) {
        diagonal[i] = diagonal[i] - ((above[i - 1] * below[i - 1]) / diagonal[i - 1]);
        constant[i] = constant[i] - ((below[i - 1] * constant[i - 1]) / diagonal[i - 1]);
    }
    tempHolder[33] = constant[33] / diagonal[33];
    for (int i = 32; i > -1; i--) {
        tempHolder[i] = (constant[i] - (above[i] * tempHolder[i + 1])) / diagonal[i];
    }
}

void ADIScheme(int xIter, int yIter, int tIter, double xDiffusionNumber, double yDiffusionNumber, double** wallTemp) {
    double* xAbove = (double*)malloc(33 * sizeof(double));
    double* xBelow = (double*)malloc(33 * sizeof(double));
    double* xConstant = (double*)malloc(34 * sizeof(double));
    double* xDiagonal = (double*)malloc(34 * sizeof(double));
    double* xSol = (double*)malloc(34 * sizeof(double));

    //"y-direction"
    double* yAbove = (double*)malloc(33 * sizeof(double));
    double* yBelow = (double*)malloc(33 * sizeof(double));
    double* yConstant = (double*)malloc(34 * sizeof(double));
    double* yDiagonal = (double*)malloc(34 * sizeof(double));
    double* ySol = (double*)malloc(34 * sizeof(double));

    double** solArray = (double**)malloc(36 * sizeof(double*));
    for (int j = 0; j < 36; j++) {
        solArray[j] = (double*)malloc(36 * sizeof(double));
    }
    for (int j = 0; j < 36; j++) {
        for (int i = 0; i < 36; i++) {
            solArray[j][i] = wallTemp[j][i];
        }
    }


    for (int n = 1; n < 51; n++) {
        for (int j = 1; j < 35; j++) {
            for (int i = 0; i < 33; i++) {
                xAbove[i] = (-0.5) * xDiffusionNumber;
                xBelow[i] = (-0.5) * xDiffusionNumber;
            }
            for (int i = 0; i < 34; i++) {
                xConstant[i] = ((0.5 * xDiffusionNumber) * wallTemp[j][i + 2]) + ((0.5 * xDiffusionNumber) * wallTemp[j][i]) + ((1 - xDiffusionNumber) * wallTemp[j][i + 1]);
                xDiagonal[i] = (1 + xDiffusionNumber); //1+{2*[(1/2)*xDiffusionNumber]}
            }
            //updating the boundaries
            xConstant[0] = xConstant[0] + ((0.5 * xDiffusionNumber) * wallTemp[j][0]);
            xConstant[33] = xConstant[33] + ((0.5 * xDiffusionNumber) * wallTemp[j][35]);

            //calling the Thomas Algorithm to solve for the tridiagonal matrix 
            ThomasAlgorithm(xAbove, xBelow, xConstant, xDiagonal, xSol);

            for (int i = 0; i < 34; i++) {
                solArray[j][i + 1] = xSol[i];
            }
        }
        for (int i = 1; i < 35; i++) {
            for (int j = 0; j < 33; j++) {
                yAbove[j] = (-0.5) * yDiffusionNumber;
                yBelow[j] = (-0.5) * yDiffusionNumber;
            }
            for (int j = 0; j < 34; j++) {
                yConstant[j] = ((0.5 * yDiffusionNumber) * solArray[j+2][i]) + ((0.5 * yDiffusionNumber) * solArray[j][i]) + ((1 - yDiffusionNumber) * solArray[j+1][i]);
                yDiagonal[j] = (1 + yDiffusionNumber);
            }
            //updating the boundaries
            yConstant[0] = yConstant[0] + ((0.5 * yDiffusionNumber) * wallTemp[0][i]);
            yConstant[33] = yConstant[33] + ((0.5 * yDiffusionNumber) * solArray[35][i]);

            //calling the Thomas Algorithm to solve for the tridiagonal matrix 
            ThomasAlgorithm(yAbove, yBelow, yConstant, yDiagonal, ySol);

            for (int j = 0; j < 34; j++) {
                wallTemp[j + 1][i] = ySol[j];
            }
        }
        if (n == 10 || n == 40) {
            printData(wallTemp, n);
        }
    }

    free(xAbove);
    free(xBelow);
    free(xDiagonal);
    free(xConstant);
    free(xSol);
}

void printData(double** wallTemp, double n) {
    if (n == 10) {
        FILE* filePtr;
        filePtr = fopen("Question2_output_10.txt", "w");
        fprintf(filePtr, "\n\nTemperature distribution at t = %.1fhr.\n\n\n", ((double)n) / 100);
        fprintf(filePtr, "\t  y     |\tx =  %.1f\t   x =  %.1f\t   x =  %.1f\t    x =  %.1f\t (ft)\n", 0.5, 1.5, 2.5, 3.5);
        fprintf(filePtr, "\t--------------------------------------------------------------\n\n");
        for (int j = 0; j < 36; j++) {
            fprintf(filePtr, "\t %3.1f\t|", ((double)j) / 10);
            for (int i = 0; i < 36; i++) {
                if (i == 5 || i == 15 || i == 25 || i == 35) {
                    if (abs(wallTemp[j][i]) < 1000) {
                        fprintf(filePtr, "   %8.3f\t", wallTemp[j][i]);
                    }
                    else if (abs(wallTemp[j][i]) > 1000) {
                        fprintf(filePtr, "   %8.3e\t", wallTemp[j][i]);
                    }
                }
            }
            fprintf(filePtr, "\n");
        }
        fclose(filePtr);
    }
    else if (n == 40) {
        FILE* filePtr;
        filePtr = fopen("Question2_output_40.txt", "w");
        fprintf(filePtr, "\n\nTemperature distribution at t = %.1fhr.\n\n\n", ((double)n) / 100);
        fprintf(filePtr, "\t  y     |\tx =  %.1f\t   x =  %.1f\t   x =  %.1f\t    x =  %.1f\t (ft)\n", 0.5, 1.5, 2.5, 3.5);
        fprintf(filePtr, "\t--------------------------------------------------------------\n\n");
        for (int j = 0; j < 36; j++) {
            fprintf(filePtr, "\t %3.1f\t|", ((double)j) / 10);
            for (int i = 0; i < 36; i++) {
                if (i == 5 || i == 15 || i == 25 || i == 35) {
                    if (abs(wallTemp[j][i]) < 1000) {
                        fprintf(filePtr, "   %8.3f\t", wallTemp[j][i]);
                    }
                    else if (abs(wallTemp[j][i]) > 1000) {
                        fprintf(filePtr, "   %8.3e\t", wallTemp[j][i]);
                    }
                }
            }
            fprintf(filePtr, "\n");
        }
        fclose(filePtr);
    }    
}
