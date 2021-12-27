//Name: Shitij Vishwakarma
// CLass: Computational Fluid Dynamics
// Project 3 Question 1
// Professor: Dr. Mark Archambault
//

//the libraries used

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//global definitions
#define Amplitude 20 //wave amplitude
#define Wavelength 20
#define startPosition 5 //initial wave position 
#define iMax 71
#define a 200 //speed of  sound (m/s)
#define totalTime 0.15 //seconds
//functions used in this program 
void FTBSExplicit(double* initData, double* currData, double courrantNumber, double tIter, double deltaT, int method);
void LaxWendroff(double* initData, double* currData, double courrantNumber, double tIter, double deltaT, int method);
void BTCSImplicit(double* initData, double* currData, double courrantNumber, double tIter, double deltaT, int method);
void ThomasAlgorithm(double* above, double* below, double* constant, double* diagonal, double* tempHolder);
void Initialize(double* initData, double deltaW);
void printData(double* currData, int method, int time, double deltaT, int j, double** printArray);
//main driving function
int main(void) {
	int tIter;
	double deltaX, deltaW, courrantNumber;
	double deltaT[4] = { 0.005, 0.0025, 0.00125 };
	double method[4] = { 1, 2, 3 };
	//1-D column arrays
	double* initData = (double*)malloc(iMax * sizeof(double));
	double* currData = (double*)malloc(iMax * sizeof(double));
	//initialization
	for (int i = 0; i < 3; i++) {
		deltaX = 1;
		deltaW = Amplitude / (Wavelength / 2); //represents the slope by which the wave is formed. 
		courrantNumber = (-a * deltaT[i]) / deltaX;
		tIter = totalTime / deltaT[i];
		for (int j = 0; j < 3; j++) {
			if (method[j] == 1) {
				//printf("\nFTBS Explicit\n");
				Initialize(currData, deltaW);
				FTBSExplicit(initData, currData, courrantNumber, tIter, deltaT[i], method[j]);
			}
			else if (method[j] == 2) {
				//printf("\nLaxWendroff\n");
				Initialize(currData, deltaW);
				LaxWendroff(initData, currData, courrantNumber, tIter, deltaT[i], method[j]);
			}
			else if (method[j] == 3) {
				Initialize(currData, deltaW);
				BTCSImplicit(initData, currData, courrantNumber, tIter, deltaT[i], method[j]);
			}
		}
	}
}
//FTBS Explicit solving function.
void FTBSExplicit(double* initData, double* currData, double courrantNumber, double tIter, double deltaT, int method) {
	//initializing the initial timescale matrix
	double** printArray = (double**)malloc(73 * sizeof(double*));
	int time = 0;
	int columns = 0;
	//the printArray is a 2D dynamic array to store all the data while solving, so that the user can print it out later. 
	//this array is also implemented in other schemes that follow.
	for (int i = 0; i < 73; i++) {
		printArray[i] = (double*)malloc(4 * sizeof(double));
	}
	for (int i = 0; i < iMax + 1; i++) {
		initData[i] = currData[i];
	}
	printData(currData, method, time, deltaT, columns, printArray);
	for (time = 1; time < tIter + 1; time++) {
		for (int i = 1; i < iMax; i++) {
			currData[i] = initData[i] + courrantNumber * (initData[i] - initData[i - 1]);
		}
		for (int i = 0; i < iMax + 1; i++) {
			initData[i] = currData[i];
		}
		if (time == tIter / 2 || time == tIter) {
			columns++;
			printData(currData, method, time, deltaT, columns, printArray);
		}
	}
}
//Lax-Wendroff scheme solving function
void LaxWendroff(double* initData, double* currData, double courrantNumber, double tIter, double deltaT, int method) {
	double** printArray = (double**)malloc(73 * sizeof(double*));
	int time = 0;
	int columns = 0;
	for (int i = 0; i < 73; i++) {
		printArray[i] = (double*)malloc(4 * sizeof(double));
	}
	for (int i = 0; i < iMax + 1; i++) {
		initData[i] = currData[i];
	}
	printData(currData, method, time, deltaT, columns, printArray);
	for (time = 1; time < tIter + 1; time++) {
		for (int i = 1; i < iMax; i++) {
			currData[i] = initData[i] + (courrantNumber / 2) * (initData[i + 1] - initData[i - 1]) + ((courrantNumber * courrantNumber) / 2) * (initData[i + 1] - (2 * initData[i]) + initData[i - 1]);
		}
		for (int i = 0; i < iMax + 1; i++) {
			initData[i] = currData[i];
		}
		if (time == tIter / 2 || time == tIter) {
			columns++;
			printData(currData, method, time, deltaT, columns, printArray);
		}
	}
}
//BTCS Implicit solving function. This function calls THomas ALgorithm function to solve the tridiagonal matrix.
void BTCSImplicit(double* initData, double* currData, double courrantNumber, double tIter, double deltaT, int method) {
	//thomas algorithm setup.
	double* above = (double*)malloc((iMax - 2) * sizeof(double));
	double* below = (double*)malloc((iMax - 2) * sizeof(double));
	double* mainDiagonal = (double*)malloc((iMax) * sizeof(double));
	double* constant = (double*)malloc((iMax) * sizeof(double));
	double* solutionArray = (double*)malloc((iMax) * sizeof(double));
	double** printArray = (double**)malloc(78 * sizeof(double*));
	int time = 0;
	int columns = 0;
	for (int i = 0; i < 74; i++) {
		printArray[i] = (double*)malloc(4 * sizeof(double));
	}
	for (int i = 0; i < iMax + 2; i++) {
		initData[i] = currData[i];
	}
	printData(currData, method, time, deltaT, columns, printArray);

	for (time = 1; time < tIter + 1; time++) {
		for (int i = 0; i < (iMax - 2); i++) {
			above[i] = -courrantNumber / 2;
			below[i] = courrantNumber / 2;
		}
		for (int i = 0; i < (iMax - 1); i++) {
			constant[i] = currData[i + 1];
			mainDiagonal[i] = 1;
		}
		constant[0] = constant[0] - currData[0];
		constant[iMax - 2] = constant[iMax - 2] - currData[iMax - 1];
		ThomasAlgorithm(above, below, constant, mainDiagonal, solutionArray);
		for (int i = 0; i < iMax; i++) {
			currData[i + 1] = solutionArray[i];
		}
		if (time == tIter / 2 || time == tIter) {
			columns++;
			printData(currData, method, time, deltaT, columns, printArray);
		}
	}
	free(above);
	free(below);
	free(mainDiagonal);
	free(constant);
	free(solutionArray);
}
//thomas algorithm called by the BTCS function to solve the problem implicitly
void ThomasAlgorithm(double* above, double* below, double* constant, double* diagonal, double* tempHolder) {
	for (int i = 1; i < iMax - 1; i++) {
		diagonal[i] = diagonal[i] - ((above[i - 1] * below[i - 1]) / diagonal[i - 1]);
		constant[i] = constant[i] - ((below[i - 1] * constant[i - 1]) / diagonal[i - 1]);
	}
	tempHolder[iMax - 2] = constant[iMax - 2] / diagonal[iMax - 2];
	for (int i = (iMax - 3); i > -1; i--) {
		tempHolder[i] = (constant[i] - (above[i] * tempHolder[i + 1])) / diagonal[i];
	}
}
//initialization function to set the matrix to the initial environment
void Initialize(double* currData, double deltaW) {
	int waveTraversal = 0;
	for (int i = 0; i < iMax + 1; i++) {
		currData[i] = 0;
	}
	for (int i = startPosition; i < (startPosition + (Wavelength / 2)) + 1; i++) {
		currData[i] = waveTraversal;
		waveTraversal += deltaW;
	}
	waveTraversal -= deltaW;
	for (int i = (startPosition + (Wavelength / 2)); i < (startPosition + Wavelength) + 1; i++) {
		currData[i] = waveTraversal;
		waveTraversal -= deltaW;
	}
}
//print data to print the data output into a file. Each output is printed onto its own file. 
void printData(double* currData, int method, int time, double deltaT, int j, double** printArray) {
	if (method == 1) {
		FILE* FTBSout;
		FTBSout = fopen("Q1_FTBSExplicit.txt", "a");
		for (int i = 0; i < iMax + 1; i++) {
			printArray[i][j] = currData[i];
		}
		if (j == 2) {
			for (int k = 1; k < iMax + 1; k++) {
				if (k == 1) {
					fprintf(FTBSout, "\n\nDelta T = %.5f\t;\tDelta X = %.3f\n\n", deltaT, 1.000);
					fprintf(FTBSout, "Time |  \t0.000  \t   0.075      0.150 |sec\nx\t |   \t  y  \t     y  \t   y\n____________________________________________\n\n");
				}
				fprintf(FTBSout, "%d \t |", k);
				for (int l = 0; l < j + 1; l++) {
					fprintf(FTBSout, "    %7.3f", printArray[k][l]);
				}
				fprintf(FTBSout, "\n");
			}
		}
		fclose(FTBSout);
	}
	else if (method == 2) {
		FILE* LAXout;
		LAXout = fopen("Q1_LaxWendroff.txt", "a");
		for (int i = 0; i < iMax + 1; i++) {
			printArray[i][j] = currData[i];
		}
		if (j == 2) {
			for (int k = 1; k < iMax + 1; k++) {
				if (k == 1) {
					fprintf(LAXout, "\n\nDelta T = %.5f\t;\tDelta X = %.3f\n\n", deltaT, 1.000);
					fprintf(LAXout, "Time |  \t0.000  \t   0.075      0.150 |sec\nx\t |   \t  y  \t     y  \t   y\n____________________________________________\n\n");
				}
				fprintf(LAXout, "%d \t |", k);
				for (int l = 0; l < j + 1; l++) {
					fprintf(LAXout, "    %7.3f", printArray[k][l]);
				}
				fprintf(LAXout, "\n");
			}
		}
		fclose(LAXout);
	}
	else if (method == 3) {
		FILE* BTCSout;
		BTCSout = fopen("Q1_BTCSImplicit.txt", "a");
		for (int i = 0; i < iMax + 2; i++) {
			printArray[i][j] = currData[i];
			if (i == 71) {
				printArray[i][j] = 0;
			}
		}
		if (j == 2) {
			for (int k = 1; k < iMax + 1; k++) {
				if (k == 1) {
					fprintf(BTCSout, "\n\nDelta T = %.5f\t;\tDelta X = %.3f\n\n", deltaT, 1.000);
					fprintf(BTCSout, "Time |  \t0.000  \t   0.075      0.150 |sec\nx\t |   \t  y  \t     y  \t   y\n____________________________________________\n\n");
				}
				fprintf(BTCSout, "%d \t |", k);
				for (int l = 0; l < j + 1; l++) {
					fprintf(BTCSout, "    %7.3f", printArray[k][l]);
				}
				fprintf(BTCSout, "\n");
			}
		}
		fclose(BTCSout);
	}
}