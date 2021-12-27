//Name: Shitij Vishwakarma
// CLass: Computational Fluid Dynamics
// Project 3 Question 2
// Professor: Dr. Mark Archambault
//

//the libraries used
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//global definitions 
#define firstHalflength 20
#define secondHalflength 40
#define firstHalfspeed 5
#define secondHalfspeed 0
#define totalTime 2.4 //seconds
//functions used 
void initialize(double* currData);
void MacCormackScheme(double* initData, double* currData, int tIter, double deltaX, double deltaT);
void printData(double* currData, int j, double** printArray, int time, double deltaX, double deltaT);
//main driving function
int main(void) {
	double deltaX;
	double deltaT[2] = { 0.1, 0.2 };
	double* initData = (double*)malloc(42 * sizeof(double));
	double* currData = (double*)malloc(42 * sizeof(double));
	int xIter1, xIter2, tIter, a;

	deltaX = 1;
	xIter1 = firstHalflength;
	xIter2 = secondHalflength;
	
	for (int i = 0; i < 2; i++) {
		tIter = (totalTime / deltaT[i]) + 1;
		initialize(currData);
		MacCormackScheme(initData, currData, tIter, deltaX, deltaT[i]);
	}
}
void initialize(double* currData) {
	for (int i = 0; i < firstHalflength; i++) {
		currData[i] = firstHalfspeed;
	}
	for (int i = firstHalflength; i < secondHalflength + 1; i++) {
		currData[i] = secondHalfspeed;
	}
}
void MacCormackScheme(double* initData, double* currData, int tIter, double deltaX, double deltaT) {
	double* arrayE = (double*)malloc(42 * sizeof(double));
	double** printArray = (double**)malloc(42 * sizeof(double*));
	int counter = 0.4 / deltaT;
	int columns;
	for (int k = 0; k < 42; k++) {
		printArray[k] = (double*)malloc(8 * sizeof(double));
	}
	for (int i = 0; i < secondHalflength; i++) {
		initData[i] = currData[i];
	}
	for (int time = 0; time < tIter + 1; time++) {
		if (time % counter == 0) {
			if (time == 0) {
				columns = 0;
			}
			printData(currData, columns, printArray, time, deltaX, deltaT);
			columns++;
		}
		for (int i = 0; i < secondHalflength; i++) {
			arrayE[i] = (currData[i] * currData[i]) / 2;
		}
		for (int i = 0 ; i < (secondHalflength - 1); i++) {
			initData[i] = currData[i] - ((deltaT / deltaX)* (arrayE[i + 1] - arrayE[i]));
		}
		for (int i = 0; i < secondHalflength; i++) {
			arrayE[i] = (initData[i] * initData[i]) / 2;
		}
		for (int i = 1; i < secondHalflength; i++) {
			currData[i] = 0.5 * (currData[i] + initData[i] - ((deltaT / deltaX) * (arrayE[i] - arrayE[i - 1])));
		}
	}
}
void printData(double* currData, int j, double** printArray, int time, double deltaX, double deltaT) {
	FILE* frmtdData;
	frmtdData = fopen("Q2_formatted_final.txt", "a");
	for (int i = 0; i < secondHalflength + 1; i++) {
		printArray[i][j] = currData[i];
	}
	if (j == 6) {
		for (int k = 0; k < secondHalflength + 1; k++) {
			if (k == 0) {
				fprintf(frmtdData, "\n\n\ndx = %.3f\t;\tdt = %.3f\n\n", deltaX, deltaT);
				fprintf(frmtdData, "time |      0.0\t      0.4\t    0.8\t      1.2\t    1.6\t      2.0\t    2.4  | sec\nx    |\n_____________________________________________________________________________\n\n");
			}
			fprintf(frmtdData, "%d \t |", k);
			for (int l = 0; l < j + 1; l++) {
				fprintf(frmtdData, "   %7.3f", printArray[k][l]);
			}
			fprintf(frmtdData, "\n");
		}
	}
	fclose(frmtdData);
}