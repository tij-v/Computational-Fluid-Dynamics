//Name: Shitij Vishwakarma
// CLass: Computational Fluid Dynamics
// Project 2
// Professor: Dr. Mark Archambault
//

//the libraries used
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//global definitions 
#define deltX 0.2
#define deltY 0.2
#define BETA (double) 1.0
#define maxPSI (double) 100.0
#define minPSI (double) 0.0
#define ERRORMAX (double) 0.01

//functions used in this code
void initializeArray(double** tempArray, int xIter, int yIter);
void PGS(double** initFlow, double** currFlow, int xIter, int yIter, int length);
void LGS(double** initFlow, double** currFlow, int xIter, int yIter, int length);
void ThomasAlgorithm(double* above, double* below, double* constant, double* diagonal, double* tempHolder);
void PSOR(double** initFlow, double** currFlow, int xIter, int yIter, int length);
void PSOR_calculator(double omegaVal, double** currFlow, double** initFlow, int i, int j);
void LSOR(double** initFlow, double** currFlow, int xIter, int yIter, int length);
void LSOR_calculator(double* above, double* below, double* mainDiagonal, double* constant, double* solutionArray, double omegaVal, double** currFlow, double** initFlow, int j, int xIter);
double errorCalc(double** initFlow, double** currFlow, int xIter, int yIter);
void printData(double** printArray, int xIter, int yIter, double length, int method, int iteration, double omegaVal);
void printRawdata(double** currFlow, int xIter, int yIter, int method);
int main(void) {
	double height = 4.0;
	double length = 6.0;

	int xIter = length / deltX;//=30
	int yIter = height / deltY; //=20

	double** initFlow = (double**)malloc(32 * sizeof(double*));
	double** currFlow = (double**)malloc(32 * sizeof(double*));

	for (int i = 0; i < xIter + 1; i++) {
		initFlow[i] = (double*)malloc(21 * sizeof(double));
		currFlow[i] = (double*)malloc(21 * sizeof(double));
	}
	initializeArray(currFlow, xIter, yIter);
	PGS(initFlow, currFlow, xIter, yIter, length);
	initializeArray(currFlow, xIter, yIter);
	LGS(initFlow, currFlow, xIter, yIter, length);
	initializeArray(currFlow, xIter, yIter);
	PSOR(initFlow, currFlow, xIter, yIter, length);
	initializeArray(currFlow, xIter, yIter);
	LSOR(initFlow, currFlow, xIter, yIter, length);
	free(initFlow);
	free(currFlow);
}

//this function sets the matrix/environment in the desired intial condition as per the question
void initializeArray(double** tempArray, int xIter, int yIter) {
	int flowInlet = 1 / deltX;
	int iter = 0;
	//set (i,0) BC: bottom-most wall
	while (iter < xIter + 1) {
		if (iter < flowInlet) {
			tempArray[iter][0] = minPSI;
		}
		else if (iter >= flowInlet) {
			tempArray[iter][0] = maxPSI;
		}
		iter++;
	}
	//set (i, y_max) BC: top-most wall
	for (int i = 0; i < xIter + 1; i++) {
		tempArray[i][yIter] = minPSI;
	}
	//set (0,j) BC: left-most wall
	for (int j = 0; j < yIter + 1; j++) {
		tempArray[0][j] = minPSI;
	}
	//complete inner array
	for (int j = 1; j < yIter; j++) {
		for (int i = 1; i < xIter + 1; i++) {
			tempArray[i][j] = 0;
		}
	}
}
//this function solves for the solution using the Point Gauss-Seidel Method scheme
void PGS(double** initFlow, double** currFlow, int xIter, int yIter, int length) {
	double constant = 1.0 / (2.0 * (1.0 + (BETA * BETA)));
	double error = 0;

	int errorCondition = 0;
	int iteration = 0;
	//set the (n-1)th array to the initial conditions
	for (int j = 0; j < yIter + 1; j++) {
		for (int i = 0; i < xIter + 1; i++) {
			initFlow[i][j] = currFlow[i][j];
		}
	}
	//error condition is 0 when the solution is not converged. When the error < 0.01, the condition changes to 1
	while (errorCondition == 0) {
		for (int j = 1; j < yIter; j++) {
			for (int i = 1; i < xIter; i++) {
				currFlow[i][j] = constant * (initFlow[i + 1][j] + currFlow[i - 1][j] + ((BETA * BETA) * (initFlow[i][j + 1] + currFlow[i][j - 1])));
			}
			currFlow[xIter][j] = currFlow[xIter - 1][j];
		}
		iteration++;
		error = 0;
		error = errorCalc(initFlow, currFlow, xIter, yIter);
		//printf("Error = %f\n", error);
		for (int j = 1; j < yIter + 1; j++) {
			for (int i = 1; i < xIter + 1; i++) {
				initFlow[i][j] = currFlow[i][j];
			}
		}
		if (error < ERRORMAX) {
			errorCondition = 1;
		}
	}
	printData(currFlow, xIter, yIter, length, 1, iteration, 1);
	printRawdata(currFlow, xIter, yIter, 1);
}
//this function sovles for the solution using the Line Gauss-Seidel Method scheme
void LGS(double** initFlow, double** currFlow, int xIter, int yIter, int length) {
	double* above = (double*)malloc(29 * sizeof(double));
	double* below = (double*)malloc(29 * sizeof(double));
	double* mainDiagonal = (double*)malloc(29 * sizeof(double));
	double* constant = (double*)malloc(29 * sizeof(double));
	double* solutionArray = (double*)malloc(29 * sizeof(double));

	int iteration = 0;
	int error = 0;
	int errorCondition = 0;

	for (int j = 0; j < yIter + 1; j++) {
		for (int i = 0; i < xIter + 1; i++) {
			initFlow[i][j] = currFlow[i][j];
		}
	}
	while (errorCondition == 0) {
		//initialiizing the tridiagonal matrix diagonals
		for (int j = 1; j < yIter; j++) {
			for (int i = 0; i < xIter - 2; i++) {
				above[i] = 1.0;
				below[i] = 1.0;
			}
			for (int i = 0; i < xIter - 1; i++) {
				mainDiagonal[i] = -(2 * (1 + (BETA * BETA)));
				constant[i] = -(BETA * BETA) * (initFlow[i + 1][j + 1] + currFlow[i + 1][j - 1]);
			}
			constant[0] = constant[0] - initFlow[0][j];
			constant[xIter - 2] = constant[xIter - 2] - initFlow[xIter][j];

			ThomasAlgorithm(above, below, constant, mainDiagonal, solutionArray);
			for (int i = 0; i < xIter - 1; i++) {
				currFlow[i + 1][j] = solutionArray[i];
			}
			currFlow[xIter][j] = currFlow[xIter - 1][j];
		}
		iteration++;
		error = 0;
		error = errorCalc(initFlow, currFlow, xIter, yIter);
		for (int j = 1; j < yIter + 1; j++) {
			for (int i = 1; i < xIter + 1; i++) {
				initFlow[i][j] = currFlow[i][j];
			}
		}
		if (error < ERRORMAX) {
			errorCondition = 1;
		}
	}
	printData(currFlow, xIter, yIter, length, 2, iteration, 1);
	printRawdata(currFlow, xIter, yIter, 2);
}
//this function solves the Thomas algorithm, utilized by the LGS and the LSOR_calculator functions
void ThomasAlgorithm(double* above, double* below, double* constant, double* diagonal, double* tempHolder) {
	for (int i = 1; i < 29; i++) {
		diagonal[i] = diagonal[i] - ((above[i - 1] * below[i - 1]) / diagonal[i - 1]);
		constant[i] = constant[i] - ((below[i - 1] * constant[i - 1]) / diagonal[i - 1]);
	}
	tempHolder[28] = constant[28] / diagonal[28];
	for (int i = 27; i > -1; i--) {
		tempHolder[i] = (constant[i] - (above[i] * tempHolder[i + 1])) / diagonal[i];
	}
}
//this function solves for the solution using the Point Successive Over Relaxation Method scheme
void PSOR(double** initFlow, double** currFlow, int xIter, int yIter, int length) {
	double* omega = (double*)malloc(2100 * sizeof(double));
	double* prevCurrflow = (double*)malloc(2100 * sizeof(double));
	double constant = 1.0 / (2.0 * (1.0 + (BETA * BETA)));
	double error = 0;

	int errorCondition = 0;
	int iteration = 0;
	int optValcatcher = 0;
	omega[0] = 0;
	//for loop to iterate the value of omega from 0 - 2.1 in 0.001 steps
	for (int k = 1; k < 2100; k++) {
		//initalizing the (n-1)th array with each new iteration
		initializeArray(currFlow, xIter, yIter);
		for (int j = 0; j < yIter + 1; j++) {
			for (int i = 0; i < xIter + 1; i++) {
				initFlow[i][j] = currFlow[i][j];
			}
		}
		//these steps below remain similar to the PGS method --the only difference being the presence of omega
		errorCondition = 0;
		iteration = 0;
		while (errorCondition == 0) {
			for (int j = 1; j < yIter; j++) {
				for (int i = 1; i < xIter; i++) {
					PSOR_calculator(omega[k - 1], currFlow, initFlow, i, j);
				}
				currFlow[xIter][j] = currFlow[xIter - 1][j];
			}
			iteration++;
			error = errorCalc(initFlow, currFlow, xIter, yIter);
			if (error < ERRORMAX) {
				errorCondition = 1;
			}
			for (int j = 1; j < yIter + 1; j++) {
				for (int i = 1; i < xIter + 1; i++) {
					initFlow[i][j] = currFlow[i][j];
					if (j == 1 && i == xIter) {
						prevCurrflow[k - 1] = currFlow[i][j];
					}
				}
			}
			
		}
		//cathing the optimal omega value:
		if (k > 1840 && (int) (prevCurrflow[k-1]*100) < (int) ((prevCurrflow[k-2])*100)) {
			if (optValcatcher == 0) {
				printData(currFlow, xIter, yIter, length, 3, iteration, omega[k-1]);
				printRawdata(currFlow, xIter, yIter, 3);
				optValcatcher = 1;
			}
		}
		omega[k] = omega[k-1] + 0.001;
	}
}
//this function uses the PSOR equation to calculate the solution
void PSOR_calculator(double omegaVal, double** currFlow, double** initFlow, int i, int j) {
	currFlow[i][j] = initFlow[i][j] + (((omegaVal) / (2.0 * (1.0 + (BETA * BETA)))) * (initFlow[i + 1][j] + currFlow[i - 1][j] + ((BETA * BETA)*(initFlow[i][j + 1] + currFlow[i][j - 1])) - (2 * (1 + (BETA * BETA)) * initFlow[i][j])));
}
//this function solves for the solution using the Line Successive Over Relaxation Method scheme
void LSOR(double** initFlow, double** currFlow, int xIter, int yIter, int length) {
	//initializing the tridiagonal diagonals
	double* above = (double*)malloc(29 * sizeof(double));
	double* below = (double*)malloc(29 * sizeof(double));
	double* mainDiagonal = (double*)malloc(29 * sizeof(double));
	double* constant = (double*)malloc(29 * sizeof(double));
	double* solutionArray = (double*)malloc(29 * sizeof(double));
	double* omega = (double*)malloc(20002 * sizeof(double));
	double* prevCurrflow = (double*)malloc(20002 * sizeof(double));

	int iteration = 0;
	int error = 0;
	int errorCondition = 0;
	int optValcatcher = 0;
	omega[0] = 0;

	for (int k = 1; k < 2100; k++) {
		//initializing the (n-1)th array
		initializeArray(currFlow, xIter, yIter);
		for (int j = 0; j < yIter + 1; j++) {
			for (int i = 0; i < xIter + 1; i++) {
				initFlow[i][j] = currFlow[i][j];
			}
		}
		errorCondition = 0;
		iteration = 0;
		while (errorCondition == 0) {
			for (int j = 1; j < yIter; j++) {
				LSOR_calculator(above, below, mainDiagonal, constant, solutionArray, omega[k - 1], currFlow, initFlow, j, xIter);
				currFlow[xIter][j] = currFlow[xIter - 1][j];
			}
			iteration++;
			error = errorCalc(initFlow, currFlow, xIter, yIter);
			for (int j = 1; j < yIter + 1; j++) {
				for (int i = 1; i < xIter + 1; i++) {
					initFlow[i][j] = currFlow[i][j];
					if (j == 1 && i == xIter) {
						prevCurrflow[k - 1] = currFlow[i][j];
					}
				}
			}
			if (error < ERRORMAX) {
				errorCondition = 1;
			}
		}
		//catching the optimal omega value
		if (k > 1100 && (int) prevCurrflow[k - 2] > (int) prevCurrflow[k - 1]){
			if (optValcatcher == 0) {
				printData(currFlow, xIter, yIter, length, 4, iteration, omega[k - 1]);
				printRawdata(currFlow, xIter, yIter, 4);
				optValcatcher = 1;
			}
		}
		omega[k] = omega[k - 1] + 0.001;
	}
} 
//this function uses the LSOR equation to calculate the solution. This function also calls the ThomasAlgorithm function
void LSOR_calculator(double* above, double* below, double* mainDiagonal, double* constant, double* solutionArray, double omegaVal, double** currFlow, double** initFlow, int j, int xIter) {
	//initializing the tridiagonal matrix
	for (int i = 0; i < xIter - 2; i++) {
		above[i] = omegaVal;
		below[i] = omegaVal;
	}
	for (int i = 0; i < xIter - 1; i++) {
		mainDiagonal[i] = -(2 * (1 + (BETA * BETA)));
		constant[i] = ((-(1 - omegaVal) * (2 * (1 + (BETA * BETA)))) * initFlow[i+1][j]) - ((omegaVal * (BETA * BETA)) * (initFlow[i + 1][j + 1] + currFlow[i + 1][j - 1]));
	}
	constant[0] = constant[0] - (omegaVal * initFlow[0][j]);
	constant[xIter - 2] = constant[xIter - 2] - (omegaVal * initFlow[xIter][j]);
	//calling the TA function to solve for the solution
	ThomasAlgorithm(above, below, constant, mainDiagonal, solutionArray);
	for (int i = 0; i < xIter - 1; i++) {
		currFlow[i + 1][j] = solutionArray[i];
	}
}
//this function  calculates the value of the error, which is needed as a convergence criteria in all the mthod schemes. it returns a double which is the value of the error
double errorCalc(double** initFlow, double** currFlow, int xIter, int yIter) {
	double sum = 0;
	for (int j = 1; j < yIter; j++) {
		for (int i = 1; i < xIter; i++) {
			sum = sum + (currFlow[i][j] - initFlow[i][j]);
		}
	}
	return(sum);
}
//printData prints the formatted data table for priniting
void printData(double** printArray, int xIter, int yIter, double length, int method, int iteration, double omegaVal) {
	FILE* frmtdData;
	
	frmtdData = fopen("Formatted_Data.txt", "a");
	if (method == 1 || method == 2 || method == 3 || method == 4) {
		if (method == 1) {
			fprintf(frmtdData, "Formatted Data\n\n\nPoint Gauss-Seidel\nIteration = %d\n\n \t x |       %d.0\t |\t   %d.0\t  |\t   %d.0\t  |\t   %d.0\n  y\n--------------------------------------------------------\n", iteration, 0, 2, 4, 6);
		}
		else if (method == 2) {
			fprintf(frmtdData, "\n\nLine Gauss-Seidel\nIteration = %d\n\n \t x |       %d.0\t |\t   %d.0\t  |\t   %d.0\t  |\t   %d.0\n  y\n--------------------------------------------------------\n", iteration, 0, 2, 4, 6);
		}
		else if (method == 3) {
			fprintf(frmtdData, "\n\nPoint SOR\nIteration = %d\nOptimal Omega = %.3f\n\n \t x |       %d.0\t |\t   %d.0\t  |\t   %d.0\t  |\t   %d.0\n  y\n--------------------------------------------------------\n", iteration, omegaVal, 0, 2, 4, 6);
		}
		else if (method == 4) {
			fprintf(frmtdData, "\n\nLine SOR\nIteration = %d\nOptimal Omega = %.3f\n\n \t x |       %d.0\t |\t   %d.0\t  |\t   %d.0\t  |\t   %d.0\n  y\n--------------------------------------------------------\n", iteration, omegaVal, 0, 2, 4, 6);
		}
		for (int j = yIter; j > -1; j--) {
			fprintf(frmtdData, "%.3f  |", (j * 0.2));
			for (int i = 0; i < xIter + 1; i++) {
				if (i % 10 == 0) {
					fprintf(frmtdData, "   %8.3f\t", printArray[i][j]);
				}
			}
			fprintf(frmtdData, "\n");
		}
	}
	fclose(frmtdData);
}
//printRawdata prints the raw data solution into the file
void printRawdata(double** currFlow, int xIter, int yIter, int method) {
	FILE* rawData;

	rawData = fopen("Raw_Data.txt", "a");
	if (method == 1 || method == 2 || method == 3 || method == 4) {
		if (method == 1) {
			fprintf(rawData, "Raw Data\n\n\nPoint Gauss-Seidel\n");
		}
		else if (method == 2) {
			fprintf(rawData, "\nLine Gauss-Seidel\n");
		}
		else if (method == 3) {
			fprintf(rawData, "\nPoint SOR\n");
		}
		else if (method == 4) {
			fprintf(rawData, "\nLine SOR\n");
		}
		for (int j = yIter; j > -1; j--) {
			for (int i = 0; i < xIter + 1; i++) {
				fprintf(rawData, "   %3.3f\t", currFlow[i][j]);
			}
			fprintf(rawData, "\n");
		}
	}
	fclose(rawData);
}

/*end of the program*/