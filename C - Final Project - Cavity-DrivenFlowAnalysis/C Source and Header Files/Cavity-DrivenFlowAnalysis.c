//Name: Shitij Vishwakarma
//Final Project
//Computational Fluid Dynamics
//Dr. Mark Archambault
//libraries used
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "ThomasAlgorithm.h"
//global definitions
double pressure[JM][IM];
double uArray[JM][IM];
double vArray[JM][IM];
double uArrayOld[JM][IM];
double vArrayOld[JM][IM];
double above[IM - 3];
double below[IM - 3];
double diagonal[IM - 2];
double constant[IM - 2];
double solution[IM - 2];
//function callings
void Initialize();
void PressureCalculation();
double pressureRHS(int i, int j);
double Dialation(int i, int j);
void MomentumCalculation();
int PrintRawData();
//main driving function
int main(void) {
	Initialize();
	for (int time = 1; time < TM; time++) {
		printf("Time step: %d\n", time);
		PressureCalculation();
		MomentumCalculation();
	}
	PrintRawData();
	return (0);
}
//this function initializes the grid into the initial condtions provided
void Initialize(){
	for (int j = 0; j < JM; j++) {
		for (int i = 0; i < IM; i++) {
			pressure[j][i] = initialP;
			uArray[j][i] = initialU;
			vArray[j][i] = initialV;
		}
	}
	for (int i = 0; i < IM; i++) {
		pressure[0][i] = BottomPlatePressure;
		pressure[JM - 1][i] = TopPlatePressure;
		uArray[JM - 1][i] = TopPlateVelocity;
	}
}
//this function calculates the pressure by solving ThomasAlgorithm, tridiagonal matrix.
void PressureCalculation(){
	double pressureSum = 0.0;
	int stopCondition = 0;
	do{
		pressureSum = 0.0;
		for (int j = 1; j < (JM - 1); j++) {
			for (int i = 0; i < (IM - 3); i++) {
				above[i] = 1.0;
				below[i] = 1.0;
			}
			for (int i = 1; i < (IM - 1); i++) {
				diagonal[i - 1] = -2.0 * (1.0 + beta);
				constant[i - 1] = (dx * dx * pressureRHS(i, j)) - (beta * (pressure[j + 1][i] + pressure[j - 1][i]));
			}
			constant[0] = constant[0] - pressure[j][0];
			constant[IM - 3] = constant[IM - 3] - pressure[j][IM - 1];

			ThomasAlgorithm(above, below, diagonal, constant, solution);

			for (int i = 0; i < (IM - 2); i++) {
				pressureSum += fabs(solution[i] - pressure[j][i + 1]);
				pressure[j][i + 1] = solution[i];
			}
			pressure[j][0] = pressure[j][1];
			pressure[j][IM - 1] = pressure[j][IM - 2];
		}
	} while (pressureSum > ERRORMAX);
}
//this function solves the right hand side of the pressure equation, called by the PressureCalculation function
double pressureRHS(int i, int j){
	double DIL1, DIL2, DIL3;
	double DIF1, DIF2, DIF3, DIF4, DIF5, DIF6;

	DIL1 = Dialation(i, j);
	DIF1 = DIL1 / dt;
	
	DIF2 = ((uArray[j][i + 1] * uArray[j][i + 1]) - (2.0 * uArray[j][i] * uArray[j][i]) + (uArray[j][i - 1] * uArray[j][i - 1])) / (pow(dx, 2));
	
	DIF3 = ((uArray[j + 1][i + 1] * vArray[j + 1][i + 1]) - (uArray[j + 1][i - 1] * vArray[j + 1][i - 1]) - (uArray[j - 1][i + 1] * vArray[j - 1][i + 1]) + (uArray[j - 1][i - 1] * vArray[j - 1][i - 1])) / (2.0 * dy * dx);
	
	DIF4 = ((vArray[j + 1][i] * vArray[j + 1][i]) - (2.0 * vArray[j][i] * vArray[j][i]) + (vArray[j - 1][i] * vArray[j - 1][i])) / (pow(dy, 2));
	
	DIL2 = Dialation((i + 1), j);
	DIL3 = Dialation((i - 1), j);
	DIF5 = (DIL2 - (2.0 * DIL1) + DIL3) / (pow(dx, 2));
	
	DIL2 = Dialation(i, (j + 1));
	DIL3 = Dialation(i, (j - 1));
	DIF6 = (DIL2 - (2.0 * DIL1) + DIL3) / (pow(dy, 2));

	return (DIF1 - DIF2 - DIF3 - DIF4 + ((1.0 / Re) * (DIF5 + DIF6)));
}
//solves for the dialation variable, called by pressureRHS
double Dialation(int i, int j){
	double du_dx, dv_dy;
	if (i == 0) {
		du_dx = (-uArray[j][i + 2] + (4.0 * uArray[j][i + 1]) - (3.0 * uArray[j][i])) / (2.0 * dx);
	}
	else if (i == (IM - 1)) {
		du_dx = ((3.0 * uArray[j][i]) - (4.0 * uArray[j][i - 1]) + uArray[j][i - 2]) / (2.0 * dx);
	}
	else {
		du_dx = (uArray[j][i + 1] - uArray[j][i - 1]) / (2.0 * dx);
	}
	if (j == 0) {
		dv_dy = (-vArray[j + 2][i] + (4.0 * vArray[j + 1][i]) - (3.0 * vArray[j][i])) / (2.0 * dy);
	}
	else if(j == (JM - 1)) {
		dv_dy = ((3.0 * vArray[j][i]) - (4.0 * vArray[j - 1][i]) + vArray[j - 2][i]) / (2.0 * dy);
	}
	else {
		dv_dy = (vArray[j + 1][i] - vArray[j - 1][i]) / (2.0 * dy);
	}
	return (du_dx + dv_dy);
}
//solves  for the momentum equation, called by main
void MomentumCalculation(){
	double DIF1, DIF2, DIF3, DIF4, DIF5;
	for (int j = 0; j < JM; j++) {
		for (int i = 0; i < IM; i++) {
			uArrayOld[j][i] = uArray[j][i];
			vArrayOld[j][i] = vArray[j][i];
		}
	}
	for (int j = 1; j < (JM - 1); j++) {
		for (int i = 1; i < (IM - 1); i++) {
			if (uArrayOld[j][i] > 0) {
				DIF1 = ((uArrayOld[j][i] * uArrayOld[j][i]) - (uArrayOld[j][i - 1] * uArrayOld[j][i - 1])) / dx;
				DIF2 = (pressure[j][i] - pressure[j][i - 1]) / dx;
			}
			else {
				DIF1 = ((uArrayOld[j][i + 1] * uArrayOld[j][i + 1]) - (uArrayOld[j][i] * uArrayOld[j][i])) / dx;
				DIF2 = (pressure[j][i + 1] - pressure[j][i]) / dx;
			}
			if (vArrayOld[j][i] > 0) {
				DIF3 = ((uArrayOld[j][i] * vArrayOld[j][i]) - (uArrayOld[j - 1][i] * vArrayOld[j - 1][i])) / dy;
			}
			else {
				DIF3 = ((uArrayOld[j + 1][i] * vArrayOld[j + 1][i]) - (uArrayOld[j][i] * vArrayOld[j][i])) / dy;
			}
			DIF4 = (uArrayOld[j][i + 1] - (2.0 * uArrayOld[j][i]) + uArrayOld[j][i - 1]) / (pow(dx, 2));
			DIF5 = (uArrayOld[j + 1][i] - (2.0 * uArrayOld[j][i]) + uArrayOld[j - 1][i]) / (pow(dy, 2));
			uArray[j][i] = uArrayOld[j][i] + (dt * (-DIF1 - DIF2 - DIF3 + ((1.0 / Re) * (DIF4 + DIF5))));

			if(uArrayOld[j][i] > 0) {
				DIF1 = ((uArrayOld[j][i] * vArrayOld[j][i]) - (uArrayOld[j][i - 1] * vArrayOld[j][i - 1])) / dx;
			}
			else {
				DIF1 = ((uArrayOld[j][i + 1] * vArrayOld[j][i + 1]) - (uArrayOld[j][i] * vArrayOld[j][i])) / dx;
			}
			if (vArrayOld[j][i] > 0) {
				DIF2 = ((vArrayOld[j][i] * vArrayOld[j][i]) - (vArrayOld[j - 1][i] * vArrayOld[j - 1][i])) / dy;
				DIF3 = (pressure[j][i] - pressure[j - 1][i]) / dy;
			}
			else {
				DIF2 = ((vArrayOld[j + 1][i] * vArrayOld[j + 1][i]) - (vArrayOld[j][i] * vArrayOld[j][i])) / dy;
				DIF3 = (pressure[j + 1][i] - pressure[j][i]) / dy;
			}
			DIF4 = (vArrayOld[j][i + 1] - (2.0 * vArrayOld[j][i]) + vArrayOld[j][i - 1]) / (pow(dx, 2));
			DIF5 = (vArrayOld[j + 1][i] - (2.0 * vArrayOld[j][i]) + vArrayOld[j - 1][i]) / (pow(dy, 2));

			vArray[j][i] = vArrayOld[j][i] + (dt * (-DIF1 - DIF2 - DIF3 + ((1.0 / Re) * (DIF4 + DIF5))));
		}
	}
}
//prints raw data to be used in matlab for plotting
int PrintRawData(){
	FILE* out;
	int iter = 0;
	out = fopen("Final_Project_Raw.txt", "w");
	if (!out) {
		perror("Failed to open the file");
		return EXIT_FAILURE;
	}
	while (iter < 3) {
		for (int j = JM; j > 0; j--) {
			if (iter == 0) {
				for (int i = 1; i < (IM + 1); i++) {
					fprintf(out, " %8.5f", pressure[j - 1][i - 1]);
				}
				fprintf(out, "\n");
			}
			else if (iter == 1) {
				if (j == JM) {
					fprintf(out, "\n\n");
				}
				for (int i = 1; i < (IM + 1); i++) {
					fprintf(out, " %8.5f", uArray[j - 1][i - 1]);
				}
				fprintf(out, "\n");
			}
			else if (iter == 2) {
				if (j == JM) {
					fprintf(out, "\n\n");
				}
				for (int i = 1; i < (IM + 1); i++) {
					fprintf(out, " %8.5f", vArray[j - 1][i - 1]);
				}
				fprintf(out, "\n");
			}
		}
		iter += 1;
	}
	fclose(out);
}