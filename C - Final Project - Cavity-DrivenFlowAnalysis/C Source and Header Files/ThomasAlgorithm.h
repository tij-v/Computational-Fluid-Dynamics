#include "constants.h"
//header file called by the source file: FinalProjectSAV.c
void ThomasAlgorithm(double* above, double* below, double* diagonal, double* constant, double* solution) {
	for (int i = 1; i < (IM - 2); i++) {
		diagonal[i] = diagonal[i] - ((above[i - 1] * below[i - 1]) / diagonal[i - 1]);
		constant[i] = constant[i] - ((constant[i - 1] * below[i - 1]) / diagonal[i - 1]);
	}
	solution[IM - 3] = constant[IM - 3] / diagonal[IM - 3];
	for (int i = (IM - 4); i > -1; i--) {
		solution[i] = (constant[i] - (above[i] * solution[i + 1])) / diagonal[i];
	}
}