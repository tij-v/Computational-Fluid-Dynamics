//header file called by the source file: FinalProjectSAV.c

#define IM 81
#define JM 81
#define TM 15001
#define dx 0.00625
#define dy 0.00625
#define dt 0.003
#define TopPlatePressure 3350.0
#define BottomPlatePressure 3350.0
#define TopPlateVelocity 1.0
#define Re 1000
#define ERRORMAX 0.001
#define LeftGradient 0.0
#define RightGradient 0.0
#define initialP ((TopPlatePressure+BottomPlatePressure) / 2)
#define initialU 0.0
#define initialV 0.0
#define beta ((dx*dx) / (dy*dy))