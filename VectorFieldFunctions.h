#include <string>

int GetNumberOfPoints(const int*);
int GetNumberOfCells(const int*);
int GetPointIndex(const int*, const int*);
int GetCellIndex(const int*, const int*);
void GetLogicalPointIndex(int*, int, const int*);
void GetLogicalCellIndex(int*, int, const int*);
void InterpolateVectors(float, float, float, float*, float*, float*);
void EvaluateVectorFieldAtLocation(const float*, const int*, const float*, const float*, const float*, float*);
void AdvectWithRK4Step(const float*, const int*, const float*, const float*, const float*, float, int, float*);
void
AdvectWithEulerStep(const float*, const int*, const float*, const float*, const float*, float, int, float*);
float* Convert3DVectorDataTo2DVectorData(const int*, const float*);
void PrintSteps(const std::string, int, float*);
