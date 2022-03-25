#include "VectorFieldFunctions.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;


// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}

void InterpolateVectors(float a, float b, float x, float* vecA, float* vecB, float* res)
{
    float t = (x-a)/(b-a);
    res[0] = vecA[0] + t*(vecB[0] - vecA[0]);
    res[1] = vecA[1] + t*(vecB[1] - vecA[1]);
}

// ****************************************************************************
//  Function: EvaluateVectorFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a vector field defined on the mesh.  Its size is 2*dims[0]*dims[1].
//        The first value in the field is the x-component for the first point.
//        The second value in the field is the y-component for the first point.
//
//     rv (output): the interpolated field value. (0,0) if the location is out of bounds.
//
// ****************************************************************************

void EvaluateVectorFieldAtLocation(const float *pt, const int *dims, const float *X, 
                              const float *Y, const float *F, float *rv)
{
    if(pt[0] < X[0] || pt[1] < Y[0]){
        rv[0] = 0;
        rv[1] = 0;
        return;
    }
    // 1a. Find the nearest X value to pt[0]
    int ind[2];
    ind[0] = 0;
    ind[1] = 0;
    while(ind[0] < dims[0] && pt[0] > X[ind[0]]){
        ind[0]++;
    }

    // 1b.if out of bounds rv = (0,0) return
    if(ind[0] >= dims[0]){
        rv[0] = 0;
        rv[1] = 0;
        return;
    }

    // 1c. Find the nearest Y value to pt[1]
    while(ind[1] < dims[1] && pt[1] > Y[ind[1]]){
        ind[1]++;
    }

    // 1d. if out of bounds rv = (0,0) return
    if(ind[1] >= dims[1]){
        rv[0] = 0;
        rv[1] = 0;
        return;
    }
    

    // 2. Save vectors, logical indices, and indices for top left and top right
    
    //Get Logical indices and indices
    int topRightInd; 
    int topRightLog[2];

    topRightLog[0] = ind[0];
    topRightLog[1] = ind[1];

    topRightInd = GetPointIndex(topRightLog, dims);

    int topLeftInd;
    int topLeftLog[2];

    topLeftLog[0] = topRightLog[0] - 1;
    topLeftLog[1] = topRightLog[1];

    topLeftInd = GetPointIndex(topLeftLog, dims);

    //Get Vectors
    float topRightVec[2];
    float topLeftVec[2];
    
    topRightVec[0] = F[2*topRightInd];
    topRightVec[1] = F[2*topRightInd+1];

    topLeftVec[0] = F[2*topLeftInd];
    topLeftVec[1] = F[2*topLeftInd+1];

    // 3. Interpolate between top left and top right vectors
    float topVecLerp[2];
    InterpolateVectors(X[topLeftLog[0]], X[topRightLog[0]], pt[0], topLeftVec, topRightVec, topVecLerp);


    // 4. Save vectors, logical indices, and indices for bottom left and bottom right
    int bottomRightLog[2]; 
    int bottomRightInd;

    bottomRightLog[0] = topRightLog[0];
    bottomRightLog[1] = topRightLog[1] - 1;

    bottomRightInd = GetPointIndex(bottomRightLog, dims);

    int bottomLeftLog[2];
    int bottomLeftInd;

    bottomLeftLog[0] = bottomRightLog[0] - 1;
    bottomLeftLog[1] = bottomRightLog[1];

    bottomLeftInd = GetPointIndex(bottomLeftLog, dims);

    float bottomRightVec[2];
    float bottomLeftVec[2];

    bottomRightVec[0] = F[2*bottomRightInd];
    bottomRightVec[1] = F[2*bottomRightInd+1];

    bottomLeftVec[0] = F[2*bottomLeftInd];
    bottomLeftVec[1] = F[2*bottomLeftInd+1];


    // 5. Interpolate between bottom left and bottom right vectors
    float bottomVecLerp[2];
    InterpolateVectors(X[bottomLeftLog[0]], X[bottomRightLog[0]], pt[0], bottomLeftVec, bottomRightVec, bottomVecLerp);

    // 6. Interpolate between topLerp and bottomLerp vectors
    InterpolateVectors(Y[bottomLeftLog[1]], Y[topLeftLog[1]], pt[1], bottomVecLerp, topVecLerp, rv);

}

void AdvectWithRK4Step(const float *pt, const int *dims, const float *X, 
                       const float *Y, const float *F, 
                       float h, int nsteps, float *output_locations)
{

    output_locations[0] = pt[0];
    output_locations[1] = pt[1];
    
    //got some odd errors when using 1/6 in the final computations so save this value up here
    //probably some integer/floating point wierdness plus division is expensive so lets save the resources
    float sixth = (1.0/6.0);

    for(int i=0; i<nsteps; ++i){

        float curPos[2], k1[2], k2[2], k3[2], k4[2];   
        
        //k1 vector is based on current position
        curPos[0] = output_locations[2*i];
        curPos[1] = output_locations[2*i+1];
        EvaluateVectorFieldAtLocation(curPos, dims, X, Y, F, k1);
       
        //k2 vector is based on a half step using the k1 vector
        curPos[0] += (h*0.5)*k1[0];
        curPos[1] += (h*0.5)*k1[1];
        EvaluateVectorFieldAtLocation(curPos, dims, X, Y, F, k2);
    
        //k3 vector is based on a half step using the k2 vector
        curPos[0] = output_locations[2*i] + (h*0.5)*k2[0];
        curPos[1] = output_locations[2*i+1] + (h*0.5)*k2[1];
        EvaluateVectorFieldAtLocation(curPos, dims, X, Y, F, k3);

        //k4 is based on a whole step with the k3 vector
        curPos[0] = output_locations[2*i] + h*k3[0];
        curPos[1] = output_locations[2*i+1] + h*k3[1];
        EvaluateVectorFieldAtLocation(curPos, dims, X, Y, F, k4);


        //I must admit that I do not fully understand the equation below but I do understand the
        //idea behind using multiple vectors to calculate a more accurate vector and why it
        //may be necessary in some cases

        //use all the vectors found to derive more accurate vector than with Euler method
        output_locations[2*i+2] = output_locations[2*i] + sixth*h*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        output_locations[2*i+3] = output_locations[2*i+1] + sixth*h*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    }
}


// ****************************************************************************
//  Function: AdvectWithEulerStep
//
//  Arguments:
//     pt: the seed location (two-dimensions)
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a vector field defined on the mesh.  Its size is 2*dims[0]*dims[1].
//     h: The size of the Euler step
//     nsteps: The number of Euler steps to take
//     output_locations (output): An array of size 2*(nsteps+1).  It's first entry
//        should be the seed location.  The second entry should be the result
//        of the first advection step.  The final entry should be the result
//        of the final advection step.
//
// ****************************************************************************


void
AdvectWithEulerStep(const float *pt, const int *dims, const float *X, 
                    const float *Y, const float *F, 
                    float h, int nsteps, float *output_locations)
{
    output_locations[0] = pt[0];
    output_locations[1] = pt[1];
    
    for(int i=0; i<nsteps; ++i){
        float curVec[2];
        float curPos[2];    
        
        //get current position
        curPos[0] = output_locations[2*i];
        curPos[1] = output_locations[2*i+1];

        //get the interpolated vector field value at the current location
        EvaluateVectorFieldAtLocation(curPos, dims, X, Y, F, curVec);
        
        //advance the position based on the interpolated vector field value and the given time step
        output_locations[2*i+2] = curPos[0] + h*curVec[0];
        output_locations[2*i+3] = curPos[1] + h*curVec[1];

    }
}

//vtkImageData *
//NewImage(int width, int height)
//{
//    vtkImageData *image = vtkImageData::New();
//    image->SetDimensions(width, height, 1);
//    //image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
//    //image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
//    //image->SetNumberOfScalarComponents(3);
//    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
//    //image->AllocateScalars();
//
//    return image;
//}

// VTK files are only 3D, so the vector data is all of the form (X,Y,0).
// Remove the 0's since it is counter-intuitive for students who are 
// thinking of this as 2D data.
float *
Convert3DVectorDataTo2DVectorData(const int *dims, const float *F)
{
    float *rv = new float[dims[0]*dims[1]*2];
    int index3D = 0;
    int index2D = 0;
    for (int i = 0 ; i < dims[0] ; i++)
       for (int j = 0 ; j < dims[1] ; j++)
       {
           rv[index2D]   = F[index3D];
           rv[index2D+1] = F[index3D+1];
           index2D += 2;
           index3D += 3;
       }

    return rv;
}

void PrintSteps(const string solver, int nsteps, float *locations)
{
   ofstream out;
   out.open(solver);
   for (int j = 0 ; j < nsteps+1 ; j++)
   {
       out << j << ": (" << locations[2*j] << ", " << locations[2*j+1] << ")" << endl;
   }
   out.close();
}
