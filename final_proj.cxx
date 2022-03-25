// See main.h for all of the includes used
#include "main.h"

using namespace std;


void usage()
{
    cout<<"Usage is ./proj4 nsteps stepSize particleCount numThreads"<<endl;
}

int main(int argc, char* argv[])
{
    if(argc != 5){
        usage();
        exit(EXIT_FAILURE);
    }

    const int nsteps = atoi(argv[1]);
    float h = atof(argv[2]);
    int numParticles = atoi(argv[3]);
    int numThreads = atoi(argv[4]);

    string fileName = to_string(numParticles)+"_Particles_Output.txt";
    fstream outFile;
    outFile.open(fileName, fstream::out | fstream::app);

    // Thread count safety check to ensure proper use of the program and protect the system being used
    if(numThreads > 10){
        cout<<"Thread count is too high, please be sure particle count and thread count were not swapped"<<endl;
        exit(EXIT_FAILURE);
    }

    int  i, j;

    //read data file
    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj4_data.vtk");
    rdr->Update();

    //create grid
    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    if (dims[0] <= 0 || dims[1] <= 0)
    { cerr << "Was not able to successfully open file \"proj4_data.vtk\"" << endl;
        exit(EXIT_FAILURE);
    }
    
    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F_3D = (float *) rgrid->GetPointData()->GetVectors()->GetVoidPointer(0);
    float *F = Convert3DVectorDataTo2DVectorData(dims, F_3D);
    
    float RK4_output_locations[2*(nsteps+1)];
    double acc_rTime = 0.0;

    // This is the serial implementation and is separated so as to avoid any uncessessary OpenMP initialization being
    // included in the final average run time
    if(numThreads == 1){
      
      outFile<<"Serial average of 20 runs is: ";
      for(i=0; i<20; ++i){

            auto startTime_serial = chrono::steady_clock::now();
            for(j=0; j<numParticles; j+=2){

                AdvectWithRK4Step(&seeds[j], dims, X, Y, F, h, nsteps, RK4_output_locations);
            }
            chrono::duration<double> totalTime_serial = chrono::steady_clock::now() - startTime_serial;
            acc_rTime += totalTime_serial.count();
      }
      outFile<<acc_rTime/20.0<<endl<<endl;
      outFile.close();

    }

    // This is the parallel implementation which runs 3 different tests a static schedule, a dynamic schedule and
    // a guided schedule. Each of these tests are run 20 times and the average of the 20 runs is what is written to
    // the output file
    else{

        omp_set_num_threads(numThreads);

        outFile<<"Static schedule and "<<numThreads<<" threads average of 20 runs is: ";
        for(i=0; i<20; ++i){

            auto startTime_parallel = chrono::steady_clock::now();
            
            #pragma omp parallel for schedule(static) private(RK4_output_locations)
            for(j=0; j<numParticles; j+=2){        

                AdvectWithRK4Step(&seeds[j], dims, X, Y, F, h, nsteps, RK4_output_locations);
            }

            chrono::duration<double> totalTime_parallel = chrono::steady_clock::now() - startTime_parallel;
            acc_rTime += totalTime_parallel.count();
        }
        outFile<<acc_rTime/20.0<<endl;

        acc_rTime = 0.0;
        outFile<<"Dynamic schedule and "<<numThreads<<" threads average of 20 runs is: ";
        for(i=0; i<20; ++i){

            auto startTime_parallel = chrono::steady_clock::now();
            
            #pragma omp parallel for schedule(dynamic) private(RK4_output_locations)
            for(j=0; j<numParticles; j+=2){        

                AdvectWithRK4Step(&seeds[j], dims, X, Y, F, h, nsteps, RK4_output_locations);
            }

            chrono::duration<double> totalTime_parallel = chrono::steady_clock::now() - startTime_parallel;
            acc_rTime += totalTime_parallel.count();
        }
        outFile<<acc_rTime/20.0<<endl;

        acc_rTime = 0.0;
        outFile<<"Guided schedule and "<<numThreads<<" threads average of 20 runs is: ";
        for(i=0; i<20; ++i){

            auto startTime_parallel = chrono::steady_clock::now();
            
            #pragma omp parallel for schedule(guided) private(RK4_output_locations)
            for(j=0; j<numParticles; j+=2){        

                AdvectWithRK4Step(&seeds[j], dims, X, Y, F, h, nsteps, RK4_output_locations);
            }

            chrono::duration<double> totalTime_parallel = chrono::steady_clock::now() - startTime_parallel;
            acc_rTime += totalTime_parallel.count();
        }
        outFile<<acc_rTime/20.0<<endl<<endl;
        outFile.close();
    }



    delete[] X;
    delete[] Y;
    delete[] F;
    delete[] F_3D;
}
