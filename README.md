
### How To Generate The Output Files

The output files that were generated using my virtual machine are included in the My_Output_Files folder and are there for reference.
In order to produce these files on your machine I have included a bash script called Generate_Output_Files.sh. Running this may take
some time especially on systems with lower core counts. It will produce 4 files for 1k 5k 10k and 15k particle counts. The output files
will contain information for 1, 2, 4, 6, and 8 thread count test runs at their respective particle counts. Alternatively you can 
./final_proj stepCount stepSize numParticles numThreads to test at individual thread counts. However this method will still run the test
at for static dynamic and guided scheduling. Each of these scheduling types will run 20 each and the average run time will be the output
to the file. 

### Notes On Setup

1. The CMakeLists.txt should be modified to reflect the location of your vtk build directory

2. Run the following command to turn the bash script into an executable:
        
        chmod u+x Generate_Output_Files.sh

3. I have generated an array of 30k random seed locations this must be in the same directory as
   the program.

4. When generating files the program will append to an existing file with the same particle count
   if it exists. So if you are not using the provided bash script I recommend running from lowest
   thread count to the highest for easier reading of the output file. 

   4a. The thread count is capped at 10 to account for accidental swapping 
       of particle count and thread count

5. Your mileage will vary. I ran my tests on a virtual machine with 4 cores allocated so the numbers
   obtained by running on a different system will almost certainly be different especially if the system
   has less than 4 cores. So long as the system has at least two cores some speed up should still be achievable.

