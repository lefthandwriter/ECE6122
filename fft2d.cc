// Distributed two-dimensional Discrete FFT transform
// ESTHER LING
// ECE8893 Project 1


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <signal.h>
#include <math.h>
#include <mpi.h>

#include "Complex.h"
#include "InputImage.h"

#include <stdlib.h> // for exit(0) on Linux GCC

using namespace std;

// Function declarations
void Transform2D(const char* inputFN);
void Transform1D(Complex* h, int w, int height, int startRow, Complex* H);
void revTransform1D(Complex* h, int w, int height, int startRow, Complex* H);


void Transform2D(const char* inputFN) 
{ // Do the 2D transform here.
  // 1) Use the InputImage object to read in the Tower.txt file and
  //    find the width/height of the input image.
  // 2) Use MPI to find how many CPUs in total, and which one
  //    this process is
  // 3) Allocate an array of Complex object of sufficient size to
  //    hold the 2d DFT results (size is width * height)
  // 4) Obtain a pointer to the Complex 1d array of input data
  // 5) Do the individual 1D transforms on the rows assigned to your CPU
  // 6) Send the resultant transformed values to the appropriate
  //    other processors for the next phase.
  // 6a) To send and receive columns, you might need a separate
  //     Complex array of the correct size.
  // 7) Receive messages from other processes to collect your columns
  // 8) When all columns received, do the 1D transforms on the columns
  // 9) Send final answers to CPU 0 (unless you are CPU 0)
  //   9a) If you are CPU 0, collect all values from other processors
  //       and print out with SaveImageData().

  // Step (1): Create the helper object for reading the image
  InputImage image(inputFN); 
  
  // Your code here, steps 2-9

  int width, height;
  width = image.GetWidth();
  height = image.GetHeight();
  Complex *h = new Complex[width*height];
  h = image.GetImageData(); /* every process (CPU) reads in the image */

  int numWorkers, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &numWorkers); // how many CPUs
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // which CPU am I
  
  // Determine how many rows each CPU is responsible for
  // Formula: Total Work to be Done / Number of Workers
  int numRows = height/numWorkers;
  // Determine the starting row for each CPU
  int startRow = numRows*rank;
  cout << "myRank is: " << rank << " and myStartRow is: " << startRow << endl;

  // Do the individual 1D transforms on the rows assigned to each CPU
  //Transform1D(h, width, numRows, startRow, rank, CPU_1D); 

  /* Now do the message passing */
  // Tag definitions: TAG = 0: 1D Transform
  //                  TAG = 1: 2D Transform
  if (rank == 0){ // CPU 0
      printf("width is: %d\n", width);
      printf("height is: %d\n", height);
      printf("numWorkers is: %d\n", numWorkers);
  
      Complex* CPU_1D  = new Complex[width*numRows];  // CPU 0 to store its 1D transform
      Complex* after1D = new Complex[width*height];   // CPU 0 to store ALL 1D transforms
      Complex* CPU_2D  = new Complex[width*numRows];  // CPU 0 to store its 2D transform
      Complex* after2D = new Complex[width*height];   // CPU 0 to Store ALL 2D Transforms
      
      int rc;
      MPI_Request* request1 = new MPI_Request[numWorkers-1];
      MPI_Request* request2 = new MPI_Request[numWorkers-1];
      MPI_Request* request3 = new MPI_Request[numWorkers-1];
      MPI_Request* request4 = new MPI_Request[numWorkers-1];
      
      MPI_Status status1[numWorkers-1];
      MPI_Status status2[numWorkers-1];
      MPI_Status status3[numWorkers-1];
      MPI_Status status4[numWorkers-1];
      //MPI_Request* request = new MPI_Request[(numWorkers-1)*2]; // 15 for DFT MPI_Irecv, 15 for DFT MPI_Isend
      
      // Set up non-blocking receives first before CPU 0 does the 1D transforms
      for (int r = 1; r < numWorkers; r++){
        // Specify starting address
        rc = MPI_Irecv(&after1D[r*width*numRows], width*numRows, MPI_C_DOUBLE_COMPLEX, r, 0,
                           MPI_COMM_WORLD, &request1[r-1]);
        //cout << "Rank 0 rc for MPI_Irecv: " << rc << endl;
        if (rc != MPI_SUCCESS){
          cout << "(DFT) Rank 0 MPI_Irecv failed, error code: " << rc << endl;
          cout << "Exiting now \n";
          MPI_Finalize();
          exit(1);
        }
      }

      // While waiting, do the individual 1D transforms on the rows assigned to each CPU
      Transform1D(h, width, numRows, startRow, CPU_1D); 
      // Store the results of CPU 0 1D transform
      for (int i = 0; i < width*numRows; i++){
          after1D[i] = CPU_1D[i];
      }
      
      // Wait for completion of MPI_Irecv in Line 102
      rc = MPI_Waitall(numWorkers-1, request1, status1);
      //cout << "Rank 0 rc for MPI_Waitall (Irecv): " << rc << endl;
      if (rc != MPI_SUCCESS){
         cout << "(DFT) Rank 0 MPI_Irecv failed, error code: " << rc << endl;
         cout << "Exiting now \n";
         MPI_Finalize();
         exit(1);
      }
      // Save 1D results to .txt file
      image.SaveImageData("MyAfter1D.txt", after1D, width, height);

      //cout << "Rank 0 TEST 0" << endl;

      // Transpose after1D to after1D_t
      Complex* after1D_t = new Complex[width*height];
      for (int i = 0; i < width; i++){
        for (int j = 0; j < width; j++)
          after1D_t[j*width + i] = after1D[i*width + j];
      }
      //cout << "Rank 0 TEST 1" << endl;

      // CPU 0 non-blocking sends [1x65536] to CPU 1-15
      for (int i = 1; i < numWorkers; i++){
        rc = MPI_Isend(after1D_t, width*height, MPI_C_DOUBLE_COMPLEX, i, 
                        1, MPI_COMM_WORLD, &request2[i-1]);
        //cout << "Rank 0 rc for MPI_Isend: " << rc << endl;
        if (rc != MPI_SUCCESS){
          cout << "(DFT) Rank 0 MPI_Isend failed, error code: " << rc << endl;
          cout << "Exiting now \n";
          MPI_Finalize();
          exit(1);
        }
      }
      //cout << "Rank 0 TEST 2" << endl;
      // CPU 0 does its 2D transforms
      Transform1D(after1D_t, width, numRows, startRow, CPU_2D); 
      // Store the results of CPU 0 2D transforms
      for (int i = 0; i < width*numRows; i++){
          after2D[i] = CPU_2D[i];
      }
      //cout << "Rank 0 TEST 3" << endl;
      
      // Wait for MPI_Isend to complete
      rc = MPI_Waitall(numWorkers-1, request2, status2);
      //cout << "Rank 0 rc for MPI_Waitall (Isend): " << rc << endl;
      if (rc != MPI_SUCCESS){
         cout << "(DFT) Rank 0 MPI_Isend failed, error code: " << rc << endl;
         cout << "Exiting now \n";
         MPI_Finalize();
         exit(1);
      }
      // Waitall effectively waits for all the Requests (Isend and Irecv)
      //MPI_Status status_wait2[numWorkers-1];
      //rc = MPI_Waitall(numWorkers-1, request, status_wait2);
      
      // Blocking Receive 2D transforms (Tag 1)
      MPI_Status status;
      for (int r = 1; r < numWorkers; r++){
        rc = MPI_Recv(&after2D[r*width*numRows], width*numRows, MPI_C_DOUBLE_COMPLEX, r, 2,
                         MPI_COMM_WORLD, &status);
        if (rc != MPI_SUCCESS){
          cout << "Rank " << rank << "Recv failed, rc " << rc << endl;
          MPI_Finalize();
          exit(1);
        }
        int count = 0;
        MPI_Get_count(&status, MPI_C_DOUBLE_COMPLEX, &count);
        cout << "(DFT) Rank 0 received " << count << " data units from " << status.MPI_SOURCE << endl;
      }
      //cout << "Rank 0 TEST 4" << endl;
      // Save 1D transpose to .txt file
      //image.SaveImageData("MyAfter1D_Transpose.txt", after1D_t, width, height);
      
      // Transpose back to original
      Complex* final2D = new Complex[width*height];
      for (int i = 0; i < width; i++){
        for (int j = 0; j < width; j++)
          final2D[j*width + i] = after2D[i*width + j];
      }
      //cout << "Rank 0 TEST 5" << endl;
      
      // Save 2D results to .txt file
      image.SaveImageData("MyAfter2D.txt", final2D, width, height);
      //cout << "Rank 0 TEST 6" << endl;

      delete [] h;
      delete [] CPU_1D;
      delete [] after1D_t;
      delete [] after2D;
      delete [] CPU_2D;

      // Compute inverse DFT
      // The DFT image is final2D
      Complex* IDFT_1D       = new Complex[width*numRows];  // CPU 0 to store its 1D IDFT
      Complex* IDFT_1D_ALL   = new Complex[width*height];   // CPU 0 to store ALL 1D IDFT
      Complex* IDFT_1D_ALL_t = new Complex[width*height];   // CPU 0 to sotre ALL 1D IDFT Tranpose
      Complex* IDFT_2D       = new Complex[width*numRows];  // CPU 0 to store its 2D IDFT
      Complex* IDFT_2D_ALL   = new Complex[width*height];   // CPU 0 to Store ALL 2D IDFT
      Complex* IDFT_2D_FINAL = new Complex[width*height];   // CPU 0 to sotre ALL 2D IDFT Tranpose

      //cout << "Rank 0 TEST 7" << endl;
      
      // CPU 0 sends DFT image to Ranks 1 -15
      for (int r = 1; r < numWorkers; r++){
        rc = MPI_Isend(final2D, width*height, MPI_C_DOUBLE_COMPLEX, r, 3, MPI_COMM_WORLD, &request3[r-1]);
        cout << "Rank 0 rc for MPI_Isend IDFT: " << rc << endl;
        if (rc != MPI_SUCCESS){
          cout << "(IDFT) Rank 0 MPI_Isend failed, error code: " << rc << endl;
          cout << "Exiting now \n";
          MPI_Finalize();
          exit(1);
        }
      }
      // CPU 0 does its IDFT row transforms
      //Transform1D(after1D_t, width, numRows, startRow, CPU_2D); 
      revTransform1D(final2D, width, numRows, startRow, IDFT_1D);
      // Store the results of CPU 0 1D IDFT transforms
      for (int i = 0; i < width*numRows; i++){
          IDFT_1D_ALL[i] = IDFT_1D[i];
      }
      //cout << "Rank 0 TEST 8" << endl;

      // Wait for the CPU 0 MPI_Isend to complete
      rc = MPI_Waitall(numWorkers-1, request3, status3);
      if (rc != MPI_SUCCESS){
         cout << "(IDFT) Rank 0 MPI_Isend failed, error code: " << rc << endl;
         cout << "Exiting now \n";
         MPI_Finalize();
         exit(1);
      }

      // Blocking Receive row IDFT from CPU 0-15
      for (int r = 1; r < numWorkers; r++){
         rc = MPI_Recv(&IDFT_1D_ALL[r*width*numRows], width*numRows, MPI_C_DOUBLE_COMPLEX, r, 4, MPI_COMM_WORLD, &status);
         if (rc != MPI_SUCCESS){
          cout << "(IDFT) Rank 0 MPI_Irecv failed, error code: " << rc << endl;
          cout << "Exiting now \n";
          MPI_Finalize();
          exit(1);
        }
        int count = 0;
        MPI_Get_count(&status, MPI_C_DOUBLE_COMPLEX, &count);
        cout << "(1D IDFT) Rank 0 received " << count << " data units from " << status.MPI_SOURCE << endl;
      }

      // Transpose the Row IDFT
      for (int i = 0; i < width; i++){
        for (int j = 0; j < width; j++)
          IDFT_1D_ALL_t[j*width + i] = IDFT_1D_ALL[i*width + j];
      }

      // Non-blocking send of tranpose to CPU 1-15
      for (int r = 1; r < numWorkers; r++){
        rc = MPI_Isend(IDFT_1D_ALL_t, width*height, MPI_C_DOUBLE_COMPLEX, r, 5, MPI_COMM_WORLD, &request4[r-1]);
        if (rc != MPI_SUCCESS){
          cout << "(IDFT) Rank 0 MPI_Isend failed, error code: " << rc << endl;
          cout << "Exiting now \n";
          MPI_Finalize();
          exit(1);
        }
      }

      // CPU 0 does its IDFT col transforms
      revTransform1D(IDFT_1D_ALL_t, width, numRows, startRow, IDFT_2D);
      for (int i = 0; i < width*numRows; i++){
        IDFT_2D_ALL[i] = IDFT_2D[i];
      }

      // Wait for Isend to finish
      rc = MPI_Waitall(numWorkers-1, request4, status4);

      // Receive 2D DFT (1x4096)
      for (int r = 1; r < numWorkers; r++){
        rc = MPI_Recv(&IDFT_2D_ALL[r*width*numRows], width*numRows, MPI_C_DOUBLE_COMPLEX, r, 6, MPI_COMM_WORLD, &status);
        if (rc != MPI_SUCCESS){
          cout << "(IDFT) Rank 0 MPI_Recv failed, error code: " << rc << endl;
          cout << "Exiting now \n";
          MPI_Finalize();
          exit(1);
        }
        int count = 0;
        MPI_Get_count(&status, MPI_C_DOUBLE_COMPLEX, &count);
        cout << "(2D IDFT) Rank 0 received " << count << " data units from " << status.MPI_SOURCE << endl;
      }

      // Tranpose 2D IDFT back
      for (int i = 0; i < width; i++){
        for (int j = 0; j < width; j++)
          IDFT_2D_FINAL[j*width + i] = IDFT_2D_ALL[i*width + j];
      }

      // Store the results
      image.SaveImageDataReal("MyAfterInverse.txt", IDFT_2D_FINAL, width, height);

      delete [] final2D;
      delete [] IDFT_1D_ALL;
      delete [] IDFT_1D;
      delete [] IDFT_2D_ALL;
      delete [] IDFT_2D;      
  }

  else{ // all other ranks:
      
      // DFT
      // 1. Blocking send of 1D DFT of Rows to CPU 0 (Tag 0)
      // 2. Blocking receive of columns from CPU 0   (Tag 1)
      // 3. Blocking send of 2D DFT of Cols to CPU 0 (Tag 2)

      // IDFT:
      // 1. Blocking receive of 2D DFT image from CPU 0 (Tag 3) ( 1x65536 )
      // 2. Blocking send of 1D IDFT to CPU 0           (Tag 4)
      // 3. Blocking receive of columns from CPU 0      (Tag 5)
      // 4. Blocking send of 2D IDFT to CPU 0           (Tag 6)
      int rc;

      // Allocate pointer array (buffer) for each CPU to store the results of:
      Complex *CPU_1D = new Complex[width*numRows]; //  1D Transform
      Complex *CPU_2D = new Complex[width*numRows]; //  2D Transform
      Complex *after1D_t = new Complex[width*height]; // incoming transpose, after1D_t

      // Do the individual 1D transforms on the rows assigned to each CPU
      Transform1D(h, width, numRows, startRow, CPU_1D); 

      //MPI_Request request; // only CPU 0 needs the array of requests
      MPI_Status status;

      // Blocking Send row transform results to CPU 0 (Tag 0)
      rc = MPI_Send(CPU_1D, width*numRows, MPI_C_DOUBLE_COMPLEX, 0, 
                      0, MPI_COMM_WORLD);
      if (rc != MPI_SUCCESS){
          cout << "(DFT) Rank" << rank << " MPI_Send failed, error code: " << rc << endl;
          cout << "Exiting now \n";
          MPI_Finalize();
          exit(1);
        }
      //cout << "Rank: " << rank << " line 1" << endl;

      // Blocking Receive transposed columns (1x65536) from CPU 0
      rc = MPI_Recv(after1D_t, width*height, MPI_C_DOUBLE_COMPLEX, 0,
                     1, MPI_COMM_WORLD, &status);
      if (rc != MPI_SUCCESS){
          cout << "(DFT) Rank " << rank << " MPI_Recv failed, error code: " << rc << endl;
          cout << "Exiting now \n";
          MPI_Finalize();
          exit(1);
        }
      // Debugging
      //int count = 0;
      //MPI_Get_count(&status, MPI_C_DOUBLE_COMPLEX, &count);
      //cout << "Rank " << rank << " received " << count << " data units from rank " 
      //    << status.MPI_SOURCE << endl;
        
      //cout << "Rank: " << rank << " line 2" << endl;

      // Do the 2D transforms
      Transform1D(after1D_t, width, numRows, startRow, CPU_2D);

      //cout << "Rank: " << rank << " line 3" << endl;
      
      // Blocking Send the column transforms to CPU 0 
      rc = MPI_Send(CPU_2D, width*numRows, MPI_C_DOUBLE_COMPLEX, 0,
                     2, MPI_COMM_WORLD);
      if (rc != MPI_SUCCESS){
          cout << "(DFT) Rank" << rank << " MPI_Send failed, error code: " << rc << endl;
          cout << "Exiting now \n";
          MPI_Finalize();
          exit(1);
        }
      //cout << "Rank: " << rank << " line 4" << endl;

      // IDFT:
      Complex* IDFT      = new Complex[width*height];
      Complex* IDFT_t    = new Complex[width*height];
      Complex* IDFT_1D   = new Complex[width*numRows];  // CPU 1-15 to store its 1D IDFT
      Complex* IDFT_2D   = new Complex[width*numRows];  // CPU 1-15 to store its 2D IDFT
      //cout << "Rank: " << rank << " line 5" << endl;
      
      // Blocking receive of DFT Image (1x65536) from CPU 0 (Tag 3)
      rc = MPI_Recv(IDFT, width*height, MPI_C_DOUBLE_COMPLEX, 0, 3, MPI_COMM_WORLD, &status);
      if (rc != MPI_SUCCESS){
          cout << "(IDFT) Rank" << rank << " MPI_Recv failed, error code: " << rc << endl;
          cout << "Exiting now \n";
          MPI_Finalize();
          exit(1);
        }
      int count = 0;
      MPI_Get_count(&status, MPI_C_DOUBLE_COMPLEX, &count);
      cout << "(IDFT) Rank " << rank << " received " << count << " data units from rank " 
          << status.MPI_SOURCE << endl;
      //cout << "Rank: " << rank << " line 6" << endl;
      
      // Do revTransform 
      revTransform1D(IDFT, width, numRows, startRow, IDFT_1D);
      
      // Send ROWS revTransform to CPU 0 (Tag 4)
      rc = MPI_Send(IDFT_1D, width*numRows, MPI_C_DOUBLE_COMPLEX, 0, 
                    4, MPI_COMM_WORLD);
      if (rc != MPI_SUCCESS){
          cout << "(IDFT) Rank " << rank << " MPI_Send failed, error code: " << rc << endl;
          cout << "Exiting now \n";
          MPI_Finalize();
          exit(1);
        }

      // Receive COLS (1x65536) from CPU 0 
      rc = MPI_Recv(IDFT_t, width*height, MPI_C_DOUBLE_COMPLEX, 0, 5, MPI_COMM_WORLD, &status);
      if (rc != MPI_SUCCESS){
          cout << "(IDFT) Rank" << rank << " MPI_Recv failed, error code: " << rc << endl;
          cout << "Exiting now \n";
          MPI_Finalize();
          exit(1);
        }
      cout << "(IDFT) Rank " << rank << " received " << count << " data units from rank " 
          << status.MPI_SOURCE << endl;
      
      // Do revTransform 
      revTransform1D(IDFT_t, width, numRows, startRow, IDFT_2D);

      // Send IDFT cols to CPU 0
      rc = MPI_Send(IDFT_2D, width*numRows, MPI_C_DOUBLE_COMPLEX, 0,
                    6, MPI_COMM_WORLD);
      if (rc != MPI_SUCCESS){
          cout << "(IDFT) Rank "<< rank << " MPI_Send failed, error code: " << rc << endl;
          cout << "Exiting now \n";
          MPI_Finalize();
          exit(1);
        }
      
      delete [] CPU_1D;
      delete [] CPU_2D;
      delete [] after1D_t;
      delete [] IDFT;
      delete [] IDFT_1D;
      delete [] IDFT_2D;
      delete [] IDFT_t;
      
  }
}

void Transform1D(Complex* h, int w, int height, int startRow, Complex* H)
{
  // Implement a simple 1-d DFT using the double summation equation
  // given in the assignment handout.  h is the time-domain input
  // data, w is the width (N), and H is the output array.

  // Edited function to receive input height (track number of rows)
  // Edited function to receive startRow
  
  // Height, height: number of rows
  // Width, w of each row is 256
  // Incoming array of pointers, h is 1 x 65536
  // Outgoing array of pointers, H is 1 x 4096   (256*16)
  
  //H[0] = 0;
  for (int row = 0; row < height; row++){
    // Computing H[n] for each element in the row
    for (int n = 0; n < w; n++){
      // Iterating h[k] over the row; from column 0-255
      for (int k = 0; k < w; k++){
        H[row*w + n] = H[row*w + n] + Complex(+cos(2.0 * M_PI * n * k / w),
                                    -sin(2.0 * M_PI * n * k / w)) 
                                    * h[(row+startRow)*w + k];
      }
      // This should be the right index now; we get the final value of an H[n] after
      // exiting the k loop 
      // Complex cc checks this so no need to put it
      if (fabs(H[row*w + n].imag) < 1e-10) H[row*w + n].imag = 0;
      if (fabs(H[row*w + n].real) < 1e-10) H[row*w + n].real = 0; 
    }
  }
}
void revTransform1D(Complex*  h, int w, int height, int startRow, Complex* H)
{
  // Reverse DFT 1D
  // Height, height: number of rows
  // Width, w of each row is 256
  // Incoming array of pointers, H is 1 x 65536
  // Outgoing array of pointers, h is 1 x 4096   (256*16)
  double N= (double) w;
  for (int row = 0; row < height; row++){
    for (int n = 0; n < w; n++){
      for (int k = 0; k < w; k++){
        H[row*w + n] = H[row*w + n] + Complex(+cos(2.0 * M_PI * n * k / w),
                                    + sin(2.0 * M_PI * n * k / w)) 
                                    * h[(row+startRow)*w + k];
      }
      H[row*w + n] = Complex(1.0/N) * H[row*w + n];
      if (fabs(H[row*w + n].imag) < 1e-10) H[row*w + n].imag = 0;
      if (fabs(H[row*w + n].real) < 1e-10) H[row*w + n].real = 0; 
    }
  }
}
int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line

  // MPI initialization here
  int rc, rank;
  rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // which CPU am I
  // Perform the 2d transform
  Transform2D(fn.c_str()); 

  // Finalize MPI here
  cout << "Rank " << rank << " Program exiting normally" << endl; 
  MPI_Finalize();
}