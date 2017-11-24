// Threaded two-dimensional Discrete FFT transform
// ESTHER LING
// ECE8893 Project 2

// Custom barrier 

#include <iostream>
#include <string>
#include <math.h>
#include <stdlib.h> // for abs and uint in Linux
#include <pthread.h>
#include <sys/time.h>

#include "Complex.h"
#include "InputImage.h"

using namespace std;

// You will likely need global variables indicating how
// many threads there are, and a Complex* that points to the
// 2d image being transformed.
Complex* h;       // Image array data
Complex* H;       // Re-ordered h array and transform in place
int height;       // height of image
int N;            // width of image
Complex* weights; // store pre-computed weights 

// Pthreads-Related Variables
int nThreads =  16; // Number of threads
pthread_mutex_t exitMutex; // Final exit mutex 
pthread_cond_t  exitCond;  // Final exit condition
pthread_mutex_t dft;

int threadCtr = 0; // Global counter for barrier
pthread_mutex_t elementCount;

// Custom Barrier related variables
bool* localSense;
bool globalSense;
pthread_cond_t barrierCond;
pthread_mutex_t barrierMutex; // cond and mutex = peas and carrots
pthread_mutex_t flagMutex;
int saveID;
bool barrierFlag = true;

int GetMilliseconds(){ // from ThreadedCount-Again.cc
  timeval tv;
  gettimeofday(&tv, 0);
  static bool first = true;
  static int startSec = 0;
  if (first){
    startSec = tv.tv_sec;
    first = false;
  }
  // Time in milliseconds
  return (tv.tv_sec - startSec) * 1000 + tv.tv_usec / 1000;
}

// Function to reverse bits in an unsigned integer
// This assumes there is a global variable N that is the
// number of points in the 1D transform.
unsigned ReverseBits(unsigned v)
{ //  Provided to students
  unsigned n = N; // Size (length) of array - (which is even - a 2 power k value)
  //unsigned n = width;
  unsigned r = 0; // Return value
   
  for (--n; n > 0; n >>= 1)
    {
      r <<= 1;        // Shift return value
      r |= (v & 0x1); // Merge in next bit
      v >>= 1;        // Shift reversal value
    }
  return r;
}

// GRAD Students implement the following 2 functions.
// Undergrads can use the built-in barriers in pthreads.

// Call MyBarrier_Init once in main and once by last thread to reset threadCtr
void MyBarrier_Init()// you will likely need some parameters)
{
    threadCtr = nThreads; // make threadCtr = 16
    barrierFlag = true;
    //pthread_mutex_init(&elementCount, 0);
    //pthread_cond_init(&barrierCond, 0);
    //pthread_mutex_init(&barrierMutex, 0);
  
}

// Each thread calls MyBarrier after completing the row-wise DFT
// All threads must enter the barrier before moving to next stage
// This is a blocking call
void MyBarrier(int myId) // Again likely need parameters
{
  //cout << "ID:" << myId << "before lock" << endl;
  pthread_mutex_lock(&barrierMutex);
  // is this the first time entering?
  /*if (barrierFlag == true){ // will two threads see 0 at the same time and deadlock?
    pthread_mutex_trylock(&elementCount); // first to enter locks elementCount
    pthread_mutex_trylock(&flagMutex); // use trylock to avoid deadlock
    barrierFlag = false;
    saveID = myId;
    pthread_mutex_unlock(&flagMutex);
    cout << "First thread: " << saveID << endl;
  }*/
  
  threadCtr--;
  //cout << "ID:" << myId << "after lock" << endl;  
  if (threadCtr == 0){
      // All threads have called the barrier
      //cout << "All threads hit the barrier." << endl;
      // time saving: last one gets to transpose the matrix
      for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
          h[j*N + i] = H[i*N + j]; // H_t is the new input to Transform1D()
        }
      }
      cout << "Last Thread: " << myId << " resetting MyBarrier_Init.." << endl;
      pthread_mutex_unlock(&barrierMutex);
      MyBarrier_Init(); // last one resets barrier parameters
          pthread_mutex_lock(&elementCount);
      pthread_cond_broadcast(&barrierCond);
          pthread_mutex_unlock(&elementCount);
      //pthread_mutex_unlock(&barrierMutex);
      cout << "Last Thread: " << myId << " unlocked" << endl;

  }
  else{
      pthread_mutex_unlock(&barrierMutex); // release the mutex, and then spin    
      cout << "Thread " << myId << " unlocked. Start spin." << endl;
      // Upon successful completion, a value of zero shall be  returned
      // pthread_cond_wait() blocks the calling thread until the specified condition is signalled
      pthread_mutex_lock(&elementCount);
      while(pthread_cond_wait(&barrierCond, &elementCount) != 0);
      pthread_mutex_unlock(&elementCount);
      /*if (myId == saveID){
        // if this was the first to enter, unlock elementCount
        pthread_mutex_unlock(&elementCount);
        cout << "First thread " << saveID << " unlocked elementCount" << endl;
      }*/
  }
}

//void Transform1D(int N, int height, int startRow, int myId)
void Transform1D(Complex* h, int N, int height, int startRow, Complex* H, int myId )
//void Transform1D(Complex* h, int N, int height)
{
  // Implement the efficient Danielson-Lanczos DFT here.
  // "h" is an input/output parameter
  // "N" is the size of the array (assume even power of 2)

  // Validate Reverse bits function
  //cout << "Reverse bit of integer 6 (00 0000 0110) -> (01 1000 0000) is: " << ReverseBits(uint(6)) << endl;
  //printf("Thread: %d startRow: %d\n", myId, startRow);

  // Step 1: Re-order h array
  for (int r = startRow; r < height+startRow; r++){
    for (int c = 0; c < N; c++){
      H[r*N + c] = h[r*N + ReverseBits(uint(c))];
    }
  }
  // Step 2: Pre-compute Wn values
  weights = new Complex[N/2];
  // // N distinct weights, all the same regardless of ROW
  for (int n = 0 ; n < N/2; n++){
    weights[n] = Complex(+cos(2*M_PI*n/N), -sin(2*M_PI*n/N));
  }
  // Step 3: Do the Transform (in-place)
  int numPtTransforms = log10(N)/log10(2);
  int w_index = 0;
  //printf("numPtTransforms: %d\n", numPtTransforms);
  int s = 0;
  for (int r = startRow; r < height+startRow; r++){
    for (int p = 0; p < numPtTransforms; p++){
      int x_point = pow(2,p+1);
      Complex* hold = new Complex[x_point/2]; // temp array of different size based on x_point
      for (int n = 0; n < N; n += x_point){ // loop over the whole row in jumps of x_point
        int save = 0;
        for (int k = 0; k < x_point; k++){ // loop over x_point
          s = n + k;
          w_index = k*N/x_point;
          if (w_index >= N/2){ // Set w_index's limit to x_point/2
            w_index = abs ((w_index % N) - N/2);
            if (w_index == N/2)
              w_index = 0;
          }
          if (k < x_point/2){ // less than halfway of the subgroup
             hold[k] = H[r*N + s]; // only save if less than x_point/2 else seg11
             H[r*N + s] = H[r*N + s] + weights[w_index] * H[r*N + s + x_point/2];
          } 
          else{ // over halfway point
             H[r*N + s] = hold[save] - weights[w_index] * H[r*N + s];
             save++;
          }
        }
      }
    }
    //cout << "row: " << r << endl;
  }
  //printf("Thread: %d transforming \n", myId);
}

void RevTransform1D(Complex* h, int N, int height, int startRow, Complex* H, int myId )
{
  // Step 1: Re-order h array
  for (int r = startRow; r < height+startRow; r++){
    for (int c = 0; c < N; c++){
      H[r*N + c] = h[r*N + ReverseBits(uint(c))];
    }
  }
  // Step 2: Pre-compute Wn values
  weights = new Complex[N/2];
  for (int n = 0 ; n < N/2; n++){
    weights[n] = Complex(+cos(2*M_PI*n/N), +sin(2*M_PI*n/N)); // change: +sin
  }
  // Step 3: Do the Transform (in-place)
  int numPtTransforms = log10(N)/log10(2);
  int w_index = 0;
  int s = 0;
  //bool status = true;
  for (int r = startRow; r < height+startRow; r++){
    for (int p = 0; p < numPtTransforms; p++){
      int x_point = pow(2,p+1);
      Complex* hold = new Complex[x_point/2]; // temp array of different size based on x_point
      for (int n = 0; n < N; n += x_point){ // loop over the whole row in jumps of x_point
        int save = 0;
        for (int k = 0; k < x_point; k++){ // loop over x_point
          s = n + k;
          w_index = k*N/x_point;
          if (w_index >= N/2){ // Set w_index's limit to x_point/2
            w_index = abs ((w_index % N) - N/2);
            if (w_index == N/2)
              w_index = 0;
          }
          if (k < x_point/2){ // less than halfway of the subgroup
             hold[k] = H[r*N + s]; // only save if less than x_point/2 else seg11
             H[r*N + s] = H[r*N + s] + (weights[w_index] * H[r*N + s + x_point/2]);
          } 
          else{ // over halfway point
             H[r*N + s] = hold[save] - (weights[w_index] * H[r*N + s]);
             save++;
          }
          if (p < 4){ // (x_point == 2 || x_point == 4 || x_point == 8 || x_point == 16)
             // 2*4*8*16 = 1024
             H[r*N + s] = Complex(1.0/x_point) * H[r*N + s];
          }
        }
      }
    }
  }
}

// Transform2D thread routine
void* Transform2DTHread(void* v)
{ // This is the thread starting point.  "v" is the thread number
  // Calculate 1d DFT for assigned rows
  // wait for all to complete
  // Calculate 1d DFT for assigned columns
  // Decrement active count and signal main if all complete

  unsigned long myId = (unsigned long)v; // which thread am I?
  int rowsPerThread = height / nThreads;
  int startRow = (int)myId * rowsPerThread;
  double my_start, my_finish;

  // Task 1: 1d DFT for rows
  pthread_mutex_lock(&dft);
  Transform1D(h, N, rowsPerThread, startRow, H, (int)myId); // H is re-ordered array
  pthread_mutex_unlock(&dft);
  
  // Call barrier - inclusive of matrix transpose
  MyBarrier((int)myId);
  //pthread_mutex_unlock(&barrierMutex);
  //cout << "Thread " << myId << " done row DFT"<< endl;
  // Task 2: 1d DFT for cols
  pthread_mutex_lock(&dft);
  //my_start = GetMilliseconds();
  Transform1D(h, N, rowsPerThread, startRow, H, (int)myId); // H_t is the transposed matrix
  //my_finish = GetMilliseconds();
  pthread_mutex_unlock(&dft);
  //cout << "Thread: " << myId << " took " << my_finish-my_start << "ms for 1D DFT\n";
  
  MyBarrier((int)myId);  
  //pthread_mutex_unlock(&barrierMutex);
  
  //cout << "Thread " << myId << " done col DFT"<< endl;
  // At the end, signal main we're done
  pthread_mutex_lock(&exitMutex);
  pthread_cond_signal(&exitCond);
  pthread_mutex_unlock(&exitMutex);
  
  return 0;
}

// Separate routine for IDFT
void* Transform2IDTHread(void* v)
{ 
  unsigned long myId = (unsigned long)v; // which thread am I?
  int rowsPerThread = height / nThreads;
  int startRow = (int)myId * rowsPerThread;
  
  // Task 3: 1d IDFT for rows
  pthread_mutex_lock(&dft);
  RevTransform1D(h, N, rowsPerThread, startRow, H, (int)myId);
  pthread_mutex_unlock(&dft);
  //cout << "test" << endl;
  MyBarrier((int)myId);
  //pthread_mutex_unlock(&barrierMutex);

  //cout << "Thread " << myId << " done row IDFT"<< endl;

  // Task 4: 1d IDFT for cols
  pthread_mutex_lock(&dft);
  RevTransform1D(h, N, rowsPerThread, startRow, H, (int)myId);
  pthread_mutex_unlock(&dft);
  
  MyBarrier((int)myId);
  //pthread_mutex_unlock(&barrierMutex);

  //cout << "Thread " << myId << " done row IDFT"<< endl;
  
  // At the end, signal main we're done
  pthread_mutex_lock(&exitMutex);
  pthread_cond_signal(&exitCond);
  pthread_mutex_unlock(&exitMutex);
  
  return 0;
}

int main(int argc, char** argv)
{  
  string fn("Tower.txt"); // default file name
  //string fn("TowerSmall.txt");
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  
  //Transform2D(fn.c_str()); // Perform the transform.
  InputImage image(fn.c_str());  // Create the helper object for reading the image
  // Get image dimensions
  N = image.GetWidth();
  height = image.GetHeight();
  printf("Image width: %d\n", N);
  printf("Image height: %d\n", height);
  h   = new Complex[N*height];
  h   = image.GetImageData(); // read in image data
  H   = new Complex[N*height];
  
  // Create 16 threads
  // Initialize all mutex and condition variables
  pthread_mutex_init(&exitMutex, 0);
  pthread_mutex_init(&dft, 0); 
  pthread_cond_init(&exitCond, 0);
  pthread_cond_init(&barrierCond, 0);
  pthread_mutex_init(&barrierMutex, 0);
  pthread_mutex_init(&elementCount, 0);
  
  // Let main hold the exit mutex while waiting for exitCond condition
  pthread_mutex_lock(&exitMutex);
  MyBarrier_Init();

  double my_start, my_finish;  
  my_start = GetMilliseconds();
  // Create 16 threads for forward DFT
  for (int i = 0; i < nThreads; i++){
    pthread_t pt;
    pthread_create(&pt, 0 ,Transform2DTHread, (void*)i);
  }

  // Wait for all threads complete
  pthread_cond_wait(&exitCond, &exitMutex);
  pthread_mutex_unlock(&exitMutex);

  // Save 1D results to .txt file
  //image.SaveImageData("MyAfter1D.txt", H, N, height);

  // Save 2D results to .txt file
  cout << "..... writing 2d results to file......" << endl;
  image.SaveImageData("MyAfter2D.txt", h, N, height);
  
  MyBarrier_Init();
  pthread_mutex_lock(&exitMutex);
  // Create 16 threads for reverse DFT
  // This works but is slower
  for (int i = 0; i < nThreads; i++){
    pthread_t pt;
    pthread_create(&pt, 0 ,Transform2IDTHread, (void*)i);
  }

  // Wait for all threads complete
  pthread_cond_wait(&exitCond, &exitMutex);
  pthread_mutex_unlock(&exitMutex);
  
  my_finish = GetMilliseconds();
  cout << "Elasped time: " << my_finish - my_start << "ms" << endl;
  
  // Write the 2D Inverse DFT data
  cout << "......writing inverse DFT results to file......." << endl;
  image.SaveImageDataReal("MyAfterInverse.txt", h, N, height);
  cout << "Done\n";
  
  //pthread_exit(NULL);
  
}