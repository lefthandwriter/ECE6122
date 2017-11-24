// Calculate and display the Mandelbrot set
// ECE4893/8893 final project
// Name: Esther
// Fall 2016

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

#include "complex.h"

#include <GLUT/glut.h>
#include <OpenGL/glext.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h> // Mac

/*#include <GL/glut.h>
#include <GL/glext.h>
#include <GL/gl.h>
#include <GL/glu.h>*/ // Linux

using namespace std;

//////   ***  Mandelbrot related variables  ***   //////
// Min and max complex plane values
//Complex minC(-2.0, -1.2);
//Complex maxC( 1.0, 1.8);
Complex* minC = new Complex[100];
Complex* maxC = new Complex[100];
//double re_dx = abs(minC.real-maxC.real)/512.0;
//double im_dx = abs(minC.imag-maxC.imag)/512.0;
//double re_dx = (2.0+1.0)/512.0;
//double im_dx = (1.2+1.8)/512.0;
double re_dx[100];
double im_dx[100];

double new_dx = 0;
int    maxIt = 8000;     // Max iterations for the set computations

int width = 512;
int height = 512;

//  Zoom history: Assume we won't / computationally can't zoom more than 100 levels
//int* mandelbrot = new int[width*height]; // store colour information
int mandelbrot[512*512][100]; // store colour information
int zoomLevel = 0;
int startX;
int startY;
int endX;
int endY;
//int PixelRange;
//double ComplexRange = 3.0; // initial value
int PixelRange[100];
double ComplexRange[100]; // in main, initialize element 0 to 3.0
//double PrevComplexRange;   // needed for back button

//////  ***   Pthreads related variables   ***  //////
pthread_mutex_t exitMutex;
pthread_cond_t exitCond;

int numThreads = 16;
int rowsPerThread = height/numThreads;

int threadCtr = 0; // Global counter for barrier
pthread_mutex_t elementCount;
bool* localSense;
bool globalSense;

//////  ***   OpenGL related variables   ***   //////
bool signalDrawSquare = false;


//////  ***   Function Declarations   ***   //////
void ComputeMandelbrot(int startRow);


//////  ***   Barrier Related Functions - Start   ***   //////
void DecrementCount()
{
    pthread_mutex_lock(&elementCount);
    threadCtr--;
    pthread_mutex_unlock(&elementCount);
}
void MyBarrier_Init()
{
    threadCtr = numThreads; // make threadCtr = 16
    //cout << "MyBarrier Initialized. threadCtr set to " << threadCtr << endl;
    pthread_mutex_init(&elementCount, 0);
    localSense = new bool[numThreads];
    for (int i = 0; i < numThreads; i++)
      localSense[i] = true;
    globalSense = true;
}
void MyBarrier(int myId)
{
    localSense[myId] = !localSense[myId]; // make all of them false
    DecrementCount();
    //cout << "Thread " << myId << " called the barrier" << endl;
    if (threadCtr == 0){
        // All threads have called the barrier
        //cout << "All threads hit the barrier." << endl;
        globalSense = localSense[myId];    
        MyBarrier_Init(); // last one resets barrier parameters
    }
    else{
        //cout << "Thread " << myId << "spinning\n";
        while (globalSense != localSense[myId]){} // Spin
    }
}
//////  ***   End - Barrier Related Functions   ***   //////


//////  ***   Thread Related Routines - Start   ***   //////
void* ComputeThread(void* v)
{
	unsigned long myId = (unsigned long)v; // which thread am I
	int startRow = (int)myId * rowsPerThread;

	ComputeMandelbrot(startRow);
	//cout << "Hello from thread ID: " << (int)myId << endl;

	// create a barrier to ensure ComputeMandelbrot finishes before glutPostRedisplay() happens
	MyBarrier((int)myId);

	// At the end, signal main we're done
  	pthread_mutex_lock(&exitMutex);
  	pthread_cond_signal(&exitCond);
  	pthread_mutex_unlock(&exitMutex);
  
	return 0;
}
void CreateThreads(){
  // Initialize all mutex and condition variables
  pthread_mutex_init(&exitMutex, 0);
  pthread_cond_init(&exitCond, 0);
  	
  // Hold the exit mutex while waiting for exitCond condition
  pthread_mutex_lock(&exitMutex);
  
  // Create 16 threads for computing Mandelbrot
  for (int i = 0; i < numThreads; i++){
    pthread_t pt;
    pthread_create(&pt, 0 ,ComputeThread, (void*)i);
  }

  // Wait for all threads complete
  pthread_cond_wait(&exitCond, &exitMutex);
  pthread_mutex_unlock(&exitMutex);
}
//////  ***   End - Thread Related Routines   ***   //////


//*    Function: Each thread calls ComputeMandelbrot to compute a fixed number of c values  *//
void ComputeMandelbrot(int startRow){
	Complex c;
	Complex z;
	Complex tmp;
	for (int y = startRow; y < rowsPerThread+startRow; y++){
		c.imag = minC[zoomLevel].imag + y*im_dx[zoomLevel];
		for (int x = 0; x < width; x++){
			c.real = minC[zoomLevel].real + x*re_dx[zoomLevel];

			// Check if c belongs to Mandelbrot set //
			z = c;
			tmp = z;
			for (int k = 0; k < maxIt; k++){
				if (z.Mag2() > 4.0){ // c is not in the mandelbrot set //
					mandelbrot[y*width + x][zoomLevel] = k;
					break;
				}
				else if(k == 1999 && z.Mag2() <= 4.0){ // c is in the mandelbrot set //
					mandelbrot[y*width + x][zoomLevel] = 0;
				}
				z.real = (tmp.real*tmp.real - tmp.imag*tmp.imag) + c.real;
				z.imag = (tmp.real*tmp.imag*2) + c.imag;
				tmp = z;
			}
		}
	}
}

void drawMandelbrot(){
	// Use mandelbrot array to draw
	glBegin(GL_POINTS);
	for (int y = 0; y < height; y++){
	    for (int x = 0; x < width; x++){
	    	float k = mandelbrot[y*width + x][zoomLevel];
			if (k == 0){} // black, do nothing
			else{ // colour that pixel
				glBegin(GL_POINTS);
				//glColor3f(k*k*sin(k)/maxIt, k*k*k/maxIt, k*k*k*sin(k)/maxIt);
				if (k < (maxIt/8 - 1)){
					glColor3f(k*k*sin(k)/maxIt, (k*k*k*k)/maxIt, k*k*sin(k)/maxIt);
				}
				else if(k > maxIt/8 && k < (maxIt/4 - 1)){
				 	glColor3f(k*k*k/maxIt, k*k*sin(k)/maxIt, k*sin(k)/maxIt);
				}
				else if(k > maxIt/4 && k < (maxIt/2 - 1)){
				 	glColor3f(k*k*k/maxIt, k*k*k/maxIt, k*sin(k)/maxIt);
				}
				else if(k > maxIt/2 && k < (maxIt*3/4 - 1)){
					glColor3f(0.5, k*k/maxIt, k*sin(k)/maxIt);
				}
				else {
					glColor3f(k*k*sin(k)/maxIt, k*k*k*sin(k)/maxIt, k*k/maxIt);
				}
				glVertex2i(x, y);
				glEnd();
			}	
		}
	}
}

void updateDX(){
	// update complex range from mouse click
	PixelRange[zoomLevel] = abs(endX - startX); // same as in drawSquare function
	cout << "PixelRange: " << PixelRange[zoomLevel] << endl;
	int dx = abs(endX - startX); // same as in drawSquare function
	
	ComplexRange[zoomLevel] = PixelRange[zoomLevel]* ComplexRange[zoomLevel-1]/PixelRange[zoomLevel-1];

	//ComplexRange[zoomLevel] = PixelRange * ComplexRange[zoomLevel-1]/512.0;
	cout << "ComplexRange: " << ComplexRange[zoomLevel] << endl;

	//PrevComplexRange = ComplexRange;
	//double newComplexRange = PixelRange * ComplexRange/512.0;
	/*minC.real = startX/512.0 * ComplexRange[zoomLevel];	
	maxC.real = endX/512.0   * ComplexRange[zoomLevel];
	minC.imag = startY/512.0 * ComplexRange[zoomLevel];
	maxC.imag = (startY+dx)/512.0  * ComplexRange[zoomLevel];*/

	// update complex coordinates
	if (startX < endX){
		minC[zoomLevel].real = minC[zoomLevel-1].real + startX * re_dx[zoomLevel-1];	
		maxC[zoomLevel].real = minC[zoomLevel-1].real + endX   * re_dx[zoomLevel-1];
		if (startY < endY){ 		// rightwards-downwards
			minC[zoomLevel].imag = minC[zoomLevel-1].imag + startY       * im_dx[zoomLevel-1];
			maxC[zoomLevel].imag = minC[zoomLevel-1].imag + (startY+dx)  * im_dx[zoomLevel-1];
		}	
		else{						// rightwards-upwards
			maxC[zoomLevel].imag = minC[zoomLevel-1].imag + startY * im_dx[zoomLevel-1];
			minC[zoomLevel].imag = minC[zoomLevel-1].imag + (startY-dx)  * im_dx[zoomLevel-1];
		}
	}
	else{
		maxC[zoomLevel].real = minC[zoomLevel-1].real + startX * re_dx[zoomLevel-1];	
		minC[zoomLevel].real = minC[zoomLevel-1].real + endX   * re_dx[zoomLevel-1];
		if (startY < endY){			// leftwards-downwards
			minC[zoomLevel].imag = minC[zoomLevel-1].imag + startY * im_dx[zoomLevel-1];
			maxC[zoomLevel].imag = minC[zoomLevel-1].imag + (startY+dx)  * im_dx[zoomLevel-1];
		}	
		else{						// leftwards-upwards
			maxC[zoomLevel].imag = minC[zoomLevel-1].imag + startY * im_dx[zoomLevel-1];
			minC[zoomLevel].imag = minC[zoomLevel-1].imag + (startY-dx) * im_dx[zoomLevel-1];
		}
	}
	// Tests that the above is correct
	// (xMaxCnew - xMinCNew)*dx = ComplexRange
	// xStart * 0.00585938 != xMinCNew
	
	// compute new dx - always divide by 512 to get effect of zoom-in
	double tmp = maxC[zoomLevel].real - minC[zoomLevel].real;
	if (tmp < 0){
		tmp = -tmp;
	}
	re_dx[zoomLevel] = tmp/512.0;
	im_dx[zoomLevel] = tmp/512.0;
	//re_dx = ComplexRange[zoomLevel]/512.0;
	//im_dx = ComplexRange[zoomLevel]/512.0;

	cout << "Updated re_dx: " << re_dx[zoomLevel] << endl;
	cout << "Updated im_dx: " << im_dx[zoomLevel] << endl;	

	cout << "New complex coordinates: "; minC[zoomLevel].Print();
	cout << "  to  "; maxC[zoomLevel].Print(); cout << endl;

	// update complex range
	//ComplexRange = newComplexRange;

}

void drawSquare(){
	glBegin(GL_LINE_LOOP);
	glColor3f(1.0, 1.0, 1.0);
	
	int dx = abs(endX - startX);

	if (startX < endX){ // draw square rightwards
		glVertex2i(startX, startY);     // top/bottom left
		glVertex2i(startX+dx, startY);  // top/bottom right

		if (startY < endY){		// rightwards-downwards	// Case 1
			glVertex2i(startX+dx, startY+dx);	// bottom-right
			glVertex2i(startX, startY+dx);      // bottom-left
		}	
		else{	// rightwards-upwards		    // Case 4
			glVertex2i(startX+dx, startY-dx);	// top-right
			glVertex2i(startX,    startY-dx);   // top-left
		}
	}
	else{	// draw square leftwards
		glVertex2i(startX, startY);     // bottom/top right
		glVertex2i(startX-dx, startY);  // bottom/top left

		if (startY < endY){  // leftwards-downwards  // Case 2
			glVertex2i(startX-dx, startY+dx);	// top-right
			glVertex2i(startX, startY+dx);      // top-left
		}	
		else{     // leftwards-upwards          // Case 3
			glVertex2i(startX-dx, startY-dx);	// bottom-right
			glVertex2i(startX, startY-dx);      // bottom-left
		}
	}
	glEnd();
	// draw square
	/*glVertex2i(startX, startY);     // start (top left)
	glVertex2i(endX, startY);       // top-right
	glVertex2i(endX, startY+dx);	// bottom-right
	glVertex2i(startX, startY+dx);    // bottom-left
	glEnd();*/
}

void mouse(int button, int state, int x, int y)
{ // Your mouse click processing here
  // state == 0 means pressed, state != 0 means released
  // Note that the x and y coordinates passed in are in 
  // PIXELS, with y = 0 at the top.
	
	//int static yPrev;

	// if mouse click
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		// toggle signalDrawSquare for use in display()
		signalDrawSquare = true; 

		startX = x;
		startY = y;
		endX = x;
		endY = y;

		//yPrev = y;
		
		cout << "left-click detected at: (x,y) " << startX << "," << startY << endl;
	}
	// if mouse release
	else if (button == GLUT_LEFT_BUTTON && state == GLUT_UP)
	{
		// save pixels
		// limit mouse movement to within the window screen
		if (x > 512){
			endX = 512;
		}
		else if (x < 0){
			endX = 0;
		}
		else{
			endX = x;
		}

		/*if ( x < startX ){
			endX = startX - abs(x-startX); 
			if (x < 0){
				endX = 0;
			}
		}
		else{
			endX = startX + abs(x-startX);
			if (x > 512){
				endX = 512;
			}
		}*/

		//endX = x;
		//endY = y;
		//endY = y + (startX-x);
		/*if (y < startY){ // mouse moving upwards
			endY = startY - abs(x-startX);	
		}
		else{	        // mouse moving downwards
			endY = startY + abs(x-startX);	
		}*/

		if (y < startY){ // mouse moving upwards
			endY = startY - abs(x-startX);	
			if (endY < 0){
				endY = 0;
			}
		}
		else{	        // mouse moving downwards
			endY = startY + abs(x-startX);
			if (endY > 512){
				endY = 512;
			}
		}

		// increment zoomLevel - will affect updateDX() and CreateThreads()
		zoomLevel++;

		// update dx value
		updateDX();

		// start threads to re-compute mandelbrot within new range
		CreateThreads();

		// after computing mandelbrot, display zoom in without square
		signalDrawSquare = false;

		cout << "mouse-release detected at: (x,y) " << endX << "," << endY << endl;

		cout << "displaying Mandelbrot at zoom level: " << zoomLevel << endl;

	}
	glutPostRedisplay();
}

void motion(int x, int y)
{ // Your mouse motion here, x and y coordinates are as above

	// Highlight drawing area (square line)
	/*int static yPrev;
	bool state = true;
	if (state == true){
		// initialize yPrev
		yPrev = y;
		state = false;
	}*/


	//yPrev = endY;

	// limit mouse movement to within the window screen
	if (x > 512){
		endX = 512;
	}
	else if (x < 0){
		endX = 0;
	}
	else{
		endX = x;
	}

	/*if (x < startX){ // mouse moving leftwards
		endX = startX - abs(x-startX); 
		if (x < 0){
				endX = 0;
		}
	}
	else{		// mouse moving rightwards
		endX = startX + abs(x-startX);
		if (x > 512){
				endX = 512;
		}
	}*/

	// even though we draw as a square, we also want to explicitly
	// constrain the mouse motion to a square as well
	if (y < startY){ // mouse moving upwards
		//endY = y - abs(x-startX);
		endY = startY - abs(x-startX);	
		if (endY < 0){
			endY = 0;
		}
	}
	else{	        // mouse moving downwards
		//endY = y + abs(x-startX);	
		endY = startY + abs(x-startX);
		if (endY > 512){
			endY = 512;
		}
	}
	
	//cout << "mouse motion detected at: (x,y) " << x << "," << y << endl;
	//cout << "mouse motion detected" << endl;
	glutPostRedisplay();

}

void reshape(int w, int h) // callback when window size is changed
{ // Your OpenGL window reshape code here
	glClearColor(0.0, 0.0, 0.0, 0.0); // create black background	
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/*GLsizei winX = glutGet(GLUT_WINDOW_WIDTH);
  	GLsizei winY = glutGet(GLUT_WINDOW_HEIGHT);

	glViewport(0,0, winX, winY); // set new dimension of viewable screen
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
  	
	glOrtho(0, winX, winY, 0, -1, 1 ); // glOrtho(left,right,bottom,top,nearVal,farVal);
  	glMatrixMode(GL_PROJECTION);
  	glLoadIdentity();*/

  	glutPostRedisplay(); // re-display the window  	
}

void keyboard(unsigned char c, int x, int y)
{ // Your keyboard processing here

	// if "b" press detected, decrement zoomLevel
	if (c == 'b'){
		if (zoomLevel != 0){
			//minC[zoomLevel].real = 0.0;
			//maxC[zoomLevel].imag = 0.0;
			zoomLevel--;
			cout << "back button pressed, decrementing zoomLevel to: " << zoomLevel << endl;
		}
		else if (zoomLevel == 0){
			// reset this in case of rounding errors
			ComplexRange[0] = 3.0;
  			PixelRange[0]   = 512;
  			minC[0].real    = -2.0;
  			minC[0].imag    = -1.2;
  			maxC[0].real    = 1.0;
  			maxC[0].imag    = 1.8;
			cout << "zoomLevel = 0, unable to decrement zoomLevel further" << endl;	
		}
		// loop back complex range so we can re-zoom from previous point
		//ComplexRange = PrevComplexRange;
	}
	glutPostRedisplay();
}

void display(void)
{ // Your OpenGL display code here
	//glClearColor(1.0, 1.0, 1.0, 0.0); // create white background	
    glClearColor(0.0, 0.0, 0.0, 0.0); // create black background	
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  	
    glShadeModel(GL_SMOOTH);
    glEnable(GL_DEPTH_TEST);

    glViewport(0,0, (GLsizei)width, (GLsizei)height); // set dimension of viewable screen
  	glMatrixMode(GL_MODELVIEW);
  	glLoadIdentity();
  	//gluOrtho2D(0, width, 0, height); // pixel-based view
  	//glOrtho(minC.real, maxC.real, minC.imag, maxC.imag, ((GLfloat)-1), (GLfloat)1); // complex-plane based view
  	//glOrtho(0, width, 0, height, ((GLfloat)-1), (GLfloat)1); // pixel-based view
    glOrtho(0, (GLsizei)width, (GLsizei)height, 0, -1, 1 ); // glOrtho(left,right,bottom,top,nearVal,farVal);
  	glMatrixMode(GL_PROJECTION);
  	glLoadIdentity();
  	
  	//ComputeMandelbrot();
  	glPushMatrix();
  	if (signalDrawSquare == true){
  		drawSquare();
  	}
	drawMandelbrot();	
	glPopMatrix();
  	glutSwapBuffers();
}


int main(int argc, char** argv)
{
  // Initialize OpenGL, but only on the "master" thread or process.
  // See the assignment writeup to determine which is "master" 
  // and which is slave.
  
  // initialize complex numbers
  ComplexRange[0] = 3.0;
  PixelRange[0]   = 512;
  minC[0].real    = -2.0;
  minC[0].imag    = -1.2;
  maxC[0].real    = 1.0;
  maxC[0].imag    = 1.8;
  re_dx[0] = (2.0+1.0)/512.0;
  im_dx[0] = (1.2+1.8)/512.0;

  cout << "re_dx: " << re_dx[0] << endl;
  cout << "im_dx: " << im_dx[0] << endl;

  // initialize barrier
  MyBarrier_Init(); 

  // compute Mandelbrot
  CreateThreads();

  // OpenGL Initialization methods
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  GLsizei winX = (glutGet(GLUT_SCREEN_WIDTH)-width)/2;
  GLsizei winY = (glutGet(GLUT_SCREEN_WIDTH)-height)/2;
  glutInitWindowSize(height, width);
  glutInitWindowPosition(winX, winY);
  glutCreateWindow("Mandelbrot Set");

  // Event handling methods
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);

  glutMouseFunc(mouse);	      // mouse click/release callback
  glutMotionFunc(motion);     // mouse motion callback
  glutKeyboardFunc(keyboard); // keyboard callback

  glutMainLoop();

  return 0;
}


// Mouse doc : https://www.opengl.org/resources/libraries/glut/spec3/node50.html
// Motion doc: https://www.opengl.org/resources/libraries/glut/spec3/node51.html
// Keyboard  : https://www.opengl.org/resources/libraries/glut/spec3/node49.html
