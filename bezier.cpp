#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <locale>
#include <stdio.h>
#include <stdlib.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include <Eigen/Dense>
#include <time.h>
#include <math.h>


#define PI 3.14159265  // Should be used from mathlib
inline float sqr(float x) { return x*x; }
using namespace Eigen;
using namespace std;

//****************************************************
// Some Classes
//****************************************************

class Point;
class Viewport;
class BezPatch;

//***************** POINT *****************//
class Point {
  public:
    Vector4f point;
    Point();
    Point(float, float, float);
    Point(Vector4f);
    //Point add(Vector);
    //Point sub(Vector);
    //Vector sub(Point); //Uses current point as the arrow side of vector
    //Point transform(Transformation); //Returns the transformed point
};

//***************** VIEWPORT *****************//
class Viewport {
  public:
    int w, h; // width and height
    Viewport();
};

//***************** BEZPATCH *****************//
class BezPatch {
  public:
    Point controlPoints[4][4];
    BezPatch(Point a[4][4]);

};

//***************** POINT METHODS *****************//
Point::Point() {
  Vector4f temp(0, 0, 0, 1);
  point = temp;
}

Point::Point(float a, float b, float c) {
  Vector4f temp(a, b, c, 1);
  point = temp;
}

Point::Point(Vector4f vec) {
  point = Vector4f(vec(0), vec(1), vec(2), 1);
}

/*Point Point::add(Vector v) {
  Vector4f temp = point + v.vector;
  return Point(temp);
  }

  Point Point::sub(Vector v) {
  Vector4f temp = point - v.vector;
  return Point(temp);
  }

  Vector Point::sub(Point p) {
  Vector4f temp = point - p.point;
  return Vector(temp);
  }

  Point Point::transform(Transformation trans) {
  Point temp;
  temp = Point(trans.matrix * point);
  return temp;
  }*/

//***************** BEZPATCH METHODS *****************//
BezPatch::BezPatch(Point cps[4][4]) {
  controlPoints = cps;
}

//***************** BEZPATCH METHODS *****************//
Viewport::Viewport() {
    w, h = 400;
}
//****************************************************
// Global Variables
//****************************************************
Viewport viewport = Viewport();
float parameter;
int patches;
bool adaptive = false;
string filename;

//****************************************************
// Simple init function
//****************************************************
void initScene(){

  // Nothing to do here for this simple example.

}


//****************************************************
// reshape viewport if the window is resized
//****************************************************
void myReshape(int w, int h) {
  viewport.w = w;
  viewport.h = h;

  glViewport (0,0,viewport.w,viewport.h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, viewport.w, 0, viewport.h);

}

//****************************************************
// A routine to set a pixel by drawing a GL point.  This is not a
// general purpose routine as it assumes a lot of stuff specific to
// this example.
//****************************************************

void setPixel(int x, int y, GLfloat r, GLfloat g, GLfloat b) {
  glColor3f(r, g, b);
  glVertex2f(x + 0.5, y + 0.5);   // The 0.5 is to target pixel
}

//****************************************************
// function that does the actual drawing of stuff
//***************************************************
void myDisplay() {
  glClear(GL_COLOR_BUFFER_BIT);				// clear the color buffer

  glMatrixMode(GL_MODELVIEW);			        // indicate we are specifying camera transformations
  glLoadIdentity();				        // make sure transformation is "zero'd"

  // Start drawing
  /*if(isTor) {
    torus(viewport.w / 2.0 , viewport.h / 2.0 , innerRad, outerRad);
    }
    else {
    circle(viewport.w / 2.0 , viewport.h / 2.0 , min(viewport.w, viewport.h) / 3.0);
    }*/
  glFlush();
  glutSwapBuffers();					// swap buffers (we earlier set double buffer)
}

//****************************************************
// Parse Scene File
//****************************************************

void loadScene(std::string file) {
  cout << "loading Scene .. "<< endl;
  ifstream inpfile(file.c_str());
  if(!inpfile.is_open()) {
    std::cout << "Unable to open file" << std::endl;
  } else {
    std::string line;
    getline(inpfile,line);
    printf("num patches: %d", atoi(line.c_str()));
    patches = atoi(line.c_str());

    int i = 0;
    while(i < patches) {
      vector<string> splitline;
      string buf;

      getline(inpfile,line);
      stringstream ss(line);

      while (ss >> buf) {
        splitline.push_back(buf);
      }
      //Ignore blank lines
      if(splitline.size() == 0) {
        continue;
      }

      //size width height
      else if(!splitline[0].compare("size")) {
        //width = atoi(splitline[1].c_str());
        //height = atoi(splitline[2].c_str());
        //printf("Outputting image of size: %d x %d\n", width, height);
      }
      //output filename
      else if(!splitline[0].compare("output")) {
        filename = splitline[1];
        printf("Writing to file: %s\n", filename.c_str());
      } else {
        std::cerr << "Unknown command: " << splitline[0] << std::endl;
      }
    }
  }

}

//****************************************************
// MAIN
//****************************************************
int main(int argc, char *argv[]) {
  loadScene(argv[1]);
  if (argc < 3) {
    cout << "Not enough arguments" << endl;
    exit(EXIT_FAILURE);
  } else {
    parameter = atoi(argv[2]);
  }

  if (argc == 3 ) {
    if (strcmp(argv[3], "-a") == 0) {
      adaptive = true;
    } else {
      cout << "Command line argument not found" << endl;
      exit(EXIT_FAILURE);
    }
  }

  //This initializes glut
  glutInit(&argc, argv);

  //This tells glut to use a double-buffered window with red, green, and blue channels 
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

  // Initalize theviewport size
  viewport.w = 400;
  viewport.h = 400;

  //The size and position of the window
  glutInitWindowSize(viewport.w, viewport.h);
  glutInitWindowPosition(0,0);
  glutCreateWindow(argv[0]);

  initScene();							// quick function to set up scene

  glutDisplayFunc(myDisplay);        // function to run when its time to draw something
  glutReshapeFunc(myReshape);        // function to run when the window gets resized

  glutMainLoop();							// infinite loop that will keep drawing and resizing
  // and whatever else

  return 0;
}






