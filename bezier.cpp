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
#include <Eigen/StdVector>


#define PI 3.14159265  // Should be used from mathlib
#define EPSILON .0001f
inline float sqr(float x) { return x*x; }
using namespace Eigen;
using namespace std;


//****************************************************
// Some Classes
//****************************************************

class Point;
class LineSeg;
class Viewport;
class BezPatch;
class BezCurve;
class Vector;
typedef std::vector<Point, Eigen::aligned_allocator<Point>> point_vector;

//***************** POINT *****************//
class Point {
public:
	Vector4f point;
	Point();
	Point(float, float, float);
	Point(Vector4f&);
	Point(Vector3f&);
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		//Point add(Vector);
		//Point sub(Vector);
		Vector sub(Point&); //Uses current point as the arrow side of vector
	//Point transform(Transformation); //Returns the transformed point
};

//***************** VECTOR *****************//
class Vector {
public:
	Vector4f vector;
	float len;
	Vector();
	Vector(float, float, float);
	Vector(Point&, Point&); //Vector(start point, end point)
	Vector(Vector4f&);
	Vector add(Vector&);
	Vector sub(Vector&);
	Vector mult(float);
	Vector div(float);
	float dot(Vector&);
	Vector cross(Vector&);
	void normalize();
	bool equals(Vector&);
	//Vector transform(Transformation); //Returns the transformed vector
};

//***************** LINESEG *****************//
class LineSeg {
public:
	Point start, end;
	LineSeg();
	LineSeg(Point&, Point&);
	Point interpolate(float);

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
	point_vector controlPointsPatch; //Yes we're gonna man mode with a 1D vector
	BezPatch(point_vector);
	Point interpolate(float, float, Vector*);
};

//***************** BEZCURVE *****************//
class BezCurve {
public:
	point_vector controlPointsCurve;  //POINTS MUST BE IN ABCD ORDER
	BezCurve(point_vector);
	BezCurve(Point&, Point&, Point&, Point&);
	Point interpolate(float u);
	Vector derivative(float);
};

//***************** POINT METHODS *****************//
Point::Point() : point(0, 0, 0, 1) {}

Point::Point(float a, float b, float c) {
	Vector4f temp(a, b, c, 1);
	point = temp;
}

Point::Point(Vector4f& vec) {
	point = Vector4f(vec(0), vec(1), vec(2), 1);
}

Point::Point(Vector3f& vec) {
	point = Vector4f(vec(0), vec(1), vec(2), 1);
}


/*Point Point::add(Vector v) {
Vector4f temp = point + v.vector;
return Point(temp);
}

Point Point::sub(Vector v) {
Vector4f temp = point - v.vector;
return Point(temp);
}*/

Vector Point::sub(Point& p) {
	Vector4f temp = point - p.point;
	return Vector(temp);
}

/*Point Point::transform(Transformation trans) {
Point temp;
temp = Point(trans.matrix * point);
return temp;
}*/

//***************** VECTOR METHODS *****************//

Vector::Vector() {
	Vector4f temp(0, 0, 0, 0);
	vector = temp;
	len = vector.norm();
}

Vector::Vector(float a, float b, float c) {
	Vector4f temp(a, b, c, 0);
	vector = temp;
	len = vector.norm();
}

Vector::Vector(Vector4f& vec) {
	vector = Vector4f(vec(0), vec(1), vec(2), 0);
	len = vector.norm();
}

Vector::Vector(Point& start, Point& end) {
	vector = end.point - start.point;
	vector = Vector4f(vector(0), vector(1), vector(2), 0);
	len = vector.norm();
}

Vector Vector::add(Vector& v) {
	Vector4f temp = vector + v.vector;
	return Vector(temp);
}

Vector Vector::sub(Vector& v) {
	Vector4f temp = vector - v.vector;
	return Vector(temp);
}

Vector Vector::mult(float k) {
	Vector4f temp = vector * k;
	return Vector(temp);
}

Vector Vector::div(float k) {
	Vector4f temp = vector / k;
	return Vector(temp);
}

float Vector::dot(Vector& v) {
	return vector.dot(v.vector);
}

Vector Vector::cross(Vector& v) {
	Vector3f temp1, temp2, temp3;
	temp1 << vector(0), vector(1), vector(2);
	temp2 << v.vector(0), v.vector(1), v.vector(2);
	temp3 = temp1.cross(temp2);
	Vector4f temp4;
	temp4 << temp3(0), temp3(1), temp3(2), 0;
	return Vector(temp4);
}

void Vector::normalize() {
	vector.normalize();
	len = 1.0f;
}

bool Vector::equals(Vector& v) {
	Vector4f temp = v.vector - vector;
	float size = temp.norm();
	return size == 0;
}

/*Vector Vector::transform(Transformation trans){
Vector temp;
temp = Vector(trans.matrix * vector);
return temp;
}*/

//***************** LINESEG METHODS *****************//
LineSeg::LineSeg() {
	start = Point();
	end = Point();
}

LineSeg::LineSeg(Point& begin, Point& finish) {
	start = begin;
	end = finish;
}

Point LineSeg::interpolate(float u) {
	Vector4f newPoint = start.point * u + end.point * (1 - u);
	return Point(newPoint);
}




//***************** BEZPATCH METHODS *****************//
BezPatch::BezPatch(point_vector cps) {
	controlPointsPatch = cps;
}

Point BezPatch::interpolate(float u, float v, Vector* norm) {
	//Direction confusion here
	point_vector vcurve, ucurve;
	Vector dPdv, dPdu;
	for(int x = 0; x < 4; x++) {
		vcurve[x] = BezCurve(controlPointsPatch[x], controlPointsPatch[4 + x], controlPointsPatch[8 + x], controlPointsPatch[12 + x]).interpolate(u);
		ucurve[x] = BezCurve(controlPointsPatch[4*x], controlPointsPatch[4*x + 1], controlPointsPatch[4*x + 2], controlPointsPatch[4*x + 3]).interpolate(v);
	}
	dPdv = BezCurve(vcurve).derivative(v);
	dPdu = BezCurve(ucurve).derivative(u);
	//DIRECTIONS WHATS GOING ON
	*norm = dPdu.cross(dPdv);
	(*norm).normalize();
	return BezCurve(vcurve).interpolate(v);
}

//***************** BEZCURVE METHODS *****************//
BezCurve::BezCurve(point_vector cps) {
	controlPointsCurve = cps;
}

BezCurve::BezCurve(Point& a, Point& b, Point& c, Point& d) : controlPointsCurve(4) {
	//cout << "asdf"<< endl;
	controlPointsCurve[0] = a;
	controlPointsCurve[1] = b;
	controlPointsCurve[2] = c;
	controlPointsCurve[3] = d;
}

Point BezCurve::interpolate(float u) {
	LineSeg AB, BC, CD;
	AB = LineSeg(controlPointsCurve.at(0), controlPointsCurve.at(1));
	BC = LineSeg(controlPointsCurve.at(1), controlPointsCurve.at(2));
	CD = LineSeg(controlPointsCurve.at(2), controlPointsCurve.at(3));
	//return Point();
	return LineSeg(LineSeg(AB.interpolate(u), BC.interpolate(u)).interpolate(u), LineSeg(BC.interpolate(u), CD.interpolate(u)).interpolate(u)).interpolate(u);
}

Vector BezCurve::derivative(float u) {
	LineSeg AB, BC, CD;
	AB = LineSeg(controlPointsCurve.at(0), controlPointsCurve.at(1));
	BC = LineSeg(controlPointsCurve.at(1), controlPointsCurve.at(2));
	CD = LineSeg(controlPointsCurve.at(2), controlPointsCurve.at(3));
	Point D = LineSeg(AB.interpolate(u), BC.interpolate(u)).interpolate(u);
	Point E = LineSeg(BC.interpolate(u), CD.interpolate(u)).interpolate(u);
	Vector out = E.sub(D);
	out = out.mult(3);
	return out;
}

//***************** SUBDIVIDEPATCH *****************//
void subdividePatch(BezPatch patch, float step) {
	float numdiv = (1 + EPSILON) / step;
	//Come confusion with the for loops here
	for(int iu = 0; iu < numdiv; iu++) {
		float u = iu * step;
		for(int iv = 0; iv < numdiv; iv++) {
			float v = iv * step;
			Vector normal = Vector();
			Point interpPoint = patch.interpolate(u, v, &normal);
			//SAVE INTERPPOINT AND NORMAL HERE
		}
	}
}

//***************** VIEWPORT METHODS *****************//
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

	/*loadScene(argv[1]);
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
	// and whatever else*/
	/*point_vector asdf(4);
	asdf[0] = Point(0, 0, 0);
	asdf[1] = Point(0, 1, 0);
	asdf[2] = Point(1, 1, 0);
	asdf[3] = Point(1, 0, 0);
	BezCurve temp = BezCurve(asdf);

	Point interpPoint = temp.interpolate(.75);

	printf("Interpolated point: %f, %f, %f\n", interpPoint.point(0), interpPoint.point(1), interpPoint.point(2));*/
	return 0;
}






