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
class Vertex;

typedef std::vector<Point, Eigen::aligned_allocator<Point> > point_vector; // DONT CHANCE THIS SPACING OTHERWISE IT DOESN'T COMPILE ON AMY'S COMPUTER

//***************** POINT *****************//
/* Class for storing 3D points using homogenous coordinates */
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
/* Class for storing 3D vectors using homogenous coordinates. Vectors are assumed
 * to start from the origin*/
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
	//bool equals(Vector&);
	//Vector transform(Transformation); //Returns the transformed vector
};

class Vertex {
public:
	Vector3f vert;
	Vertex();
	Vertex(Point&);
};

//***************** LINESEG *****************//
/* Class to store a line segment. Represented by a beginning and end Poing */
class LineSeg {
public:
	Point start, end;
	LineSeg();
	LineSeg(Point&, Point&);
	void interpolate(float, Point*);

};

//***************** VIEWPORT *****************//
/* Class to represent the viewing window */
class Viewport {
public:
	int w, h; // width and height
	Viewport();
};

//***************** BEZPATCH *****************//
/* Class to represent a cubic Bezier Patch with 16 control points */
class BezPatch {
public:
	point_vector controlPointsPatch;
	BezPatch(point_vector);
	void interpolate(float, float, Vector*, Point*);
};

//***************** BEZCURVE *****************//
/* Class to represent a cubic Bezier Curve with 4 control points */
class BezCurve {
public:
	point_vector controlPointsCurve;  //POINTS MUST BE IN ABCD ORDER
	BezCurve(point_vector);
	BezCurve(Point&, Point&, Point&, Point&);
	void interpolate(float u, Point*);
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

Vector Point::sub(Point& p) {
	Vector4f temp = point - p.point;
	return Vector(temp);
}

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

//***************** VERTEX METHODS *****************//
Vertex::Vertex() {
	vert = Vector3f(0, 0, 0);
}
Vertex::Vertex(Point& p) {
	vert = Vector3f(p.point(0), p.point(1), p.point(2));
}

//***************** LINESEG METHODS *****************//
LineSeg::LineSeg() {
	start = Point();
	end = Point();
}

LineSeg::LineSeg(Point& begin, Point& finish) {
	start = begin;
	end = finish;
}

void LineSeg::interpolate(float u, Point* interp) {
	/* Parametric interpolation of a line segment, assuming start is u = 0 end is u = 1*/
	Vector4f newPoint = start.point * (1 - u) + end.point * (u);
	(*interp) = Point(newPoint);
}


//***************** BEZPATCH METHODS *****************//
BezPatch::BezPatch(point_vector cps) {
	controlPointsPatch = cps;
}

void BezPatch::interpolate(float u, float v, Vector* norm, Point* pt) {
	//Direction confusion here
	point_vector vcurve = point_vector(4);
	point_vector ucurve = point_vector(4);
	Vector dPdv, dPdu;
	for(int x = 0; x < 4; x++) {
		Point upoint = Point();
		Point vpoint = Point();
		
		BezCurve(controlPointsPatch[x], controlPointsPatch[4 + x], controlPointsPatch[8 + x], controlPointsPatch[12 + x]).interpolate(u, &upoint);
		BezCurve(controlPointsPatch[4*x], controlPointsPatch[4*x + 1], controlPointsPatch[4*x + 2], controlPointsPatch[4*x + 3]).interpolate(v, &vpoint);
		ucurve[x] = vpoint;
		vcurve[x] = upoint;
	}
	dPdv = BezCurve(vcurve).derivative(v);
	dPdu = BezCurve(ucurve).derivative(u);
	//DIRECTIONS WHATS GOING ON
	(*norm) = dPdu.cross(dPdv);
	(*norm).normalize();
	Point p = Point();
	BezCurve(ucurve).interpolate(u, &p);
	*pt = p;
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

void BezCurve::interpolate(float u, Point* pt) {
	LineSeg AB, BC, CD, EF, FG, HI;
	Point E, F, G, H, I, J;
	AB = LineSeg(controlPointsCurve.at(0), controlPointsCurve.at(1));
	BC = LineSeg(controlPointsCurve.at(1), controlPointsCurve.at(2));
	CD = LineSeg(controlPointsCurve.at(2), controlPointsCurve.at(3));
	AB.interpolate(u, &E); //A
	BC.interpolate(u, &F); //B
	CD.interpolate(u, &G); //C
	EF = LineSeg(E, F); //AB
	FG = LineSeg(F, G); //BC
	EF.interpolate(u, &H); //D
	FG.interpolate(u, &I); //E
	HI = LineSeg(H, I);
	HI.interpolate(u, &J);
	*pt = J;
}

Vector BezCurve::derivative(float u) {
	LineSeg AB, BC, CD, EF, FG, HI;
	Point E, F, G, H, I, J;
	AB = LineSeg(controlPointsCurve.at(0), controlPointsCurve.at(1));
	BC = LineSeg(controlPointsCurve.at(1), controlPointsCurve.at(2));
	CD = LineSeg(controlPointsCurve.at(2), controlPointsCurve.at(3));
	AB.interpolate(u, &E);
	BC.interpolate(u, &F);
	CD.interpolate(u, &G);
	EF = LineSeg(E, F);
	FG = LineSeg(F, G);
	EF.interpolate(u, &H);
	FG.interpolate(u, &I);
	Vector out = H.sub(I);
	out = out.mult(3);
	return out;
}

//***************** SUBDIVIDEPATCH *****************//
void subdividePatch(BezPatch patch, float step, point_vector* VertexArray) {
	int x = 0;
	int numdiv = 1 / step;
	printf("numdiv = %f, step = %f \n", numdiv, step);
	//Come confusion with the for loops here
	for(int iu = 0; iu < numdiv; iu++) {
		float u = iu * step;
		for(int iv = 0; iv < numdiv; iv++) {
			float v = iv * step;
			Vector normal = Vector();
            printf("u = %f, v = %f \n", u, v);
            Point interpPoint0 = Point();
            Point interpPoint1 = Point();
            Point interpPoint2 = Point();
            Point interpPoint3 = Point();

            patch.interpolate(u, v, &normal, &interpPoint0);
            patch.interpolate(u, min(v+step, 1.0f), &normal, &interpPoint1);
            patch.interpolate(min(u+step, 1.0f), min(v+step, 1.0f), &normal, &interpPoint2);
            patch.interpolate(min(u+step, 1.0f), v, &normal, &interpPoint3);

            VertexArray->push_back(interpPoint0);
            VertexArray->push_back(interpPoint1);
            VertexArray->push_back(interpPoint2);
            VertexArray->push_back(interpPoint3);
            x++;
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
bool isFlat = false;
bool isWireframe = false;
float tipangle = 0.0f;
float turnangle = 0.0f;
vector<BezPatch> patchList;

//****************************************************
// DRAW FUNCTION
//****************************************************
//Currently draws a single bezier patch given a step
void drawBezPatch(BezPatch patch, float step) {
	float numdiv = 1 / step;
	//vertexArray size needs to be related to numdiv
	int vertexArraySize =  (int) numdiv * numdiv * 4;
	printf("vertexArraySize = %d\n", vertexArraySize);
	point_vector vertexArray;
	subdividePatch(patch, step, &vertexArray);
	glBegin(GL_QUADS);
	printf("vertexArray.size() = %d\n", vertexArray.size());
    for(int i = 0; i < vertexArray.size(); i +=4) {
    	printf("drawing point (%f, %f, %f) \n", vertexArray[i].point(0), vertexArray[i].point(1), vertexArray[i].point(2));
        printf("drawing point (%f, %f, %f) \n", vertexArray[i+1].point(0), vertexArray[i+1].point(1), vertexArray[i+1].point(2));
        printf("drawing point (%f, %f, %f) \n", vertexArray[i+2].point(0), vertexArray[i+2].point(1), vertexArray[i+2].point(2));
        printf("drawing point (%f, %f, %f) \n", vertexArray[i+3].point(0), vertexArray[i+3].point(1), vertexArray[i+3].point(2));
        glVertex3f(vertexArray[i].point(0), vertexArray[i].point(1), vertexArray[i].point(2));
        glVertex3f(vertexArray[i+1].point(0), vertexArray[i+1].point(1), vertexArray[i+1].point(2));
        glVertex3f(vertexArray[i+2].point(0), vertexArray[i+2].point(1), vertexArray[i+2].point(2));
        glVertex3f(vertexArray[i+3].point(0), vertexArray[i+3].point(1), vertexArray[i+3].point(2));
    }
    glEnd();
}

//****************************************************
// Simple init function
//****************************************************
void initScene(){
	glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
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
	//gluOrtho2D(0, viewport.w, 0, viewport.h);

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
	glRotatef (tipangle, 1,0,0);  // Up and down arrow keys 'tip' view.
    glRotatef (turnangle, 0,0,1);  // Right/left arrow keys 'turn' view.
	for(int i = 0; i < patchList.size(); i++) {
		drawBezPatch(patchList[i], parameter);
	}
	glEnd();                            // Finished Drawing The Triangle

	glFlush();
	glutSwapBuffers();					// swap buffers (we earlier set double buffer)
}

//****************************************************
// Deal with Keyboard input
//***************************************************
void keyboard( unsigned char key, int x, int y )
{
	switch(key) {
		case 's':
			if(isFlat) {
				glShadeModel(GL_SMOOTH);
			} else {
				glShadeModel(GL_FLAT);
			}
			break;
		case 'w':
			if(isWireframe) {
				glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
			} else {
				glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
			}
			break;
	}
	glutPostRedisplay();
}

void arrowkeys( int key, int x, int y )
{
	switch(key) {
		case GLUT_KEY_LEFT :
			printf("left keyboard \n");
			turnangle += 1;
			break;
       	case GLUT_KEY_RIGHT: 
       		printf("right keyboard \n");
       		turnangle -= 1;
       		break;
       	case GLUT_KEY_UP   :
       		printf("up keyboard \n");
       		tipangle += 1;
       		break;  
       	case GLUT_KEY_DOWN :
       		printf("down keyboard \n");  
       		tipangle -= 1;
       		break;
	}
	glutPostRedisplay();
}
//****************************************************
// Parse Scene File
//****************************************************

void loadScene(std::string file) {
	cout << "loading Scene .. "<< endl;
	ifstream inpfile(file.c_str());
    printf("blah\n");
	if(!inpfile.is_open()) {
		cout << "Unable to open file" << endl;
	} else {
		string line;
		getline(inpfile,line);
		printf("num patches: %d\n", atoi(line.c_str()));
		patches = atoi(line.c_str());
		int i = 0;
		while(i < patches) {
            printf("i = %d\n", i);
			point_vector points;
			for(int j = 0; j < 4; j++) {
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
				} else {
					float a1 = atof(splitline[0].c_str());
					float a2 = atof(splitline[1].c_str());
					float a3 = atof(splitline[2].c_str());
					float b1 = atof(splitline[3].c_str());
					float b2 = atof(splitline[4].c_str());
					float b3 = atof(splitline[5].c_str());
					float c1 = atof(splitline[6].c_str());
					float c2 = atof(splitline[7].c_str());
					float c3 = atof(splitline[8].c_str());
					float d1 = atof(splitline[9].c_str());
					float d2 = atof(splitline[10].c_str());
					float d3 = atof(splitline[11].c_str());
					points.push_back(Point(a1, a2, a3));
					points.push_back(Point(b1, b2, b3));
					points.push_back(Point(c1, c2, c3));
					points.push_back(Point(d1, d2, d3));
                	printf("%f, %f, %f,  %f,  %f,  %f,  %f,  %f,  %f,  %f,  %f,  %f\n", a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3);
                }
			}
			patchList.push_back(BezPatch(points));
			printf("added new patch\n");
			i++;

		}
	}

}

//****************************************************
// MAIN
//****************************************************
int main(int argc, char *argv[]) {
	loadScene(argv[1]);
	printf("after loadScene()\n");
	if (argc < 3) {
		cout << "Not enough arguments" << endl;
		exit(EXIT_FAILURE);
	} else {
		printf("argv[2]: %s\n", argv[2]);
		parameter = atof(argv[2]);
		printf("subdivision parameter: %f\n", parameter);
	}
	printf("argc = %d\n", argc);
	if (argc > 3 ) {
		if (strcmp(argv[3], "-a") == 0) {
			adaptive = true;
		} else {
			cout << "Command line argument not found" << endl;
			exit(EXIT_FAILURE);
		}
	}
	printf("after reading in parameters\n");
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
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(arrowkeys);
	

	/*point_vector asdf(16);

	asdf[0] = Point(0, 0, 0);
	asdf[1] = Point(0, .333, 0);
	asdf[2] = Point(0, .667, 0);
	asdf[3] = Point(0, 1.00, 0);
	asdf[4] = Point(.333, 0, 0);
	asdf[5] = Point(.333, .333, 0);
	asdf[6] = Point(.333, .667, 0);
	asdf[7] = Point(.333, 1.00, 0);
	asdf[8] = Point(.667, 0, 0);
	asdf[9] = Point(.667, .333, 0);
	asdf[10] = Point(.667, .667, 0);
	asdf[11] = Point(.667, 1.00, 0);
	asdf[12] = Point(1.00, 0, 0);
	asdf[13] = Point(1.00, .333, 0);
	asdf[14] = Point(1.00, .667, 0);
	asdf[15] = Point(1.00, 1.00, 0);
	
	BezPatch temp = BezPatch(asdf);
	Vector norm = Vector();
	int size = 16;
	point_vector vertexList = point_vector(size);
	Point interpPoint = temp.interpolate(.33, .33, &norm);
	//drawBezPatch(temp, .33);
	//subdividePatch(temp, .33, &vertexList);

	//printf("Interpolated point: %f, %f, %f\n", interpPoint.point(0), interpPoint.point(1), interpPoint.point(2));
	printf("vertexList point: %f, %f, %f\n", vertexList[0].point(0), vertexList[0].point(1), vertexList[0].point(2));
	printf("vertexList point: %f, %f, %f\n", vertexList[1].point(0), vertexList[1].point(1), vertexList[1].point(2));
	printf("vertexList point: %f, %f, %f\n", vertexList[2].point(0), vertexList[2].point(1), vertexList[2].point(2));
	printf("vertexList point: %f, %f, %f\n", vertexList[3].point(0), vertexList[3].point(1), vertexList[3].point(2));
	printf("vertexList point: %f, %f, %f\n", vertexList[4].point(0), vertexList[4].point(1), vertexList[4].point(2));
	printf("vertexList point: %f, %f, %f\n", vertexList[5].point(0), vertexList[5].point(1), vertexList[5].point(2));
	printf("vertexList point: %f, %f, %f\n", vertexList[6].point(0), vertexList[6].point(1), vertexList[6].point(2));
	printf("vertexList point: %f, %f, %f\n", vertexList[7].point(0), vertexList[7].point(1), vertexList[7].point(2));
	printf("vertexList point: %f, %f, %f\n", vertexList[8].point(0), vertexList[8].point(1), vertexList[8].point(2));*/


	glutMainLoop();							// infinite loop that will keep drawing and resizing
	// and whatever else*/

	
	return 0;
}






