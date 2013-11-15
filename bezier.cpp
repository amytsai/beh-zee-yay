#include <vector>
#include <queue>
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

//Adaptive Subdivision Classes
class Vertex;
class TriangleSide;
class Triangle;

typedef std::vector<Point, Eigen::aligned_allocator<Point> > point_vector; // DONT CHANGE THIS SPACING OTHERWISE IT DOESN'T COMPILE ON AMY'S COMPUTER
typedef std::vector<Vector, Eigen::aligned_allocator<Vector> > normal_vector; 
typedef std::vector<Triangle> triangle_vector;
typedef std::queue<Triangle> triangle_queue;

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
};

//***************** VERTEX *****************//
class Vertex {
public:
	Vector3f worldCoord;
	Vector2f bezierCoord;
	Vertex();
	Vertex(Point&, Vector2f&);
	Vertex(Vector3f&, Vector2f&);
	float u();
	float v();
};

//***************** TRIANGLE SIDE *****************//

class TriangleSide {
public:
	Vertex start, end;
	TriangleSide(Vertex&, Vertex&);
	TriangleSide();
	void midpoint(Vertex*);
};

//***************** TRIANGLE *****************//
class Triangle {
public:
	Vertex a, b, c;
	TriangleSide ab, bc, ca;
	Triangle(Vertex&, Vertex&, Vertex&);
	bool subdivide(BezPatch, float, triangle_vector*);
	void draw(BezPatch);

private:
	bool checkAB(BezPatch, float);
	bool checkBC(BezPatch, float);
	bool checkCA(BezPatch, float);
};

//***************** LINESEG *****************//
/* Class to store a line segment. Represented by a beginning and end Poing */
class LineSeg {
public:
	Point start, end;
	LineSeg();
	LineSeg(Point&, Point&);
	LineSeg(Vertex&, Vertex&);
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
	void interpolate(float, float, Point*);
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
	//printf("NEW POINT : ( %f, %f, %f )\n", vec(0), vec(1), vec(2));
	point = Vector4f(vec(0), vec(1), vec(2), 1);
}

Point::Point(Vector3f& vec) {
	point = Vector4f(vec(0), vec(1), vec(2), 1);
}

Vector Point::sub(Point& p) {
	Vector4f temp = point - p.point;
	return Vector(temp);
}

float dist(Point& a, Point& b) {
	return sqrt(pow(b.point(0) - a.point(0), 2) + pow(b.point(1) - a.point(1), 2) + pow(b.point(2) - a.point(2), 2));
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
	worldCoord = Vector3f(0, 0, 0);
}

Vertex::Vertex(Point& p, Vector2f& uv) {
	worldCoord = Vector3f(p.point(0), p.point(1), p.point(2));
	bezierCoord = uv;
}

Vertex::Vertex(Vector3f& p, Vector2f& uv) {
	worldCoord = p;
	bezierCoord = uv;
}

float Vertex::u() {
	return bezierCoord(0);
}

float Vertex::v() {
	return bezierCoord(1);
}
//***************** TRIANGLE SIDE *****************//
TriangleSide::TriangleSide(Vertex& a, Vertex& b) {
	start = a;
	end = b;
}

TriangleSide::TriangleSide() {
	start = Vertex();
	end = Vertex();
}

void TriangleSide::midpoint(Vertex *mdpt) {
	Point worldmdpt;
	Vector2f beziermdpt;
	Vertex midpoint;
	float midU, midV;
	midU = (start.u() + end.u()) / 2;
	midV = (start.v() + end.v()) / 2;

	beziermdpt = Vector2f(midU, midV);
	LineSeg(start,end).interpolate(0.5, &worldmdpt);
	//printf("Triangle Side Midpoint (%f, %f, %f) \n", worldmdpt.point(0), worldmdpt.point(1), worldmdpt.point(2));

	midpoint = Vertex(worldmdpt, beziermdpt);
	*mdpt = midpoint;
}

//***************** TRIANGLE METHODS *****************//
Triangle::Triangle(Vertex& p1, Vertex& p2, Vertex& p3) {
	a = p1;
	b = p2;
	c = p3;
	ab = TriangleSide(a,b);
	bc = TriangleSide(b,c);
	ca = TriangleSide(c,a);
}

bool Triangle::subdivide(BezPatch patch, float error, triangle_vector* triangles) {
	bool e1, e2, e3;
	e1 = checkAB(patch, error); //e1 = ab
	e2 = checkBC(patch, error); //e2 = bc
	e3 = checkCA(patch, error); //e3 = ca
	//printf("A : (%f, %f, %f) (%f, %f) \n", a.worldCoord(0), a.worldCoord(1), a.worldCoord(2), a.bezierCoord(0), a.bezierCoord(1));

	if(e1 && e2 && e3) {
		return false;
	}
	if(!e1 && e2 && e3) {
		Point beziermdpt;
		Vertex mid = Vertex();
		ab.midpoint(&mid);
		patch.interpolate(mid.u(), mid.v(), &beziermdpt);
		mid = Vertex(beziermdpt, mid.bezierCoord);
		triangles->push_back(Triangle(mid, b, c));
		triangles->push_back(Triangle(mid, c, a));
	}
	else if(e1 && !e2 && e3) {
		Point beziermdpt;
		Vertex mid = Vertex();
		bc.midpoint(&mid);
		patch.interpolate(mid.u(), mid.v(), &beziermdpt);
		mid = Vertex(beziermdpt, mid.bezierCoord);
		triangles->push_back(Triangle(mid, c, a));
		triangles->push_back(Triangle(mid, a, b));
	}
	else if(e1 && e2 && !e3) {
		Point beziermdpt;
		Vertex mid = Vertex();
		ca.midpoint(&mid);
		patch.interpolate(mid.u(), mid.v(), &beziermdpt);
		mid = Vertex(beziermdpt, mid.bezierCoord);
		triangles->push_back(Triangle(mid, a, b));
		triangles->push_back(Triangle(mid, b, c));
	}
	else if(!e1 && !e2 && e3) {
		Point beziermdpt1, beziermdpt2;
		Vertex mid1 = Vertex();
		Vertex mid2 = Vertex();
		ab.midpoint(&mid1);
		bc.midpoint(&mid2);
		patch.interpolate(mid1.u(), mid1.v(), &beziermdpt1);
		mid1 = Vertex(beziermdpt1, mid1.bezierCoord);
		patch.interpolate(mid2.u(), mid2.v(), &beziermdpt2);
		mid2 = Vertex(beziermdpt2, mid2.bezierCoord);
		triangles->push_back(Triangle(mid2, mid1, b));
		triangles->push_back(Triangle(mid2, a, mid1));
		triangles->push_back(Triangle(mid2, c, a));
	}
	else if(e1 && !e2 && !e3) {
		Point beziermdpt1, beziermdpt2;
		Vertex mid1 = Vertex();
		Vertex mid2 = Vertex();
		bc.midpoint(&mid1);
		ca.midpoint(&mid2);
		patch.interpolate(mid1.u(), mid1.v(), &beziermdpt1);
		mid1 = Vertex(beziermdpt1, mid1.bezierCoord);
		patch.interpolate(mid2.u(), mid2.v(), &beziermdpt2);
		mid2 = Vertex(beziermdpt2, mid2.bezierCoord);
		triangles->push_back(Triangle(mid2, mid1, c));
		triangles->push_back(Triangle(mid2, b, mid1));
		triangles->push_back(Triangle(mid2, a, b));
	}
	else if(!e1 && e2 && !e3) {
		Point beziermdpt1, beziermdpt2;
		Vertex mid1 = Vertex();
		Vertex mid2 = Vertex();
		ca.midpoint(&mid1);
		ab.midpoint(&mid2);
		patch.interpolate(mid1.u(), mid1.v(), &beziermdpt1);
		mid1 = Vertex(beziermdpt1, mid1.bezierCoord);
		patch.interpolate(mid2.u(), mid2.v(), &beziermdpt2);
		mid2 = Vertex(beziermdpt2, mid2.bezierCoord);
		triangles->push_back(Triangle(mid2, mid1, a));
		triangles->push_back(Triangle(mid2, c, mid1));
		triangles->push_back(Triangle(mid2, b, c));
	}
	else if(!e1 && !e2 && !e3) {
		Point beziermdpt1, beziermdpt2, beziermdpt3;
		Vertex mid1 = Vertex();
		Vertex mid2 = Vertex();
		Vertex mid3 = Vertex();
		ab.midpoint(&mid1);
		bc.midpoint(&mid2);
		ca.midpoint(&mid3);
		patch.interpolate(mid1.u(), mid1.v(), &beziermdpt1);
		mid1 = Vertex(beziermdpt1, mid1.bezierCoord);
		patch.interpolate(mid2.u(), mid2.v(), &beziermdpt2);
		mid2 = Vertex(beziermdpt2, mid2.bezierCoord);
		patch.interpolate(mid3.u(), mid3.v(), &beziermdpt3);
		mid3 = Vertex(beziermdpt3, mid3.bezierCoord);
		triangles->push_back(Triangle(mid1, b, mid2));
		triangles->push_back(Triangle(mid2, c, mid3));
		triangles->push_back(Triangle(mid3, a, mid1));
		triangles->push_back(Triangle(mid1, mid2, mid3));
	}
	return true;
}

void Triangle::draw(BezPatch patch) {
	Vector normA, normB, normC;
	Point pointA, pointB, pointC;
	patch.interpolate(a.u(), a.v(), &normA, &pointA);
	patch.interpolate(b.u(), b.v(), &normB, &pointB);
	patch.interpolate(c.u(), c.v(), &normC, &pointC);
	glNormal3f(normA.vector(0), normA.vector(1), normA.vector(2));
	glVertex3f(pointA.point(0), pointA.point(1), pointA.point(2));
	glNormal3f(normB.vector(0), normB.vector(1), normB.vector(2));
	glVertex3f(pointB.point(0), pointB.point(1), pointB.point(2));
	glNormal3f(normC.vector(0), normC.vector(1), normC.vector(2));
	glVertex3f(pointC.point(0), pointC.point(1), pointC.point(2));
}

bool Triangle::checkAB(BezPatch patch, float error) {

	Point worldmdpt, beziermdpt;
	Vertex midpoint;
	ab.midpoint(&midpoint);
	worldmdpt = Point(midpoint.worldCoord);
	patch.interpolate(midpoint.u(), midpoint.v(), &beziermdpt);
	
	if(dist(worldmdpt, beziermdpt) < error) {
		return true;
	} else {
		//printf("triangle midpoint \t (%f, %f, %f)\n", midpoint.worldCoord(0), midpoint.worldCoord(1), midpoint.worldCoord(2));
		//printf("surface midpoint \t (%f, %f, %f)\n", beziermdpt.point(0), beziermdpt.point(1), beziermdpt.point(2));
		//printf("current error = %f\n", dist(worldmdpt, beziermdpt));

		return false;
	}	
}

bool Triangle::checkBC(BezPatch patch, float error) {
	Point worldmdpt, beziermdpt;
	Vertex midpoint;
	bc.midpoint(&midpoint);
	worldmdpt = Point(midpoint.worldCoord);
	patch.interpolate(midpoint.u(), midpoint.v(), &beziermdpt);
	

	if(dist(worldmdpt, beziermdpt) < error) {
		return true;
	} else {

		//printf("triangle midpoint \t (%f, %f, %f)\n", midpoint.worldCoord(0), midpoint.worldCoord(1), midpoint.worldCoord(2));
		//printf("surface midpoint \t (%f, %f, %f)\n", beziermdpt.point(0), beziermdpt.point(1), beziermdpt.point(2));
		//printf("current error = %f\n", dist(worldmdpt, beziermdpt));

		return false;
	}	
}

bool Triangle::checkCA(BezPatch patch, float error) {
	Point worldmdpt, beziermdpt;
	Vertex midpoint;
	ca.midpoint(&midpoint);
	worldmdpt = Point(midpoint.worldCoord);
	patch.interpolate(midpoint.u(), midpoint.v(), &beziermdpt);
	

	if(dist(worldmdpt, beziermdpt) < error) {
		return true;
	} else {
		//printf("triangle midpoint \t (%f, %f, %f)\n", midpoint.worldCoord(0), midpoint.worldCoord(1), midpoint.worldCoord(2));
		//printf("surface midpoint \t (%f, %f, %f)\n", beziermdpt.point(0), beziermdpt.point(1), beziermdpt.point(2));
		//printf("current error = %f\n", dist(worldmdpt, beziermdpt));
		return false;
	}	
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

LineSeg::LineSeg(Vertex& begin, Vertex& finish) {
	start = Point(begin.worldCoord(0), begin.worldCoord(1), begin.worldCoord(2));
	end = Point(finish.worldCoord(0), finish.worldCoord(1), finish.worldCoord(2));
}

void LineSeg::interpolate(float u, Point* interp) {
	/* Parametric interpolation of a line segment, assuming start is u = 0 end is u = 1*/
	Vector4f newPoint = start.point * (1.0f - u) + end.point * (u);
	*interp = Point(newPoint);
	//printf("interpolated point (%f, %f, %f) \n", newPoint(0), newPoint(1), newPoint(2));
}


//***************** BEZPATCH METHODS *****************//
BezPatch::BezPatch(point_vector cps) {
	controlPointsPatch = cps;
}

void BezPatch::interpolate(float u, float v, Vector* norm, Point* pt) {
	point_vector vcurve;
	point_vector ucurve;
	Vector dPdv, dPdu;
	for(int x = 0; x < 4; x++) {
		/*if(controlPointsPatch.size() != 16) {
			printf("OH NO controlPointsPatch.size() = %d\n", controlPointsPatch.size());
		}*/
		Point upoint = Point();
		Point vpoint = Point();		
		BezCurve(controlPointsPatch[x], controlPointsPatch[4 + x], controlPointsPatch[8 + x], controlPointsPatch[12 + x]).interpolate(u, &upoint);
		BezCurve(controlPointsPatch[4*x], controlPointsPatch[4*x + 1], controlPointsPatch[4*x + 2], controlPointsPatch[4*x + 3]).interpolate(v, &vpoint);
		ucurve.push_back(vpoint);
		vcurve.push_back(upoint);
	}
	//printf("ucurve.size = %d, vcurve.size = %d\n", ucurve.size(), vcurve.size());
	dPdv = BezCurve(vcurve).derivative(v);
	dPdu = BezCurve(ucurve).derivative(u);
	(*norm) = dPdu.cross(dPdv);
	(*norm).normalize();
	Point p = Point();
	BezCurve(ucurve).interpolate(u, &p);
	(*pt) = p;

}

void BezPatch::interpolate(float u, float v, Point* pt) {
	point_vector vcurve;
	point_vector ucurve;
	Vector dPdv, dPdu;
	for(int x = 0; x < 4; x++) {
		Point upoint = Point();
		Point vpoint = Point();		
		BezCurve(controlPointsPatch[x], controlPointsPatch[4 + x], controlPointsPatch[8 + x], controlPointsPatch[12 + x]).interpolate(u, &upoint);
		BezCurve(controlPointsPatch[4*x], controlPointsPatch[4*x + 1], controlPointsPatch[4*x + 2], controlPointsPatch[4*x + 3]).interpolate(v, &vpoint);
		ucurve.push_back(vpoint);
		vcurve.push_back(upoint);
	}
	Point p = Point();
	BezCurve(ucurve).interpolate(u, &p);
	(*pt) = p;

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
void subdividePatch(BezPatch patch, float step, point_vector* VertexArray, normal_vector* NormalArray) {
	int x = 0;
	int numdiv = 1 / step;
	//printf("numdiv = %d, step = %f \n", numdiv, step);

	for(int iu = 0; iu <= numdiv; iu++) {
		float u = iu * step;
		for(int iv = 0; iv <= numdiv; iv++) {
			float v = iv * step;
			Vector normal0, normal1, normal2, normal3;
            Point interpPoint0, interpPoint1, interpPoint2, interpPoint3;

            patch.interpolate(u, v, &normal0, &interpPoint0);
            patch.interpolate(u, min(v+step, 1.0f), &normal1, &interpPoint1);
            patch.interpolate(min(u+step, 1.0f), min(v+step, 1.0f), &normal2, &interpPoint2);
            patch.interpolate(min(u+step, 1.0f), v, &normal3, &interpPoint3);

            VertexArray->push_back(interpPoint3);
            VertexArray->push_back(interpPoint2);
            VertexArray->push_back(interpPoint1);
            VertexArray->push_back(interpPoint0);

            NormalArray->push_back(normal3);
            NormalArray->push_back(normal2);
            NormalArray->push_back(normal1);
            NormalArray->push_back(normal0);
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
float dx = 0.0f;
float dz = 0.0f;
float scale = 1.0f;
vector<BezPatch> patchList;
float leftB  = 0.0;
float rightB = 0.0;
float bottomB = 0.0;
float topB = 0.0; 
float nearamt = 1;
float faramt = 1;

//****************************************************
// DRAW BEZ PATCH
//****************************************************
//Currently draws a single bezier patch given a step
void drawBezPatch(BezPatch patch, float step) {
	float numdiv = 1 / step;
	point_vector vertexArray;
	normal_vector normalArray;

	if(!adaptive) {
		subdividePatch(patch, step, &vertexArray, &normalArray);
		glBegin(GL_QUADS);

		printf("vertexArray.size() = %d\n", vertexArray.size());
	    for(int i = 0; i < vertexArray.size(); i +=4) {
			//asdfasdf
	    	glNormal3f(normalArray[i].vector(0), normalArray[i].vector(1), normalArray[i].vector(2));
	        glVertex3f(vertexArray[i].point(0), vertexArray[i].point(1), vertexArray[i].point(2));
	        glNormal3f(normalArray[i+1].vector(0), normalArray[i+1].vector(1), normalArray[i+1].vector(2));
	        glVertex3f(vertexArray[i+1].point(0), vertexArray[i+1].point(1), vertexArray[i+1].point(2));
	        glNormal3f(normalArray[i+2].vector(0), normalArray[i+2].vector(1), normalArray[i+2].vector(2));
	        glVertex3f(vertexArray[i+2].point(0), vertexArray[i+2].point(1), vertexArray[i+2].point(2));
	        glNormal3f(normalArray[i+3].vector(0), normalArray[i+3].vector(1), normalArray[i+3].vector(2));
	        glVertex3f(vertexArray[i+3].point(0), vertexArray[i+3].point(1), vertexArray[i+3].point(2));
	                
	    }
	    glEnd();
	} else {
		triangle_queue triangles;
		Point bl = patch.controlPointsPatch[0].point;
		Point br = patch.controlPointsPatch[3].point;
		Point tl = patch.controlPointsPatch[12].point;
		Point tr = patch.controlPointsPatch[15].point;

		Vector2f bottomL = Vector2f(0.0, 0.0);
		Vector2f bottomR = Vector2f(0.0, 1.0);
		Vector2f topL = Vector2f(1.0, 0.0);
		Vector2f topR = Vector2f(1.0, 1.0);

		Vertex bottomleft = Vertex(bl, bottomL);
		Vertex bottomright = Vertex(br, bottomR);
		Vertex topleft = Vertex(tl, topL);
		Vertex topright = Vertex(tr, topR);

		triangles.push(Triangle(bottomleft, bottomright, topright));
		triangles.push(Triangle(bottomleft, topright, topleft));
		glBegin(GL_TRIANGLES);
		while(!triangles.empty()) {
			triangle_vector subtriangles;
			Triangle curtriangle = triangles.front();
			/*printf("curtriangle: A(%f, %f, %f) (%f, %f)\n", curtriangle.a.worldCoord(0), curtriangle.a.worldCoord(1), curtriangle.a.worldCoord(2), curtriangle.a.bezierCoord(0), curtriangle.a.bezierCoord(1));
			printf("curtriangle: B(%f, %f, %f) (%f, %f)\n", curtriangle.b.worldCoord(0), curtriangle.b.worldCoord(1), curtriangle.b.worldCoord(2), curtriangle.b.bezierCoord(0), curtriangle.b.bezierCoord(1));
			printf("curtriangle: C(%f, %f, %f) (%f, %f)\n", curtriangle.c.worldCoord(0), curtriangle.c.worldCoord(1), curtriangle.c.worldCoord(2), curtriangle.c.bezierCoord(0), curtriangle.c.bezierCoord(1));*/
			triangles.pop();
			if(curtriangle.subdivide(patch, parameter, &subtriangles)) {
				while(!subtriangles.empty()) {
					//printf("subtriangles.size() = %d\n", subtriangles.size());
					triangles.push(subtriangles.back());
					subtriangles.pop_back();

			 }

			} else {
				curtriangle.draw(patch);

			}

		}
		glEnd();
	}
}

//****************************************************
// Simple init function
//****************************************************
void initScene(){
    glShadeModel(GL_SMOOTH);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    //glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHT2);

	//glCullFace(GL_BACK);

    GLfloat light_position[] = { 1.0, 1.0, 0.0, 0.0 };
    GLfloat light_position1[] = {0.0, 1.0, 0.0, 0.0};
    GLfloat light_position2[] = {-1.0, 0.0, -1.0, 0.0};
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glLightfv(GL_LIGHT1, GL_POSITION, light_position1);
    glLightfv(GL_LIGHT2, GL_POSITION, light_position2);
    glEnable(GL_NORMALIZE);
	GLfloat kd[] = {.1f, .74f, .61f, 1.0f};
    GLfloat ka[] = {0.1f, 0.74f, 0.61f, 1.0f};
    GLfloat ks[] = {0.1f, 0.1f, 0.1f, 1.0f};
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, kd);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ka);
	glMaterialfv(GL_FRONT, GL_SPECULAR, ks);

	// Nothing to do here for this simple example.

}


//****************************************************
// reshape viewport if the window is resized
//****************************************************
void myReshape(int w, int h) {
	viewport.w = w;
	viewport.h = h;
	glClear(GL_COLOR_BUFFER_BIT);				// clear the color buffer
	glClear(GL_DEPTH_BUFFER_BIT);
	glViewport (0,0,viewport.w,viewport.h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//gluOrtho2D(0, viewport.w, 0, viewport.h);

}

//****************************************************
// function that does the actual drawing of stuff
//***************************************************
void myDisplay() {
	bool test = true;
	glClear(GL_COLOR_BUFFER_BIT);				// clear the color buffer
	glClear(GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);		        // indicate we are specifying camera transformation
	glLoadIdentity();
	float xBound = max(abs(leftB), rightB);
	float yBound = max(abs(bottomB),topB);
	float zBound = max(abs(nearamt), faramt);
	float Bound = max(max(xBound,yBound), zBound);
	glOrtho(-Bound - .5,  Bound + .5, -Bound - .5, Bound + .5, -Bound - .5, Bound + .5);
	//glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
	glTranslatef(dx, dz, 0);

	// make sure transformation is "zero'd"
	glRotatef (tipangle, 1,0,0);  // Up and down arrow keys 'tip' view.
    glRotatef (turnangle, 0,0,1);  // Right/left arrow keys 'turn' view.
    glScalef(scale, scale, scale);
	//glColor3f(1.0, 0, 0);

	for(int i = 0; i < patchList.size(); i++) {
		drawBezPatch(patchList[i], parameter);
		//glColor3f(1.0, 1.0, 0);
		//printf("Finished bez patch \n");
	}

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
				isFlat = false;
				glShadeModel(GL_SMOOTH);
			} else {
				isFlat = true;
				glShadeModel(GL_FLAT);
			}
			break;
		case 'w':
			if(isWireframe) {
				isWireframe = false;
				glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
			} else {
				isWireframe = true;
				glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
			}
			break;
		case '=':
			if(scale >= 10.0) {
				scale = 10.0;
			} else {
				scale += .1;
			}
			break;
		case '-':
			if(scale <= .1f) {
				scale = .1f;
			} else {
				scale -= .1f;
			}
			break;
	}
	glutPostRedisplay();

}

void arrowkeys( int key, int x, int y )
{
	int mod = glutGetModifiers();
	if(mod == 0) {
		switch(key) {
			case GLUT_KEY_LEFT :
				//printf("left keyboard \n");
				turnangle += 2;
				break;
	       	case GLUT_KEY_RIGHT: 
	       		//printf("right keyboard \n");
	       		turnangle -= 2;
	       		break;
	       	case GLUT_KEY_UP   :
	       		//printf("up keyboard \n");
	       		tipangle -= 2;
	       		break;  
	       	case GLUT_KEY_DOWN :
	       		//printf("down keyboard \n");  
	       		tipangle += 2;
	       		break;
		}
	} else if (mod > 0) {
		if(mod == GLUT_ACTIVE_SHIFT) {
			switch(key) {
				case GLUT_KEY_LEFT:
					dx -= .1;
					break;
				case GLUT_KEY_RIGHT:
					dx += .1;
					break;
				case GLUT_KEY_UP:
					dz += .1;
					break;
				case GLUT_KEY_DOWN:
					dz -= .1;
					break;
			}
		} else if(mod == GLUT_ACTIVE_ALT) {
			switch(key) {
				case GLUT_KEY_RIGHT:
					parameter *= 1.5;
					break;
				case GLUT_KEY_LEFT:
					parameter *= .66;
					break;
				case GLUT_KEY_UP:
					parameter *= 2.0;
					break;
				case GLUT_KEY_DOWN:
					parameter *= .5;
					break;
			}
		}
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
			for(int j = 0; j < 5; j++) {
				vector<string> splitline;
				string buf;
				getline(inpfile,line);
				cout << line << endl;
				stringstream ss(line);
				while (ss >> buf) {
					splitline.push_back(buf);
				}
				//Ignore blank lines
				if(splitline.size() == 0) {
					continue;
				} else {
					float a1 = atof(splitline[0].c_str());
					leftB = min(a1, leftB);
					rightB = max(a1, rightB);
					float a2 = atof(splitline[1].c_str());
					bottomB = min(a2, bottomB);
					topB = max(a2, topB);
					float a3 = atof(splitline[2].c_str());
					nearamt = min(a3, nearamt);
					faramt = max(a3, faramt);

					float b1 = atof(splitline[3].c_str());
					leftB = min(b1, leftB);
					rightB = max(b1, rightB);
					float b2 = atof(splitline[4].c_str());
					bottomB = min(b2, bottomB);
					topB = max(b2, topB);
					float b3 = atof(splitline[5].c_str());
					nearamt = min(b3, nearamt);
					faramt = max(b3, faramt);

					float c1 = atof(splitline[6].c_str());
					leftB = min(c1, leftB);
					rightB = max(c1, rightB);
					float c2 = atof(splitline[7].c_str());
					bottomB = min(c2, bottomB);
					topB = max(c2, topB);
					float c3 = atof(splitline[8].c_str());
					nearamt = min(c3, nearamt);
					faramt = max(c3, faramt);

					float d1 = atof(splitline[9].c_str());
					leftB = min(d1, leftB);
					rightB = max(d1, rightB);
					float d2 = atof(splitline[10].c_str());
					bottomB = min(c2, bottomB);
					topB = max(c2, topB);
					float d3 = atof(splitline[11].c_str());
					nearamt = min(d3, nearamt);
					faramt = max(d3, faramt);

					points.push_back(Point(a1, a2, a3));
					points.push_back(Point(b1, b2, b3));
					points.push_back(Point(c1, c2, c3));
					points.push_back(Point(d1, d2, d3));
                	printf("%f, %f, %f,  %f,  %f,  %f,  %f,  %f,  %f,  %f,  %f,  %f\n", a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3);
                	printf("left: %f, right: %f, top: %f, bottom: %f near: %f, far: %f \n", leftB, rightB, topB, bottomB, nearamt, faramt);
                }
			}
			if (points.size() != 16) {
				printf("OH NO points.size() = %d \n", points.size());
			}
			BezPatch newP = BezPatch(points);
			patchList.push_back(newP);
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
			printf("ADAPTIVE IS TRUE\n");
		} else {
			cout << "Command line argument not found" << endl;
			exit(EXIT_FAILURE);
		}
	}
	printf("after reading in parameters\n");
	//This initializes glut
	glutInit(&argc, argv);

	//This tells glut to use a double-buffered window with red, green, and blue channels 
	glutInitDisplayMode(GLUT_DEPTH | GLUT_RGB);

	// Initalize theviewport size
	viewport.w = 600;
	viewport.h = 600;

	//The size and position of the window
	glutInitWindowSize(viewport.w, viewport.h);
	glutInitWindowPosition(0,0);
	glutCreateWindow(argv[0]);
	initScene();							// quick function to set up scene

	glutDisplayFunc(myDisplay);        // function to run when its time to draw something
	glutReshapeFunc(myReshape);        // function to run when the window gets resized
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(arrowkeys);
	
	glutMainLoop();							// infinite loop that will keep drawing and resizing
	// and whatever else*/

	
	return 0;
}






