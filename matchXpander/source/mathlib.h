//////////////////////////////////////////////////////////////////////
////Header File for Math Library
////
///////////////////////////////


#if !defined(_MATHLIB_H_)
#define _MATHLIB_H_

#define PI 3.1415926535897932384626433832795



double * transformVector(double * matrix, double * vector);
int * sort3(int first, int second, int third);				///Make this better
double * transpose(double * mat);
double * transpose4x4(double * mat);
double * invertMatrix(double * mat);	
double * normalizeVector(double * vector);
double vectorSize(double x, double y, double z);
double * crossProd(double * first, double * second);
double * getCentroid(double * vertices, int size);
double * matrixMult(double * first, double * second);
double determinant(double * matrix);
double * scalarMultiplication(double * vector, double scale);
double dotProd(double * first, double * second);
char * stringCopy(char* str, int start, int end);
char * firstChar(char * str);
bool makeBool(int * i);
double getRand(double bound);
double SSStriangleAngle(double a, double b, double c);
void printMatrix(double * mat, int size);
double * transformVector3x4(double * mat, double * vector);		////transforms a 3-vector with a 4matrix, by padding a 1 for scale.
double * matrixMult4x4(double * first, double * second);		///4x4 matrix times 4x4 matrix multiplication
bool testMatrix(double * matrix);
bool testMatrix3x3(double * matrix);
//void generateSideChainRMSD(int * matchArray, int numRes, Motif * motif, caList * target, double * RMSD, double * transform);
/*----------------Macros----------------*/
#define sq(a) (a)*(a)


#endif


