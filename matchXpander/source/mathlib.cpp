//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//Stand in Math Library by Brian Chen/////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
/*********************************************************

	Contents:

	3x3 matrix Inversion
	3x3 matrix X 3d vector transformation
	3 element sort
	3 dimensional Cross Product
	3x3 matrix determinant
	3x3 times 3x3 matrix multiplication
	3d vector normalization
	3d vector magnitude
	3 3d vertex centroid

	3 dimensional scalar multiplication
	3 dimensional dot product

*/

#include "StdAfx.h"


double * transformVector(double * matrix, double * vector)
{
	//Matrix Arrangement			Vector Arrangement
	//[0, 3. 6]						[0]
	//[1, 4. 7]				X		[1]
	//[2, 5. 8]						[2]

	double * result = new double[3];

	////apply matrix multiplication to the vector
		result[0] = (matrix[0]*vector[0]) + (matrix[3]*vector[1]) + (matrix[6]*vector[2]);
		result[1] = (matrix[1]*vector[0]) + (matrix[4]*vector[1]) + (matrix[7]*vector[2]);
		result[2] = (matrix[2]*vector[0]) + (matrix[5]*vector[1]) + (matrix[8]*vector[2]);

	return result;

}



////////////////////////////
///3 element sorter
////FIXME: make me more efficient.
int * sort3(int first, int second, int third)
{

	int * result = new int[3];
	int state = 0;

	/////find the first element
		if		( ( first >= second) && (first >= third) )	{result[0] = first; state = 0;}
		else if ( ( second >= first) && (second >= third) )	{result[0] = second; state = 1;}
		else if ( ( third >= first) && (third >= second) )	{result[0] = third; state = 2;}


	/////find the second element, and set the third element
		switch (state){
		case 0 :	if (second <= third)	{result[1] = third; result[2] = second; break;}
					else					{result[1] = second; result[2] = third; break;}
		
		case 1 :	if (first <= third)		{result[1] = third; result[2] = first; break;}
					else					{result[1] = first; result[2] = third; break;}

		case 2 :	if (second <= first)	{result[1] = first; result[2] = second; break;}
					else					{result[1] = second; result[2] = first; break;}
		}


	return result;

}

////simple determinant for 3x3 matrices
double determinant(double * matrix){

	double det = 0;

	////calculate determinant
		det = matrix[0]*( (matrix[4]*matrix[8]) - (matrix[7]*matrix[5]) );
		det = det - matrix[3]*( (matrix[1]*matrix[8]) - (matrix[7]*matrix[2]) );
		det = det + matrix[6]*( (matrix[1]*matrix[5]) - (matrix[4]*matrix[2]) );


	return det;

}

double * transpose(double * mat)
{
	//Matrix Arrangement
	//[0, 3. 6]
	//[1, 4. 7]
	//[2, 5. 8]

	double * result = new double[9];

	result[0] = mat[0];
	result[1] = mat[3];
	result[2] = mat[6];
	result[3] = mat[1];
	result[4] = mat[4];
	result[5] = mat[7];
	result[6] = mat[2];
	result[7] = mat[5];
	result[8] = mat[8];

	return result;

}

double * transpose4x4(double * mat)
{
	//Matrix Arrangement
	//[0, 4. 8  12]
	//[1, 5. 9  13]
	//[2, 6. 10 14]
	//[3, 7, 11,15] 

	double * result = new double[16];

	result[0]  = mat[0];
	result[1]  = mat[4];
	result[2]  = mat[8];
	result[3]  = mat[12];
	result[4]  = mat[1];
	result[5]  = mat[5];
	result[6]  = mat[9];
	result[7]  = mat[13];
	result[8]  = mat[2];
	result[9]  = mat[6];
	result[10] = mat[10];
	result[11] = mat[14];
	result[12] = mat[3];
	result[13] = mat[7];
	result[14] = mat[11];
	result[15] = mat[15];

	return result;

}

double * invertMatrix(double * mat)
{

	//Matrix Arrangement
	//[0, 3. 6]
	//[1, 4. 7]
	//[2, 5. 8]

	double * matrix = new double[9];
	double det;
	double * adjoint = new double[9];
	double * inverse = new double[9];

	////copy out the values
		matrix[0] = mat[0];	
		matrix[1] = mat[1];	
		matrix[2] = mat[2];	
		matrix[3] = mat[3];	
		matrix[4] = mat[4];	
		matrix[5] = mat[5];	
		matrix[6] = mat[6];	
		matrix[7] = mat[7];	
		matrix[8] = mat[8];	

	///calculate determinant
		det = determinant(matrix);

	///calculate adjoint - it's the Transpose of the matrix of component determinants.
		adjoint[0] =		( (matrix[4]*matrix[8]) - (matrix[7]*matrix[5]) );
		adjoint[3] =	-	( (matrix[3]*matrix[8]) - (matrix[6]*matrix[5]) );
		adjoint[6] =		( (matrix[3]*matrix[7]) - (matrix[6]*matrix[4]) );
		adjoint[1] =	-	( (matrix[1]*matrix[8]) - (matrix[7]*matrix[2]) );
		adjoint[4] =		( (matrix[0]*matrix[8]) - (matrix[6]*matrix[2]) );
		adjoint[7] =	-	( (matrix[0]*matrix[7]) - (matrix[6]*matrix[1]) );
		adjoint[2] =		( (matrix[1]*matrix[5]) - (matrix[4]*matrix[2]) );
		adjoint[5] =	-	( (matrix[0]*matrix[5]) - (matrix[3]*matrix[2]) );
		adjoint[8] =		( (matrix[0]*matrix[4]) - (matrix[3]*matrix[1]) );

	///calculate inverse
		inverse[0] = adjoint[0]/det;
		inverse[1] = adjoint[1]/det;
		inverse[2] = adjoint[2]/det;
		inverse[3] = adjoint[3]/det;
		inverse[4] = adjoint[4]/det;
		inverse[5] = adjoint[5]/det;
		inverse[6] = adjoint[6]/det;
		inverse[7] = adjoint[7]/det;
		inverse[8] = adjoint[8]/det;

	delete(matrix);
	delete(adjoint);

	return inverse;
}



////////////////////
//Returns normalized vector in input direction
double * normalizeVector(double * vector)
{
	double * result = new double[3];

	////calculate the magnitude
		double mag = sqrt(	(vector[0]*vector[0]) + 
							(vector[1]*vector[1]) + 
							(vector[2]*vector[2]) );

	////set the normalized vector
		result[0] = vector[0]/mag;
		result[1] = vector[1]/mag;
		result[2] = vector[2]/mag;

	return result;

}


///finds the matrix magnitude
double vectorSize(double x, double y, double z)
{
	return sqrt( (x*x) + (y*y) + (z*z) );
}



double * crossProd(double * first, double * second)
{
	///Fix size of array to 3 doubles
	double * result = new double[3];	

	/////calculate cross product
		result[0] = (first[1] * second[2]) - (first[2] * second[1]);
		result[1] = -( (first[0] * second[2]) - (first[2] * second[0]) );
		result[2] = (first[0] * second[1]) - (first[1] * second[0]);

	return result;

}

double * getCentroid(double * vertices, int size)
{
	double * result = new double[3];
	double tempSum = 0.0;
	int i = 0;
	
	for (i = 0; i<size; i++){tempSum = tempSum + vertices[3*i+0];}
	result[0] = tempSum/size;

	tempSum = 0.0;
	for (i = 0; i<size; i++){tempSum = tempSum + vertices[3*i+1];}
	result[1] = tempSum/size;

	tempSum = 0.0;
	for (i = 0; i<size; i++){tempSum = tempSum + vertices[3*i+2];}
	result[2] = tempSum/size;

	return result;
}



///////////////////////////
///3x3 matrix times 3x3 matrix multiplication
double * matrixMult(double * first, double * second)
{
	double * result = new double[9];
	
		result[0] = (first[0]*second[0]) + (first[3]*second[1]) + (first[6]*second[2]);
		result[1] = (first[1]*second[0]) + (first[4]*second[1]) + (first[7]*second[2]);
		result[2] = (first[2]*second[0]) + (first[5]*second[1]) + (first[8]*second[2]);
		result[3] = (first[0]*second[3]) + (first[3]*second[4]) + (first[6]*second[5]);
		result[4] = (first[1]*second[3]) + (first[4]*second[4]) + (first[7]*second[5]);
		result[5] = (first[2]*second[3]) + (first[5]*second[4]) + (first[8]*second[5]);
		result[6] = (first[0]*second[6]) + (first[3]*second[7]) + (first[6]*second[8]);
		result[7] = (first[1]*second[6]) + (first[4]*second[7]) + (first[7]*second[8]);
		result[8] = (first[2]*second[6]) + (first[5]*second[7]) + (first[8]*second[8]);

	return result;
}


double * scalarMultiplication(double * vector, double scale)
{
	double * result = new double[3];
	
	/////calculate the scaled vector
		result[0] = vector[0]*scale;
		result[1] = vector[1]*scale;
		result[2] = vector[2]*scale;

	return result;

}



double dotProd(double * first, double * second)
{
	double result = 0;

	/////calculate the dot product
		result =	(first[0]*second[0]) +
					(first[1]*second[1]) +
					(first[2]*second[2]);

	return result;

}


/////////////////////////////////////////////////////////////////////////////////////////////////
/////This function is for getting values out of the file when the position is fixed//////////////
/////Unlike strcpy which copies everything; I want to just copy specifics////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
char * stringCopy(char* str, int start, int end)
{
	int size = end - start + 1;
	char * result = new char[size + 1];

	for(int i = 0; i<size; i++){
		result[i] = str[start+i];	
	}
	result[size] = '\0';

	return result;

}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////This function is for getting the first character out of a string
/////////////////////////////////////////////////////////////////////////////////////////////////
char * firstChar(char * str)
{
	char * result = new char[2];
	if (str != NULL){
		result[0] = str[0];
		result[1] = '\0';
	}
	return result;

}


/////////////////////////////////////////////////////////////////////////////////////////////////
/////This function is for getting a bool fast
/////////////////////////////////////////////////////////////////////////////////////////////////
bool makeBool(int * i)
{
	if(i[0] == 1){
		return true;
	}
	else{
		return false;
	}


}



///returns a random double between 0 and the bound
///MAKE SURE TO SEED THE RANDOM NUMBER GENERATOR FIRST
double getRand(double bound)
{
	if(bound == 0.0){
		return 0.0;
	}
	double stuff = (double) rand();
	double stuff2 = (double)(stuff/RAND_MAX);
	double stuff3 = stuff2 * bound;
	return stuff3;

}

/*
double sq(double val)
{
	return pow(val, 2.0);
}
*/

///Given side lengths a, b, c, returns the angle opposite A, in degrees
double SSStriangleAngle(double a, double b, double c)
{
	if( (a+b-c > 0) && (a+b-c < .0001 )){			//process extreme triangles to prevent strange vals
		return 0.0;
	}
	if( (c+a-b > 0) && (c+a-b < .0001 )){			//process extreme triangles to prevent strange vals
		return 0.0;
	}
	if( (c+b-a > 0) && (c+b-a < .0001 )){			//process extreme triangles to prevent strange vals
		return 180.0;
	}

	double param, result;
	param = ((sq(a)-sq(b)-sq(c))/(2*b*c));
	result = 180 - (acos(param) * 180 / PI);
	return result;
}


////prints out a matrix for debugging
void printMatrix(double * mat, int size)
{
	switch(size){
	case 9:
		printf("%06lf %06lf %06lf \n", mat[0], mat[3], mat[6]);
		printf("%06lf %06lf %06lf \n", mat[1], mat[4], mat[7]);
		printf("%06lf %06lf %06lf \n", mat[2], mat[5], mat[8]);
		break;
	
	case 16:
		printf("%06lf %06lf %06lf %06lf\n", mat[0], mat[4], mat[8], mat[12]);
		printf("%06lf %06lf %06lf %06lf\n", mat[1], mat[5], mat[9], mat[13]);
		printf("%06lf %06lf %06lf %06lf\n", mat[2], mat[6], mat[10], mat[14]);
		printf("%06lf %06lf %06lf %06lf\n\n", mat[3], mat[7], mat[11], mat[15]);
		break;
	
	}
}

////transforms a 3-vector with a 4matrix, by padding a 1 for scale.
double * transformVector3x4(double * mat, double * vector)
{
	double * result = new double[3];

	result[0] = mat[0]*vector[0] + mat[4]*vector[1] + mat[8]*vector[2] + mat[12];
	result[1] = mat[1]*vector[0] + mat[5]*vector[1] + mat[9]*vector[2] + mat[13];
	result[2] = mat[2]*vector[0] + mat[6]*vector[1] + mat[10]*vector[2] + mat[14];

	return result;
}


///////////////////////////
///4x4 matrix times 4x4 matrix multiplication
double * matrixMult4x4(double * first, double * second)
{
	double * result = new double[16];

	result[0]  = first[ 0]*second[ 0] + first[ 4]*second[ 1] + first[ 8]*second[ 2] + first[12]*second[ 3];
	result[1]  = first[ 1]*second[ 0] + first[ 5]*second[ 1] + first[ 9]*second[ 2] + first[13]*second[ 3];
	result[2]  = first[ 2]*second[ 0] + first[ 6]*second[ 1] + first[10]*second[ 2] + first[14]*second[ 3];
	result[3]  = first[ 3]*second[ 0] + first[ 7]*second[ 1] + first[11]*second[ 2] + first[15]*second[ 3];

	result[4]  = first[ 0]*second[ 4] + first[ 4]*second[ 5] + first[ 8]*second[ 6] + first[12]*second[ 7];
	result[5]  = first[ 1]*second[ 4] + first[ 5]*second[ 5] + first[ 9]*second[ 6] + first[13]*second[ 7];
	result[6]  = first[ 2]*second[ 4] + first[ 6]*second[ 5] + first[10]*second[ 6] + first[14]*second[ 7];
	result[7]  = first[ 3]*second[ 4] + first[ 7]*second[ 5] + first[11]*second[ 6] + first[15]*second[ 7];

	result[8]  = first[ 0]*second[ 8] + first[ 4]*second[ 9] + first[ 8]*second[10] + first[12]*second[11];
	result[9]  = first[ 1]*second[ 8] + first[ 5]*second[ 9] + first[ 9]*second[10] + first[13]*second[11];
	result[10] = first[ 2]*second[ 8] + first[ 6]*second[ 9] + first[10]*second[10] + first[14]*second[11];
	result[11] = first[ 3]*second[ 8] + first[ 7]*second[ 9] + first[11]*second[10] + first[15]*second[11];

	result[12] = first[ 0]*second[12] + first[ 4]*second[13] + first[ 8]*second[14] + first[12]*second[15];
	result[13] = first[ 1]*second[12] + first[ 5]*second[13] + first[ 9]*second[14] + first[13]*second[15];
	result[14] = first[ 2]*second[12] + first[ 6]*second[13] + first[10]*second[14] + first[14]*second[15];
	result[15] = first[ 3]*second[12] + first[ 7]*second[13] + first[11]*second[14] + first[15]*second[15];

	return result;
}

bool testMatrix(double * matrix)
{
	double * matrixa = new double[16];
	matrixa[0] = matrix[0];
	matrixa[1] = matrix[1];
	matrixa[2] = matrix[2];
	matrixa[3] = matrix[3];
	matrixa[4] = matrix[4];
	matrixa[5] = matrix[5];
	matrixa[6] = matrix[6];
	matrixa[7] = matrix[7];
	matrixa[8] = matrix[8];
	matrixa[9] = matrix[9];
	matrixa[10] = matrix[10];
	matrixa[11] = matrix[11];
	matrixa[12] = 0;
	matrixa[13] = 0;
	matrixa[14] = 0;
	matrixa[15] = matrix[15];
//	printMatrix(matrixa, 16);
	double * matrix1 = transpose4x4(matrixa);
//	printMatrix(matrix1, 16);

	double * matrix2 = matrixMult4x4(matrix1, matrixa);
//	printMatrix(matrix2, 16);
//
	double * matrixI = new double[16];
	matrixI[0]  = 1 - matrix2[0];
	matrixI[1]  = 0 - matrix2[1];
	matrixI[2]  = 0 - matrix2[2];
	matrixI[3]  = 0 - matrix2[3];
	matrixI[4]  = 0 - matrix2[4];
	matrixI[5]  = 1 - matrix2[5];
	matrixI[6]  = 0 - matrix2[6];
	matrixI[7]  = 0 - matrix2[7];
	matrixI[8]  = 0 - matrix2[8];
	matrixI[9]  = 0 - matrix2[9];
	matrixI[10] = 1 - matrix2[10];
	matrixI[11] = 0 - matrix2[11];
	matrixI[12] = 0 - matrix2[12];
	matrixI[13] = 0 - matrix2[13];
	matrixI[14] = 0 - matrix2[14];
	matrixI[15] = 1 - matrix2[15];

//	printMatrix(matrixI, 16);

	int i = 0;
	for(i = 0; i<16; i++){
		if(matrixI[i] > 0.01 || matrixI[i] < -0.01 ){
			printf("WHOA\n");

			delete[](matrixa);
			delete[](matrix1);
			delete[](matrix2);
			delete[](matrixI);
			return false;
		}
	}

	delete[](matrixa);
	delete[](matrix1);
	delete[](matrix2);
	delete[](matrixI);
	return true;

}



bool testMatrix3x3(double * matrix)
{
	double * matrixa = new double[9];
	matrixa[0] = matrix[0];
	matrixa[1] = matrix[1];
	matrixa[2] = matrix[2];
	matrixa[3] = matrix[3];
	matrixa[4] = matrix[4];
	matrixa[5] = matrix[5];
	matrixa[6] = matrix[6];
	matrixa[7] = matrix[7];
	matrixa[8] = matrix[8];

	printMatrix(matrixa, 9);
	double * matrix1 = transpose(matrixa);
	printMatrix(matrix1, 9);

	double * matrix2 = matrixMult(matrix1, matrixa);
	printMatrix(matrix2, 9);

	double * matrixI = new double[9];
	matrixI[0]  = 1 - matrix2[0];
	matrixI[1]  = 0 - matrix2[1];
	matrixI[2]  = 0 - matrix2[2];
	matrixI[3]  = 0 - matrix2[3];
	matrixI[4]  = 1 - matrix2[4];
	matrixI[5]  = 0 - matrix2[5];
	matrixI[6]  = 0 - matrix2[6];
	matrixI[7]  = 0 - matrix2[7];
	matrixI[8]  = 1 - matrix2[8];

	printMatrix(matrixI, 9);

	int i = 0;
	for(i = 0; i<9; i++){
		if(matrixI[i] > 0.01 || matrixI[i] < -0.01 ){
			printf("WHOA\n");

			delete[](matrixa);
			delete[](matrix1);
			delete[](matrix2);
			delete[](matrixI);
			return false;
		}
	}

	delete[](matrixa);
	delete[](matrix1);
	delete[](matrix2);
	delete[](matrixI);

	return true;

}


/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void generateSideChainRMSD(int * matchArray, int numRes, Motif * motif, caList * target, double * RMSD, double * transform)
{
	int i = 0, j = 0, k = 0;
	int numPts = 0;
	PdbAtom * tempAtom;
	aminoAcid * tempAcid;

	///sort the match array by source indices
	int smallest = INT_MAX;
	int smallestInd = 0;
	for(i = 0; i<numRes; i++){
		for(j = i; j<numRes; j++){
			if(i == j){continue;}
			if(matchArray[2*j] <= smallest){
				smallest = matchArray[2*j];
				smallestInd = j;
			}
		}
		int tempSource = matchArray[2*i];
		int tempTarg = matchArray[2*i+1];
		matchArray[2*i] = matchArray[2*smallestInd];
		matchArray[2*i+1] = matchArray[2*smallestInd+1];
		matchArray[2*smallestInd] = tempSource;
		matchArray[2*smallestInd+1] = tempTarg;
		smallest = INT_MAX;
		smallestInd = 0;
	}

	///first count the number of point pairs
	for(i = 0; i<numRes; i++){
		for(j = 0; j<motif->numAtoms[matchArray[2*i]]; j++){
			tempAtom = motif->atoms[matchArray[2*i]][j];
			tempAcid = target->acids[matchArray[2*i+1]];
			if(translateAtomCode(tempAcid->Calpha->atom_name) == translateAtomCode(tempAtom->atom_name)){
				numPts++;
				continue;
			}
			else{
				for(k = 0; k<tempAcid->numAtoms; k++){
					if(translateAtomCode(tempAcid->atoms[k]->atom_name) == translateAtomCode(tempAtom->atom_name)){
						numPts++;
						break;
					}
				}	
			}
		}
	}

	double * * targetCoords = new double*[numPts];
	double * * sourceCoords = new double*[numPts];

	numPts = 0;

	///set all the coordinates
	for(i = 0; i<numRes; i++){
		for(j = 0; j<motif->numAtoms[matchArray[2*i]]; j++){
			tempAtom = motif->atoms[matchArray[2*i]][j];
			tempAcid = target->acids[matchArray[2*i+1]];
			if(translateAtomCode(tempAcid->Calpha->atom_name) == translateAtomCode(tempAtom->atom_name)){
				targetCoords[numPts] = tempAcid->Calpha->coords;
				sourceCoords[numPts] = tempAtom->coords;
				numPts++;
				continue;
			}
			else{
				for(k = 0; k<tempAcid->numAtoms; k++){
					if(translateAtomCode(tempAcid->atoms[k]->atom_name) == translateAtomCode(tempAtom->atom_name)){
						targetCoords[numPts] = tempAcid->atoms[k]->coords;
						sourceCoords[numPts] = tempAtom->coords;
						numPts++;
						break;
					}
				}	
			}
		}
	}

	double * resultArray;
	if(transform == NULL){
		resultArray = min_rmsd (targetCoords, sourceCoords, numPts, NULL);
	}
	else{
		resultArray = min_rmsd (targetCoords, sourceCoords, numPts, transform);
	}

	if(RMSD != NULL){
		RMSD[0] = resultArray[0];
	}

	////Do not delete to preserve the data in the lists.
//	for(i = 0; i<numPts; i++){
//		delete[](targetCoords[i]);
//		delete[](sourceCoords[i]);
//	}
	delete[](targetCoords);
	delete[](sourceCoords);
	delete[](resultArray);

}


  */




