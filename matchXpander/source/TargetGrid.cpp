// TargetGrid.cpp: implementation of the TargetGrid class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

///the reason ids is used this way is bc ids MUST be the same size as the # of atoms in teh AtomList
TargetGrid::TargetGrid( AtomList * a, int * ids )
{
	int i = 0;
	list = a;

	xneg = INT_MAX;
	xpos = INT_MIN;
	yneg = INT_MAX;
	ypos = INT_MIN;
	zneg = INT_MAX;
	zpos = INT_MIN;


	num = list->size;
	xcentroid = 0;
	ycentroid = 0;
	zcentroid = 0;

	for (i = 0; i<list->size; i++){
		xcentroid = xcentroid + list->atomList[i]->coords[0];
		ycentroid = ycentroid + list->atomList[i]->coords[1];
		zcentroid = zcentroid + list->atomList[i]->coords[2];
	}
	xcentroid = xcentroid/list->size;
	ycentroid = ycentroid/list->size;
	zcentroid = zcentroid/list->size;

	testProbeCoordx = xcentroid;
	testProbeCoordy = ycentroid;
	testProbeCoordz = zcentroid;

	testProbeAssocAA = new char[27];
	testProbeAssocAA[ 0] = 'A';
	testProbeAssocAA[ 1] = 'B';
	testProbeAssocAA[ 2] = 'C';
	testProbeAssocAA[ 3] = 'D';
	testProbeAssocAA[ 4] = 'E';
	testProbeAssocAA[ 5] = 'F';
	testProbeAssocAA[ 6] = 'G';
	testProbeAssocAA[ 7] = 'H';
	testProbeAssocAA[ 8] = 'I';
	testProbeAssocAA[ 9] = 'J';
	testProbeAssocAA[10] = 'K';
	testProbeAssocAA[11] = 'L';
	testProbeAssocAA[12] = 'M';
	testProbeAssocAA[13] = 'N';
	testProbeAssocAA[14] = 'O';
	testProbeAssocAA[15] = 'P';
	testProbeAssocAA[16] = 'Q';
	testProbeAssocAA[17] = 'R';
	testProbeAssocAA[18] = 'S';
	testProbeAssocAA[19] = 'T';
	testProbeAssocAA[20] = 'U';
	testProbeAssocAA[21] = 'V';
	testProbeAssocAA[22] = 'W';
	testProbeAssocAA[23] = 'X';
	testProbeAssocAA[24] = 'Y';
	testProbeAssocAA[25] = 'Z';
	testProbeAssocAA[26] = '\0';

	currentProbeFind = -1;

	setupTable(ids);
}



TargetGrid::~TargetGrid()
{
	///Do not delete list - it is provided and used elsewhere.
	///caList * list;
	int i = 0; 
	for( i = 0; i<xdim*ydim*zdim; i++){
		delete[](grid[i]);
	}
	if(xdim*ydim*zdim > 0){
		delete[](grid);
	}

	delete[](testProbeAssocAA);
}

///the reason ids is used this way is bc ids MUST be the same size as the # of atoms in teh AtomList
void TargetGrid::setupTable(int * ids)
{
	int i = 0;
	int j = 0;
	int k = 0;
	double tempx;
	double tempy;
	double tempz;
	PdbAtom * tempAtom;
	int x;
	int y;
	int z;

	///Find int boundaries
	for (i = 0; i<list->size; i++){
		tempAtom = list->atomList[i];
		tempx = tempAtom->coords[0];
		tempy = tempAtom->coords[1];
		tempz = tempAtom->coords[2];
		if (tempx < 0 && tempx != (int) tempx){tempx--;}
		if (tempy < 0 && tempy != (int) tempy){tempy--;}
		if (tempz < 0 && tempz != (int) tempz){tempz--;}

		if(tempx < xneg) {xneg = (int) tempx;}			/////inclusive boundaries
		if(tempx > xpos) {xpos = (int) tempx;}
		if(tempy < yneg) {yneg = (int) tempy;}
		if(tempy > ypos) {ypos = (int) tempy;}
		if(tempz < zneg) {zneg = (int) tempz;}
		if(tempz > zpos) {zpos = (int) tempz;}
	}

	////create a table of ints of the proper size
	xdim = xpos-xneg+1;
	ydim = ypos-yneg+1;
	zdim = zpos-zneg+1;

	///Set the dimensions to the number of Bins per side
	if (xdim % TARGETGRIDSIZE != 0){ xdim = (int) ( (xdim/TARGETGRIDSIZE) + 1 ); }
	else{ xdim = (int) xdim/TARGETGRIDSIZE; }
	if (ydim % TARGETGRIDSIZE != 0){ ydim = (int) ( (ydim/TARGETGRIDSIZE) + 1 ); }
	else{ ydim = (int) ydim/TARGETGRIDSIZE; }
	if (zdim % TARGETGRIDSIZE != 0){ zdim = (int) ( (zdim/TARGETGRIDSIZE) + 1 ); }
	else{ zdim = (int) zdim/TARGETGRIDSIZE; }

	///Allocate the grid
	grid = new int*[xdim*ydim*zdim];

	///Fill the grid with empties
	for (i = 0; i<xdim*ydim*zdim; i++){
		grid[i] = new int[1];
		grid[i][0] = -1;
	}

	///fill the grid.
	for (i = 0; i<list->size; i++){
		j = 0;
		tempAtom = list->atomList[i];

		tempx = tempAtom->coords[0];
		if (tempx < 0 && tempx != (int) tempx){tempx--;}
		x = ( ((int) tempx)-xneg )/TARGETGRIDSIZE;
		tempy = tempAtom->coords[1];
		if (tempy < 0 && tempy != (int) tempy){tempy--;}
		y = ( ((int) tempy)-yneg )/TARGETGRIDSIZE;
		tempz = tempAtom->coords[2];
		if (tempz < 0 && tempz != (int) tempz){tempz--;}
		z = ( ((int) tempz)-zneg )/TARGETGRIDSIZE;
	
		///x, y, and z are now the coordinates in which to place this point. So now we add it.
		///first measure the size of the bin:
		while (grid[(z*(xdim*ydim)) + (y*(xdim)) + x][j] != -1){
			j++;
		}

		j = j+2;	///j + 2 is the size of the array with 1 more element (j++ doesnt terminate at the end, so need +2
					///for the -1 terminator and one for the new value                          )

		///fill the array with the correct elements
		int * tempArray = new int[j];
		for (k = 0; k<j-2; k++){
			tempArray[k] = grid[(z*(xdim*ydim)) + (y*(xdim)) + x][k];
		}
		tempArray[j-2] = ids[i];
		tempArray[j-1] = -1;

		///swap the new arrya with the old one.
		delete[](grid[(z*(xdim*ydim)) + (y*(xdim)) + x]);
		grid[(z*(xdim*ydim)) + (y*(xdim)) + x] = tempArray;
	}
}




void TargetGrid::currentBoxPlus()
{
	if (currentBox + 1 < xdim*ydim*zdim){
		currentBox++;
	}
}

void TargetGrid::currentBoxMinus()
{
	if (currentBox - 1 >= 0){
		currentBox--;
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
///Shell VERSION (Seed Finding Use)
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////THe result is packed with the size in the 0th position
int * TargetGrid::queryShell(double x, double y, double z, set_t set, double radius)
{
	int * result = NULL;

	double maxVal = radius + EXTMATCHTHRESHOLD;
	double minVal = radius - EXTMATCHTHRESHOLD;
	if(minVal < 0){ minVal = 0.0; }

	///elimiate things which are way off
	if (x<( (double)xneg - maxVal) || 
		x>( (double)xpos + maxVal + 1) ||
		y<( (double)yneg - maxVal) ||
		y>( (double)ypos + maxVal + 1) ||
		z<( (double)zneg - maxVal) ||
		z>( (double)zpos + maxVal + 1) ) {
		result = new int[1];
		result[0] = 0;
		return result;
	}

	///Initialize variables
	int xmax, xmin, ymax, ymin, zmax, zmin;		///these define the boundary of the cube bounding the threshold radius sphere	int boxX, boxY, boxZ;
	int * tempInd;
	int tempIndex;
	int boxX, boxY, boxZ;
	double tempClosest;
//	double tempx, tempy, tempz;
	double tempVal;
	int * *testArrays;
	int i = 0;
	int j = 0;
	int k = 0;
	int counter = 0;
	PdbAtom * tempAtom;

	///find out which box the point is in.
	if (x < 0 && x != ((int) x) ){tempVal = x-1;}
	else {tempVal = x;}
	boxX = ( ((int) tempVal)-xneg )/TARGETGRIDSIZE;

	if (y < 0 && y != ((int) y) ){tempVal = y-1;}
	else {tempVal = y;}
	boxY = ( ((int) tempVal)-yneg )/TARGETGRIDSIZE;

	if (z < 0 && z != ((int) z) ){tempVal = z-1;}
	else {tempVal = z;}
	boxZ = ( ((int) tempVal)-zneg )/TARGETGRIDSIZE;

	///Clamp outer values
	if (boxX < 0){boxX = 0;}
	if (boxX >= xdim){boxX = xdim-1;}
	if (boxY < 0){boxY = 0;}
	if (boxY >= ydim){boxY = ydim-1;}
	if (boxZ < 0){boxZ = 0;}
	if (boxZ >= zdim){boxZ = zdim-1;}

	///find out the boundaries of the bounding cube before considering the boundaries of the table
	if ( (maxVal == (int) maxVal) && ( ((int) maxVal) % TARGETGRIDSIZE == 0) ){ 
		xmax = boxX + ((int) (maxVal/TARGETGRIDSIZE) );
		xmin = boxX - ((int) (maxVal/TARGETGRIDSIZE) );
		ymax = boxY + ((int) (maxVal/TARGETGRIDSIZE) );
		ymin = boxY - ((int) (maxVal/TARGETGRIDSIZE) );
		zmax = boxZ + ((int) (maxVal/TARGETGRIDSIZE) );
		zmin = boxZ - ((int) (maxVal/TARGETGRIDSIZE) );
	}
	else{ 
		xmax = boxX + ((int) (maxVal/TARGETGRIDSIZE) ) + 1;
		xmin = boxX - ((int) (maxVal/TARGETGRIDSIZE) ) - 1;
		ymax = boxY + ((int) (maxVal/TARGETGRIDSIZE) ) + 1;
		ymin = boxY - ((int) (maxVal/TARGETGRIDSIZE) ) - 1;
		zmax = boxZ + ((int) (maxVal/TARGETGRIDSIZE) ) + 1;
		zmin = boxZ - ((int) (maxVal/TARGETGRIDSIZE) ) - 1;
	}

	///now correct xmax, xmin, etc so that they are not outside table coordinates
	if (xmin < 0){xmin = 0;}
	if (xmax >= xdim){xmax = xdim-1;}
	if (ymin < 0){ymin = 0;}
	if (ymax >= ydim){ymax = ydim-1;}
	if (zmin < 0){zmin = 0;}
	if (zmax >= zdim){zmax = zdim-1;}

	///now we construct the list of points to test
	testArrays = new int*[(xmax-xmin+1)*(ymax-ymin+1)*(zmax-zmin+1)];

	for(i = zmin; i <= zmax; i++){
		for(j = ymin; j <= ymax; j++){
			for(k = xmin; k <= xmax; k++){
				testArrays[counter] = grid[( i*(ydim)*(xdim) ) + (j*(xdim)) + k];
				counter++;
			}
		}
	}

	///run through all arrays that need to be considered and find the point 
	///fulfilling the property that it
	///1) under the threshold (must)
	j = 0;
	counter = 0;
	for (i = 0; i<(xmax-xmin+1)*(ymax-ymin+1)*(zmax-zmin+1); i++){
		while (testArrays[i][j] != -1){
			////use the set mapping to get the index of the associated PdbAtom in this targetGrid's list
			tempInd = (int *) mapsto_set(set, testArrays[i][j]);
			tempIndex = tempInd[0];
			tempAtom = list->atomList[tempIndex];

			tempClosest = vectorSize(tempAtom->coords[0] - x, tempAtom->coords[1] - y, tempAtom->coords[2] - z);

			if(tempClosest <= maxVal && tempClosest >= minVal){
				counter++;
			}
			j++;
		}
		j = 0;
	}

	result = new int[counter+1];
	result[0] = counter;
	counter = 1;					////offset by 1 to fit the size in first.

	for (i = 0; i<(xmax-xmin+1)*(ymax-ymin+1)*(zmax-zmin+1); i++){
		while (testArrays[i][j] != -1){
			////use the set mapping to get the index of the associated PdbAtom in this targetGrid's list
			tempInd = (int *) mapsto_set(set, testArrays[i][j]);
			tempIndex = tempInd[0];
			tempAtom = list->atomList[tempIndex];

			tempClosest = vectorSize(tempAtom->coords[0] - x, tempAtom->coords[1] - y, tempAtom->coords[2] - z);

			if(tempClosest <= maxVal && tempClosest >= minVal){
				result[counter] = testArrays[i][j];		/////output the residue # of the matching atom
				counter++;
			}
			j++;
		}
		j = 0;
	}


	delete[](testArrays);

	return result;
}





////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
///Sphere VERSION (Match Augmentation use)
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////THe result is packed with the size in the 0th position
set_t TargetGrid::querySphere(double x, double y, double z, set_t set, double radius, set_t output)
{
	double maxVal = radius;

	///elimiate things which are way off
	if (x<( (double)xneg - maxVal) || 
		x>( (double)xpos + maxVal + 1) ||
		y<( (double)yneg - maxVal) ||
		y>( (double)ypos + maxVal + 1) ||
		z<( (double)zneg - maxVal) ||
		z>( (double)zpos + maxVal + 1) ) {
		return output;
	}

	///Initialize variables
	int xmax, xmin, ymax, ymin, zmax, zmin;		///these define the boundary of the cube bounding the threshold radius sphere	int boxX, boxY, boxZ;
	int * tempInd;
	int tempIndex;
	int boxX, boxY, boxZ;
	double tempClosest;
//	double tempx, tempy, tempz;
	double tempVal;
	int * *testArrays;
	int i = 0;
	int j = 0;
	int k = 0;
	int counter = 0;
	PdbAtom * tempAtom;

	///find out which box the point is in.
	if (x < 0 && x != ((int) x) ){tempVal = x-1;}
	else {tempVal = x;}
	boxX = ( ((int) tempVal)-xneg )/TARGETGRIDSIZE;

	if (y < 0 && y != ((int) y) ){tempVal = y-1;}
	else {tempVal = y;}
	boxY = ( ((int) tempVal)-yneg )/TARGETGRIDSIZE;

	if (z < 0 && z != ((int) z) ){tempVal = z-1;}
	else {tempVal = z;}
	boxZ = ( ((int) tempVal)-zneg )/TARGETGRIDSIZE;

	///Clamp outer values
	if (boxX < 0){boxX = 0;}
	if (boxX >= xdim){boxX = xdim-1;}
	if (boxY < 0){boxY = 0;}
	if (boxY >= ydim){boxY = ydim-1;}
	if (boxZ < 0){boxZ = 0;}
	if (boxZ >= zdim){boxZ = zdim-1;}

	///find out the boundaries of the bounding cube before considering the boundaries of the table
	if ( (maxVal == (int) maxVal) && ( ((int) maxVal) % TARGETGRIDSIZE == 0) ){ 
		xmax = boxX + ((int) (maxVal/TARGETGRIDSIZE) );
		xmin = boxX - ((int) (maxVal/TARGETGRIDSIZE) );
		ymax = boxY + ((int) (maxVal/TARGETGRIDSIZE) );
		ymin = boxY - ((int) (maxVal/TARGETGRIDSIZE) );
		zmax = boxZ + ((int) (maxVal/TARGETGRIDSIZE) );
		zmin = boxZ - ((int) (maxVal/TARGETGRIDSIZE) );
	}
	else{ 
		xmax = boxX + ((int) (maxVal/TARGETGRIDSIZE) ) + 1;
		xmin = boxX - ((int) (maxVal/TARGETGRIDSIZE) ) - 1;
		ymax = boxY + ((int) (maxVal/TARGETGRIDSIZE) ) + 1;
		ymin = boxY - ((int) (maxVal/TARGETGRIDSIZE) ) - 1;
		zmax = boxZ + ((int) (maxVal/TARGETGRIDSIZE) ) + 1;
		zmin = boxZ - ((int) (maxVal/TARGETGRIDSIZE) ) - 1;
	}

	///now correct xmax, xmin, etc so that they are not outside table coordinates
	if (xmin < 0){xmin = 0;}
	if (xmax >= xdim){xmax = xdim-1;}
	if (ymin < 0){ymin = 0;}
	if (ymax >= ydim){ymax = ydim-1;}
	if (zmin < 0){zmin = 0;}
	if (zmax >= zdim){zmax = zdim-1;}

	///now we construct the list of points to test
	testArrays = new int*[(xmax-xmin+1)*(ymax-ymin+1)*(zmax-zmin+1)];

	for(i = zmin; i <= zmax; i++){
		for(j = ymin; j <= ymax; j++){
			for(k = xmin; k <= xmax; k++){
				testArrays[counter] = grid[( i*(ydim)*(xdim) ) + (j*(xdim)) + k];
				counter++;
			}
		}
	}

	///run through all arrays that need to be considered and find the point 
	///fulfilling the property that it
	///1) under the threshold (must)
	j = 0;
	for (i = 0; i<(xmax-xmin+1)*(ymax-ymin+1)*(zmax-zmin+1); i++){
		while (testArrays[i][j] != -1){
			////use the set mapping to get the index of the associated PdbAtom in this targetGrid's list
			tempInd = (int *) mapsto_set(set, testArrays[i][j]);
			tempIndex = tempInd[0];
			tempAtom = list->atomList[tempIndex];

 			tempClosest = vectorSize(tempAtom->coords[0] - x, tempAtom->coords[1] - y, tempAtom->coords[2] - z);

			if(tempClosest <= maxVal){
				output = put_set(output, testArrays[i][j]);
			}
			j++;
		}
		j = 0;
	}

	delete[](testArrays);

	return output;
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////RENDERING SECTION////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef INTERACTIVEMODE
//////////////////////////////////////////


void TargetGrid::renderMe()
{
	int i = 0;
	int j = 0;
	int k = 0;
	int tempx = 0;
	int tempy = 0;
	int tempz = 0;
	int currentInd = 0;
	int tempVal = 0;
	float x = 0.0;
	float y = 0.0;
	float z = 0.0;

	///translate to the centroid:
	glTranslatef(-xcentroid, -ycentroid, -zcentroid);
	glPushMatrix();

////////////////////////////////////////DEBUG
	//render the test probe
	glColor3f(1.0, 1.0, 1.0);

	glPushMatrix();
	glTranslatef( (float) testProbeCoordx, (float) testProbeCoordy, (float) testProbeCoordz );

	glutWireSphere(EXTMATCHTHRESHOLD, 2*ATOMVERTSUBS, 2*ATOMHORISUBS);

	glPopMatrix();

	//render the found point
	if (currentProbeFind != -1){
		glColor3f(1.0, 0.0, 1.0);
		glPushMatrix();

		glTranslatef( 
			(float) motif->constellation[3*currentProbeFind+0], 
			(float) motif->constellation[3*currentProbeFind+1], 
			(float) motif->constellation[3*currentProbeFind+2]
		);
		glutWireSphere(1.0, 8*ATOMVERTSUBS, 8*ATOMHORISUBS);

		glPopMatrix();
	}

////////////////////////////////////////DEBUG

	//render the points
	glColor3f(0.0, 1.0, 0.0);
	for (i = 0; i<num; i++){

		glPushMatrix();
		glTranslatef( (float) motif->constellation[3*i+0], (float) motif->constellation[3*i+1], (float) motif->constellation[3*i+2] );

		glutSolidSphere(.3, ATOMVERTSUBS, ATOMHORISUBS);

		glPopMatrix();
	}

	//render the grid
	tempx = xneg;
	tempy = yneg;
	tempz = zneg;

	glBegin(GL_LINES);
	for (i = 0; i<zdim; i++){
		for (j = 0; j<ydim; j++){
			for (k = 0; k<xdim; k++){
				tempx = xneg + (TARGETGRIDSIZE * k);
				tempy = yneg + (TARGETGRIDSIZE * j);
				tempz = zneg + (TARGETGRIDSIZE * i);
				//make the cube here.
				glColor3f(1.0, 1.0, 1.0);

				///highlight the "currentBox"
				if (i*(xdim*ydim) + j*(xdim) + k == currentBox){
					glColor3f(0.0, 0.0, 1.0);
				}
				glVertex3f(tempx,	tempy,	tempz);		glVertex3f(tempx + TARGETGRIDSIZE,	tempy,	tempz);
				glVertex3f(tempx,	tempy,	tempz);		glVertex3f(tempx,	tempy + TARGETGRIDSIZE,	tempz);
				glVertex3f(tempx,	tempy,	tempz);		glVertex3f(tempx,	tempy,	tempz + TARGETGRIDSIZE);
				glVertex3f(tempx + TARGETGRIDSIZE,	tempy,	tempz);		glVertex3f(tempx + TARGETGRIDSIZE,	tempy + TARGETGRIDSIZE,	tempz);
				glVertex3f(tempx,	tempy + TARGETGRIDSIZE,	tempz);		glVertex3f(tempx + TARGETGRIDSIZE,	tempy + TARGETGRIDSIZE,	tempz);
				glVertex3f(tempx + TARGETGRIDSIZE,	tempy,	tempz);		glVertex3f(tempx + TARGETGRIDSIZE,	tempy,	tempz + TARGETGRIDSIZE);
				glVertex3f(tempx,	tempy + TARGETGRIDSIZE,	tempz);		glVertex3f(tempx,	tempy + TARGETGRIDSIZE,	tempz + TARGETGRIDSIZE);
				glVertex3f(tempx + TARGETGRIDSIZE,	tempy + TARGETGRIDSIZE,	tempz);		glVertex3f(tempx + TARGETGRIDSIZE,	tempy + TARGETGRIDSIZE,	tempz + TARGETGRIDSIZE);
				glVertex3f(tempx,	tempy,	tempz + TARGETGRIDSIZE);		glVertex3f(tempx + TARGETGRIDSIZE,	tempy,	tempz + TARGETGRIDSIZE);
				glVertex3f(tempx,	tempy,	tempz + TARGETGRIDSIZE);		glVertex3f(tempx,	tempy + TARGETGRIDSIZE,	tempz + TARGETGRIDSIZE);
				glVertex3f(tempx + TARGETGRIDSIZE,	tempy,	tempz + TARGETGRIDSIZE);	glVertex3f(tempx + TARGETGRIDSIZE,	tempy + TARGETGRIDSIZE,	tempz + TARGETGRIDSIZE);
				glVertex3f(tempx,	tempy + TARGETGRIDSIZE,	tempz + TARGETGRIDSIZE);	glVertex3f(tempx + TARGETGRIDSIZE,	tempy + TARGETGRIDSIZE,	tempz + TARGETGRIDSIZE);
			}
		}
	}
	glEnd();

	//render the highlighted points
	i = 0;
	tempVal = grid[currentBox][i];
	while (tempVal != -1){
		x = (float) motif->constellation[3*(tempVal)+0];
		y = (float) motif->constellation[3*(tempVal)+1];
		z = (float) motif->constellation[3*(tempVal)+2];

		glColor3f(1.0, 0.0, 0.0);

		glPushMatrix();
		glTranslatef( x, y, z );
		glutSolidSphere(.3, ATOMVERTSUBS, ATOMHORISUBS);
		glPopMatrix();

		i++;
		tempVal = grid[currentBox][i];
	}

	///pop the centroid translation
	glPopMatrix();

}


//////////////////////////////////////////
#endif
//////////////////////////////////////////



