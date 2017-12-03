// AtomList.cpp: implementation of the AtomList class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

AtomList::AtomList(int s){
	
	size = s;
	atomList = (PdbAtom **) new PdbAtom*[s];
	maxSeq = 0;
	highlights = new bool[s];

	int i = 0;
	for(i = 0; i<s; i++)
	{
		atomList[i] = new PdbAtom();
		atomList[i]->significance = -1;			//set all sigs to -1.
		highlights[i] = false;
	}

	name = new char[80];
	selectx = 0;
	selecty = 0;
	selectz = 0;
	radius = 5;
	renderSelectBox = true;

}

AtomList::AtomList(int s, PdbAtom * *atoms)
{
	size = s;
	maxSeq = 0;
	name = NULL;
	highlights = NULL;
	selectx = 0;
	selecty = 0;
	selectz = 0;
	radius = 5;
	renderSelectBox = true;

	atomList = (PdbAtom **) new PdbAtom*[s];

	int i = 0;
	for(i = 0; i<s; i++){
		///copy all atoms
		atomList[i] = atoms[i]->copy();
	}
}


AtomList * AtomList::copy(void)
{
	AtomList * result = new AtomList(size);
	result->maxSeq = maxSeq;
	strcpy(result->name, this->name);

	int i = 0;
	for(i = 0; i<size; i++){
		///copy all highlights
		result->highlights[i] = this->highlights[i];
		///copy all atoms
		result->atomList[i]->assocAA = new char[NUMCOLUMNS];
		strcpy(result->atomList[i]->assocAA, this->atomList[i]->assocAA);
		strcpy(result->atomList[i]->atom_name, this->atomList[i]->atom_name);
		strcpy(result->atomList[i]->residue_name, this->atomList[i]->residue_name);
		strcpy(result->atomList[i]->fileLine, this->atomList[i]->fileLine);
		result->atomList[i]->atom_number = this->atomList[i]->atom_number;
		result->atomList[i]->residue_number = this->atomList[i]->residue_number;
		result->atomList[i]->significance = this->atomList[i]->significance;
		result->atomList[i]->coords[0] = this->atomList[i]->coords[0];
		result->atomList[i]->coords[1] = this->atomList[i]->coords[1];
		result->atomList[i]->coords[2] = this->atomList[i]->coords[2];
	}

	return result;
}




AtomList::~AtomList()
{
	int i = 0;
	for (i = 0; i<size; i++){
		delete(atomList[i]);
	}
	delete[](atomList);
	delete[](highlights);
	delete[](name);				///??
}

////////////////////////////////
/////accessor method
PdbAtom * AtomList::atomNumber(int i)
{

	return atomList[i];

}


/////////////////////////////////////
//////debugging output and class info
void AtomList::toString()
{
	for (int i = 0; i<size; i++)
	{
		atomList[i]->toString();
	}
}


///////////////////////
///Set the highlight bit if it's not highlighted
///Set back to zero if highlighted
/// - make a pass and find the atoms of the same AA and highlight. (we cant ensure data locality)
void AtomList::highlightAA(int index)
{
	int aaNum = atomNumber(index)->residue_number;

	int i = 0;
	for(i = 0; i<size; i++){
		if (atomNumber(i)->residue_number == aaNum){
			if (highlights[i] == false){ 
				highlights[i] = true; 
			}
			else if (highlights[i] == true){ highlights[i] = false; }
		}
	}

}


/////////////////////////////////////
///Highlights all Atoms under a given threshold
void AtomList::highlightUnderThresh(int thresh)
{
	int i = 0;
	for(i = 0; i<size; i++){
		if (atomNumber(i)->significance < thresh){
			highlights[i] = true;
		}
	}

}


void AtomList::highlightAllAAs(void){
	int i = 0;
	for (i = 0; i<size; i++){
		highlights[i] = true;	
	}
}


void AtomList::highlightNoAAs(void){
	int i = 0;
	for (i = 0; i<size; i++){
		highlights[i] = false;	
	}
}


void AtomList::moveSelectionBox(float dx, float dy, float dz)
{
	selectx = selectx + dx;
	selecty = selecty + dy;
	selectz = selectz + dz;
}


void AtomList::resizeSelectBox(float dr)
{
	if ( (radius + dr) >= 1 ){
		radius = radius + dr;
	}
}

void AtomList::highlightInsideBox(bool yesNo)
{
	///If yesNo = true, we select
	///if yesNo = false we deselect
	int i = 0;
	bool xVal;
	bool yVal;
	bool zVal;
	double * pts = new double[3*size];

	///Get the Centroid
	for (i = 0; i<size; i++){
		pts[3*i+0] = atomList[i]->coords[0];
		pts[3*i+1] = atomList[i]->coords[1];
		pts[3*i+2] = atomList[i]->coords[2];
	}
	double * cen = getCentroid(pts, size);

	for(i = 0; i<size; i++){
		xVal = ( ( (atomNumber(i)->coords[0]-cen[0]) >= selectx-(radius/2) ) && ((atomNumber(i)->coords[0]-cen[0]) <= selectx+(radius/2)) );
		yVal = ( ( (atomNumber(i)->coords[1]-cen[1]) >= selecty-(radius/2) ) && ((atomNumber(i)->coords[1]-cen[1]) <= selecty+(radius/2)) );
		zVal = ( ( (atomNumber(i)->coords[2]-cen[2]) >= selectz-(radius/2) ) && ((atomNumber(i)->coords[2]-cen[2]) <= selectz+(radius/2)) );
		///check to see if it's inside the box
		if ( xVal & yVal & zVal ){
			if (yesNo){ highlights[i] = true; }
			else { highlights[i] = false; }
		}
	}

	highlightAssociatedAAs(yesNo);
}


int * AtomList::getHighlightedAtoms()
{
	int i = 0;
	int numHighlighted = 0;
	for(i = 0; i<size; i++){
		if(highlights[i]){
			numHighlighted++;
		}
	}

	int * result = new int[numHighlighted+1];
	result[0] = numHighlighted;

	int counter = 1;
	for(i = 0; i<size; i++){
		if(highlights[i]){
			result[counter] = i;
			counter++;
		}	
	}

	return result;
}


int * AtomList::getHighlightedAtomNumbers()
{
	int i = 0;
	int numHighlighted = 0;
	for(i = 0; i<size; i++){
		if(highlights[i]){
			numHighlighted++;
		}
	}

	int * result = new int[numHighlighted+1];
	result[0] = numHighlighted;

	int counter = 1;
	for(i = 0; i<size; i++){
		if(highlights[i]){
			result[counter] = atomList[i]->atom_number;
			counter++;
		}	
	}

	return result;
}

///first entry packed with size
void AtomList::highlightAssociatedAAs(int * atoms)
{
	int i = 0, j = 0;
	int numInts = atoms[0];
	set_t temp = alloc_set(0);
	set_t temp2 = alloc_set(0);

	for(i = 0; i<size; i++){
		highlights[i] = false;
	}

	for(i = 0; i<numInts; i++){
		int residueNum = atomList[atoms[i+1]]->residue_number;
		temp2 = put_set(temp2, residueNum);
		for(j = 0; j<size; j++){
			if(atomList[j]->residue_number == residueNum){
				highlights[j] = true;
				temp  = put_set(temp, atomList[j]->residue_number);
			} 
		}
	}

	free_set(temp);
	free_set(temp2);


}

///////////////////////////
///input true: If any atom of an AA is highlighted, highlight every other atom of it
///input false: If any atom of an AA is not highlighted, de-highlight every other atom of it
///Performance: Very slow n^2 algorithm.  Executed in user time, so its ok.
void AtomList::highlightAssociatedAAs(bool info)
{
	int i = 0;
	int j = 0;
	int tempVal1;
	int tempVal2;

	///Process postive Case
	if(info)
	{	///Loop outside
		for(i = 0; i<size; i++){
			///ONly process if we will make assignment
			if(!highlights[i]){
				///Loop through all atoms
				for (j = 0; j<size; j++){
					tempVal1 = atomNumber(i)->residue_number;
					tempVal2 = atomNumber(j)->residue_number;
					///Set the value only if matching residues and other is highlighted
					if( highlights[j] && (tempVal1 == tempVal2) ){
						highlights[i] = true;
						break;					
					}
				}
			}
		}
	}
	///Process negative Case
	else
	{	///Loop outside
		for(i = 0; i<size; i++){
			///ONly process if we will make assignment
			if(highlights[i]){
				///Loop through all atoms
				for(j = 0; j<size; j++){
					tempVal1 = atomNumber(i)->residue_number;
					tempVal2 = atomNumber(j)->residue_number;
					///Set the value only if matching residues and other is highlighted
					if(!highlights[j] && (tempVal1 == tempVal2) ){
						highlights[i] = false;
						break;
					}
				}
			}
		}		
	}
}


///////////////////////////
//Gives a mapping of all MX-qualified (top tier, or top 3) residues to their positions in the motif as they ought to be,
//So that a proper FAT file can be reconstructed even though GeoHash really only operates on the qualified residues.
//The mapping as as follows:
//	the first int is the size of the array not including the first.  this is always equal to the size of the original motif specified
//  the value of each remaining position is either -1 if that member was not in the seed motif, some number if it was.
int * AtomList::getMXmapping(bool AAorAtomMode)
{
	int i = 0;
	int j = 0;
	PdbAtom * tempAtom;
	int counter = 0;		//this is the number of highlighted residues
	int counter2 = 1;
	int * tempArray1 = NULL;
	int * tempArray2 = NULL;
	int tempVal;
	int numSignificant;
	int * result = NULL;

	switch(AAorAtomMode){
		///Amino Acids
		case true: 
			////count the selected Amino Acids
			for(i = 0; i<size; i++){
				tempAtom = atomList[i];
				if(strcmp(tempAtom->atom_name, "CA") == 0 && highlights[i]){
					counter++;
				}
			}

			tempArray1 = new int[counter];

			////record the mapping of highlighted residues by their signifiance
			counter = 0;
			for(i = 0; i<size; i++){
				tempAtom = atomList[i];
				if(strcmp(tempAtom->atom_name, "CA") == 0 && highlights[i]){
					tempArray1[counter] = tempAtom->significance;
					counter++;	///this is the number of highlighted residues
				}
			}

			////sort these significances
			insertionSortAscending(tempArray1, counter);

			//find the significance level we care about
			tempVal = tempArray1[2];	///we care about the significance of the third atom.

			//count how many are significant
			numSignificant = 0;
			for(i = 0; i<counter; i++){
				if(tempArray1[i] <= tempVal){
					numSignificant++;
				}
			}

			///get indices of every significant atom
			counter2 = 0;
			tempArray2 = new int[numSignificant];
			for(i = 0; i<size; i++){
				tempAtom = atomList[i];
				if(strcmp(tempAtom->atom_name, "CA") == 0 && highlights[i] && tempAtom->significance <= tempVal){
					tempArray2[counter2] = tempAtom->residue_number;
					counter2++;
				}
			
			}

			//define the new array.  indices in this array are the indices of the pts in the big motif
			result = new int[counter+1];
			result[0] = counter;

			///fill the mapping with -1s where nothign is mapped.
			for(i = 1; i<counter+1; i++){
				result[i] = -1;
			}

			///now insert the maps.
			counter2 = 1;
			for(i = 0; i<size; i++){
				tempAtom = atomList[i];
				if(strcmp(tempAtom->atom_name, "CA") == 0 && highlights[i]){
					for(j = 0; j<numSignificant; j++){
						if(tempAtom->residue_number == tempArray2[j]){
							result[counter2] = j;
							break;
						}
					}
					counter2++;
				}
			}

		break;

		//////////////////////////////////////////////////////////////////////////
		///Atoms//////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////
		case false:
			////count the selected Amino Acids
			for(i = 0; i<size; i++){
				tempAtom = atomList[i];
				if(highlights[i]){
					counter++;
				}
			}

			tempArray1 = new int[counter];

			////record the mapping of highlighted residues by their signifiance
			counter = 0;
			for(i = 0; i<size; i++){
				tempAtom = atomList[i];
				if(highlights[i]){
					tempArray1[counter] = tempAtom->significance;
					counter++;	///this is the number of highlighted residues
				}
			}

			////sort these significances
			insertionSortAscending(tempArray1, counter);

			//find the significance level we care about
			tempVal = tempArray1[2];	///we care about the significance of the third atom.

			//count how many are significant
			numSignificant = 0;
			for(i = 0; i<counter; i++){
				if(tempArray1[i] <= tempVal){
					numSignificant++;
				}
			}

			///get indices of every significant atom
			counter2 = 0;
			tempArray2 = new int[numSignificant];
			for(i = 0; i<size; i++){
				tempAtom = atomList[i];
				if(highlights[i] && tempAtom->significance <= tempVal){
					tempArray2[counter2] = i;
					counter2++;
				}
			
			}

			//define the new array.  indices in this array are the indices of the pts in the big motif
			result = new int[counter+1];
			result[0] = counter;

			///fill the mapping with -1s where nothign is mapped.
			for(i = 1; i<counter+1; i++){
				result[i] = -1;
			}

			///now insert the maps.
			counter2 = 1;
			for(i = 0; i<size; i++){
				tempAtom = atomList[i];
				if(highlights[i]){
					for(j = 0; j<numSignificant; j++){
						if(i == tempArray2[j]){
							result[counter2] = j;
							break;
						}
					}
					counter2++;
				}
			}
	
		break;
	
	}

	delete[](tempArray1);
	delete[](tempArray2);

	return result;

}


///returns the MX threshold we wish to use, and also de-highlights all atoms below this threshold.
int AtomList::setMXthresh(bool AAorAtomMode)
{
	int i = 0;
	PdbAtom * tempAtom;
	int counter = 0;
	int * tempArray1 = NULL;
	int tempVal;

	/////Amino Acid Processing
	if(AAorAtomMode){
		for(i = 0; i<size; i++){
			tempAtom = atomList[i];
			if(strcmp(tempAtom->atom_name, "CA") == 0 && highlights[i]){
				counter++;
			}
		}

		tempArray1 = new int[counter];
		counter = 0;

		for(i = 0; i<size; i++){
			tempAtom = atomList[i];
			if(strcmp(tempAtom->atom_name, "CA") == 0 && highlights[i]){
				tempArray1[counter] = tempAtom->significance;
				counter++;
			}
		}

		insertionSortAscending(tempArray1, counter);

		tempVal = tempArray1[2];	///we care about the significance of the third atom.

		for(i = 0; i<size; i++){
			tempAtom = atomList[i];
			if(tempAtom->significance > tempVal && highlights[i]){
				highlights[i] = false;
			}
		}
	}

	/////Atom processing
	else{
		for(i = 0; i<size; i++){
			tempAtom = atomList[i];
			if(highlights[i]){
				counter++;
			}
		}

		tempArray1 = new int[counter];
		counter = 0;

		for(i = 0; i<size; i++){
			tempAtom = atomList[i];
			if(highlights[i]){
				tempArray1[counter] = tempAtom->significance;
				counter++;
			}
		}

		insertionSortAscending(tempArray1, counter);

		tempVal = tempArray1[2];	///we care about the significance of the third atom.

		for(i = 0; i<size; i++){
			tempAtom = atomList[i];
			if(tempAtom->significance > tempVal && highlights[i]){
				highlights[i] = false;
			}
		}
	}

	delete[](tempArray1);
	return tempVal;

}





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////RENDERING SECTION////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef INTERACTIVEMODE
//////////////////////////////////////////

//////////////////////////////////////////////////////
///////////Drawing method for molecule visualization
void AtomList::renderMolecule(int etDisplay, int pickMode)
{
	double xpos = 0.0;
	double ypos = 0.0;
	double zpos = 0.0;
	int sig = 0;
	int rendermode = etDisplay; 
	int i = 0;
	double * pts = new double[3*size];
	float rad;
	int horiz = ATOMHORISUBS;
	int vert = ATOMVERTSUBS;

	for (i = 0; i<size; i++){
		pts[3*i+0] = atomList[i]->coords[0];
		pts[3*i+1] = atomList[i]->coords[1];
		pts[3*i+2] = atomList[i]->coords[2];
	}

	////FIXME: REMOVE THIS IN FUTURE VERSIONS
	double * centroid = getCentroid(pts, size);

	///bounds checking
	if (rendermode < 0) rendermode = 0;
	if (rendermode > maxSeq) rendermode = maxSeq;

	for (i = 0; i<size; i++){

		///set the temp vals from the atomList and render a sphere in that location
		xpos = pts[3*i+0] - centroid[0];
		ypos = pts[3*i+1] - centroid[1];
		zpos = pts[3*i+2] - centroid[2];

		if(pickMode == GL_RENDER){
			if (etDisplay == -1){
				glClearColor( 0.0, 0.0, 0.0, 0.0 );
				switch(atomNumber(i)->atom_name[0]){
					case 'N' : glColor3f( .4, .4, .6); break;
					case 'C' : glColor3f( .5, .5, .5); break;
					case 'O' : glColor3f( .8, .0, .0); break;
					case 'S' : glColor3f( .6, .5, .2); break;
					case 'H' : glColor3f( 1, 1, 1); break;
				}
			}
			else{
				sig = atomNumber(i)->significance;
				if (rendermode == 0){			////render the spectrum of ET vals simultaneously
					glClearColor( .5, .5, .5, 0.0 );
					glColor3f( 1.0 - .8*((float) sig/(float)maxSeq), 0, 0); ////Find out appropriate color
				}
				else {							///if not, render them independantly
					glClearColor( 0.0, 0.0, 0.0, 0.0 );
					if (rendermode <= sig){ glColor3f( .8, 0, 0); }
					if (rendermode > sig) { glColor3f( .7, .7, .7); }
				}
			}
		}

		glPushMatrix();

		if(pickMode == GL_SELECT){
			glPushName(i);
		}

		switch(atomNumber(i)->atom_name[0]){
			case 'N' : rad = 1.6; horiz = ATOMHORISUBS-2; vert = ATOMVERTSUBS+1; break;
			case 'C' : rad = 2.0; horiz = ATOMHORISUBS-1; vert = ATOMVERTSUBS+2; break;
			case 'O' : rad = 1.5; horiz = ATOMHORISUBS-2; vert = ATOMVERTSUBS+1; break;
			case 'S' : rad = 2.0; horiz = ATOMHORISUBS-2; vert = ATOMVERTSUBS; break;
			case 'H' : rad = 1.2; horiz = ATOMHORISUBS; vert = ATOMVERTSUBS; break;
			default : rad = ATOMRADIUS; break;
		}

		glTranslatef( (float) xpos, (float) ypos, (float) zpos);
		glutSolidSphere(rad, vert, horiz);

		if(pickMode == GL_SELECT){
			glPopName();
		}

		///Render Highlights if applicable for this atom
		if(pickMode == GL_RENDER && highlights[i] == true){
//			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
				glColor3f( 0, 1, 0);
				glutWireSphere( (rad+.3), HIGHLIGHTEDATOMVERTSUBS, HIGHLIGHTEDATOMHORISUBS);
//			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}

		glPopMatrix();
	}

	delete(centroid);
	delete(pts);

	if(pickMode == GL_RENDER && renderSelectBox){
		glColor3f( 0, 1, 0);
		glPushMatrix();
		glTranslatef(selectx, selecty, selectz);
		glutWireCube(radius);
		glPopMatrix();
	}

}


///////////////////
///Processes Selection Hits
///
void AtomList::processHits(GLint hits, GLuint buffer[], bool aminoAcidsOrAtoms, float * rot, float * pos){
	
	float closestZ = -HUGE_VAL;
	int closest = -1;
	unsigned int i, j;
	GLuint names, *ptr;
	ptr = (GLuint *) buffer;

	int * hitList = new int[hits];
	double * tempVec = new double[3];
	double * tempVec2;
	double bestZ = -HUGE_VAL;

	////FIXME: REMOVE THIS IN FUTURE VERSIONS// - this should not be calculated on line
	double * pts = new double[3*size];
	for (i = 0; (int)i<size; i++){
		pts[3*i+0] = atomList[i]->x;
		pts[3*i+1] = atomList[i]->y;
		pts[3*i+2] = atomList[i]->z;
	}
	////FIXME: REMOVE THIS IN FUTURE VERSIONS// - this should not be calculated on line
	double * centroid = getCentroid(pts, size);
	///////////////////////////////////////////

	double * rotation = new double[9];
	rotation[0] = (double) rot[0];	rotation[1] = (double) rot[1];	rotation[2] = (double) rot[2];
	rotation[3] = (double) rot[4];	rotation[4] = (double) rot[5];	rotation[5] = (double) rot[6];
	rotation[6] = (double) rot[8];	rotation[7] = (double) rot[9];	rotation[8] = (double) rot[10];

	double * translation = new double[3];
	translation[0] = ((double) pos[0]) - centroid[0];
	translation[1] = ((double) pos[1]) - centroid[1];
	translation[2] = ((double) pos[2]) - centroid[2] + ZDISPLACEMENT ;

	///Record the list of atoms which were hit
	for (i = 0; (int) i < hits; i++) { /*  for each hit  */
		names = *ptr;
		ptr++;		///skip past names
		ptr++;		///skip past first threshold
		ptr++;		///skip past second threshold

		///DEFENSIVE
		if (names != 1){printf("Un-named Object detected in GL SELECTION");}
		///DEFENSIVE

		////Store the records in HitList
		for (j = 0; j < names; j++) { /*  for each name */
			hitList[i] = *ptr;
			ptr++;
		}

////////////////////////////////////////////////////////////////////DEBUG
//		printf("%i:\n ", closest);
//		printf("%s\n---\n", atomList[closest]->fileLine);
//		highlightAA(closest);
////////////////////////////////////////////////////////////////////DEBUG
	}
////////////////////////////////////////////////////////////////////DEBUG
//	printf("\nhitList: ");
//
//	for (i = 0; (int)i<hits; i++){
//		printf("%i ", hitList[i]);
//	}
//	printf("\n");
////////////////////////////////////////////////////////////////////DEBUG

	////test each member of the hitList to find the closest one
	for (i = 0; (int)i<hits; i++){
		tempVec[0] = translation[0] + atomList[hitList[i]]->x;
		tempVec[1] = translation[1] + atomList[hitList[i]]->y;
		tempVec[2] = translation[2] + atomList[hitList[i]]->z;
		tempVec2 = transformVector(rotation, tempVec);

		if (bestZ < tempVec2[2]){
			bestZ = tempVec2[2];
			closest = hitList[i];
		}
		delete[](tempVec2);
	}

	if (hits > 0 ){
		if(!aminoAcidsOrAtoms){
			if(highlights[closest] == true){highlights[closest] = false;}
			else if(highlights[closest] == false){highlights[closest] = true;}
		}
		else if (closest != -1){ highlightAA(closest); }
	}
	
	delete[](hitList);
	delete[](rotation);
	delete[](translation);
	delete[](tempVec);
//	delete[](tempVec2);		//this was already deleted
	delete[](centroid);
	delete[](pts);

}


//////////////////////////////////////////
#endif //#ifdef INTERACTIVEMODE
//////////////////////////////////////////


