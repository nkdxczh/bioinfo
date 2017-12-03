// caList.cpp: implementation of the caList class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

caList::caList()
{
	numAcids = 0;
	name = new char[NUMCOLUMNS];
	bag = new atomBag();
	
	maxSeq = 1;
	highlights = NULL;
	acids = NULL;

}

///makes a caList out of an atomList
caList::caList(AtomList * aList)
{
	numAcids = 0;
	maxSeq = 1;
	highlights = NULL;

	int i = 0;
//	int j = 0;
//	int k = 0;
	int currentResNumber = 0;	///encodes the current res number
	int lastResNumber = -1;		///encodes the resNumber of the prev atom
	PdbAtom * tempAtom;
	aminoAcid * tempAcid = NULL;

	bag = new atomBag();
	
	for(i = 0; i<aList->size; i++){
		tempAtom = aList->atomList[i]->copy();

		currentResNumber = tempAtom->residue_number;

		if(!aList->highlights[i] ){
			bag->addAtom(tempAtom);

			if(currentResNumber != lastResNumber){
				lastResNumber = currentResNumber;
			}
			continue;
		}

		if(currentResNumber != lastResNumber){
			tempAcid = new aminoAcid();
			if(tempAtom->assocAA != NULL){
				strcpy(tempAcid->assocArray, tempAtom->assocAA);
			}
			else{
				delete[](tempAcid->assocArray);
				tempAcid->assocArray = NULL;		
			}
			addAcid(tempAcid);						//////ok
		}
	
		tempAcid->addAtom(tempAtom);

		lastResNumber = currentResNumber;
	}

	///set the globals
	name = new char[NUMCOLUMNS];
	strcpy(name, aList->name);

}

caList * caList::subSet(int size, int * vals)
{
	caList * result = new caList();
	aminoAcid * tempAcid;

	int i = 0;
	for(i = 0; i<size; i++){
		tempAcid = acids[vals[i]]->copy();
		result->addAcid(tempAcid);
	}

	return result;
}



caList::~caList()
{
	int i = 0;
	for(i = 0; i<numAcids; i++){
		delete(acids[i]);
	}
	if(acids != NULL){
		delete[](acids);
	}

	delete(bag);
	delete[](name);
	if(highlights != NULL){
		delete[](highlights);
	}

}

void caList::addAcid(aminoAcid * a)
{
	int i = 0;
	aminoAcid * *tempAcids = new aminoAcid*[numAcids+1];
	for(i = 0; i<numAcids; i++){
		tempAcids[i] = acids[i];
	}

	tempAcids[numAcids] = a;

	if (numAcids > 0){
		delete[](acids);
	}
	acids = tempAcids;

	numAcids++;
}

void caList::setName(char * n){
	strcpy(name, n);
}

void caList::calcMaxSeq()
{
	int i = 0;
	maxSeq = 0;
	for (i = 0; i<numAcids; i++){
		if (maxSeq < acids[i]->significance){
			maxSeq = acids[i]->significance;
		}
	}
	for (i = 0; i<bag->numAtoms; i++){
		if (maxSeq < bag->atoms[i]->significance){
			maxSeq = bag->atoms[i]->significance;
		}
	}

	///Instantiate the Highlight Bin now that the caList is complete
	highlights = new bool[numAcids];
	for (i = 0; i<numAcids; i++){
	  highlights[i] = false;

	}

}

void caList::setSourceHighlightBin(int num, int * vals)
{
	int i = 0;
	///clear the bin
	for(i = 0; i<numAcids; i++){
		highlights[i] = false;
	}

	///set the vals
	for(i = 0; i<num; i++){
		highlights[vals[2*i]] = true;
		//printf("%i\n", vals[2*i]);
	}

}

void caList::setTargetHighlightBin(int num, int * vals)
{
	int i = 0;
	///clear the bin
	for(i = 0; i<numAcids; i++){
		highlights[i] = false;
	}

	///set the vals
	for(i = 0; i<num; i++){
		highlights[vals[2*i+1]] = true;
		//printf("%i\n", vals[2*i+1]);
	}
//	for(i = 0; i<numAcids; i++){
//		printf("highlights[%i] = %i\n", i, highlights[i]);
//	}

}


void caList::clearHighlightBin()
{
	int i = 0;

	///clear the bin
	for(i = 0; i<numAcids; i++){
		highlights[i] = false;
	}

}

void caList::setHighlightsForTarget(){
	int i = 0;
	highlights = new bool[numAcids];
	for(i = 0; i<numAcids; i++){
		highlights[i] = false;
	}

}

void caList::toString()
{
  int i = 0;
  for (i = 0; i<bag->numAtoms; i++){
    PdbAtom * tempAtom = bag->atoms[i];
    tempAtom->toString();

  }
}




/*

///////////////////////////
///Switches for render Sources:
///	Background yes/no (bool background)
///	Background colored/not colored (bool colored)
/// render ETlevel (bool colormode, int et)
///	Atoms or Motifs (bool atoms)
/// rendering matching atoms only (matches)
///	render no hydrogens (nohydro)
///
///Rendering is done by first rendering the c-alpha, then the sidechain atoms
///Then by rendering all the unselected atoms
void caList::renderMe(bool background, bool colored, bool atoms, bool colormode, bool searchArea, int et, bool matches, bool nohydro)
{
	double xpos = 0.0;
	double ypos = 0.0;
	double zpos = 0.0;
	int i = 0;
	int j = 0;
	double * pts = new double[3*numAcids];
	float rad;
	int horiz, vert, sig = 0;
	aminoAcid * tempAcid;
	PdbAtom * tempAtom;
	double * centroid;

	////////////////////////////////////////////////////
	///Begin by rendering the Selected Amino Acids//////
	////////////////////////////////////////////////////

	///We operate on the center of mass of the selected C-Alphas
	for (i = 0; i<numAcids; i++){
		pts[3*i+0] = acids[i]->Calpha->x;
		pts[3*i+1] = acids[i]->Calpha->y;
		pts[3*i+2] = acids[i]->Calpha->z;
	}
	centroid = getCentroid(pts, numAcids);
	////////////////////////////////////////////////////////////
	glClearColor( .5, .5, .5, 0.0 );

	///Render the C-alpha
	for	(i = 0; i<numAcids; i++){
		tempAcid = acids[i];
		tempAtom = tempAcid->Calpha;
		if (!highlights[i] && matches){continue;}

		xpos = pts[3*i+0] - centroid[0];
		ypos = pts[3*i+1] - centroid[1];
		zpos = pts[3*i+2] - centroid[2];
		if (colormode){
			if (!colored){
				switch(tempAtom->atom_name[0]){		///These are all for CA's, so none should matter except the "C" case
					case 'N' : glColor3f( .4, .4, .6); break;
					case 'C' : glColor3f( .5, .5, .5); break;
					case 'O' : glColor3f( .8, .0, .0); break;
					case 'S' : glColor3f( .6, .5, .2); break;
					case 'H' : glColor3f( 1, 1, 1); break;
				}
			}
			else{glColor3f( .7, 0, 0);}///otherwise color it red
		}
		else{
			if (!colored){
				sig = tempAtom->significance;
				if (et == 0){			////render the spectrum of ET vals simultaneously
					glColor3f( 1.0 - .8*((float) sig/(float)maxSeq), 0, 0); ////Find out appropriate color
				}
				else {							///if not, render them independantly
					if (et < sig){ glColor3f( .8, 0, 0); }
					if (et >= sig) { glColor3f( .7, .7, .7); }
				}
			}
			else{glColor3f( .7, 0, 0);}///otherwise color it red
		}
		glPushMatrix();

		switch(tempAtom->atom_name[0]){				///These are all for CA's, so none should matter except the "C" case
			case 'N' : rad = 1.6; horiz = ATOMHORISUBS-2; vert = ATOMVERTSUBS+1; break;
			case 'C' : rad = 2.0; horiz = ATOMHORISUBS-1; vert = ATOMVERTSUBS+2; break;
			case 'O' : rad = 1.5; horiz = ATOMHORISUBS-2; vert = ATOMVERTSUBS+1; break;
			case 'S' : rad = 2.0; horiz = ATOMHORISUBS-2; vert = ATOMVERTSUBS; break;
			case 'H' : rad = 1.2; horiz = ATOMHORISUBS; vert = ATOMVERTSUBS; break;
			default : rad = ATOMRADIUS; break;
		}

		if(searchArea){glColor3f( .6, .6, 0);}

		glTranslatef( (float) xpos, (float) ypos, (float) zpos);
		glutSolidSphere(rad, vert, horiz);

		if(highlights[i]){
			glColor3f( 0, 1, 0);
			glutWireSphere( (rad+.3), HIGHLIGHTEDATOMVERTSUBS, HIGHLIGHTEDATOMHORISUBS);
		}

		glPopMatrix();

		///////////////////////////
		///Render the rest/////////
		///////////////////////////
		if(atoms){
			for (j = 0; j<tempAcid->numAtoms; j++){

				tempAtom = tempAcid->atoms[j];
				if (nohydro && (tempAtom->atom_name[0] == 'H')){continue;}

				///set the temp vals from the atomList and render a sphere in that location
				xpos = tempAtom->x - centroid[0];
				ypos = tempAtom->y - centroid[1];
				zpos = tempAtom->z - centroid[2];

				if (colormode){
					if (!colored){
						switch(tempAtom->atom_name[0]){
							case 'N' : glColor3f( .4, .4, .6); break;
							case 'C' : glColor3f( .5, .5, .5); break;
							case 'O' : glColor3f( .8, .0, .0); break;
							case 'S' : glColor3f( .6, .5, .2); break;
							case 'H' : glColor3f( 1, 1, 1); break;
						}
					}
					else{glColor3f( .7, 0, 0);}///otherwise color it grey
				}
				else{
					if (!colored){
						sig = tempAtom->significance;
						if (et == 0){			////render the spectrum of ET vals simultaneously
							glColor3f( 1.0 - .8*((float) sig/(float)maxSeq), 0, 0); ////Find out appropriate color
						}
						else {							///if not, render them independantly
							if (et < sig){ glColor3f( .8, 0, 0); }
							if (et >= sig) { glColor3f( .7, .7, .7); }
						}
					}
					else{glColor3f( .7, 0, 0);}///otherwise color it red
				}
				glPushMatrix();

				switch(tempAtom->atom_name[0]){
					case 'N' : rad = 1.6; horiz = ATOMHORISUBS-2; vert = ATOMVERTSUBS+1; break;
					case 'C' : rad = 2.0; horiz = ATOMHORISUBS-1; vert = ATOMVERTSUBS+2; break;
					case 'O' : rad = 1.5; horiz = ATOMHORISUBS-2; vert = ATOMVERTSUBS+1; break;
					case 'S' : rad = 2.0; horiz = ATOMHORISUBS-2; vert = ATOMVERTSUBS; break;
					case 'H' : rad = 1.2; horiz = ATOMHORISUBS; vert = ATOMVERTSUBS; break;
					default : rad = ATOMRADIUS; break;
				}

				if(searchArea){glColor3f( .6, .6, 0);}

				glTranslatef( (float) xpos, (float) ypos, (float) zpos);
				glutSolidSphere(rad, vert, horiz);

				if(highlights[i]){
					glColor3f( 0, 1, 0);
					glutWireSphere( (rad+.3), HIGHLIGHTEDATOMVERTSUBS, HIGHLIGHTEDATOMHORISUBS);
				}

				glPopMatrix();

			}
		}
	}

	////////////////////////////////////////////////////
	////Now Render the Extraneous unselected atoms//////
	////////////////////////////////////////////////////
	if (background && atoms && !matches){
		for (i = 0; i<bag->numAtoms; i++){
			tempAtom = bag->atoms[i];
			if (nohydro && (tempAtom->atom_name[0] == 'H')){continue;}

			///set the temp vals from the atomList and render a sphere in that location
			xpos = tempAtom->x - centroid[0];
			ypos = tempAtom->y - centroid[1];
			zpos = tempAtom->z - centroid[2];

			if (colormode){
				if (!colored){
					switch(tempAtom->atom_name[0]){
						case 'N' : glColor3f( .4, .4, .6); break;
						case 'C' : glColor3f( .5, .5, .5); break;
						case 'O' : glColor3f( .8, .0, .0); break;
						case 'S' : glColor3f( .6, .5, .2); break;
						case 'H' : glColor3f( 1, 1, 1); break;
					}
				}
				else{glColor3f( .7, 0, 0);}///otherwise color it grey
			}
			else{
				if (!colored){
					sig = tempAtom->significance;
					if (et == 0){			////render the spectrum of ET vals simultaneously
						glColor3f( 1.0 - .8*((float) sig/(float)maxSeq), 0, 0); ////Find out appropriate color
					}
					else {							///if not, render them independantly
						if (et < sig){ glColor3f( .8, 0, 0); }
						if (et >= sig) { glColor3f( .7, .7, .7); }
					}
				}
				else{glColor3f( .7, 0, 0);}///otherwise color it red
			}
			glPushMatrix();

			switch(tempAtom->atom_name[0]){
				case 'N' : rad = 1.6; horiz = ATOMHORISUBS-2; vert = ATOMVERTSUBS+1; break;
				case 'C' : rad = 2.0; horiz = ATOMHORISUBS-1; vert = ATOMVERTSUBS+2; break;
				case 'O' : rad = 1.5; horiz = ATOMHORISUBS-2; vert = ATOMVERTSUBS+1; break;
				case 'S' : rad = 2.0; horiz = ATOMHORISUBS-2; vert = ATOMVERTSUBS; break;
				case 'H' : rad = 1.2; horiz = ATOMHORISUBS; vert = ATOMVERTSUBS; break;
				default : rad = ATOMRADIUS; break;
			}

			glTranslatef( (float) xpos, (float) ypos, (float) zpos);
			glutSolidSphere(rad, vert, horiz);

			glPopMatrix();
		}
	}

	delete(centroid);
	delete(pts);
}

*/

