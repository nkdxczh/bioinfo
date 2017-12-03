// pdbParser.cpp: implementation of the pdbParser class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////



pdbParser::pdbParser(char *filename)
{

	pdbFile = fopen(filename, "r");					//Open File
	//////////////////////////////////////////////////////////////////
	//printf("Reading PDB File....\n");///////////////////////////////
	//////////////////////////////////////////////////////////////////
	if(pdbFile == NULL) {
		fprintf(stderr, "ERROR: Unable to Find PDB File: %s\n", filename);
		exit(1);
	}

	int atomNumber = atomCount();					//Count the # of atoms so we can size the AtomList

	atomList = new AtomList(atomNumber);		//Allocate the atomlist
//	setAtomList(alist);								//set the list
	readPdbFile();									//get all coords/info from Pdb file
//	alist->recenterAtoms();							//recenter the atoms to Centroid; to make rotation intuitive

	strcpy(atomList->name, filename);
	fclose(pdbFile);
}

pdbParser::~pdbParser()
{

	delete(atomList);

}

/*
void pdbParser::setAtomList(AtomList *alist)
{

	alist;

}
*/

//////////////////////////
///Retursn a copy of it's AtomList.  (pass by value, baby)
////
AtomList * pdbParser::getAtomList(void)
{

	int s = atomList->size;
	int i = 0;

	AtomList * result = new AtomList(s);

	///check string lengths/space allocated./

	for (i = 0; i<s; i++){
		result->atomList[i]->atom_number	= atomList->atomNumber(i)->atom_number;

		strcpy(result->atomList[i]->atom_name, atomList->atomNumber(i)->atom_name);
		strcpy(result->atomList[i]->residue_name, atomList->atomNumber(i)->residue_name);

		result->atomList[i]->residue_number	= atomList->atomNumber(i)->residue_number;
		result->atomList[i]->coords[0]		= atomList->atomNumber(i)->coords[0];
		result->atomList[i]->coords[1]		= atomList->atomNumber(i)->coords[1];
		result->atomList[i]->coords[2]		= atomList->atomNumber(i)->coords[2];
		result->atomList[i]->significance	= atomList->atomNumber(i)->significance;
		strcpy(result->atomList[i]->fileLine, atomList->atomNumber(i)->fileLine);
	}

	strcpy(result->name, atomList->name);

	return result;

}


int pdbParser::atomCount(void){

	char *line = new char[NUMCOLUMNS];
	int listLength = 0;

	//////////////////////////////////////////////////////////////////
	//printf("Counting Atoms in PDB File....\n");/////////////////////
	//////////////////////////////////////////////////////////////////
	while( fgets(line, NUMCOLUMNS, pdbFile) != NULL ){
		if( line[0]=='A' && line[1]=='T' && line[2]=='O' && line[3]=='M' ){
			listLength++;
		}
	}

	//////////////////////////////////////////////////////////////////
	//printf("Evaluating PDB File Integrity....\n");//////////////////
	//////////////////////////////////////////////////////////////////
	if (listLength < 1){
		printf("ERROR : no atoms found in input file! \n");
	}

	rewind(pdbFile);

	delete[](line);

	return listLength;

}


//////////////////////////////////////////////////////////////////////////////
////////This function reads the file and stores the info//////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////Note: This Function goes out of date if they decide to change/////////////
//////////The PDB File Format.  IF they do, you must adjust the //////////////
//////////In the smallest If brackets to match the new file format setting////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///
///
///
void pdbParser::readPdbFile()
{
	
	char *line = new char[NUMCOLUMNS];
	int currentAtom = 0;
	char *lineBuffer;
	char *lineBuffer2;

	while(fgets(line, NUMCOLUMNS, pdbFile) != NULL){	////need to check for null lines and such to catch the case of a blank line

		lineBuffer = &line[0];					/////mirror line as linebuffer
		lineBuffer2 = &line[0];					/////mirror line as linebuffer2
		char * lineBuffer3 = new char[NUMCOLUMNS];
		strcpy(lineBuffer3, line);
		char * tempString;

		if( strcmp(strtok(lineBuffer, " "), "ATOM") == 0 ){	///only process a line with ATOM info

			tempString = stringCopy(lineBuffer2, 6,10);
			atomList->atomNumber(currentAtom)->atom_number	  = atoi(strtok(tempString, " "));
			delete[](tempString);

			tempString = stringCopy(lineBuffer2,12,15);
			strcpy(atomList->atomNumber(currentAtom)->atom_name, strtok(tempString, " ") );
			delete[](tempString);

			tempString = stringCopy(lineBuffer2,17,19);
			strcpy(atomList->atomNumber(currentAtom)->residue_name, strtok(tempString, " ") );
			delete[](tempString);

			tempString = stringCopy(lineBuffer2,22,25);
			atomList->atomNumber(currentAtom)->residue_number = atoi(strtok(tempString, " "));
			delete[](tempString);

			tempString = stringCopy(lineBuffer2,30,37);
			atomList->atomNumber(currentAtom)->coords[0]			  = atof(strtok(tempString, " "));
			delete[](tempString);
			
			tempString = stringCopy(lineBuffer2,38,45);
			atomList->atomNumber(currentAtom)->coords[1]			  = atof(strtok(tempString, " "));
			delete[](tempString);
			
			tempString = stringCopy(lineBuffer2,46,53);
			atomList->atomNumber(currentAtom)->coords[2]			  = atof(strtok(tempString, " "));
			delete[](tempString);

			strcpy(atomList->atomNumber(currentAtom)->fileLine, lineBuffer3);
			
			currentAtom++;
		}
		delete[](lineBuffer3);
	}

	//////////////////////////////////////////////////////////////////
	//printf("PDB File Parsing Complete. \n");//////////////////////////
	//////////////////////////////////////////////////////////////////

	delete[](line);

}



/*
char * pdbParser::stringCopy(char* str, int start, int end)
{
	int size = end - start + 1;

	char * result = new char[size + 1];

	for(int i = 0; i<size; i++){
		result[i] = str[start+i];	
	}
	result[size] = NULL;

	return result;

}
*/



