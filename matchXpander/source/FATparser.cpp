// FATparser.cpp: implementation of the FATparser class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FATparser::FATparser(char * filename)
{
	numSources = 0;
	sources = NULL;
	target = NULL;
	sourceFileName = new char[NUMCOLUMNS];
	targetFileName = new char[NUMCOLUMNS];

	name = new char[NUMCOLUMNS];
	strcpy(name, filename);

	FATfile = fopen(filename, "r");
	AAMode = true;

}

FATparser::~FATparser()
{
	delete[](sourceFileName);
	delete[](targetFileName);
	delete[](name);
	int i = 0;
	for(i = 0; i<numSources; i++){
		delete(sources[i]);
	}
	delete[](sources);
	delete(target);
}


bool FATparser::readFile()
{
	if (FATfile == NULL){ return false; }

	char * line = new char[NUMCOLUMNS+1];
	int i = 0;
	////read the global information
	if (readGlobals() == false){
		///waste the blank line
		fgets(line, NUMCOLUMNS, FATfile);
	}

	sources = new AtomList*[numSources];
	for (i = 0; i<numSources; i++){
		sources[i] = readSource();
	}
	target = readTarget();

	delete[](line);
	return true;
}


bool FATparser::readGlobals(){
	char *line = new char[NUMCOLUMNS];
	char *lineBuffer;
	char *temp;
	bool result = false;
	char *temp2;

////skip past the # lines

	fgets(line, NUMCOLUMNS, FATfile);
	temp2 = firstChar(line);
	while( strcmp( temp2, "#") == 0 ){
		fgets(line, NUMCOLUMNS, FATfile);
		delete[](temp2);
		temp2 = firstChar(line);
	}
	delete[](temp2);

	lineBuffer = &line[0];					/////mirror line as linebuffer

////Get the number of sources

	if( strcmp(strtok(line, " "), "NUMSOURCES=") == 0){
		temp = stringCopy(lineBuffer, 12,17);
		numSources = atoi(strtok(temp, " \r\n"));
		delete[](temp);
	}

////Get the number of arrays

	fgets(line, NUMCOLUMNS, FATfile);
	lineBuffer = &line[0];					/////mirror line as linebuffer

	if( strcmp(strtok(lineBuffer, " "), "NUMARRAYS=") == 0){
		//Currently ignore the NUMARRAYS info
	}

////Find out if we are using Amino ACids or Atoms
	fgets(line, NUMCOLUMNS, FATfile);
	lineBuffer = &line[0];					/////mirror line as linebuffer

	if( strcmp(strtok(lineBuffer, " "), "AMINOACIDS=") == 0){
		temp = strtok(stringCopy(lineBuffer, 12, 18), " \r\n");
		if (strcmp(temp, "true") == 0){
			AAMode = true;
			printf("PROCEEDING WITH AMINO ACID (C-alpha) EXPERIMENT FILE\n");
		}
		else if (strcmp(temp, "false") == 0){
			AAMode = false;
			printf("PROCEEDING WITH ATOM-LEVEL EXPERIMENT FILE\n");
		}
		else{
			printf("PARSE ERROR: INDETERMINANT TOKEN IN AMINOACIDS= LINE\n"); 
		}
		delete[](temp);
	}
	else{

#ifndef NO_OUTPUT//////////////////////////////////////////////////////////////////////////////////////////////////
			printf("LEGACY MODE: AMINOACIDS= NOT SPECIFIED. PROCEEDING USING C-Alpha MOTIFS\n"); 
#endif/////////////////////////////////////////////////////////////////////////////////////////////////////////////
			AAMode = true;
			result = true;
	}

	delete[](line);

	return result;
}


AtomList * FATparser::readSource(){
	char *line = new char[NUMCOLUMNS];
	char *lineBuffer = new char[NUMCOLUMNS];
	char *tempString = new char[NUMCOLUMNS];
	int numAtoms = 0;
	int i = 0;
	PdbAtom * tempAtom;
	int maxSequence = 0;


	////Detect the Source Header Info
	fgets(line, NUMCOLUMNS, FATfile);
	strcpy(lineBuffer, line);
	if( strcmp(strtok(lineBuffer, " "), "SOURCEPATTERN") == 0 ){
		///process the SOURCE header info
		strtok(NULL, " ");
		strcpy(	sourceFileName, strtok(NULL, " \r\n") );
		numAtoms = atoi(strtok(NULL, " \r\n"));
		maxSequence = atoi(strtok(NULL, " \r\n"));
	}
	else{
		printf("ERROR: %s has bad spacing in it's SOURCE definition", name);
	}

	delete[](line); line = NULL;
	delete[](lineBuffer); lineBuffer = NULL;
	line = new char[numAtoms*NUMCOLUMNS];
	lineBuffer = new char[numAtoms*NUMCOLUMNS];

#ifndef NO_OUTPUT//////////////////////////////////////////////////////////////////////////////////////////////////
	printf("Source File name :              %s\n", sourceFileName);
	printf("Num Atoms:                      %i\n", numAtoms);
#endif ////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///Generate the AtomList to fill
	AtomList * result = new AtomList(numAtoms);
	result->maxSeq = maxSequence;
	strcpy(result->name, sourceFileName);
	///////////////////////////

	fgets(line, NUMCOLUMNS, FATfile);
	strcpy(lineBuffer, line);
	for (i = 0; i<numAtoms; i++){
		if ( strcmp( strtok(lineBuffer, " \r\n"), "SOURCEATOM") == 0 ){
			tempAtom = result->atomNumber(i);

			////process the Atom header info
			tempAtom->coords[0] = atof( strtok(NULL, " \r\n"));
			tempAtom->coords[1] = atof( strtok(NULL, " \r\n"));
			tempAtom->coords[2] = atof( strtok(NULL, " \r\n"));
			tempAtom->residue_number = atoi( strtok(NULL, " \r\n"));
			tempAtom->significance = atoi( strtok(NULL, " \r\n"));
			strcpy( tempAtom->atom_name, strtok(NULL, " \r\n") );
			strcpy( tempAtom->residue_name, strtok(NULL, " \r\n") );

			strcpy(tempString, strtok(NULL, " \r\n"));
			tempAtom->assocAA = new char[strlen(tempString)+1];
			strcpy(tempAtom->assocAA, tempString);
			tempAtom->atom_number = atoi(strtok(NULL, " \r\n"));
//			printf("TESTING: %s\n", tempAtom->assocAA);

			////Load the fileLine for this amino acid, to be used for PDB output
			fgets(line, numAtoms*NUMCOLUMNS, FATfile);
			strcpy(tempAtom->fileLine, line);

			////Load the next line to be read
			fgets(line, numAtoms*NUMCOLUMNS, FATfile);
			strcpy(lineBuffer, line);
		}
		else{
			printf("%s has bad formatting in Source Description\n", name);
		}
	}

	char * linePointer = lineBuffer;

	if ( strcmp( strtok(lineBuffer, " \r\n"), "MOTIF:") == 0 ){
		lineBuffer = strtok(NULL, " \r\n");
		while (lineBuffer != NULL){
			result->highlights[atoi(lineBuffer)] = true;
			lineBuffer = strtok(NULL, " \r\n");
		}
	}

	set_t stuff = alloc_set(0);
	for(i= 0; i<result->size; i++){
		if(result->highlights[i]){
			//printf("%i %i\n", i, result->atomList[i]->residue_number);
			stuff = put_set(stuff, result->atomList[i]->residue_number); 
		}
	}

	delete[](tempString);

	///waste the last line
	fgets(line, numAtoms*NUMCOLUMNS, FATfile);

	delete[](line);
	delete[](linePointer);
	free_set(stuff);
//	delete[](lineBuffer);

	return result;
}


AtomList * FATparser::readTarget(){
	char *line = new char[NUMCOLUMNS];
	char *lineBuffer = new char[NUMCOLUMNS];
	int numAtoms = 0;
	int i = 0;
	PdbAtom * tempAtom;


	////Detect the Source Header Info
	fgets(line, NUMCOLUMNS, FATfile);
	strcpy(lineBuffer, line);
	if( strcmp(strtok(lineBuffer, " "), "TARGETPATTERN") == 0 ){
		///process the TARGET header info
		strcpy(	targetFileName, strtok(NULL, " \r\n") );
		numAtoms = atoi(strtok(NULL, " \r\n"));
	}
	else{
		printf("ERROR: %s has bad spacing in it's TARGET definition\n", name);
		printf("%s", line);
	}

	delete[](line); line = NULL;
	delete[](lineBuffer); lineBuffer = NULL;
	line = new char[numAtoms*NUMCOLUMNS];
	lineBuffer = new char[numAtoms*NUMCOLUMNS];

#ifndef NO_OUTPUT//////////////////////////////////////////////////////////////////////////////////////////////////
	printf("Target File name :              %s\n", targetFileName);
	printf("Num Atoms:                      %i\n", numAtoms);
#endif ////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///Generate the AtomList to fill
	AtomList * result = new AtomList(numAtoms);
	strcpy(result->name, targetFileName);
	///////////////////////////

	fgets(line, numAtoms*NUMCOLUMNS, FATfile);

	strcpy(lineBuffer, line);
	for (i = 0; i<numAtoms; i++){
		if ( strcmp( strtok(lineBuffer, " \r\n"), "TARGETATOM") == 0 ){
			tempAtom = result->atomNumber(i);

			////process the Atom header info
			tempAtom->coords[0] = atof( strtok(NULL, " \r\n"));
			tempAtom->coords[1] = atof( strtok(NULL, " \r\n"));
			tempAtom->coords[2] = atof( strtok(NULL, " \r\n"));
			tempAtom->residue_number = atoi( strtok(NULL, " \r\n"));
			tempAtom->significance = atoi( strtok(NULL, " \r\n"));
			strcpy( tempAtom->atom_name, strtok(NULL, " \r\n") );
			strcpy( tempAtom->residue_name, strtok(NULL, " \r\n") );
			tempAtom->atom_number = atoi(strtok(NULL, " \r\n"));

//			printf("r#: %i residue_name: %s\n", tempAtom->residue_number, tempAtom->residue_name);

			////Load the fileLine for this amino acid, to be used for PDB output
			fgets(line, numAtoms*NUMCOLUMNS, FATfile);
			strcpy(tempAtom->fileLine, line);

			////Load the next line to be read
			fgets(line, numAtoms*NUMCOLUMNS, FATfile);
			strcpy(lineBuffer, line);
		}
		else{
			printf("%s has bad formatting in Target Description\n", name);
		}
	}

	char * linePointer = lineBuffer;

	if ( strcmp( strtok(lineBuffer, " \r\n"), "MOTIF:") == 0 ){
		lineBuffer = strtok(NULL, " \r\n");
		while (lineBuffer != NULL){
			result->highlights[atoi(lineBuffer)] = true;
			lineBuffer = strtok(NULL, " \r\n");
		}
	}

	delete[](line);
	delete[](linePointer);
//	delete[](lineBuffer);

	return result;
}




