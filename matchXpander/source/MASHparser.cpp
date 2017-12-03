

#include "MASHparser.h"


///this determines if the file is in the MASH format by testing the first line.
bool testMASHformat(char * fileName)
{
	char * line = new char[NUMCOLUMNS];
	char * ptr;
	bool test1 = false;
	bool test2 = false;
	bool result = false;

	///test file	
	FILE * file = fopen(fileName, "r");
	if(file == NULL){ printf("ERROR: Filename [%s] Not Found.  Exitting\n", fileName); exit(1); } 

	///parse first line of file
//	printf("LINE: [%s]\n", line);
	fgets(line, NUMCOLUMNS, file);
	ptr = strtok(line, " \r\t\n");
//	printf("PTR: [%s]\n", ptr);
	if(strcmp(ptr, "MASH")==0){ test1 = true; }
	ptr = strtok(NULL, " \r\t\n");
//	printf("PTR: [%s]\n", ptr);
	if(strcmp(ptr, "format")==0){ test2 = true; }

	///test results
	if(test1 & test2){ result = true; }

	fclose(file);
	
	return result;
}



///This parses the MASH file for the points, and returns an atomList
AtomList * parseMASHfilePoints(char * fileName)
{
	char * line = new char[NUMCOLUMNS];
	char * lineCopy = new char[NUMCOLUMNS];
	char *lineBuffer2;
	char * ptr;
	int i = 0;

	///test file	
	FILE * file = fopen(fileName, "r");
	if(file == NULL){ printf("ERROR: Filename [%s] Not Found.  Exitting\n", fileName); exit(1); } 

	///parse the name of PDB file
	fgets(line, NUMCOLUMNS, file);
	fgets(line, NUMCOLUMNS, file);
	ptr = strtok(line, " []\n\r\t");
	ptr = strtok(NULL, " []\n\r\t");
	ptr = strtok(NULL, " []\n\r\t");
	char * thisName = new char[80];
	strcpy(thisName, ptr);	

	///parse until we get to the points
	fgets(line, NUMCOLUMNS, file);
	strcpy(lineCopy, line);
	ptr = strtok(lineCopy, " \n\r\t");
	while( strcmp(ptr, "POINTNUM") != 0 ){
		fgets(line, NUMCOLUMNS, file);
		strcpy(lineCopy, line);
		ptr = strtok(lineCopy, " \n\r\t");
	}
	
	///get the number of atoms
	ptr = strtok(NULL, " \n\r\t");
	int numAtoms = atoi(ptr);
	
	///create the result
	AtomList * result = new AtomList(numAtoms);
	strcpy(result->name, thisName);
	
	///now parse and fill
	for(i = 0; i<numAtoms; i++){
		///read in the ATOM card from the PDB format
		fgets(line, NUMCOLUMNS, file);
		lineBuffer2 = &line[0];					/////mirror line as linebuffer2
		char * lineBuffer3 = new char[NUMCOLUMNS];
		strcpy(lineBuffer3, line);
		char * tempString;

		///parse Atom Number
		tempString = stringCopy(lineBuffer2, 6,10);
		result->atomNumber(i)->atom_number = atoi(strtok(tempString, " \t\r\n")); 		delete[](tempString);
		///parse Atom Name
		tempString = stringCopy(lineBuffer2,12,15);
		strcpy(result->atomNumber(i)->atom_name, strtok(tempString, " \t\r\n") ); 		delete[](tempString);
		///parse Residue Name
		tempString = stringCopy(lineBuffer2,17,19);
		strcpy(result->atomNumber(i)->residue_name, strtok(tempString, " \t\r\n") ); 	delete[](tempString);
		///parse Residue Number
		tempString = stringCopy(lineBuffer2,22,25);
		result->atomNumber(i)->residue_number = atoi(strtok(tempString, " \t\r\n")); 	delete[](tempString);
		///parse X pos
		tempString = stringCopy(lineBuffer2,30,37);
		result->atomNumber(i)->coords[0] = atof(strtok(tempString, " \t\r\n")); 		delete[](tempString);
		///parse y pos
		tempString = stringCopy(lineBuffer2,38,45);
		result->atomNumber(i)->coords[1] = atof(strtok(tempString, " \t\r\n")); 		delete[](tempString);
		///parse z pos
		tempString = stringCopy(lineBuffer2,46,53);
		result->atomNumber(i)->coords[2] = atof(strtok(tempString, " \t\r\n")); 		delete[](tempString);
		///parse file line
		strcpy(result->atomNumber(i)->fileLine, lineBuffer3);
		delete[](lineBuffer3);

		///Now read the Evolutionary significance and alternative amino acids
		fgets(line, NUMCOLUMNS, file);
		lineBuffer2 = strtok(line, " \r\n\t");
		lineBuffer2 = strtok(NULL, " \r\n\t");
		result->atomNumber(i)->significance = atoi(lineBuffer2);

		result->atomNumber(i)->assocAA = new char[21];
		lineBuffer2 = strtok(NULL, " \r\n\t");
		strcpy(result->atomNumber(i)->assocAA, lineBuffer2);
		
		///now highlight the atom
		result->highlights[i] = true;
	}
	delete[](line);
	delete[](thisName);

	return result;
}
