
#include "StdAfx.h"


caList * generateCaListFromBinary(char * binary, int startingByte)
{
	int i = 0;
	int j = 0;

	caList * result = new caList();

	char * thisProt = &binary[startingByte];
	int currentOffset = 0;
	int numAA = 0;

        ////THIS INDEX (aaCounter=1) ENSURES THAT BINARY PDB FILES HAVE THE SAME RESIDUE NUMBERS
        ////AS PURIFIED PDB FILES, WHICH START AA NUMBERING AT 1, BYC 07/07/07
	int aaCounter = 1;
        ////THIS INDEX (aaCounter=1) ENSURES THAT BINARY PDB FILES HAVE THE SAME RESIDUE NUMBERS
        ////AS PURIFIED PDB FILES, WHICH START AA NUMBERING AT 1, BYC 07/07/07

	int atomCounter = 0;

	////first read the Protein Header
	currentOffset = readProteinHeader(thisProt, currentOffset, &numAA, result->name);

	for(i = 0; i<numAA; i++){
		aminoAcid * newAcid = new aminoAcid();
		int numAtoms;
		int aaType = 0;
		int atomType = 0;

		////next write hte Amino Acid headers
		currentOffset = readAminoAcid(thisProt, currentOffset, &numAtoms, &aaType);
		char * resName = decodeAAcode(aaType);
	
		////read C-alpha entry
		PdbAtom * newAtom = new PdbAtom();
		strcpy(newAtom->fileLine, "Artificial File Line");

		currentOffset = readAtom(thisProt, currentOffset, &newAtom->coords[0], &newAtom->coords[1], &newAtom->coords[2], &atomType);

		char * atomName = encodeAtomCode(atomType);
		strcpy(newAtom->atom_name, atomName);
		delete[](atomName);
		strcpy(newAtom->residue_name, resName);
		newAtom->residue_number = aaCounter;
		newAtom->atom_number = atomCounter;
		atomCounter++;
		newAtom->significance = -1;

		newAcid->addAtom(newAtom);

		////read entries for each atom.
		for(j = 0; j<numAtoms-1; j++){
			newAtom = new PdbAtom();
			strcpy(newAtom->fileLine, "Artificial File Line");
			currentOffset = readAtom(thisProt, currentOffset, &newAtom->coords[0], &newAtom->coords[1], &newAtom->coords[2], &atomType);
			atomName = encodeAtomCode(atomType);
			strcpy(newAtom->atom_name, atomName);
			delete[](atomName);
			strcpy(newAtom->residue_name, resName);
			newAtom->residue_number = aaCounter;
			newAtom->atom_number = atomCounter;
			atomCounter++;
			newAtom->significance = -1;

			newAcid->addAtom(newAtom);
		}

		result->addAcid(newAcid);

		aaCounter++;
		delete[](resName);
	}

	return result;
}


AtomList * generateAtomListFromBinary(char * binary, int startingByte)
{
	caList * tempList = generateCaListFromBinary(binary, startingByte);

	int i = 0;
	int j = 0;
	int numAtoms = 0;
	int counter = 0;

	for(i = 0; i<tempList->numAcids; i++){
		numAtoms += tempList->acids[i]->numAtoms+1;
	}

	PdbAtom * *aList = new PdbAtom*[numAtoms];

	for(i = 0; i<tempList->numAcids; i++){
		PdbAtom * tempAtom;
		for(j = 0; j<tempList->acids[i]->numAtoms; j++){
			tempAtom = tempList->acids[i]->atoms[j]->copy();
			aList[counter] = tempAtom;
			counter++;
		}
	}

	AtomList * result = new AtomList(numAtoms, aList);

	for(i = 0; i<numAtoms; i++){
		delete(aList[i]);
	}
	delete[](aList);

	return result;
}


int * parseTOC(char * toc, int * binFileSize, int * nFiles)
{
	int i = 0;
	FILE * tocFile = fopen(toc, "r");
	char * line = new char[NUMCOLUMNS];
	char * line2 = new char[NUMCOLUMNS];
	
	/////Skip past the comments				
	fgets(line, NUMCOLUMNS, tocFile);			
	strcpy(line2, line);					
	while(line[0] == '#'){					
		fgets(line, NUMCOLUMNS, tocFile);		
		strcpy(line2, line);				
	}										

	///read the fileSize
	char * ptr = strtok(line2, " \n\t");
	ptr = strtok(NULL, " \n\t");
	int numFiles = atoi(ptr);

	///output the fileSize
	nFiles[0] = numFiles;

	///alloc table of contents
	int * tableOfContents = new int[numFiles];

	//read each file size
	for(i = 0; i<numFiles; i++){
		fgets(line, NUMCOLUMNS, tocFile);		
		strcpy(line2, line);				
		ptr = strtok(line2, " \n\t");
		tableOfContents[i] = atoi(ptr);
	}

	///read the overall filesize
	fgets(line, NUMCOLUMNS, tocFile);		
	strcpy(line2, line);				
	ptr = strtok(line2, " \n\t");
	ptr = strtok(NULL, " \n\t");

	int fileSize = atoi(ptr);
	binFileSize[0] = fileSize;

	fclose(tocFile);
	return tableOfContents;
}

char * loadBinFile(char * binFile, int fileSize)
{
	///now read the bitFile
	char * result = new char[fileSize];
	
	if(result == NULL){
		printf("ERROR: OS refuses to allocate %i bytes for binaryStorage! Quitting!\n", fileSize);
		exit(1);
	}

	FILE * binaryFile = fopen(binFile, "rb");
	fread(result, fileSize, 1, binaryFile);
	fclose(binaryFile);

	return result;
}


void parsePdbListAndOutputBinFile(char * list, char * outputFile, char * tableOfContentsName)
{
	FILE * pdbList = fopen(list, "r");
	char * line = new char[NUMCOLUMNS];
	char * line2 = new char[NUMCOLUMNS];

	/////Skip past the comments				//##	Parseing for listFile
	fgets(line, NUMCOLUMNS, pdbList);		//##	Parseing for listFile
	strcpy(line2, line);					//##	Parseing for listFile
	while(line[0] == '#'){					//##	Parseing for listFile
		fgets(line, NUMCOLUMNS, pdbList);	//##	Parseing for listFile
		strcpy(line2, line);				//##	Parseing for listFile
	}										//##	Parseing for listFile
											//##	Parseing for listFile
	///now we read the size					//##	Parseing for listFile
	char * ptr = strtok(line2, " \t\n");	//##	Parseing for listFile
	ptr = strtok(NULL, " \t\n");			//##	Parseing for listFile

	int listSize = atoi(ptr);

	///output Table of Contents File											//$$$$$ Output for Table of Contents
	FILE * tocFile = fopen(tableOfContentsName, "w");						//$$$$$ Output for Table of Contents
	char * tempString = new char[NUMCOLUMNS];									//$$$$$ Output for Table of Contents
	sprintf(tempString, "## Byte-wise table of contents for %s\n", outputFile);	//$$$$$ Output for Table of Contents
	fprintf(tocFile, "##\n");														//$$$$$ Output for Table of Contents
	fprintf(tocFile, tempString);													//$$$$$ Output for Table of Contents
	fprintf(tocFile, "##\n");														//$$$$$ Output for Table of Contents
	fprintf(tocFile, "##\n");														//$$$$$ Output for Table of Contents
	fprintf(tocFile, "SIZE %i\n", listSize);									//$$$$$ Output for Table of Contents
	int runningTotal = 0;														//$$$$$ Output for Table of Contents

	///Open the binary file for reading								//@@@	Output for Binary File
	FILE * myFile = fopen(outputFile, "wb");						//@@@	Output for Binary File
																	
	///now we read the files and their given names					
	int i = 0;														
	for(i = 0; i<listSize; i++){									
		fgets(line, NUMCOLUMNS, pdbList);							//##	Parseing for listFile
		strcpy(line2, line);										//##	Parseing for listFile
		ptr = strtok(line2, "=\t\n");	///thius is the filename	//##	Parseing for listFile

		////immidiately get the file and parse it					//====  Parsing for PDB file
		pdbParser * pparse = new pdbParser(ptr);					//====  Parsing for PDB file
		AtomList * alist = pparse->getAtomList();					//====  Parsing for PDB file
		alist->highlightAllAAs();									//====  Parsing for PDB file
		caList * clist = new caList(alist);							//====  Parsing for PDB file

		///continue parsing the list input							//##	Parseing for listFile
		ptr = strtok(line, "=\t\n");	///thius is the filename	//##	Parseing for listFile
		ptr = strtok(NULL, "=\t\n");///this is the given name		//##	Parseing for listFile

		/////write the file contents								//@@@	Output for Binary File
		int ending = writeBinaryPdbFile(clist, ptr, myFile);		//@@@	Output for Binary File

		sprintf(tempString, "%i\t%s\n",runningTotal, ptr);			//$$$$$ Output for Table of Contents
		fprintf(tocFile, tempString);									//$$$$$ Output for Table of Contents
		runningTotal += ending;

		///clear the temp data
		delete(clist);
		delete(alist);
		delete(pparse);
	}

	///close the binaryFile 
	fclose(myFile);

	///close the list file
	fclose(pdbList);

	///close the table of contents
	sprintf(tempString, "FILESIZE %i\n", runningTotal);							//$$$$$ Output for Table of Contents
	fprintf( tocFile, tempString);													//$$$$$ Output for Table of Contents
	fclose( tocFile);

	delete[](tempString);
	delete[](line);
	delete[](line2);
}


void readJoinedBinFile(char * filename)
{
	FILE * binFile = fopen(filename, "rb");

	bool stuff = true;
	while(stuff){
		stuff = readBinaryPdbFile(binFile);	
	}
}



void writePdbFile(char * pdbFile, char * outputFile)
{
	pdbParser * pparse = new pdbParser(pdbFile);
	AtomList * alist = pparse->getAtomList();
	alist->highlightAllAAs();
	caList * clist = new caList(alist);

	FILE * myFile = fopen(outputFile, "wb");
	writeBinaryPdbFile(clist, "1BQK1 ", myFile);

	fclose(myFile);
}


void readBinFile(char * binFile)
{
	FILE * myFile = fopen(binFile, "rb");
	readBinaryPdbFile(myFile);
	fclose(myFile);

}

int writeBinaryPdbFile(caList * clist, char * gname, FILE * binFile)
{
	int i = 0;
	int j = 0;
	int k = 0;
	char * protBuffer = new char[7];
	char * aaBuffer = new char[2];
	char * atomBuffer = new char[11];
	int lastByte = 0;

	////first write the Protein Header
	for(i = 0; i<7*8; i++){ CLEAR(protBuffer, i); }
	writeProteinHeader(protBuffer, 0, gname, clist->numAcids);
	fwrite(protBuffer, 7, 1, binFile);
	lastByte += 7;

	for(i = 0; i<clist->numAcids; i++){
		////next write hte Amino Acid headers
		for(k = 0; k<2*8; k++){ CLEAR(aaBuffer, k); }
		writeAminoAcid(aaBuffer, 0, encodeAAcode(clist->acids[i]->residue_name), clist->acids[i]->numAtoms);
		fwrite(aaBuffer, 2, 1, binFile);
		lastByte += 2;

		////last write entries for each atom.
		for(j = 0; j<clist->acids[i]->numAtoms; j++){
			////last write entries for each atom.
			for(k = 0; k<11*8; k++){ CLEAR(atomBuffer, k); }
			writeAtom(atomBuffer, 0, clist->acids[i]->atoms[j]->coords, translateAtomCode(clist->acids[i]->atoms[j]->atom_name) );
			fwrite(atomBuffer, 11, 1, binFile);
			lastByte +=11;
		}
	}

	delete[](protBuffer);
	delete[](aaBuffer);
	delete[](atomBuffer);

	return lastByte;
}

bool readBinaryPdbFile(FILE * binFile)
{
	int i = 0;
	int j = 0;
	int k = 0;
	char * protBuffer = new char[7];
	char * aaBuffer = new char[2];
	char * atomBuffer = new char[11];

	int numAA = 0;
	int numAtoms = 0;

	bool numRead = true;

	////first write the Protein Header
	for(i = 0; i<7*8; i++){ CLEAR(protBuffer, i); }
	numRead = fread(protBuffer, 7, 1, binFile);	
	if(numRead == 0){ 
		delete[](protBuffer);
		delete[](aaBuffer);
		delete[](atomBuffer);
		return false; 
	}
	readProteinHeader(protBuffer, 0, &numAA);

	for(i = 0; i<numAA; i++){
		////next write hte Amino Acid headers
		for(k = 0; k<2*8; k++){ CLEAR(aaBuffer, k); }
		fread(aaBuffer, 2, 1, binFile);
		readAminoAcid(aaBuffer, 0, &numAtoms);

		////read entries for each atom.
		for(j = 0; j<numAtoms; j++){
			////last write entries for each atom.
			for(k = 0; k<11*8; k++){ CLEAR(atomBuffer, k); }
			fread(atomBuffer, 11, 1, binFile);
			readAtom(atomBuffer, 0);
		}
	}

	delete[](protBuffer);
	delete[](aaBuffer);
	delete[](atomBuffer);
	return true;
}


////56 bits, or 7 byte header
//////////////////////////////////////////////////////////////////////
//Encodes a header for a PDB file into a binary array/////////////////
//////////////////////////////////////////////////////////////////////
int writeProteinHeader(char * array, int offset, char * name, int size)
{
	int currentOffset = offset;
	if(strlen(name) != 6){ printf("Protein Name %s incorrect number of Alphanumeric characters\n", name); }
	currentOffset = encodeString(name, array, currentOffset);		///36 bits
	currentOffset = encodeInt(size, 14, array, currentOffset);		///14 bits
	return currentOffset+6;											//// Total: 50, + 6 spacer bits to make 56, or 7 bytes.
}


//////////////////////////////////////////////////////////////////////
//Decodes a header for a PDB file from a binary array, for debugging//
//////////////////////////////////////////////////////////////////////
void readProteinHeader(char * array, int offset, int * numAA)
{
	char * theString = new char[NUMCOLUMNS];
	int * theInt= new int[1];
	int currentOffset = offset;
	currentOffset = decodeString(array, currentOffset, 6, theString);
	currentOffset = decodeInt(array, currentOffset, 14, theInt);

	theString[6] = '\0';
	printf("Given PDB Code: %s\n", theString);
	printf("Number of Amino Acids: %i\n", theInt[0]);

	numAA[0] = theInt[0];

	delete[](theString);
	delete[](theInt);
}


//////////////////////////////////////////////////////////////////////
//Decodes a header for a PDB file from a binary array, and outputs////
//////////////////////////////////////////////////////////////////////
int readProteinHeader(char * array, int offset, int * numAA, char * givenName)
{
	char * theString = new char[NUMCOLUMNS];
	int * theInt= new int[1];
	int currentOffset = offset;
	currentOffset = decodeString(array, currentOffset, 6, theString);
	currentOffset = decodeInt(array, currentOffset, 14, theInt);

	theString[6] = '\0';
	numAA[0] = theInt[0];
	strcpy(givenName, theString);

	delete[](theString);
	delete[](theInt);

	return offset + 56;
}



////16 bits, or 2 byte AA header
//////////////////////////////////////////////////////////////////////
//Encodes a header for amino acid into a binary array/////////////////
//////////////////////////////////////////////////////////////////////
int writeAminoAcid(char * array, int offset, int aaType, int numAtoms)
{
	int currentOffset = offset;
	currentOffset = encodeInt(aaType, 5, array, currentOffset);			///encode Type.			5 bits
	currentOffset = encodeInt(numAtoms, 5, array, currentOffset);		///encode num Atoms		5 bits
	return currentOffset+6;												///Total 10 + 6 bits spacer, for 16 bits, or 2 bytes
}


//////////////////////////////////////////////////////////////////////
//Decodes amino acid information, from a binary array, for debugging//
//////////////////////////////////////////////////////////////////////
void readAminoAcid(char * array, int offset, int * numAtoms)
{
	int int1 = 0;
	int int2 = 0;
	int currentOffset = offset;
	currentOffset = decodeInt(array, currentOffset, 5, &int1);
	currentOffset = decodeInt(array, currentOffset, 5, &int2);

	char * thisAA = decodeAAcode(int1);
	printf("AA type: %s\n", thisAA);
	printf("numAtoms: %i\n", int2);

	numAtoms[0] = int2;

	delete[](thisAA);
}


//////////////////////////////////////////////////////////////////////
//Decodes amino acid information, from a binary array, and outputs////
//////////////////////////////////////////////////////////////////////
int readAminoAcid(char * array, int offset, int * numAtoms, int * aaType)
{
	int int1 = 0;
	int int2 = 0;
	int currentOffset = offset;
	currentOffset = decodeInt(array, currentOffset, 5, &int1);
	currentOffset = decodeInt(array, currentOffset, 5, &int2);

	numAtoms[0] = int2;
	aaType[0] = int1;

	return offset + 16;
}



////88 bits, or 11 byte atom header
//////////////////////////////////////////////////////////////////////
//Encodes an atom from a PDB file into a binary array/////////////////
//////////////////////////////////////////////////////////////////////
int writeAtom(char * array, int offset, double * coords, int atomType)
{
	int currentOffset = offset;
	currentOffset = encodeInt(atomType, 6, array, currentOffset);		///encode atomType		6 bits
	currentOffset = encodeDouble(coords[0], array, currentOffset);				///encode x				25 bits
	currentOffset = encodeDouble(coords[1], array, currentOffset);				///						25 bits
	currentOffset = encodeDouble(coords[2], array, currentOffset);				///						25 bits
	return currentOffset+7;												///Total 81 bits + 7 bits spacer for 88 bits, or 11 bytes.
}


//////////////////////////////////////////////////////////////////////
//Decodes atom information, from a binary array, for debugging////////
//////////////////////////////////////////////////////////////////////
void readAtom(char * array, int offset)
{
	int currentOffset = offset;
	int atomType = 0;
	double x = 0;
	double y = 0;
	double z = 0;
	currentOffset = decodeInt(array, currentOffset, 6, &atomType);
	currentOffset = decodeDouble(array, currentOffset, &x);
	currentOffset = decodeDouble(array, currentOffset, &y);
	currentOffset = decodeDouble(array, currentOffset, &z);

	printf("Atom: %s %lf, %lf, %lf\n", encodeAtomCode(atomType), x, y, z);
}



//////////////////////////////////////////////////////////////////////
//Decodes atom information, from a binary array, for debugging////////
//////////////////////////////////////////////////////////////////////
int readAtom(char * array, int offset, double * x, double * y, double * z, int * atomType)
{
	int currentOffset = offset;
	currentOffset = decodeInt(array, currentOffset, 6, atomType);
	currentOffset = decodeDouble(array, currentOffset, x);
	currentOffset = decodeDouble(array, currentOffset, y);
	currentOffset = decodeDouble(array, currentOffset, z);

	return offset + 88;
}





/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
////Helper functions for encoding and decoding binary information////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

//Encodes chars into appropriate integers for binary storage.  Alphanumeric only.
int encodeChar(char input, char * array, int offset)
{

	int i = 0;
	int result = 0;

	switch(input){
		case 'A': result =  0; break;
		case 'B': result =  1; break;
		case 'C': result =  2; break;
		case 'D': result =  3; break;
		case 'E': result =  4; break;
		case 'F': result =  5; break;
		case 'G': result =  6; break;
		case 'H': result =  7; break;
		case 'I': result =  8; break;
		case 'J': result =  9; break;
		case 'K': result = 10; break;
		case 'L': result = 11; break;
		case 'M': result = 12; break;
		case 'N': result = 13; break;
		case 'O': result = 14; break;
		case 'P': result = 15; break;
		case 'Q': result = 16; break;
		case 'R': result = 17; break;
		case 'S': result = 18; break;
		case 'T': result = 19; break;
		case 'U': result = 20; break;
		case 'V': result = 21; break;
		case 'W': result = 22; break;
		case 'X': result = 23; break;
		case 'Y': result = 24; break;
		case 'Z': result = 25; break;
		case '0': result = 26; break;
		case '1': result = 27; break;
		case '2': result = 28; break;
		case '3': result = 29; break;
		case '4': result = 30; break;
		case '5': result = 31; break;
		case '6': result = 32; break;
		case '7': result = 33; break;
		case '8': result = 34; break;
		case '9': result = 35; break;
		case ' ': result = 36; break;
	}

	for(i = 0; i<6; i++){
		if(result & 1){
			SET( array, (offset+i) ); 
		}
		else{
			CLEAR( array, (offset+i) );
		}
		result = result >> 1;
	}

	return offset + 6;
}


///Decodes chars, alphanumeric only
char decodeChar(char * array, int offset, int numBits)
{
	int i;
	int target = 0;
	for(i = 0; i<numBits; i++){
		if( TEST( array, (offset+i) ) ){
			target += (int) pow(2, i);
		}
	}

	char result;
	switch(target){
		case  0: result = 'A'; break;
		case  1: result = 'B'; break;
		case  2: result = 'C'; break;
		case  3: result = 'D'; break;
		case  4: result = 'E'; break;
		case  5: result = 'F'; break;
		case  6: result = 'G'; break;
		case  7: result = 'H'; break;
		case  8: result = 'I'; break;
		case  9: result = 'J'; break;
		case 10: result = 'K'; break;
		case 11: result = 'L'; break;
		case 12: result = 'M'; break;
		case 13: result = 'N'; break;
		case 14: result = 'O'; break;
		case 15: result = 'P'; break;
		case 16: result = 'Q'; break;
		case 17: result = 'R'; break;
		case 18: result = 'S'; break;
		case 19: result = 'T'; break;
		case 20: result = 'U'; break;
		case 21: result = 'V'; break;
		case 22: result = 'W'; break;
		case 23: result = 'X'; break;
		case 24: result = 'Y'; break;
		case 25: result = 'Z'; break;
		case 26: result = '0'; break;
		case 27: result = '1'; break;
		case 28: result = '2'; break;
		case 29: result = '3'; break;
		case 30: result = '4'; break;
		case 31: result = '5'; break;
		case 32: result = '6'; break;
		case 33: result = '7'; break;
		case 34: result = '8'; break;
		case 35: result = '9'; break;
		case 36: result = ' '; break;
	}

	return result;
}


//////encodes ints
int encodeInt(int input, int charSize, char * array, int offset)
{
	int i = 0;
	int copy = input;

	for(i = 0; i<charSize; i++){
		if(copy & 1){
			SET( array, (offset+i) ); 
		}
		else{
			CLEAR( array, (offset+i) );
		}
		copy = copy >> 1;
	}

	return offset + charSize;
}


//////decodes ints
int decodeInt(char * array, int offset, int numBits, int * myInt)		///serialize int to bitarray
{
	int i;
	int target = 0;
	for(i = 0; i<numBits; i++){
		if( TEST( array, (offset+i) ) ){
			target += (int) pow(2, i);
		}
	}
	
	myInt[0] = target;
	return offset+numBits;
}


////encodes doubles for Fortran Real(8.3)
int encodeDouble(double input, char * array, int offset)
{
	int myInt = (int) (input*1000);	///PDB doubles have max precision to 3 decimal points.  we multiply by 1000 and store the integer.
	int currentOffset = offset;

	if(myInt > 0){
		SET( array, currentOffset);
		currentOffset++;
		currentOffset = encodeInt(myInt, 24, array, currentOffset);
	}
	else{
		CLEAR(array, currentOffset);
		currentOffset++;
		currentOffset = encodeInt(-myInt, 24, array, currentOffset);
	}
	return currentOffset;
}


////decodes doubles for Fortran Real(8.3)
int decodeDouble(char * array, int offset, double * myDouble)
{
	bool isPos = (bool) TEST(array, offset);
	int currentOffset = offset+1;
	int myInt = 0;

	currentOffset = decodeInt(array, currentOffset, 24, &myInt);

	if(isPos){
		myDouble[0] = (((double) myInt) / ((double) 1000) );
	}
	else{
		myDouble[0] = -(((double) myInt) / ((double) 1000) );
	}

	return currentOffset;
}


///MyChars use size 6 bits
int encodeString(char * string, char * array, int offset)
{
	int length = strlen(string);
	int i =0;
	int currentOffset = offset;
	for(i = 0; i<length; i++){
		currentOffset = encodeChar(string[i], array, currentOffset);
	}

	return currentOffset;
}


///MyChars use size 6 bits
int decodeString(char * array, int offset, int numChars, char * output)
{
	int i = 0;
	for(i = 0; i<numChars; i++){
		output[i] = decodeChar(array, offset+i*6, 6);
	}

	return offset+(i*6);
}


int encodeAAcode(char * code)
{
	int result;
	
	if		(strcmp(code, "GLY") == 0) result = 0;		//	G - Glycine (Gly) 
	else if (strcmp(code, "PRO") == 0) result = 1; 	//	P - Proline (Pro) 
	else if (strcmp(code, "ALA") == 0) result = 2;		//	A - Alanine (Ala) 
	else if (strcmp(code, "VAL") == 0) result = 3;		//	V - Valine (Val) 
	else if (strcmp(code, "LEU") == 0) result = 4;		//	L - Leucine (Leu) 
	else if (strcmp(code, "ILE") == 0) result = 5;		//	I - Isoleucine (Ile) 
	else if (strcmp(code, "MET") == 0) result = 6;		//	M - Methionine (Met) 
	else if (strcmp(code, "CYS") == 0) result = 7;		//	C - Cysteine (Cys) 
	else if (strcmp(code, "PHE") == 0) result = 8;		//	F - Phenylalanine (Phe) 
	else if (strcmp(code, "TYR") == 0) result = 9;		//	Y - Tyrosine (Tyr) 
	else if (strcmp(code, "TRP") == 0) result = 10;		//	W - Tryptophan (Trp) 
	else if (strcmp(code, "HIS") == 0) result = 11;		//	H - Histidine (His) 
	else if (strcmp(code, "LYS") == 0) result = 12;		//	K - Lysine (Lys) 
	else if (strcmp(code, "ARG") == 0) result = 13;		//	R - Arginine (Arg) 
	else if (strcmp(code, "GLN") == 0) result = 14;		//	Q - Glutamine (Gln) 
	else if (strcmp(code, "ASN") == 0) result = 15;		//	N - Asparagine (Asn) 
	else if (strcmp(code, "GLU") == 0) result = 16;		//	E - Glutamic Acid (Glu) 
	else if (strcmp(code, "ASP") == 0) result = 17;		//	D - Aspartic Acid (Asp) 
	else if (strcmp(code, "SER") == 0) result = 18;		//	S - Serine (Ser) 
	else if (strcmp(code, "THR") == 0) result = 19;		//	T - Threonine (Thr) 
	else {
		printf("ERROR: Unidentified Amino Acid Tag Specified\n");
		result = 'X';
		exit(1);
	}

	return result;

}



char * decodeAAcode(int code)
{
	char * result = new char[4];
	
	switch(code){
		case 0:			strcpy(result, "GLY"); break;
		case 1:         strcpy(result, "PRO"); break;
		case 2:         strcpy(result, "ALA"); break;
		case 3:         strcpy(result, "VAL"); break;
		case 4:         strcpy(result, "LEU"); break;
		case 5:         strcpy(result, "ILE"); break;
		case 6:         strcpy(result, "MET"); break;
		case 7:         strcpy(result, "CYS"); break;
		case 8:         strcpy(result, "PHE"); break;
		case 9:         strcpy(result, "TYR"); break;
		case 10:        strcpy(result, "TRP"); break;
		case 11:        strcpy(result, "HIS"); break;
		case 12:        strcpy(result, "LYS"); break;
		case 13:        strcpy(result, "ARG"); break;
		case 14:        strcpy(result, "GLN"); break;
		case 15:        strcpy(result, "ASN"); break;
		case 16:        strcpy(result, "GLU"); break;
		case 17:        strcpy(result, "ASP"); break;
		case 18:        strcpy(result, "SER"); break;
		case 19:        strcpy(result, "THR"); break;
		default:
			printf("ERROR: Unidentified Amino Acid Tag Specified\n");
			strcpy(result, "XXX");
			exit(1);
			break;
	}
	
	return result;

}
