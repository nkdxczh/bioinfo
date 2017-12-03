#define CHARSIZE 8

#define MASK(y) ( 1 << y % CHARSIZE )

#define BITSLOT(y) ( y / CHARSIZE )

#define SET( x, y) ( x[BITSLOT(y)] |= MASK(y) )

#define CLEAR(x, y) ( x[BITSLOT(y)] &= ~MASK(y) )

#define TEST(x, y) ( x[BITSLOT(y)] & MASK(y) )

#define NUMSLOTS(n) ((n + CHARSIZE - 1) / CHARSIZE) 

	caList * generateCaListFromBinary(char * binary, int startingByte);
	AtomList * generateAtomListFromBinary(char * binary, int startingByte);

	int * parseTOC(char * toc, int * binFileSize, int * nFiles);
	char * loadBinFile(char * binFile, int fileSize);

	void parsePdbListAndOutputBinFile(char * list, char * outputFile, char * tableOfContentsName);
	void readJoinedBinFile(char * filename);

	void writePdbFile(char * pdbFile, char * outputFile);
	void readBinFile(char * binFile);

	///returns the last BYTE (NOT the last BIT) of the PDB file in the current File)
	int writeBinaryPdbFile(caList * clist, char * gname, FILE * binFile);
	bool readBinaryPdbFile(FILE * binFile);


	///File Level Binary Parsing and Writing Functions
	int writeProteinHeader(char * array, int offset, char * name, int size);
	void readProteinHeader(char * array, int offset, int * numAA);
	int readProteinHeader(char * array, int offset, int * numAA, char * givenName);

	int writeAminoAcid(char * array, int offset, int aaType, int numAtoms);
	void readAminoAcid(char * array, int offset, int * numAtoms);
	int readAminoAcid(char * array, int offset, int * numAtoms, int * aaType);

	int writeAtom(char * array, int offset, double * coords, int atomType);
	void readAtom(char * array, int offset);
	int readAtom(char * array, int offset, double * x, double * y, double * z, int * atomType);


	///Binary Encoding Functions
	int encodeChar(char input, char * array, int offset);			//serialize char to bitArray
	char decodeChar(char * array, int offset, int numBits);			//deserialize bitArray to char

	int encodeInt(int input, int charSize, char * array, int offset);		///serialize int to bitarray
	int decodeInt(char * array, int offset, int numBits, int * myInt);		///deserialize int to bitarray

	int encodeDouble(double input, char * array, int offset);
	int decodeDouble(char * array, int offset, double * myDouble);			///deserialize binary to double

	int encodeString(char * string, char * array, int offset);
	int decodeString(char * array, int offset, int numChars, char * output);

	int encodeAAcode(char * code);
	char * decodeAAcode(int code);


