// pdbParser.h: interface for the pdbParser class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_GEOHASH_PDBPARSER_)
#define _GEOHASH_PDBPARSER_

class pdbParser  
{
public:
	///**you must remember to set the atomList after calling the constructor
	///**The order is as follows:
	///**
	///**1) Construct pdbParser
	///**2) call atomCount to find out how many atoms
	///**3) construct appropriate sized AtomList
	///**4) call setAtomList() when you have the new AtomList
	///**5) call readPdbFile() to fill the AtomList



	pdbParser(char * filename);			
	virtual ~pdbParser();
	void readPdbFile();									//parse tthe file
	int atomCount(void);								//count the # of atoms
//	void setAtomList(AtomList * alist);
//	char * stringCopy(char * str, int start, int end);	//utility method
	AtomList * getAtomList(void);

	AtomList * atomList;

	FILE * pdbFile;

};

#endif


