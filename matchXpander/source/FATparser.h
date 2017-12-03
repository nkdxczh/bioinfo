// FATparser.h: interface for the FATparser class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FATPARSER_H_)
#define _FATPARSER_H_


class FATparser  
{
public:
	FATparser(char *filename);
	virtual ~FATparser();

	bool readFile();
	bool readGlobals();
	AtomList * readSource();
	AtomList * readTarget();

	char * name;	///the name of the GAT file
	int numSources;
	AtomList ** sources;
	AtomList * target;
	FILE * FATfile;
	char * sourceFileName;
	char * targetFileName;
	bool AAMode;	///true for Amino Acids, false for atoms

};



#endif



