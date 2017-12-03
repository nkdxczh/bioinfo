// StatOutput.h: interface for the StatOutput class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_GEOHASH_STATOUTPUT_)
#define _GEOHASH_STATOUTPUT_

class StatOutput  
{
public:
	StatOutput(GATenv * e, HierarchialController * h);
	virtual ~StatOutput();

	FILE * currentOutputFile;
	HierarchialController * controller;
	GATenv * env;

	//////MetaCall for GAT file output
	void outputStatFile(char * idHeader);

	//////Call for printing a hashed alignment of the residues
	void printAlign(void);

	//////Call for printing statistical information about the execution
	void printStats(void);

	//////Call for printing the match pairs in residue numbers + residue names
	void printMatches(void);

};

#endif


