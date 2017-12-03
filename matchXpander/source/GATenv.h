// GATenv.h: interface for the GATenv class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_GAT_ENVIRONMENT_)
#define _GAT_ENVIRONMENT_

class GATenv  
{
public:
	GATenv();
	virtual ~GATenv();
	void addSource(AtomList * alist);
	void setTarget(caList * targ);
	void addMatchArray(int index, int * array);
	void setUpMatchArray();
	void setUpMatchArray(int num);
	void setUpMatchArraySizes(int i, int val);
	void setUpRMSDArraySizes();	///call this only after numArrays has been set in env.
	void sortMatchArrays();
	void filterMatchArrays(double rmsd);
	void filterMatchArrays(int numBest);
	void resetEnvironment(void);
	set_t query(int news, double * accumulationMatrix);

//	///debug
//	void resetMatchArrays(void);
//	///debug

	/////Tests RMSD alignments from GAT file.
//	void testRMSDs();	DEPRECATED

//	int numSources;
	int numArrays;
//	caList * *lists;
	caList * sourceList;
	caList * targetList;
	int * *matchArray;
	int * matchArraySizes;
	LatticeHash ***hashSets;
	Motif * motif;						////stores sidechain info accessibly
	bool hasAssocArrays;

	////RMSD information stored from target
	double * RMSDarray;
	double ** transArray;
	double ** rotArray;
};


#endif
