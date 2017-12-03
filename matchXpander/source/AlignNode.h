// AlignNode.h: interface for the AlignNode class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_ALIGNNODE_H_)
#define _ALIGNNODE_H_

class AlignNode  
{
public:
	AlignNode(int ident, int s, GATenv * e);
	virtual ~AlignNode();

	int id;
	int size;
	int * sourceIds;			///These are the residues we need to match onto the target
	int * targetIds;			///these are the residues the source pts match to
	int currentTier;
	bool isdead;
	double * accumulationMatrix; ////4x4 matrix which contains all transformations from global to target.
	GATenv * env;
//	double * pointBuffer;
	double RMSD;
	int startingIndex;			////This keeps track of the position that this will search from on the tier
								////Multimatch-copies will start from something positive, rest start at 0.
	double * tempMatrix;		////This is used for holding the transformation of the last match, for future use.

//	bool sane;					///debugging flag

	////////////////////////////
	////Geometric Set/gets
	////////////////////////////
	void initializeAccumulationMatrix(double * trans, double * rot);///only call this after both trans and rot are set.
	void initializeAccumulationMatrix(double * accum);///only call this after both trans and rot are set.
	double * getTotalRotation(void);						///outputs accumulated translation/rotation info
	double * getTotalTranslation(void);
	AlignNode * copy(int id);
	bool element(int val);			///This is negative if the source pt is not in the Node.  Pos if so.

	/////////////////////////////
	////Set the Point info
	/////////////////////////////
	void setGeometricPts(int * matchArray);		/////Call this only after Init. AccumMAtrix() has been called
	void addGeometricPt(int peg, int hole);		/////Adds the index of a peg (source) and a hole (target), incremets size, 

	/////////////////////////////
	////Actual Geometric Calls
	/////////////////////////////
	double testPt(int source, int target);		///Set newAlign to true if you want rot/trans applied to pts before 
												///Returns -1 if the match has a poitn outside of EXTMATCHTHRESHOLD.  returns RMSD otherwise.
//	bool sanityCheck(void);
	bool rmsdCheck(void);
	double getBestFitRMSD(void);

	/////////////////////////////
	////Helper Functions
	/////////////////////////////
//	double * *genCoords(int queryPt);

	/////////////////////////////
	////Output Calls
	/////////////////////////////
	void sortMatchIds(void);
	void printResults(void);
};

#endif // !defined(_ALIGNNODE_H_)



