// HierarchialController.h: interface for the HierarchialController class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_HIERARCHIALCONTROLLER_H_)
#define _HIERARCHIALCONTROLLER_H_

class HierarchialController  
{
public:
	HierarchialController(GATenv * e);
	virtual ~HierarchialController();

	GATenv * env;

	PriorityStack * stack;
	ActiveList * active;
	ActiveList * newList;			////holeds multimatch AlignNodes for processing
	AlignNode * *tierNodes;			////holds one for each tier
	TermList * term;
	int startTime;
	int endTime;

	/////Various levels of testing implemented
	void beginTesting();											///Global test management
	bool testIndex(int ind, AlignNode * node);						///tests a particular Index target, with a particular node
	bool testArray(int * array, int length, AlignNode * tempNode);	///tests a single tier
	void multiLevelTest(AlignNode * n, bool initial);				///tests all levels for a single node

	/////Output
	void printResults();



};

#endif // !defined(_HIERARCHIALCONTROLLER_H_)



