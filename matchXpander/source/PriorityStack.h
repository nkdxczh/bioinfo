// PriorityStack.h: interface for the PriorityStack class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_PRIORITYSTACK_H_)
#define _PRIORITYSTACK_H_

class PriorityStack  
{
public:
	PriorityStack(GATenv * e);
	virtual ~PriorityStack();

	int ** flags;
	int * sizes;
	int height;
	int currentLevel;
	GATenv * env;

	void nextLevel();
	int * currentFlags();
	int currentSize();
	void resetStack();

	bool changeLevel(int i);
	void printStatus();


};

#endif // !defined(_PRIORITYSTACK_H_)


