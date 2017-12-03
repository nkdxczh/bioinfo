// TermList.h: interface for the TermList class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_TERMLIST_H_)
#define _TERMLIST_H_

class TermList  
{
public:
	TermList(int c, GATenv * e);
	virtual ~TermList();

	GATenv * env;

	int size;
	int capacity;
	AlignNode * *nodes;
	bool hasZeroElement;

	void addNode(AlignNode * n);
	bool checkDuplicate(AlignNode * node1, AlignNode * node2);
	void sortEntries(void);

	int * idList;
	int idListSize;
	int idListCap;
	void printIdList();

};

#endif // !defined(_TERMLIST_H_)


