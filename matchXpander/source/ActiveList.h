// ActiveList.h: interface for the ActiveList class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_ACTIVELIST_H_)
#define _ACTIVELIST_H_


class ActiveList  
{
public:
	ActiveList(GATenv * e);
	ActiveList(GATenv * e, int c);
	virtual ~ActiveList();

	GATenv * env;
	int size;
	int capacity;
	AlignNode * *list;
	int count;

	void expand();
	void addNode(AlignNode * n);
	AlignNode * remove(int ind);
	AlignNode * remove();
	int getNewId();



};

#endif // !defined(_ACTIVELIST_H_)


