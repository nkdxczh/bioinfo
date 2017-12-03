// EdgeSet.h: interface for the EdgeSet class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_EDGESET_H__A80B745E_1486_4F49_98A4_40F4E752028A__INCLUDED_)
#define AFX_EDGESET_H__A80B745E_1486_4F49_98A4_40F4E752028A__INCLUDED_

class EdgeSet  
{
public:
	EdgeSet();
	virtual ~EdgeSet();

	set_t edgeRoot;
	int size;

	int addEdge(int s, int e);
	bool contains(int a, int b); 
	set_t getNeighbors(int s);		///packed with size first

};

#endif // !defined(AFX_EDGESTACK_H__A80B745E_1486_4F49_98A4_40F4E752028A__INCLUDED_)



