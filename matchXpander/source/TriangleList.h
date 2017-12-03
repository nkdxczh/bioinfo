// TriangleList.h: interface for the TriangleList class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TRIANGLELIST_H__AB4B96E2_A499_46D7_9C9C_017A43AACA4C__INCLUDED_)
#define AFX_TRIANGLELIST_H__AB4B96E2_A499_46D7_9C9C_017A43AACA4C__INCLUDED_

class TriangleList  
{
public:
	TriangleList();
	virtual ~TriangleList();

	int * *list;
	int size;
	int capacity;

	void expand();
	void addTriangle(int s1, int t1, int s2, int t2, int s3, int t3);
	void print();
};

#endif 


