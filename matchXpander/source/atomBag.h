// atomBag.h: interface for the atomBag class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_ATOMBAG_H_)
#define _ATOMBAG_H_

class atomBag  
{
public:
	atomBag();
	virtual ~atomBag();
	void addAtom(PdbAtom * a);
	
	PdbAtom * *atoms;
	int numAtoms;



};

#endif

