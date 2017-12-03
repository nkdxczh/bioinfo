// LatticeHash.h: interface for the LatticeHash class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_LATTICEHASH_H__27606015_B6EF_4DD7_997A_F9486ECFA2FB__INCLUDED_)
#define AFX_LATTICEHASH_H__27606015_B6EF_4DD7_997A_F9486ECFA2FB__INCLUDED_

class LatticeHash  
{
public:
	LatticeHash(caList * target, PdbAtom * templateAtom);
	virtual ~LatticeHash();

	int size;
	int * ids;
	int * backIds;
	TargetGrid * grid;
	set_t set;
	set_t backSet;
	AtomList * list;

	int * queryShell(double x, double y, double z, double range);
	set_t querySphere(double x, double y, double z, double range, set_t output);
	int mapBack(int k);

};

#endif // 



