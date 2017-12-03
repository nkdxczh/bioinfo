// Motif.h: interface for the Motif class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(__MX_MOTIF__)
#define __MX_MOTIF__

class Motif  
{
public:
	Motif(caList * clist, int * selectedAtoms);
	virtual ~Motif();

	int numResidues;
	int * numAtoms;
	PdbAtom ***atoms;

};

#endif // !defined(AFX_MOTIF_H__73C104F3_0F15_48EE_8991_91B11BDE925E__INCLUDED_)
