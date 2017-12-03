// aminoAcid.h: interface for the aminoAcid class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_AMINOACID_H_)
#define _AMINOACID_H_

class aminoAcid  
{
public:
	aminoAcid();
	virtual ~aminoAcid();
	void addAtom(PdbAtom * a);
//	void setCalpha(PdbAtom * a);
	aminoAcid * copy();

//	PdbAtom * Calpha;
	PdbAtom * *atoms;
	int numAtoms;
	char * assocArray;
	char * residue_name;
	int residue_number;
	int significance;

	////debug
	int clearAA();/////DEBUG

};


#endif



