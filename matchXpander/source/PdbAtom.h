// PdbAtom.h: interface for the PdbAtom class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_PDBATOM_H_)
#define _PDBATOM_H_

class PdbAtom  
{
public:
	PdbAtom();
	virtual ~PdbAtom();
	void toString(void);					///output variable vals
	PdbAtom * copy();

	int atom_number;
	char * atom_name;
	char * residue_name;
	int residue_number;
	double * coords;

//	double x;			DEPRECATED
//	double y;			DEPRECATED
//	double z;			DEPRECATED

	int significance;
	float spread_color;

	char * fileLine;
	char * assocAA;

};


#endif


