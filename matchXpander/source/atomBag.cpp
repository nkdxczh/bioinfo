// atomBag.cpp: implementation of the atomBag class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

atomBag::atomBag()
{
	numAtoms = 0;
	atoms = NULL;
}

atomBag::~atomBag()
{
	int i = 0;
	for(i = 0; i<numAtoms; i++){
		delete(atoms[i]);
	}
	if(atoms != NULL){
		delete[](atoms);
	}
//	if(numAtoms > 0){
//		delete[](atoms[i]);
//	}

}

void atomBag::addAtom(PdbAtom * a){
	int i = 0;
	PdbAtom * *tempAtoms = new PdbAtom*[numAtoms +1];

	for (i = 0; i<numAtoms; i++){
		tempAtoms[i] = atoms[i];
	}
	tempAtoms[numAtoms] = a;

	if(numAtoms > 0){
		delete[](atoms);	
	}

	atoms = tempAtoms;

	numAtoms++;
}



