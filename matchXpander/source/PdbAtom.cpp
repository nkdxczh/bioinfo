// PdbAtom.cpp: implementation of the PdbAtom class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

PdbAtom::PdbAtom()
{

	atom_number = 0;
	atom_name = new char[15];
	residue_name = new char[15];
	residue_number = 0;
	coords = new double[3];

	significance = 0;
	spread_color = 0;

	fileLine = new char[NUMCOLUMNS];
	assocAA = NULL;							///in MatchXpander this is vestigial.  it is used as a convenience for parsing.

}

PdbAtom::~PdbAtom()
{

	delete[](atom_name);
	delete[](residue_name);
	delete[](fileLine);
	delete[](coords);
	if(assocAA != NULL){
		delete[](assocAA);
	}

}

PdbAtom * PdbAtom::copy()
{
	PdbAtom * result = new PdbAtom();
	strcpy(result->atom_name, atom_name);
	strcpy(result->residue_name, residue_name);
	strcpy(result->fileLine, fileLine);
	if(assocAA != NULL){
		result->assocAA = new char[strlen(assocAA)+1];
		strcpy(result->assocAA, assocAA);
	}else{
		result->assocAA = NULL;
	}

	result->atom_number = atom_number;
	result->residue_number = residue_number;
	result->significance = significance;
	result->spread_color = spread_color;
	result->coords[0] = coords[0];
	result->coords[1] = coords[1];
	result->coords[2] = coords[2];

	return result;
}

void PdbAtom::toString(void)
{
	printf("%i ", atom_number);
	printf("%s ", atom_name);
	printf("%s ", residue_name);
	printf("%i ", residue_number);
	printf("%lf ", coords[0]);
	printf("%lf ", coords[1]);
	printf("%lf ", coords[2]);
	printf("%i \n", significance);

}



