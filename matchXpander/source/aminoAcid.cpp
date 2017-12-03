// aminoAcid.cpp: implementation of the aminoAcid class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

aminoAcid::aminoAcid()
{
	numAtoms = 0;
	assocArray = new char[25];
	residue_name = NULL;
	residue_number = -1;
	significance = -1;

}

aminoAcid::~aminoAcid()
{
	if(assocArray != NULL){
		delete[](assocArray);
	}
	if(residue_name != NULL){
		delete[](residue_name);
	}

	int i = 0;
	for(i = 0; i<numAtoms; i++){
		delete(atoms[i]);
	}
	if(numAtoms > 0){
		delete[](atoms);
	}
}

/////DEBUG
int aminoAcid::clearAA()
{
	int temp = numAtoms;

	if(assocArray != NULL){
		delete(assocArray);
	}
	if(residue_name != NULL){
		delete[](residue_name);
	}

	int i = 0;
	for(i = 0; i<numAtoms; i++){
		delete(atoms[i]);
	}
	if(numAtoms > 0 || atoms != NULL){
		delete[](atoms);
	}

	numAtoms = 0;
	return temp;
}
/////DEBUG

void aminoAcid::addAtom(PdbAtom * a)
{
	///Doesnt need to instantiate atoms even if atoms is null.

	PdbAtom * *tempAtoms = new PdbAtom*[numAtoms+1];
	int i = 0;

	for (i = 0; i<numAtoms; i++){
		/////copy the references over
		tempAtoms[i] = atoms[i];
	}

	tempAtoms[numAtoms] = a;

	if (numAtoms > 0){
		delete[](atoms);
	}
	atoms = tempAtoms;
	numAtoms++;

	if(numAtoms == 1){
		if(residue_name != NULL){ delete[](residue_name); }
		residue_name = new char[4];
		sprintf(residue_name, a->residue_name);
		residue_number = a->residue_number;
		significance = a->significance;
	}
	else{
		if( (strcmp(residue_name, a->residue_name) != 0) || (residue_number != a->residue_number) || (significance != a->significance) ){
			printf("ERROR: Inconsistant Residue name (%s vs %s) or Inconsistant residue # (%i vs %i) or Inconsistant significance (%i vs %i)\n", residue_name, a->residue_name, residue_number, a->residue_number, significance, a->significance);
		}
	}

}

/* DEPRECATED
void aminoAcid::setCalpha(PdbAtom * a)
{
	Calpha = a;

}
*/


aminoAcid * aminoAcid::copy()
{
	aminoAcid * result = new aminoAcid();
	if(assocArray == NULL){
		delete[](result->assocArray);
		result->assocArray = NULL;
	}
	else{
		strcpy(result->assocArray, assocArray);	///copy the assoc array into the new amino acid
	}
	if(residue_name == NULL){
		result->residue_name = NULL;
	}
	else{
		result->residue_name = new char[4]; 
		strcpy(result->residue_name, residue_name);	///copy the assoc array into the new amino acid
	}

	int i = 0;
	PdbAtom * tempAtom;
	for(i = 0; i<numAtoms; i++){
		tempAtom = atoms[i]->copy();
		result->addAtom(tempAtom);
	}
	return result;
}

