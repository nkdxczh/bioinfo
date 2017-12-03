// Motif.cpp: implementation of the Motif class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Motif::Motif(caList * clist, int * selectedAtoms)
{
	int i = 0;
	int j = 0;
	int k = 0;

	numResidues = clist->numAcids;
	int atomsUsed = selectedAtoms[0];
	numAtoms = new int[numResidues];
	atoms = new PdbAtom**[numResidues];

	printf("Motif Residues: ");
	for(i = 0; i<clist->numAcids; i++){
		printf("%i ", clist->acids[i]->residue_number);
	}
	printf("\n");


	int tempCount = 0;
	for(i = 0; i<numResidues; i++){
		for(j = 0; j<atomsUsed; j++){
			for(k = 0; k<clist->acids[i]->numAtoms; k++){
				if(clist->acids[i]->atoms[k]->atom_number == selectedAtoms[j+1]){
					tempCount++;
					break;
				}
			}
		}
		numAtoms[i] = tempCount;
		tempCount = 0;

		atoms[i] = new PdbAtom*[numAtoms[i]];
		for(j = 0; j<atomsUsed; j++){
			for(k = 0; k<clist->acids[i]->numAtoms; k++){
				if(clist->acids[i]->atoms[k]->atom_number == selectedAtoms[j+1]){
					atoms[i][tempCount] = clist->acids[i]->atoms[k];
					tempCount++;
					break;
				}
			}
		}

		tempCount = 0;
	}

	/////ORder the C-alphas in the first position.  Seed Matches should prefer C-alpha if it is available.
	for(i = 0; i<numResidues; i++){
		for(j = 0; j<numAtoms[i]; j++){
			if(strcmp(atoms[i][j]->atom_name, "CA") == 0){
				PdbAtom * tempAtom = atoms[i][j];
				atoms[i][j] = atoms[i][0];
				atoms[i][0] = tempAtom;
				break;
			}
		}
		for(j = 0; j<numAtoms[i]; j++){
			printf("%s", atoms[i][j]->fileLine);
		}
	}
}

Motif::~Motif()
{
	int i = 0;
	for(i = 0; i<numResidues; i++){
		delete[](atoms[i]);
	}
	delete[](atoms);
	delete[](numAtoms);

}
