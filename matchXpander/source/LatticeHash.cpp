// LatticeHash.cpp: implementation of the LatticeHash class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
////Lattice Hashes assemble atom Lists to generate the grid, but each grid is stored
////with residue index information.  This way it can be looked up in a CaList again.
LatticeHash::LatticeHash(caList * target, PdbAtom * templateAtom)
{
	int i = 0, j = 0;
	size = 0;

	///find the number (size) of matching atoms
	for(i = 0; i<target->numAcids; i++){
		aminoAcid * tempAcid = target->acids[i];
		if( subChar(templateAtom->assocAA, translateAAcode(tempAcid->residue_name)) ){
			for(j = 0; j<tempAcid->numAtoms; j++){
				if( translateAtomCode(tempAcid->atoms[j]->atom_name) == translateAtomCode(templateAtom->atom_name) ){
					size++;
					break;
				}
			}
		}
	}

	////get pointers to all matching atoms
	PdbAtom * *newList = new PdbAtom*[size];
	ids = new int[size];/////////////////////////////////////ids[i] takes local index returns global residue index
	backIds = new int[target->numAcids];/////////////////////backIds[i] takes global residue index returns local index
	for(i = 0; i<target->numAcids; i++){ backIds[i] = -1; }

	int tempCount = 0;
	for(i = 0; i<target->numAcids; i++){
		aminoAcid * tempAcid = target->acids[i];
		if( subChar(templateAtom->assocAA, translateAAcode(tempAcid->residue_name)) ){
			for(j = 0; j<tempAcid->numAtoms; j++){
				if( translateAtomCode(tempAcid->atoms[j]->atom_name) == translateAtomCode(templateAtom->atom_name) ){
					newList[tempCount] = tempAcid->atoms[j]->copy();
					ids[tempCount] = i;
					backIds[i] = tempCount;
					tempCount++;
					break;
				}
			}
		}
	}

	///store a copy of all matching atoms.
	list = new AtomList(size, newList);

	set = alloc_set(SP_MAP);
	backSet = alloc_set(SP_MAP);
	
	for(i = 0; i<size; i++){
		set = associate_set(set, i, &ids[i]);
	}
	for(i = 0; i<target->numAcids; i++){
		backSet = associate_set(backSet, i, &backIds[i]);
	}

	grid = new TargetGrid(list, ids);

	for(i = 0; i<size; i++){
		delete(newList[i]);
	}
	delete[](newList);

}

LatticeHash::~LatticeHash()
{
	delete[](ids);
	delete[](backIds);
	delete(grid);
	free_set(set);
	free_set(backSet);
	delete(list);
}

int * LatticeHash::queryShell(double x, double y, double z, double range)
{
	int * result = grid->queryShell(x, y, z, backSet, range);

	return result;
}


set_t LatticeHash::querySphere(double x, double y, double z, double range, set_t output)
{
	output = grid->querySphere(x, y, z, backSet, range, output);

	return output;
}



int LatticeHash::mapBack(int k)
{
	int * info = (int *) mapsto_set(backSet, k);
	int result = -1;

	if(info == NIL){
		return result;
	}
	else{
		result = info[0];
	}

	return result;
}
