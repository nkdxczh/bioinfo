// GATenv.cpp: implementation of the GATenv class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//Implementation Note: Multiple source patterns are not supported in matchXpander
///////As a result, only the first source is ever considered.
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

GATenv::GATenv()
{
	hasAssocArrays = false;
	rotArray = NULL;
	transArray = NULL;
	hashSets = NULL;
	sourceList = NULL;
	matchArray = NULL;
	matchArraySizes = NULL;
	RMSDarray = NULL;
	targetList = NULL;
	motif = NULL;
}

GATenv::~GATenv()
{
	delete(sourceList);
	delete(motif);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////Builds the motif, builds the CaList.////////////////////////////////////////////////////////////////////////////
void GATenv::addSource(AtomList * alist)
{
	///first get the selected Atoms
	int * selectedAtoms = alist->getHighlightedAtomNumbers();

	///then set all atoms in those residues highlighted (for the CA construction
	alist->highlightAssociatedAAs(true);

	///Build the CAList and the motif
	sourceList = new caList(alist);
	motif = new Motif(sourceList, selectedAtoms);

	delete[](selectedAtoms);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///this assumes the source Motif has been set.//////////////////////////////////////////////////////////////////////
void GATenv::setTarget(caList * targ)
{
	int i = 0, j = 0;

	targetList = targ;
	hashSets = new LatticeHash**[motif->numResidues];

	for(i = 0; i<motif->numResidues; i++){
		hashSets[i] = new LatticeHash*[motif->numAtoms[i]];
		for(j = 0; j<motif->numAtoms[i]; j++){
			hashSets[i][j] = new LatticeHash(targetList, motif->atoms[i][j]);
		}
	}
}

void GATenv::setUpMatchArray()
{
	matchArray = new int*[numArrays];
	matchArraySizes = new int[numArrays];

}

void GATenv::setUpMatchArray(int num)
{
	numArrays = num;

	if(num == 0){
		matchArray = NULL;
		matchArraySizes = NULL;
		return;
	}

	matchArray = new int*[num];
	matchArraySizes = new int[num];

	int i = 0;
	for(i = 0; i<num; i++){
		matchArray[i] = NULL;
		matchArraySizes[i] = 0;
	}

}


///requires that setUpMatchArraySizes() was called first.
void GATenv::addMatchArray(int index, int * array)
{
	if(matchArray[index] != NULL){
		delete[](matchArray[index]);
	}
	int i = 0;

	matchArray[index] = new int[2*matchArraySizes[index]];

	for(i = 0; i<2*matchArraySizes[index]; i++){
		matchArray[index][i] = array[i];
	}
}

void GATenv::setUpMatchArraySizes(int i, int val)
{
	if (i<numArrays){
		matchArraySizes[i] = val;
	} 
	else{printf("Attempted to set MatchArray out of bounds");}
}

void GATenv::setUpRMSDArraySizes()
{
	if(numArrays == 0){
		RMSDarray = NULL;
		transArray = NULL;
		rotArray = NULL;
		return;
	}

	int i = 0;
	RMSDarray = new double[numArrays];
	transArray = new double*[numArrays];
	rotArray = new double*[numArrays];

	for (i = 0; i<numArrays; i++)
	{
		RMSDarray[i] = 0.0;
		transArray[i] = new double[3];
			transArray[i][0] = 0.0;
			transArray[i][1] = 0.0;
			transArray[i][2] = 0.0;
		rotArray[i] = new double[9];
			rotArray[i][0] = 0.0;
			rotArray[i][1] = 0.0;
			rotArray[i][2] = 0.0;
			rotArray[i][3] = 0.0;
			rotArray[i][4] = 0.0;
			rotArray[i][5] = 0.0;
			rotArray[i][6] = 0.0;
			rotArray[i][7] = 0.0;
			rotArray[i][8] = 0.0;
	}
}

///this function clears memory in the env so that processing can restart without leaking.
void GATenv::resetEnvironment(void)
{
	///Not Modified
	///int numSources;
	///caList * *lists;
	///bool hasAssocArrays;
	int i = 0, j = 0;

	if(numArrays > 0){
		///stuff that needs to get cleared
		delete(targetList);
//		delete(hashSets);
//		delete(grid);
		delete[](matchArraySizes);
		delete[](RMSDarray);
		matchArraySizes = NULL;
		RMSDarray = NULL;

		///stuff that needs to get emptied
		///FIXME: if we ever deal with multiple sources, seedFinder resets numArrays to be smaller than
		/// its actual allocation.  This prevents null pointers, but if multiple 
		/// sources exist more careful deletion must occur.
		///  this is stymied by checking to see if they are null first, but will not catch everything if there are multiple sources
		int i = 0;
		for(i = 0; i<numArrays; i++){
			if(matchArray[i] != NULL){delete[](matchArray[i]);}
			if(transArray[i] != NULL){delete[](transArray[i]);}
			if(rotArray[i] != NULL){delete[](rotArray[i]);}
		}
		for(i = 0; i<motif->numResidues; i++){
			for(j = 0; j<motif->numAtoms[i]; j++){
				delete(hashSets[i][j]);
			}
			delete[](hashSets[i]);
		}
		delete[](hashSets);
		delete[](matchArray);
		matchArray = NULL;
		delete[](transArray);
		transArray = NULL;
		delete[](rotArray);
		rotArray = NULL;
	}
	else{
		delete(targetList);
//		delete(hashSets);
		for(i = 0; i<motif->numResidues; i++){
			for(j = 0; j<motif->numAtoms[i]; j++){
				delete(hashSets[i][j]);
			}
			delete[](hashSets[i]);
		}
		delete[](hashSets);
//		delete(grid);
		matchArray = NULL;
		transArray = NULL;
		rotArray = NULL;	
		matchArraySizes = NULL;
		RMSDarray = NULL;
	}

	////back to zero matching arrays
	numArrays = 0;

}




void GATenv::sortMatchArrays()
{
	int * indices = new int[numArrays];
	double * rmsds = new double[numArrays];

	int i = 0;
	for(i = 0; i<numArrays; i++){
		indices[i] = i;
		rmsds[i] = RMSDarray[i];
	}

	///now sort by RMSD.
	mergeSort(rmsds, indices, numArrays);

	///now remap the rmsds.
	delete[](RMSDarray);
	RMSDarray = rmsds;
	
	///remap and repoint the other arrays
	double ** transArrayTEMP = new double*[numArrays];
	double ** rotArrayTEMP = new double*[numArrays];
	int * *matchArrayTEMP = new int*[numArrays];
	int * matchArraySizesTEMP = new int[numArrays];

	for(i = 0; i<numArrays; i++){
		transArrayTEMP[i] = transArray[indices[i]];
		rotArrayTEMP[i] = rotArray[indices[i]];
		matchArrayTEMP[i] = matchArray[indices[i]];
		matchArraySizesTEMP[i] = matchArraySizes[indices[i]];
	}

	///clear the old arrays;
	delete[](transArray);
	delete[](rotArray);
	delete[](matchArray);
	delete[](matchArraySizes);

	////replace
	transArray = transArrayTEMP;
	rotArray = rotArrayTEMP;
	matchArray = matchArrayTEMP;
	matchArraySizes = matchArraySizesTEMP;

	delete[](indices);

}


///takes everything below an RMSD cutoff
void GATenv::filterMatchArrays(double rmsd)
{
	int * indices = new int[numArrays];
	double * rmsds = new double[numArrays];
	int stoppingThresh = numArrays;

	int i = 0;
	for(i = 0; i<numArrays; i++){
		indices[i] = i;
		rmsds[i] = RMSDarray[i];
	}

	///now sort by RMSD.
	mergeSort(rmsds, indices, numArrays);

	for(i = 0; i<numArrays; i++){
		if(rmsds[i] > rmsd){
			stoppingThresh = i;
			break;
		}
	}

	///now remap the rmsds.
	delete[](RMSDarray);
	RMSDarray = new double[stoppingThresh];
	for(i = 0; i<stoppingThresh; i++){
		RMSDarray[i] = rmsds[i];
	}
	
	///remap and repoint the other arrays
	double ** transArrayTEMP = new double*[stoppingThresh];
	double ** rotArrayTEMP = new double*[stoppingThresh];
	int * *matchArrayTEMP = new int*[stoppingThresh];
	int * matchArraySizesTEMP = new int[stoppingThresh];

	for(i = 0; i<stoppingThresh; i++){
		transArrayTEMP[i] = transArray[indices[i]];
		rotArrayTEMP[i] = rotArray[indices[i]];
		matchArrayTEMP[i] = matchArray[indices[i]];
		matchArraySizesTEMP[i] = matchArraySizes[indices[i]];
	}
	for(i = stoppingThresh; i<numArrays; i++){
		delete[](transArray[indices[i]]);
		delete[](rotArray[indices[i]]);
		delete[](matchArray[indices[i]]);
	}

	numArrays = stoppingThresh;

	///clear the old arrays;
	delete[](transArray);
	delete[](rotArray);
	delete[](matchArray);
	delete[](matchArraySizes);

	////replace
	transArray = transArrayTEMP;
	rotArray = rotArrayTEMP;
	matchArray = matchArrayTEMP;
	matchArraySizes = matchArraySizesTEMP;

	delete[](rmsds);
	delete[](indices);
}


///takes everything below an absolute cutoff
void GATenv::filterMatchArrays(int numBest)
{
	int * indices = new int[numArrays];
	double * rmsds = new double[numArrays];
	int stoppingThresh = numBest;
	if(stoppingThresh > numArrays){
		stoppingThresh = numArrays;
	}

	int i = 0;
	for(i = 0; i<numArrays; i++){
		indices[i] = i;
		rmsds[i] = RMSDarray[i];
	}

	///now sort by RMSD.
	mergeSort(rmsds, indices, numArrays);

	///now remap the rmsds.
	delete[](RMSDarray);
	RMSDarray = new double[stoppingThresh];
	for(i = 0; i<stoppingThresh; i++){
		RMSDarray[i] = rmsds[i];
	}

//	for(i = 0; i<numArrays; i++){
//		printf("%i\n", indices[i]);
//	}
	
	///remap and repoint the other arrays
	double ** transArrayTEMP = new double*[stoppingThresh];
	double ** rotArrayTEMP = new double*[stoppingThresh];
	int * *matchArrayTEMP = new int*[stoppingThresh];
	int * matchArraySizesTEMP = new int[stoppingThresh];

	for(i = 0; i<stoppingThresh; i++){
		transArrayTEMP[i] = transArray[indices[i]];
		rotArrayTEMP[i] = rotArray[indices[i]];
		matchArrayTEMP[i] = matchArray[indices[i]];
		matchArraySizesTEMP[i] = matchArraySizes[indices[i]];
//
//		printf("\t%i\n", indices[i]);
	}
	for(i = stoppingThresh; i<numArrays; i++){
		delete[](transArray[indices[i]]);
		delete[](rotArray[indices[i]]);
		delete[](matchArray[indices[i]]);
//		printf("\t\t%i\n", indices[i]);
	}

	numArrays = stoppingThresh;

	///clear the old arrays;
	delete[](transArray);
	delete[](rotArray);
	delete[](matchArray);
	delete[](matchArraySizes);

	////replace
	transArray = transArrayTEMP;
	rotArray = rotArrayTEMP;
	matchArray = matchArrayTEMP;
	matchArraySizes = matchArraySizesTEMP;

	delete[](indices);
	delete[](rmsds);
}


set_t GATenv::query(int news, double * accumulationMatrix)
{
	int i = 0;
	set_t result = alloc_set(0);

	///find the union of all residues with potential matches to any motif point in this amino acid.
	int numAtoms = motif->numAtoms[news];

	for(i = 0; i<numAtoms; i++){
		double * queryPt = motif->atoms[news][i]->coords;
		double * newQPt =  transformVector3x4(accumulationMatrix, queryPt);

		result = hashSets[news][i]->querySphere(newQPt[0], newQPt[1], newQPt[2], RANGESEARCHSIZE, result); 
	
//		delete[](queryPt);	////do not delete this to preserve coords info
		delete[](newQPt);
	}

	return result;
}




/*
///DEBUG
void GATenv::resetMatchArrays(void)
{
	///Not Modified
	///int numSources;
	///caList * *lists;
	///bool hasAssocArrays;

	if(numArrays > 0){
		///stuff that needs to get cleared
		delete[](matchArraySizes);
		delete[](RMSDarray);
		matchArraySizes = NULL;
		RMSDarray = NULL;

		///stuff that needs to get emptied
		///FIXME: if we ever deal with multiple sources, seedFinder resets numArrays to be smaller than
		/// its actual allocation.  This prevents null pointers, but if multiple 
		/// sources exist more careful deletion must occur.
		///  this is stymied by checking to see if they are null first, but will not catch everything if there are multiple sources
		int i = 0;
		for(i = 0; i<numArrays; i++){
			if(matchArray[i] != NULL){delete[](matchArray[i]);}
			if(transArray[i] != NULL){delete[](transArray[i]);}
			if(rotArray[i] != NULL){delete[](rotArray[i]);}
		}
		delete[](matchArray);
		matchArray = NULL;
		delete[](transArray);
		transArray = NULL;
		delete[](rotArray);
		rotArray = NULL;
	}
	else{
		matchArray = NULL;
		transArray = NULL;
		rotArray = NULL;	
		matchArraySizes = NULL;
		RMSDarray = NULL;
	}

	////back to zero matching arrays
	numArrays = 0;
}
///DEBUG
*/




