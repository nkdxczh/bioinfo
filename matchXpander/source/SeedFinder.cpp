
#include "StdAfx.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void findSeedMatches(PriorityStack * stack, caList * target, GATenv * env)
{
	int counter = 0;
	int i = 0, j = 0, k = 0, l = 0, m = 0, n = 0;
	int * tempArray;
	
	////force seed to 3 residues.
	int numResInSeed = 3;

	////describe the source
	int * seedResidueIndices = new int[numResInSeed];	///the source indices of the seed.
	counter = 0;
	for(i = 0; i<stack->height; i++){
		tempArray = stack->flags[i];
		for(j = 0; j<stack->sizes[i]; j++){
			if(counter < numResInSeed){
				seedResidueIndices[counter] = tempArray[j];
				counter++;
			}
			else{ break; }
		}
		if(counter == numResInSeed){
			break;
		}
	}

	///Now generate a lattice based geometric lookup structure for each point in the seed.
	counter = 0;
	LatticeHash * *hashBins = new LatticeHash*[numResInSeed];
	for(i = 0; i<numResInSeed; i++){
		hashBins[i] = new LatticeHash(target, env->motif->atoms[seedResidueIndices[i]][0]);
	}

	///Now generate all pairwise pair-hashes:  should be in O((m^3)(n^2))
	EdgeSet * *edgeHashes = new EdgeSet*[numResInSeed*numResInSeed];
	for(i = 0; i<numResInSeed; i++){
		for(j = 0; j<numResInSeed; j++){
			if(i <= j){continue;}
			else{
				edgeHashes[numResInSeed*i+j] = new EdgeSet();
			}
		}
	}
	
	int * residuesInShell;
	set_t residuesEndingFirst = NULL;
	set_t residuesEndingSecond = NULL;

	PdbAtom * *seedAtoms = new PdbAtom*[3];
	for(i = 0; i<3; i++){
		///we use the first atom from each amino acid in the seed
		seedAtoms[i] = env->motif->atoms[seedResidueIndices[i]][0];
	}

	TriangleList * triList = new TriangleList();
	double tempDistance;
	for(i = 0; i<numResInSeed; i++){
		for(j = 0; j<numResInSeed; j++){
			if(i <= j){
				///process only nonredundantly.
				continue;
			}
			tempDistance = vectorSize(seedAtoms[i]->coords[0] - seedAtoms[j]->coords[0],
									seedAtoms[i]->coords[1] - seedAtoms[j]->coords[1],
									seedAtoms[i]->coords[2] - seedAtoms[j]->coords[2]);

			////as long a we are on the top triangle, then calculate the pair-hash between latticeHashes i and j
			///list each hash
			for(k = 0; k<hashBins[i]->size; k++){
				residuesInShell = hashBins[j]->queryShell(hashBins[i]->list->atomList[k]->coords[0],	///here we can use k bc its all we care about
													hashBins[i]->list->atomList[k]->coords[1],
													hashBins[i]->list->atomList[k]->coords[2],
													tempDistance);
				///for every point that is within the shell, add it to the proper pairHash
				///And check to see if any triangles were generated
				///remember to clear mem for residuesInShell
				for(l = 1; l<residuesInShell[0]+1; l++){
					edgeHashes[numResInSeed*i+j]->addEdge(hashBins[i]->ids[k], residuesInShell[l]);	///add to the proper pairHash

					///now query to see if a triangle was generated
					for(m = 0; m<numResInSeed; m++){
						if(m == i || m == j){ continue; }

						////get a list of the edges connected to the mth hashList for i
						if(m < i){ residuesEndingFirst = edgeHashes[numResInSeed*i+m]->getNeighbors( hashBins[i]->ids[k] ); }
						else{ residuesEndingFirst = edgeHashes[numResInSeed*m+i]->getNeighbors(hashBins[i]->ids[k]); }
						
						////get a list of the edges connected to the mth hashList for j
						if(m < j){ residuesEndingSecond = edgeHashes[numResInSeed*j+m]->getNeighbors(residuesInShell[l]); }
						else{ residuesEndingSecond = edgeHashes[numResInSeed*m+j]->getNeighbors(residuesInShell[l]); }
						
						///if either is empty then theres nothing to do.
						if( (residuesEndingFirst == NULL) || ( residuesEndingSecond == NULL) ){
							if(residuesEndingFirst != NULL){
								free_set(residuesEndingFirst);
							}
							if(residuesEndingSecond != NULL){
								free_set(residuesEndingSecond);
							}
							continue;
						}
						
						///IF both lists are non-empty, copy out all triples by finding points in both lists
						if(size_set(residuesEndingFirst) < size_set(residuesEndingSecond)){
							for(n = 0; n<size_set(residuesEndingFirst); n++){
								if(contains_set(residuesEndingSecond, residuesEndingFirst[n])){
									///only if there are no residue number doubles (applies in the multi-atom-per-residue case)
									if(!( (target->acids[hashBins[i]->ids[k]]->residue_number == target->acids[residuesInShell[l]]->residue_number) || 
										(target->acids[hashBins[i]->ids[k]]->residue_number == target->acids[residuesEndingFirst[n]]->residue_number) || 
										(target->acids[residuesInShell[l]]->residue_number == target->acids[residuesEndingFirst[n]]->residue_number) ))
									{
									//add new triangle
									triList->addTriangle(seedResidueIndices[i], hashBins[i]->ids[k], 
														seedResidueIndices[j], residuesInShell[l], 
														seedResidueIndices[m], residuesEndingFirst[n]);
									}
// DEBUG							triList->addTriangle(seedResidueIndices[i], target->acids[hashBins[i]->mapBack(k)]->Calpha->residue_number, 
// THIS CODE SUBSTITUTES								seedResidueIndices[j], target->acids[residuesInShell[l]]->Calpha->residue_number, 
// INDICES FOR RESIDUE #								seedResidueIndices[m], target->acids[residuesEndingFirst[n]]->Calpha->residue_number);
								}
							}
						}else{
							for(n = 0; n<size_set(residuesEndingSecond); n++){
								if(contains_set(residuesEndingFirst, residuesEndingSecond[n])){
									///only if there are no residue number doubles (applies in the multi-atom-per-residue case)
									if(!( (target->acids[hashBins[i]->ids[k]]->residue_number == target->acids[residuesInShell[l]]->residue_number) || 
										(target->acids[hashBins[i]->ids[k]]->residue_number == target->acids[residuesEndingSecond[n]]->residue_number) || 
										(target->acids[residuesInShell[l]]->residue_number == target->acids[residuesEndingSecond[n]]->residue_number) ))
									{
									//add new triangle
									triList->addTriangle(seedResidueIndices[i], hashBins[i]->ids[k], 
														seedResidueIndices[j], residuesInShell[l], 
														seedResidueIndices[m], residuesEndingSecond[n]);
									}
// DEBUG							triList->addTriangle(seedResidueIndices[i], target->acids[hashBins[i]->mapBack(k)]->Calpha->residue_number, 
// THIS CODE SUBSTITUTES								seedResidueIndices[j], target->acids[residuesInShell[l]]->Calpha->residue_number, 
// INDICES FOR RESIDUE #								seedResidueIndices[m], target->acids[residuesEndingSecond[n]]->Calpha->residue_number);
								}
							}						
						}
						free_set(residuesEndingFirst);
						residuesEndingFirst = NULL;
						free_set(residuesEndingSecond);
						residuesEndingSecond = NULL;
					}
				}
				///Clear residuesINShell now that the edges are added.
				delete[](residuesInShell);
			}
		}
	}

	///debug
	///triList->print();

	env->setUpMatchArray(triList->size);
	env->setUpRMSDArraySizes();
	
	/*////////////////////////////////////////////////
	///////////////////////////// DEBUG: do not delete
	printf("FOUND: %i Seed matches\n", triList->size);
	///////////////////////////// DEBUG: do not delete
	*/////////////////////////////////////////////////

	double * transform = new double[16];
	counter = 0;
	for(i = 0; i<triList->size; i++){
	///First check to see if the triangle is acceptable within bounds
		///copy the match Array into a usable array format (tempArray will be used later)
		double RMSD;
		generateSideChainRMSD(triList->list[i], 3, env->motif, target, &RMSD, transform);

		/*////////////////////////////////////////////////
		///////////////////////////// DEBUG: do not delete
		printf("Match %i, (%i, %i), (%i, %i), (%i, %i)\n", i, triList->list[i][0], target->acids[triList->list[i][1]]->Calpha->residue_number, 
									triList->list[i][2], target->acids[triList->list[i][3]]->Calpha->residue_number, 
									triList->list[i][4], target->acids[triList->list[i][5]]->Calpha->residue_number);
		printf("RMSD: %lf\n", RMSD);
		///////////////////////////// DEBUG: do not delete
		*/////////////////////////////////////////////////

		///set up the match Array Size
		env->setUpMatchArraySizes(counter, 3);						///covered (memory wise)

		///add the match array.& env data
		env->addMatchArray(counter, triList->list[i]);	//tempArray gets absorbed	///covered (memory wise)
		env->RMSDarray[counter] = RMSD;

		env->transArray[counter][0] = transform[12];
		env->transArray[counter][1] = transform[13];
		env->transArray[counter][2] = transform[14];
		env->rotArray[counter][0] = transform[0];	env->rotArray[counter][3] = transform[4];	env->rotArray[counter][6] = transform[8];
		env->rotArray[counter][1] = transform[1];	env->rotArray[counter][4] = transform[5];	env->rotArray[counter][7] = transform[9];
		env->rotArray[counter][2] = transform[2];	env->rotArray[counter][5] = transform[6];	env->rotArray[counter][8] = transform[10];

		counter++;
	}

	///Clear all rot and trans arrays above counter up to triList->size
	for(i = counter; i<triList->size; i++){
		delete[](env->rotArray[i]);
		delete[](env->transArray[i]);
	}

	///reset the size of the matchArrays so that no null entries are deleted.
	env->numArrays = counter;

	////clear memory
	///**do not delete tempArray bc it gets absorbed in every case**
	for(i = 0; i<numResInSeed; i++){
		delete(hashBins[i]);
	}
	delete[](hashBins);
	for(i = 0; i<numResInSeed; i++){
		for(j = 0; j<numResInSeed; j++){
			if(i <= j){continue;}
			delete(edgeHashes[numResInSeed*i+j]);
		}
	}
	delete[](edgeHashes);
	delete(triList);
	delete[](seedResidueIndices);
	delete[](seedAtoms);				//	PdbAtom * *seedAtoms = new PdbAtom*[3];
	delete[](transform);

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
double * generateSideChainRMSD(int * sourceArray, int * targArray, int numRes, Motif * motif, caList * target, double * RMSD, double * transform)
{
	int i = 0, j = 0, k = 0;
	int numPts = 0;
	PdbAtom * tempAtom;
	aminoAcid * tempAcid;

	///sort the match array by source indices
	smallSort(sourceArray, targArray, numRes);

	///first count the number of point pairs
	for(i = 0; i<numRes; i++){
		for(j = 0; j<motif->numAtoms[sourceArray[i]]; j++){
			tempAtom = motif->atoms[sourceArray[i]][j];
			tempAcid = target->acids[targArray[i]];
			for(k = 0; k<tempAcid->numAtoms; k++){
				if(translateAtomCode(tempAcid->atoms[k]->atom_name) == translateAtomCode(tempAtom->atom_name)){
					numPts++;
					break;
				}
			}	
		}
	}

	double * * targetCoords = new double*[numPts];
	double * * sourceCoords = new double*[numPts];

	numPts = 0;

	///set all the coordinates
	for(i = 0; i<numRes; i++){
		for(j = 0; j<motif->numAtoms[sourceArray[i]]; j++){
			tempAtom = motif->atoms[sourceArray[i]][j];
			tempAcid = target->acids[targArray[i]];
			for(k = 0; k<tempAcid->numAtoms; k++){
				if(translateAtomCode(tempAcid->atoms[k]->atom_name) == translateAtomCode(tempAtom->atom_name)){
					targetCoords[numPts] = tempAcid->atoms[k]->coords;
					sourceCoords[numPts] = tempAtom->coords;
					numPts++;
					break;
				}
			}	
		}
	}

	double * resultArray;
	if(transform == NULL){
		resultArray = min_rmsd (targetCoords, sourceCoords, numPts, NULL);
	}
	else{
		resultArray = min_rmsd (targetCoords, sourceCoords, numPts, transform);
	}

	if(RMSD != NULL){
		RMSD[0] = resultArray[0];
	}

//	for(i = 0; i<numPts; i++){				/// do not delete these topreserve PdbAtom information.
//		delete[](targetCoords[i]);
//		delete[](sourceCoords[i]);
//	}
	delete[](targetCoords);
	delete[](sourceCoords);
//	delete[](resultArray);
	return resultArray;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void generateSideChainRMSD(int * matchArray, int numRes, Motif * motif, caList * target, double * RMSD, double * transform)
{
	int i = 0, j = 0, k = 0;
	int numPts = 0;
	PdbAtom * tempAtom;
	aminoAcid * tempAcid;

	///sort the match array by source indices
	int smallest = INT_MAX;
	int smallestInd = 0;
	for(i = 0; i<numRes; i++){
		for(j = i; j<numRes; j++){
			if(i == j){continue;}
			if(matchArray[2*j] <= smallest){
				smallest = matchArray[2*j];
				smallestInd = j;
			}
		}
		int tempSource = matchArray[2*i];
		int tempTarg = matchArray[2*i+1];
		matchArray[2*i] = matchArray[2*smallestInd];
		matchArray[2*i+1] = matchArray[2*smallestInd+1];
		matchArray[2*smallestInd] = tempSource;
		matchArray[2*smallestInd+1] = tempTarg;
		smallest = INT_MAX;
		smallestInd = 0;
	}

	///first count the number of point pairs
	for(i = 0; i<numRes; i++){
		for(j = 0; j<motif->numAtoms[matchArray[2*i]]; j++){
			tempAtom = motif->atoms[matchArray[2*i]][j];
			tempAcid = target->acids[matchArray[2*i+1]];
			for(k = 0; k<tempAcid->numAtoms; k++){
				if(translateAtomCode(tempAcid->atoms[k]->atom_name) == translateAtomCode(tempAtom->atom_name)){
					numPts++;
					break;
				}
			}	
		}
	}

	double * * targetCoords = new double*[numPts];
	double * * sourceCoords = new double*[numPts];

	numPts = 0;

	///set all the coordinates
	for(i = 0; i<numRes; i++){
		for(j = 0; j<motif->numAtoms[matchArray[2*i]]; j++){
			tempAtom = motif->atoms[matchArray[2*i]][j];
			tempAcid = target->acids[matchArray[2*i+1]];
			for(k = 0; k<tempAcid->numAtoms; k++){
				if(translateAtomCode(tempAcid->atoms[k]->atom_name) == translateAtomCode(tempAtom->atom_name)){
					targetCoords[numPts] = tempAcid->atoms[k]->coords;
					sourceCoords[numPts] = tempAtom->coords;
					numPts++;
					break;
				}
			}	
		}
	}

	double * resultArray;
	if(transform == NULL){
		resultArray = min_rmsd (targetCoords, sourceCoords, numPts, NULL);
	}
	else{
		resultArray = min_rmsd (targetCoords, sourceCoords, numPts, transform);
	}

	if(RMSD != NULL){
		RMSD[0] = resultArray[0];
	}

	////Do not delete to preserve the data in the lists.
//	for(i = 0; i<numPts; i++){
//		delete[](targetCoords[i]);
//		delete[](sourceCoords[i]);
//	}
	delete[](targetCoords);
	delete[](sourceCoords);
	delete[](resultArray);

}

int * generateAtomMatchListings( int numSourcePts, PdbAtom * *source, aminoAcid * target)
{
	int i = 0, j = 0;

	int numPts = 0;

	for(i = 0; i<numSourcePts; i++){
		for(j = 0; j<target->numAtoms; j++){
			if( translateAtomCode(source[i]->atom_name) == translateAtomCode(target->atoms[j]->atom_name) ){
				numPts++;
				break;
			}
		}	
	}

	int * results = new int[numPts+1];
	results[0] = numPts;
	numPts = 0;

	for(i = 0; i<numSourcePts; i++){
		for(j = 0; j<target->numAtoms; j++){
			if( translateAtomCode(source[i]->atom_name) == translateAtomCode(target->atoms[j]->atom_name) ){
				results[numPts+1] = translateAtomCode(source[i]->atom_name);
				numPts++;
				break;
			}
		}	
	}

	return results;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//This function aligns three pairs of points, which are index correlated in sourceCoords and targCoords.
//E.g., sourceCoords[0] is to align with targCoords[0], etc.  the function first finds translations to
// equate their centroids, then rotates the patterns so that they are coplanar.  Then, the amount of rotation
// necessary to move each set of points within EXTMATCHTHRESHOLD of each other is determined   THe overall
// rotation matrix (with translatioon) is returned.
//
// if alignment of some pair within EXTMATCHTHRESHOLD is impossible, then a rotation matrix is returned with 
// -1 in the scale position (i.e. result[15])

double alignThreePoints(double ** sourceCoords, double ** targCoords, double * transform)
{
	int i = 0;
	double * tempMatrix = new double[16];
	for(i = 0; i<16; i++){
		tempMatrix[i] = 0;
	}
	tempMatrix[0] = 1;
	tempMatrix[5] = 1;
	tempMatrix[10] = 1;
	tempMatrix[15] = 1;

	double * targCentroid = new double[3];
	targCentroid[0] = (targCoords[0][0] + targCoords[1][0] + targCoords[2][0]) / 3;
	targCentroid[1] = (targCoords[0][1] + targCoords[1][1] + targCoords[2][1]) / 3;
	targCentroid[2] = (targCoords[0][2] + targCoords[1][2] + targCoords[2][2]) / 3;

	double * sourceCentroid = new double[3];
	sourceCentroid[0] = (sourceCoords[0][0] + sourceCoords[1][0] + sourceCoords[2][0]) / 3;
	sourceCentroid[1] = (sourceCoords[0][1] + sourceCoords[1][1] + sourceCoords[2][1]) / 3;
	sourceCentroid[2] = (sourceCoords[0][2] + sourceCoords[1][2] + sourceCoords[2][2]) / 3;
	
	///Generate translations to the centroids.
	double * sourceCenteringTranslation = new double[16];
	double * targCenteringTranslation = new double[16];
	for(i = 0; i<16; i++){ sourceCenteringTranslation[i] = 0.0; targCenteringTranslation[i] = 0.0; }
	sourceCenteringTranslation[0] = 1.0;	targCenteringTranslation[0] = 1.0;
	sourceCenteringTranslation[5] = 1.0;	targCenteringTranslation[5] = 1.0;
	sourceCenteringTranslation[10] = 1.0;	targCenteringTranslation[10] = 1.0;
	sourceCenteringTranslation[15] = 1.0;	targCenteringTranslation[15] = 1.0;

	sourceCenteringTranslation[12] = -sourceCentroid[0];	/////global to source-local
	sourceCenteringTranslation[13] = -sourceCentroid[1];
	sourceCenteringTranslation[14] = -sourceCentroid[2];
	targCenteringTranslation[12] = -targCentroid[0];  /////global to target-local
	targCenteringTranslation[13] = -targCentroid[1];
	targCenteringTranslation[14] = -targCentroid[2];

	///make the patterns coplanar.
	double * sourceAngle1 = new double[3];
	sourceAngle1[0] = sourceCoords[1][0] - sourceCoords[0][0];
	sourceAngle1[1] = sourceCoords[1][1] - sourceCoords[0][1];
	sourceAngle1[2] = sourceCoords[1][2] - sourceCoords[0][2];

	double * sourceAngle2 = new double[3];
	sourceAngle2[0] = sourceCoords[2][0] - sourceCoords[0][0];
	sourceAngle2[1] = sourceCoords[2][1] - sourceCoords[0][1];
	sourceAngle2[2] = sourceCoords[2][2] - sourceCoords[0][2];

	double * sourcePerp = crossProd(sourceAngle1, sourceAngle2);
	double * sourcePerpN = normalizeVector(sourcePerp);
	double * sourceSecond = crossProd(sourcePerp, sourceAngle1);
	double * sourceSecondN = normalizeVector(sourceSecond);
	double * sourceFirstN = normalizeVector(sourceAngle1);

	double * sourceRotArray = new double[9];
	sourceRotArray[0] = sourcePerpN[0]; 	sourceRotArray[1] = sourcePerpN[1]; 	sourceRotArray[2] = sourcePerpN[2]; 
	sourceRotArray[3] = sourceFirstN[0]; 	sourceRotArray[4] = sourceFirstN[1]; 	sourceRotArray[5] = sourceFirstN[2]; 
	sourceRotArray[6] = sourceSecondN[0]; 	sourceRotArray[7] = sourceSecondN[1]; 	sourceRotArray[8] = sourceSecondN[2]; 

	double * invSourceRotArray = transpose(sourceRotArray);
	double * invSourceRotArray4x4 = new double[16];				///inverse moves centroid-centered points to "flat" points.
	for(i = 0; i<16; i++){invSourceRotArray4x4[i] = 0.0;}
	invSourceRotArray4x4[0] = invSourceRotArray[0]; invSourceRotArray4x4[1] = invSourceRotArray[1]; invSourceRotArray4x4[2] = invSourceRotArray[2];
	invSourceRotArray4x4[4] = invSourceRotArray[3]; invSourceRotArray4x4[5] = invSourceRotArray[4]; invSourceRotArray4x4[6] = invSourceRotArray[5];
	invSourceRotArray4x4[8] = invSourceRotArray[6]; invSourceRotArray4x4[9] = invSourceRotArray[7]; invSourceRotArray4x4[10] = invSourceRotArray[8];
	invSourceRotArray4x4[15] = 1.0;

	double * targAngle1 = new double[3];
	targAngle1[0] = targCoords[1][0] - targCoords[0][0];
	targAngle1[1] = targCoords[1][1] - targCoords[0][1];
	targAngle1[2] = targCoords[1][2] - targCoords[0][2];

	double * targAngle2 = new double[3];
	targAngle2[0] = targCoords[2][0] - targCoords[0][0];
	targAngle2[1] = targCoords[2][1] - targCoords[0][1];
	targAngle2[2] = targCoords[2][2] - targCoords[0][2];

	double * targPerp = crossProd(targAngle1, targAngle2);
	double * targPerpN = normalizeVector(targPerp);
	double * targSecond = crossProd(targPerp, targAngle1);
	double * targSecondN = normalizeVector(targSecond);
	double * targFirstN = normalizeVector(targAngle1);

	double * targRotArray = new double[9];
	targRotArray[0] = targPerpN[0]; 	targRotArray[1] = targPerpN[1]; 	targRotArray[2] = targPerpN[2]; 
	targRotArray[3] = targFirstN[0]; 	targRotArray[4] = targFirstN[1]; 	targRotArray[5] = targFirstN[2]; 
	targRotArray[6] = targSecondN[0];	targRotArray[7] = targSecondN[1];	targRotArray[8] = targSecondN[2]; 

	double * invTargRotArray = transpose(targRotArray);
	double * invTargRotArray4x4 = new double[16];				///inverse moves centroid-centered points to "flat" points.
	for(i = 0; i<16; i++){invTargRotArray4x4[i] = 0.0;}
	invTargRotArray4x4[0] = invTargRotArray[0]; invTargRotArray4x4[1] = invTargRotArray[1]; invTargRotArray4x4[2] = invTargRotArray[2];
	invTargRotArray4x4[4] = invTargRotArray[3]; invTargRotArray4x4[5] = invTargRotArray[4]; invTargRotArray4x4[6] = invTargRotArray[5];
	invTargRotArray4x4[8] = invTargRotArray[6]; invTargRotArray4x4[9] = invTargRotArray[7]; invTargRotArray4x4[10] = invTargRotArray[8];
	invTargRotArray4x4[15] = 1.0;
	
	///cleanup everything we dont need anymore
	delete[](targCentroid);
	delete[](sourceCentroid);
	delete[](sourceAngle1);
	delete[](sourceAngle2);
	delete[](sourcePerp);
	delete[](sourcePerpN);
	delete[](sourceSecond);
	delete[](sourceSecondN);
	delete[](sourceFirstN);
	delete[](sourceRotArray);
	delete[](invSourceRotArray);
	delete[](targAngle1);
	delete[](targAngle2);
	delete[](targPerp);
	delete[](targPerpN);
	delete[](targSecond);
	delete[](targSecondN);
	delete[](targFirstN);
	delete[](targRotArray);
	delete[](invTargRotArray);
//	delete[](sourceCenteringTranslation);		///leave these alone for future calculation
//	delete[](targCenteringTranslation);
//	delete[](invSourceRotArray4x4);
//	delete[](invTargRotArray4x4);

	//now we need a rotation about the x-axis that will put the flattened points in one close to the flattened points in the other.
	double * flatteningSourceTrans = matrixMult4x4(invSourceRotArray4x4, sourceCenteringTranslation);
	double * flatteningTargTrans = matrixMult4x4(invTargRotArray4x4, targCenteringTranslation);

	bool sanityCheck = true;
	double ** flatSourceCoords = new double*[3];
	for(i = 0; i<3; i++){
		flatSourceCoords[i] = transformVector3x4(flatteningSourceTrans, sourceCoords[i]);
		if(flatSourceCoords[i][0] < .0000001 && flatSourceCoords[i][0] > -.0000001 ){flatSourceCoords[i][0] = 0;}
		if(flatSourceCoords[i][0] != 0.0){
			sanityCheck = false;
//			printf("weird\n");
		}
	}
	double ** flatTargCoords = new double*[3];
	for(i = 0; i<3; i++){
		flatTargCoords[i] = transformVector3x4(flatteningTargTrans, targCoords[i]);
		if(flatTargCoords[i][0] < .0000001 && flatTargCoords[i][0] > -.0000001 ){flatTargCoords[i][0] = 0;}
		if(flatTargCoords[i][0] != 0.0){
			sanityCheck = false;
//			printf("weird\n");
		}
	}

	if(!sanityCheck){
		for(i = 0; i<3; i++){
			delete[](flatSourceCoords[i]);
			delete[](flatTargCoords[i]);
		}
		delete[](flatSourceCoords);
		delete[](flatTargCoords);

		delete[](flatteningSourceTrans);
		delete[](flatteningTargTrans);
		delete[](invSourceRotArray4x4);
		delete[](invTargRotArray4x4);
		delete[](sourceCenteringTranslation);
		delete[](targCenteringTranslation);
		delete[](tempMatrix);

		for(i = 0; i<16; i++){
			transform[i] = 0.0;
		}
		return -1;
	}

	///generate an array of angle displacements
	double ** angleDisplacements = new double*[3];
	angleDisplacements[0] = new double[2];
	angleDisplacements[1] = new double[2];
	angleDisplacements[2] = new double[2];

	double * tempPerp;
	for(i = 0; i<3; i++){
		tempPerp = crossProd(flatSourceCoords[i], flatTargCoords[i]);
		double dist1 = vectorSize(0, flatSourceCoords[i][1] - flatTargCoords[i][1], flatSourceCoords[i][2] - flatTargCoords[i][2]);
		double dist2 = vectorSize(0, flatTargCoords[i][1], flatTargCoords[i][2]);			////target to origin
		double dist3 = vectorSize(0, flatSourceCoords[i][1], flatSourceCoords[i][2]);		////source to origin
	
		if( (dist2 + dist3) < EXTMATCHTHRESHOLD ){
			angleDisplacements[i][0] = +HUGE_VAL;
			angleDisplacements[i][1] = -HUGE_VAL;
			delete[](tempPerp);
			continue;
		}
		if( dist3 > (dist2 + EXTMATCHTHRESHOLD) || dist3 < (dist2 - EXTMATCHTHRESHOLD) ){
			angleDisplacements[i][0] = -HUGE_VAL;
			angleDisplacements[i][1] = +HUGE_VAL;
			delete[](tempPerp);
			continue;
		}

		if(tempPerp[0] > 0 && tempPerp[1] == 0 && tempPerp[2] == 0){	///if the target is a positive rotation
			if( dist1 < EXTMATCHTHRESHOLD ){		///if inside the circle
				double tempAngle1 = SSStriangleAngle(EXTMATCHTHRESHOLD, dist2, dist3);
				double tempAngle2 = SSStriangleAngle(dist1, dist2, dist3);
				angleDisplacements[i][0] = tempAngle1 + tempAngle2;
				angleDisplacements[i][1] = -(tempAngle1 - tempAngle2);
			}
			else{///if outside the circle
				double tempAngle1 = SSStriangleAngle(EXTMATCHTHRESHOLD, dist2, dist3);
				double tempAngle2 = SSStriangleAngle(dist1, dist2, dist3);
				angleDisplacements[i][0] = tempAngle1 + tempAngle2;
				angleDisplacements[i][1] = tempAngle2 - tempAngle1;
			}
			delete[](tempPerp);
			continue;
		}
		if(tempPerp[0] < 0 && tempPerp[1] == 0 && tempPerp[2] == 0){	///if the target is a positive rotation
			if( dist1 < EXTMATCHTHRESHOLD ){		///if inside the circle
				double tempAngle1 = SSStriangleAngle(EXTMATCHTHRESHOLD, dist2, dist3);
				double tempAngle2 = SSStriangleAngle(dist1, dist2, dist3);
//				angleDisplacements[i][0] = tempAngle1 + tempAngle2;			////on the flip side, the angles must be reversed, and negated	
//				angleDisplacements[i][1] = -(tempAngle1 - tempAngle2);		
				angleDisplacements[i][0] = tempAngle1 - tempAngle2;			
				angleDisplacements[i][1] = -(tempAngle1 + tempAngle2);
			}
			else{///if outside the circle
				double tempAngle1 = SSStriangleAngle(EXTMATCHTHRESHOLD, dist2, dist3);
				double tempAngle2 = SSStriangleAngle(dist1, dist2, dist3);
//				angleDisplacements[i][0] = tempAngle1 + tempAngle2;
//				angleDisplacements[i][1] = tempAngle2 - tempAngle1;
				angleDisplacements[i][0] = -(tempAngle2 - tempAngle1);
				angleDisplacements[i][1] = -(tempAngle1 + tempAngle2);
			}
			delete[](tempPerp);
			continue;
		}
		///if the target is an identical point (happens only for identical structures)
		if(tempPerp[0] == 0 && tempPerp[1] == 0 && tempPerp[2] == 0){	
			double tempAngle = SSStriangleAngle(EXTMATCHTHRESHOLD, dist2, dist3);	///dist2 = dist3 in this case
			angleDisplacements[i][0] = tempAngle;
			angleDisplacements[i][1] = -tempAngle;
			delete[](tempPerp);
			continue;
		}
//		printf("weird\n");
		delete[](tempPerp);
	}

	////now that we have the angle displacements, generate the ideal angle of rotation about the first axis.
	double minRotation = -HUGE_VAL; 
	double maxRotation = HUGE_VAL;
	for(i = 0; i<3; i++){
		if(maxRotation > angleDisplacements[i][0]){ maxRotation = angleDisplacements[i][0];}
		if(minRotation < angleDisplacements[i][1]){ minRotation = angleDisplacements[i][1];}
	}

	if(minRotation > maxRotation){
		//then there is no rotation possible.
		//we return tempMatrix with -1 in the scale position
//		tempMatrix[15] = -1;
		transform[15] = -1;

		///clean up again.
		for(i = 0; i<3; i++){
			delete[](flatSourceCoords[i]);
			delete[](flatTargCoords[i]);
			delete[](angleDisplacements[i]);
		}
		//delete[](tempPerp);	///deleted in loop above
		delete[](flatSourceCoords);
		delete[](flatTargCoords);
		delete[](angleDisplacements);
		delete[](flatteningSourceTrans);
		delete[](flatteningTargTrans);
//		delete[](tempMatrix);
		delete[](sourceCenteringTranslation);		///do not delete these.
		delete[](targCenteringTranslation);
		delete[](invSourceRotArray4x4);
		delete[](invTargRotArray4x4);

		delete[](tempMatrix);
		return -1.0;
	}
	
	////IN RADIANS
	double angleToRotate = (((maxRotation - minRotation)/2) + minRotation) * (PI/180);

	////Generate a Rotation matrix about the x axis to rotate the points into alignment
	double sin_a = sin( angleToRotate );
    double cos_a = cos( angleToRotate );

	double * centerRotation = new double[16];
	centerRotation[0]  = 1;		centerRotation[4]  = 0;			centerRotation[8]  = 0;			centerRotation[12] = 0;
	centerRotation[1]  = 0;		centerRotation[5]  = cos_a;		centerRotation[9]  = -sin_a;		centerRotation[13] = 0;
	centerRotation[2]  = 0;		centerRotation[6]  = sin_a;		centerRotation[10] = cos_a;		centerRotation[14] = 0;
	centerRotation[3]  = 0;		centerRotation[7]  = 0;			centerRotation[11] = 0;			centerRotation[15] = 1;
	
	///Invert (transpose) the Target-side translation matrix flatteningTargTrans to get a full map from source to target
//	double * invFlatteningTargTrans = transpose4x4(flatteningTargTrans);

	///Multiply the Rotation matrix flatteningSourceTrans, the perpendicular rotation, and the inverse flatteningTargTrans 
	double * tempMatrix1 = matrixMult4x4(centerRotation, flatteningSourceTrans);
	double * inv = transpose4x4(invTargRotArray4x4);	///get the inverse rotation to target-local
	double * prod = matrixMult4x4(inv, tempMatrix1);	///multiply it in
	tempMatrix[12] = -targCenteringTranslation[12];			///get the inverse translation to target-local
	tempMatrix[13] = -targCenteringTranslation[13];
	tempMatrix[14] = -targCenteringTranslation[14];
	double * result = matrixMult4x4(tempMatrix, prod);		///multiply it in - this is the answer

	///now copy the values of result over to transform
	for(i = 0; i<16; i++){
		if(result[i] < .00000001 && result[i] > -.00000001){
			result[i] = 0;
		}
		transform[i] = result[i];
	}
	delete[](result);	///this matrix is no longer necessary, as it has been copied out into transform.

	///Now calculate the RMSD of the current alignment
	double tempVal = 0;
	for(i = 0; i<3; i++){
		double * temp = transformVector3x4(transform, sourceCoords[i]);
		double tempDistance = vectorSize(temp[0] - targCoords[i][0], temp[1] - targCoords[i][1], temp[2] - targCoords[i][2] );
		tempVal = tempVal + (tempDistance*tempDistance);
		delete[](temp);
	}

	double RMSD = sqrt(tempVal/3);

	delete[](inv);
	delete[](prod);

	///clean up again.
	for(i = 0; i<3; i++){
		delete[](flatSourceCoords[i]);
		delete[](flatTargCoords[i]);
		delete[](angleDisplacements[i]);
	}
	//delete[](tempPerp);	///deleted in loop above
	delete[](flatSourceCoords);
	delete[](flatTargCoords);
	delete[](angleDisplacements);
	delete[](flatteningSourceTrans);
	delete[](flatteningTargTrans);
//	delete[](invFlatteningTargTrans);
	delete[](tempMatrix1);
	delete[](tempMatrix);
	delete[](centerRotation);
	
	delete[](sourceCenteringTranslation);		///do not delete these.
	delete[](targCenteringTranslation);
	delete[](invSourceRotArray4x4);
	delete[](invTargRotArray4x4);

	return RMSD;

}



