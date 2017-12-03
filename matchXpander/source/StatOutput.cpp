// StatOutput.cpp: implementation of the StatOutput class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
////THESE HAVE BEEN MODIFIED FOR OUTPUTTING GATenv STATS "XSTATS".  DO NOT USE FOR GEOHASH OUTPUT

StatOutput::StatOutput(GATenv * e, HierarchialController * h)
{
	env = e;
	controller = h;
}

StatOutput::~StatOutput()
{

}



///////////////
///Returns the alignment match of the 0th source's biggest (0th) match.
void StatOutput::outputStatFile(char * idHeader)
{
	///Compose the name of the file
	char * filename = new char[NUMCOLUMNS]; 
	sprintf(filename, "%s-output.XSTAT", idHeader);

	///Define the output file.
	currentOutputFile = fopen(filename, "w");

	///print out the hashed "alignment"
	printAlign();

	///print out the statistics
	printStats();

	///print out the match arrays with residue numbers
	printMatches();
	
	///Close the file
	fclose(currentOutputFile);
}




void StatOutput::printAlign()
{
	int i = 0;
	int j = 0;
	int k = 0;
	int tempVal = 0;
	char * outputString;		//this string gets outputted to file
	PdbAtom * tempAtom;
	PdbAtom * tempAtom2;
	int * residueLoc;
	int maxRes;
	int currentResidue;
	PdbAtom * *sourceAtomList;
	PdbAtom * *targetAtomList;
	bool tempFlag;

	residueLoc = new int[env->sourceList->numAcids];

	///print out the name of the source file
	fprintf(currentOutputFile, "Source: %s\n", env->sourceList->name);
	

	//find the highest residue Number
	tempVal = 0;
	for(i = 0; i<env->sourceList->numAcids; i++){
		if(tempVal < env->sourceList->acids[i]->residue_number){
			tempVal = env->sourceList->acids[i]->residue_number;
		}
	}
	for(i = 0; i<env->sourceList->bag->numAtoms; i++){
		tempAtom = env->sourceList->bag->atoms[i];
		if(tempVal < tempAtom->residue_number){
			tempVal = tempAtom->residue_number;
		}
	}
	maxRes = tempVal;

	/////Go through all possible residues as you fill the array
	sourceAtomList = new PdbAtom*[maxRes];
	residueLoc = new int[maxRes];
	currentResidue = 0;
	tempAtom = NULL;
	tempAtom2 = NULL;
	for(i = 0; i<maxRes; i++){
		///search for it here
		for(j = 0; j<env->sourceList->numAcids; j++){
			if(env->sourceList->acids[j]->residue_number == i){
				tempAtom = env->sourceList->acids[j]->atoms[0];
			}
		}
		///search for it here
		for(j = 0; j<env->sourceList->bag->numAtoms; j++){
			if(strcmp(env->sourceList->bag->atoms[j]->atom_name, "CA") == 0
				&& env->sourceList->bag->atoms[j]->residue_number == i){
				tempAtom2 = env->sourceList->bag->atoms[j];
			}
		}
		////should never happen
		if(tempAtom2 != NULL && tempAtom != NULL){
			printf("ERROR: INCONSISTANT SOURCE FILE DETECTED");
			exit(0);
		}
		///put the real one in the array
		if(tempAtom != NULL && tempAtom2 == NULL){
			sourceAtomList[currentResidue] = tempAtom;
			residueLoc[currentResidue] = 1;
			currentResidue++;
		}
		if(tempAtom == NULL && tempAtom2 != NULL){
			sourceAtomList[currentResidue] = tempAtom2;
			residueLoc[currentResidue] = 0;
			currentResidue++;
		}
		if(tempAtom == NULL && tempAtom2 == NULL){
			///do nothing
		}
		tempAtom = NULL;
		tempAtom2 = NULL;
	}
	///now this becomes our total number of residues in sourceAtomList
	maxRes = currentResidue;

	///print 1st source letters
	outputString = new char[10*maxRes+1];
	outputString[0] = '\0';
	fprintf(currentOutputFile, "[");
	for(i = 0; i<maxRes; i++){
		outputString[i] = translateAAcode(sourceAtomList[i]->residue_name);
	}
	outputString[maxRes] = '\0';
	fprintf(currentOutputFile, outputString);
	fprintf(currentOutputFile, "]");
	///now print source hashes
	fprintf(currentOutputFile, "[");
	for(i = 0; i<maxRes; i++){
		if(residueLoc[i] == 1){
			outputString[i] = '|';
		}else{
			outputString[i] = ' ';
		}
	}
	outputString[maxRes] = '\0';
	fprintf(currentOutputFile, outputString);
	fprintf(currentOutputFile, "]\n\n)");
	delete[](outputString);

	delete[](residueLoc);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///now we do the target name
	fprintf(currentOutputFile, "Target: %s\n", env->targetList->name);

	///size the target

	//find the highest residue Number
	tempVal = 0;
	for(i = 0; i<env->targetList->numAcids; i++){
		if(tempVal < env->targetList->acids[i]->residue_number){
			tempVal = env->targetList->acids[i]->residue_number;
		}
	}
	for(i = 0; i<env->targetList->bag->numAtoms; i++){
		tempAtom = env->targetList->bag->atoms[i];
		if(tempVal < tempAtom->residue_number){
			tempVal = tempAtom->residue_number;
		}
	}
	maxRes = tempVal;

	/////Go through all possible residues as you fill the array
	targetAtomList = new PdbAtom*[maxRes];
	residueLoc = new int[maxRes];
	currentResidue = 0;
	tempAtom = NULL;
	tempAtom2 = NULL;
	for(i = 0; i<maxRes; i++){
		///search for it here
		for(j = 0; j<env->targetList->numAcids; j++){
			if(env->targetList->acids[j]->residue_number == i){
				tempAtom = env->targetList->acids[j]->atoms[0];
			}
		}
		///search for it here
		for(j = 0; j<env->targetList->bag->numAtoms; j++){
			if(strcmp(env->targetList->bag->atoms[j]->atom_name, "CA") == 0
				&& env->targetList->bag->atoms[j]->residue_number == i){
				tempAtom2 = env->targetList->bag->atoms[j];
			}
		}
		////should never happen
		if(tempAtom2 != NULL && tempAtom != NULL){
			printf("ERROR: INCONSISTANT TARGET FILE DETECTED");
			exit(0);
		}
		///put the real one in the array
		if(tempAtom != NULL && tempAtom2 == NULL){
			targetAtomList[currentResidue] = tempAtom;
			residueLoc[currentResidue] = 1;
			currentResidue++;
		}
		if(tempAtom == NULL && tempAtom2 != NULL){
			targetAtomList[currentResidue] = tempAtom2;
			residueLoc[currentResidue] = 0;
			currentResidue++;
		}
		if(tempAtom == NULL && tempAtom2 == NULL){
			///do nothing
		}
		tempAtom = NULL;
		tempAtom2 = NULL;
	}
	///now this becomes our total number of residues in sourceAtomList
	maxRes = currentResidue;

	tempFlag = false;
	////each target is a match
	for(k = 0; k<controller->term->size; k++){
		///first print the match #
		fprintf(currentOutputFile , "Match: %i\n", k);
		
		///now print the target residues
		outputString = new char[10*maxRes+10];
		fprintf(currentOutputFile, "[");
		for (i = 0; i<maxRes; i++){
			outputString[i] = translateAAcode(targetAtomList[i]->residue_name);
		}
		outputString[maxRes] = '\0';
		fprintf(currentOutputFile, outputString);
		fprintf(currentOutputFile, "]\n");
	
		///now print the target hashes
		fprintf(currentOutputFile, "[");
		for (i = 0; i<maxRes; i++){
			///find this residue number in the termList's current AlignNode
			for(j = 0; j<controller->term->nodes[k]->size; j++){
				if(targetAtomList[i]->residue_number == 
					env->targetList->acids[ controller->term->nodes[k]->targetIds[j] ]->residue_number ){
					fprintf(currentOutputFile, "|");
					tempFlag = true;
					break;
				}
			}
			////try to find it in the selected target
			for(j = 0; j<env->targetList->numAcids; j++){
				if(targetAtomList[i]->residue_number == env->targetList->acids[j]->residue_number){
					if(!tempFlag){
						fprintf(currentOutputFile, "+");
					}
					tempFlag = true;
					break;
				}
			}
			if(!tempFlag){
				fprintf(currentOutputFile, " ");
			}
			tempFlag = false;			
		}	
		fprintf(currentOutputFile, "]\n\n");
		
	}


}

///prints out the statistics
void StatOutput::printStats()
{
	int i = 0;
	int bestIndex = 0;
	int tempVal = 0;
	double tempRMSD = HUGE_VAL;
	for(i = 0; i<controller->term->size; i++){
		if(tempVal < controller->term->nodes[i]->size){
			tempVal = controller->term->nodes[i]->size;
			bestIndex = i;
			tempRMSD = controller->term->nodes[i]->RMSD;
		}
		if(tempVal == controller->term->nodes[i]->size){
			if(tempRMSD > controller->term->nodes[i]->RMSD){
				bestIndex = i;
				tempRMSD = controller->term->nodes[i]->RMSD;
			}
		}
	}

	//print best match precentage
	double tempVal1 = (double) tempVal;
	double tempVal2 = (double) env->sourceList->numAcids;
	char * tempString = new char[NUMCOLUMNS];
	sprintf(tempString, "Best Structural Match: %lf%% of source pattern\n", 100*(tempVal1/tempVal2));
	fprintf( currentOutputFile, tempString);
	sprintf(tempString, "Best number of structural source matches: %i\n", tempVal);
	fprintf( currentOutputFile, tempString);
	sprintf(tempString, "Best number of structural target matches: %i\n", tempVal);
	fprintf( currentOutputFile, tempString);
	sprintf(tempString, "Source size: %i\n", env->sourceList->numAcids);
	fprintf( currentOutputFile, tempString);
	sprintf(tempString, "Target size: %i\n", env->targetList->numAcids);
	fprintf( currentOutputFile, tempString);
	sprintf(tempString, "Number of unmatched Source residues in best match: %i\n", (env->sourceList->numAcids - tempVal) );
	fprintf( currentOutputFile, tempString);
	/////
	///output RMSD
	if(	controller->term->size > 0 && controller->term->nodes[bestIndex]->rmsdCheck() ){	///this sets RMSD
		sprintf(tempString, "RMSD of best match: %lf\n", controller->term->nodes[bestIndex]->RMSD);
		fprintf( currentOutputFile, tempString);
	}
	else{
		if(controller->term->size == 0){
			sprintf(tempString, "RMSD of best match: N/A (no matches)\n");
			fprintf( currentOutputFile, tempString);
	}
		else{
			sprintf(tempString, "RMSD of best match: %lf (error)\n", controller->term->nodes[bestIndex]->RMSD);
			fprintf( currentOutputFile, tempString);
		}
	}

	/////
	///output execution time
	sprintf(tempString, "Execution time: %lf seconds\n", ( (double)(controller->endTime-controller->startTime) )/CLOCKS_PER_SEC);
	fprintf( currentOutputFile, tempString);

	fprintf( currentOutputFile, "\n\n");
	delete[](tempString);

}


/////
///output match arrays, with residue numbers
void StatOutput::printMatches()
{

	int i = 0;
	int j = 0;
	char * tempString;

	int temp0;
	char temp1;
	int temp2;
	char temp3;

	for (i = 0; i<controller->term->size; i++){
		tempString = new char[NUMCOLUMNS];
		sprintf(tempString, "Source %i, Match %i: ", 0, i );
		fprintf( currentOutputFile, tempString);
		delete[](tempString);

		for (j = 0; j<controller->term->nodes[i]->size; j++){
			tempString = new char[NUMCOLUMNS];

			temp0 = env->sourceList->acids[ controller->term->nodes[i]->sourceIds[j] ]->residue_number;
			temp1 = translateAAcode(env->sourceList->acids[ controller->term->nodes[i]->sourceIds[j] ]->residue_name);
			temp2 = env->targetList->acids[ controller->term->nodes[i]->targetIds[j] ]->residue_number;
			temp3 = translateAAcode(env->targetList->acids[ controller->term->nodes[i]->targetIds[j] ]->residue_name);

			sprintf(tempString, "(%i%c, %i%c) ", 
				env->sourceList->acids[ controller->term->nodes[i]->sourceIds[j] ]->residue_number,
				translateAAcode(env->sourceList->acids[ controller->term->nodes[i]->sourceIds[j] ]->residue_name),
				env->targetList->acids[ controller->term->nodes[i]->targetIds[j] ]->residue_number,
				translateAAcode(env->targetList->acids[ controller->term->nodes[i]->targetIds[j] ]->residue_name)
			);
			fprintf( currentOutputFile, tempString);
			delete[](tempString);
		}
		fprintf(currentOutputFile, "\n");
	}


}




