
////////////////////////////////////////////////////////
//Execution_CPP.  This file defines functions which operate MatchXpander
//				Which are called from the main method.
////////////////////////////////////////////////////////
#include "StdAfx.h"

///in this file:
////////void initBatchMXexperiment(char * inputFile, char * batchListFile, bool ind, bool out, bool checkPoint)
////////void initMXExperiment(char * input, char * output, char * targ)


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void initMXExperiment(AtomList * aMotif, char * output, char * targ)
{
	GATenv * env = NULL;
	HierarchialController * controller = NULL;
	pdbParser * pparse = NULL;
	AtomList * tempList = NULL;
	caList * tempCaList = NULL;
	GatOutput * printer = NULL;
	StatOutput * statPrinter = NULL;
	char *pdbFile = NULL, *outputFile = NULL;
	int startTime = 0;
	int endTime = 0;

	////Copy the file names/////////////////////////////////////////////////////////////
//	inputFile = new char[NUMCOLUMNS];
	outputFile = new char[NUMCOLUMNS];
	strcpy(outputFile, output);

	/////If a new target is specified, take down the filename///////////////////////////
	if(targ == NULL){ printf("ERROR: Target file [%s] Not Found.  Exitting.\n", targ); exit(1); }
	else{
		pdbFile = new char[NUMCOLUMNS];                
		strcpy(pdbFile, targ);                      
	}

	////Generate the Operating Environment
	env = new GATenv();

	////If a new Target was specified, parse it now
	pparse = new pdbParser(pdbFile);

	////Start Performance Timing
	startTime = clock();

	////Take parsed data structures
	env->addSource(aMotif);

	////Take Parsed Target data
	tempList = pparse->getAtomList();
	tempList->highlightAllAAs();				///assume all atoms are to be used
	tempCaList = new caList(tempList);
	env->setTarget(tempCaList);
	delete(tempList);

	////Generate the Priority Stack for processing
	PriorityStack * stack = new PriorityStack(env);

	///find the seed matches
	findSeedMatches(stack, env->targetList, env);
	env->hasAssocArrays = true;

	///restrict triangles
	env->filterMatchArrays(NUMBESTTRIANGLES);
//	env->filterMatchArrays(2.0);

	///Start Match Augmentation
	controller = new HierarchialController(env);
	controller->beginTesting();
	controller->term->sortEntries();

	///Output Results
	endTime = clock();
	printf("RUNTIME: %lf\n", ( (double)(endTime-startTime) )/CLOCKS_PER_SEC);
	controller->startTime = startTime;
	controller->endTime = endTime;
	printer = new GatOutput(env, controller);
	printer->outputXGatFiles(outputFile);
	statPrinter = new StatOutput(env, controller);
	statPrinter->outputStatFile(outputFile);				
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////Processing Function for database scanning////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initBatchMXexperiment(AtomList * aMotif, char * batchListFile, bool ind, bool out, bool checkPoint)
{
	GATenv * env;
	char * outputFile;
	char * tempFileName;
	char * tempFile;
	HierarchialController * controller;
	GatOutput * printer;
	StatOutput * statPrinter;

	AtomList * tempList;
	caList * tempCaList;
	int startTime = 0;
	int endTime = 0;
	int i = 0, j = 0, k = 0;

	env = new GATenv();
	env->addSource(aMotif);

	////Throw an error if less than 3 amino acids are used.
	if(env->motif->numResidues < 3){
		printf("ERROR: less than 3 residues used.  Exitting\n");
		exit(1);	
	}

	////If the METAFile already exists and -checkpoint was enabled, quit.
	if(checkPoint == true){
		outputFile = new char[50];
		tempFileName = new char[50];
		strcpy(tempFileName, env->sourceList->name);
		strcpy(outputFile, strtok(tempFileName, "."));
		strcat(outputFile, "-query.META");

		if(fopen(outputFile, "r") != NULL){
			printf("Existing Meta File %s found.  Quitting\n", outputFile);
			exit(1);
		}
		else{
			delete[](outputFile);
			delete[](tempFileName);
		}
	}

	/// now loop through the target lists provided by the BatchList
	pdbParser * parse;
	PriorityStack * stack = new PriorityStack(env);
	BatchList * bList = new BatchList(batchListFile, env->motif);
	printf("##############################################################\n");
	printf("##############################################################\n");
	for(i = 0; i<bList->size; i++){
//		startTime = clock();
		bList->beginTime(clock());

		parse = new pdbParser(bList->fileList[i]);
		tempList = parse->getAtomList();
		tempList->highlightAllAAs();
		delete(parse);

		tempCaList = new caList(tempList);
		env->setTarget(tempCaList);
		delete(tempList);

		printf("Now Processing Target: %s\n", env->targetList->name);

		findSeedMatches(stack, env->targetList, env);
		env->hasAssocArrays = true;

		if(env->numArrays > 0){
		///restrict triangles
		//  env->filterMatchArrays(2.0);
			env->filterMatchArrays(NUMBESTTRIANGLES);

			controller = new HierarchialController(env);
			controller->beginTesting();
			controller->term->sortEntries();

			////prepare output name for results
			outputFile = new char[50];
			tempFileName = new char[50];
			tempFile = new char[50];
			strcpy(tempFileName, env->sourceList->name);
			strcpy(tempFile, env->targetList->name);
			sprintf(outputFile, "%s-%s", strtok(tempFileName, "."), strtok(tempFile, " .,\n\r") );
			///output the individual info
			if(ind){
				endTime = clock();
				printf("RUNTIME: %lf\n", ( (double)(endTime-startTime) )/CLOCKS_PER_SEC);
				controller->startTime = startTime;
				controller->endTime = endTime;

			#ifndef NO_OUTPUT////////////////////////////////////////////////////////////////////////////////////////////////
				printf("serializing results to disk...");
				printer = new GatOutput(env, controller);
				printer->outputXGatFiles(outputFile);
				delete(printer);
				printf("\n...done\n\n");

				printf("outputting statistics to disk...");
				statPrinter = new StatOutput(env, controller);
				statPrinter->outputStatFile(outputFile);
				delete(statPrinter);						
				printf("done\n\n");

			#else /////this version has no printfs
				printer = new GatOutput(env, controller);
				printer->outputXGatFiles(outputFile);
				delete(printer);
				statPrinter = new StatOutput(env, controller);
				statPrinter->outputStatFile(outputFile);								
				delete(statPrinter);
			#endif///////////////////////////////////////////////////////////////////////////////////////////////////////////
			}

			///Save the results in the bList
			if(env->matchArray[0][0] == env->sourceList->numAcids){
				bList->numMatchedPerfect++;
			}
			bList->avgPointsMatched += (double) ((double)controller->term->nodes[0]->size/(double)bList->size);
			bList->avgMatchRMSD += (double) ((double)controller->term->nodes[0]->RMSD/(double)bList->size);

			///Save results in BatchList arrays, for stat output
			bList->expList[i] = new char[strlen(outputFile)+1];
			sprintf(bList->expList[i], "%s", outputFile);
			bList->matchSizes[i] = new char[15];
			sprintf(bList->matchSizes[i], "%i/%i(%i)", controller->term->nodes[0]->size, env->sourceList->numAcids, env->targetList->numAcids);
			bList->matchRMSDs[i] = new char[15];
			sprintf(bList->matchRMSDs[i], "%lf A", controller->term->nodes[0]->RMSD);

			bList->matchListings[i] = new char[NUMCOLUMNS];
			for (j = 0; j<controller->term->nodes[0]->size; j++){
				char * tempString = new char[NUMCOLUMNS];
				sprintf(tempString, "(%i%c, %i%c)", 
					env->sourceList->acids[ controller->term->nodes[0]->sourceIds[j] ]->residue_number, 
					translateAAcode(env->sourceList->acids[ controller->term->nodes[0]->sourceIds[j] ]->residue_name),
					env->targetList->acids[ controller->term->nodes[0]->targetIds[j] ]->residue_number, 
					translateAAcode(env->targetList->acids[ controller->term->nodes[0]->targetIds[j] ]->residue_name)
				);
				/////state matching Atoms
				int * atoms = generateAtomMatchListings(env->motif->numAtoms[ controller->term->nodes[0]->sourceIds[j] ], 
					env->motif->atoms[ controller->term->nodes[0]->sourceIds[j] ],
					env->targetList->acids[ controller->term->nodes[0]->targetIds[j] ]);
				for(k = 0; k<atoms[0]; k++){
					char * extraTempString = encodeAtomCode(atoms[k+1]);
					sprintf(tempString, "%s %s,", tempString, extraTempString);
					delete[](extraTempString);
				}
				delete[](atoms);

				if(j == 0){ strcpy(bList->matchListings[i], tempString); }
				else{ strcat(bList->matchListings[i], tempString); }
				delete[](tempString);
			}

			env->resetEnvironment();
			delete(controller);
			delete[](outputFile);
			delete[](tempFileName);
			delete[](tempFile);
		}
		else{	///else there was nothing found at all.
			///generate outputFile name
			outputFile = new char[50];
			tempFileName = new char[50];
			tempFile = new char[50];
			strcpy(tempFileName, env->sourceList->name);
			strcpy(tempFile, env->targetList->name);
			strcpy(outputFile, strtok(tempFileName, "."));
			strcat(outputFile, "-");
			strcat(outputFile, strtok(tempFile, " .,\n\r"));

			///Save results in BatchList arrays, for stat output
			bList->expList[i] = new char[strlen(outputFile)+1];
			strcpy(bList->expList[i], outputFile);
			bList->matchSizes[i] = new char[15];
			sprintf(bList->matchSizes[i], "       ");
			bList->matchRMSDs[i] = new char[10];
			sprintf(bList->matchRMSDs[i], "N/A");
			bList->matchListings[i] = new char[NUMCOLUMNS];
			sprintf(bList->matchListings[i], "        ");
		
			env->resetEnvironment();
			delete[](outputFile);
			delete[](tempFileName);
			delete[](tempFile);
		}

		bList->updateTime(clock());
	}

	///output the statistics file
	if(out){
		outputFile = new char[50];
		tempFileName = new char[50];
		strcpy(tempFileName, env->sourceList->name);
		strcpy(outputFile, strtok(tempFileName, "."));
		strcat(outputFile, "-query.META");

		printf("Printing Global Statistics\n");
		bList->printOutMetaStats(outputFile);

		delete[](outputFile);
		delete[](tempFileName);
	}

	printf("##############################################################\n");
	printf("##############################################################\n");
	printf("Processing Complete.\n");
//	printf("Total Runtime: %i\n", fullEndTime-fullStartTime);

	bList->printTime();
	delete(bList);
	delete(stack);
	delete(env);

//	fullEndTime = clock();

}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///Bin File Batch IMPLEMENTATION of MATCHXPANDER//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initBinBatchExperiment(AtomList * aMotif, char * binFile, char * binToc, char * outputFileName)
{
	GATenv * env;
	char * outputFile;
	char * tempFileName;
	char * tempFile;
	HierarchialController * controller;

//	AtomList * tempList;
	caList * tempCaList;
//	int startTime = 0;
//	int endTime = 0;
//	int fullStartTime = clock();
//	int fullEndTime = 0;
	int i = 0, j = 0, k = 0;

	env = new GATenv();
	env->addSource(aMotif);

	////Read the TOC for the Bin File
	int binFileSize = 0;	///number of bytes the binFile is
	int numFiles = 0;
	int * TOC = parseTOC(binToc, &binFileSize, &numFiles);

	////Read the BinFile
	char * pdbArray = loadBinFile(binFile, binFileSize);

	////Generate the BatchList
	BatchList * bList = new BatchList(env->motif, numFiles);

	/// now loop through the target lists provided by the BatchList
	PriorityStack * stack = new PriorityStack(env);
	printf("##############################################################\n");
	printf("##############################################################\n");
	for(i = 0; i<numFiles; i++){
//		startTime = clock();
		bList->beginTime(clock());

		tempCaList = generateCaListFromBinary(pdbArray, TOC[i]);
		env->setTarget(tempCaList);

		printf("Now Processing Target: %s\n", env->targetList->name);

		findSeedMatches(stack, env->targetList, env);
		env->hasAssocArrays = true;

		if(env->numArrays > 0){
		///restrict triangles
		//  env->filterMatchArrays(2.0);
			env->filterMatchArrays(NUMBESTTRIANGLES);

			controller = new HierarchialController(env);
			controller->beginTesting();
			controller->term->sortEntries();

			////prepare output name for results
			outputFile = new char[50];
			tempFileName = new char[50];
			tempFile = new char[50];
			strcpy(tempFileName, env->sourceList->name);
			strcpy(tempFile, env->targetList->name);
			sprintf(outputFile, "%s-%s", strtok(tempFileName, "."), strtok(tempFile, " .,\n\r") );
//			///output the individual info
//			if(ind){
//				endTime = clock();
//				printf("RUNTIME: %lf\n", ( (double)(endTime-startTime) )/CLOCKS_PER_SEC);
//				controller->startTime = startTime;
//				controller->endTime = endTime;
//
//			#ifndef NO_OUTPUT////////////////////////////////////////////////////////////////////////////////////////////////
//				printf("serializing results to disk...");
//				printer = new GatOutput(env, controller);
//				printer->outputXGatFiles(outputFile);
//				delete(printer);
//				printf("\n...done\n\n");
//
//				printf("outputting statistics to disk...");
//				statPrinter = new StatOutput(env, controller);
//				statPrinter->outputStatFile(outputFile);
//				delete(statPrinter);						
//				printf("done\n\n");
//
//			#else /////this version has no printfs
//				printer = new GatOutput(env, controller);
//				printer->outputXGatFiles(outputFile);
//				delete(printer);
//				statPrinter = new StatOutput(env, controller);
//				statPrinter->outputStatFile(outputFile);								
//				delete(statPrinter);
//			#endif///////////////////////////////////////////////////////////////////////////////////////////////////////////
//			}

			///Save the results in the bList
			if(env->matchArray[0][0] == env->sourceList->numAcids){ bList->numMatchedPerfect++; }
			bList->avgPointsMatched += (double) ((double)controller->term->nodes[0]->size/(double)bList->size);
			bList->avgMatchRMSD += (double) ((double)controller->term->nodes[0]->RMSD/(double)bList->size);

			///Save results in BatchList arrays, for stat output
			bList->expList[i] = new char[strlen(outputFile)+1];
			sprintf(bList->expList[i], "%s", outputFile);
			bList->matchSizes[i] = new char[15];
			sprintf(bList->matchSizes[i], "%i/%i(%i)", controller->term->nodes[0]->size, env->sourceList->numAcids, env->targetList->numAcids);
			bList->matchRMSDs[i] = new char[15];
			sprintf(bList->matchRMSDs[i], "%lf A", controller->term->nodes[0]->RMSD);

			bList->matchListings[i] = new char[NUMCOLUMNS];
			for (j = 0; j<controller->term->nodes[0]->size; j++){
				char * tempString = new char[NUMCOLUMNS];
				sprintf(tempString, "(%i%c, %i%c)", 
					env->sourceList->acids[ controller->term->nodes[0]->sourceIds[j] ]->residue_number, 
					translateAAcode(env->sourceList->acids[ controller->term->nodes[0]->sourceIds[j] ]->residue_name),
					env->targetList->acids[ controller->term->nodes[0]->targetIds[j] ]->residue_number, 
					translateAAcode(env->targetList->acids[ controller->term->nodes[0]->targetIds[j] ]->residue_name)
				);
				/////state matching Atoms
				int * atoms = generateAtomMatchListings(env->motif->numAtoms[ controller->term->nodes[0]->sourceIds[j] ], 
					env->motif->atoms[ controller->term->nodes[0]->sourceIds[j] ],
					env->targetList->acids[ controller->term->nodes[0]->targetIds[j] ]);
				for(k = 0; k<atoms[0]; k++){
					char * extraTempString = encodeAtomCode(atoms[k+1]);
					sprintf(tempString, "%s %s,", tempString, extraTempString);
					delete[](extraTempString);
				}
				delete[](atoms);

				if(j == 0){ strcpy(bList->matchListings[i], tempString); }
				else{ strcat(bList->matchListings[i], tempString); }
				delete[](tempString);
			}

			env->resetEnvironment();
			delete(controller);
			delete[](outputFile);
			delete[](tempFileName);
			delete[](tempFile);
		}
		else{	///else there was nothing found at all.
			///generate outputFile name
			outputFile = new char[50];
			tempFileName = new char[50];
			tempFile = new char[50];
			strcpy(tempFileName, env->sourceList->name);
			strcpy(tempFile, env->targetList->name);
			strcpy(outputFile, strtok(tempFileName, "."));
			strcat(outputFile, "-");
			strcat(outputFile, strtok(tempFile, " .,\n\r"));

			///Save results in BatchList arrays, for stat output
			bList->expList[i] = new char[strlen(outputFile)+1];
			strcpy(bList->expList[i], outputFile);
			bList->matchSizes[i] = new char[15];
			sprintf(bList->matchSizes[i], "       ");
			bList->matchRMSDs[i] = new char[10];
			sprintf(bList->matchRMSDs[i], "N/A");
			bList->matchListings[i] = new char[NUMCOLUMNS];
			sprintf(bList->matchListings[i], "        ");
		
			env->resetEnvironment();
			delete[](outputFile);
			delete[](tempFileName);
			delete[](tempFile);
		}

		bList->updateTime(clock());
	}

	///output the statistics file
	printf("Printing Global Statistics\n");
	bList->printOutMetaStats(outputFileName);

	printf("##############################################################\n");
	printf("##############################################################\n");
	printf("Processing Complete.\n");
//	printf("Total Runtime: %i\n", fullEndTime-fullStartTime);

	bList->printTime();
	delete(bList);
	delete(stack);
	delete(env);

//	fullEndTime = clock();

}









//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///MPI CODE IMPLEMENTATION AND CONNECTION TO MATCHXPANDER/////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
////calculator for this work
double initMXExperimentMPI(char * target, AtomList * motif)
{
	int startTime = 0;
	int endTime = 0;
	int i = 0;
	AtomList * tempList;
	caList * tempCaList;
	HierarchialController * controller;
	GATenv * env = new GATenv();;
	pdbParser * pparse = new pdbParser(target);

	////If a new Target was specified, parse it now
	AtomList * newAtomList = pparse->getAtomList();
	tempList = motif;
	for (i = 0;i<tempList->size;i++){
		sprintf(tempList->atomList[i]->fileLine,"artificial fileLine inserted");
	}
	env->addSource(tempList);

	delete(tempList);
	////Set the target.  if targ is not NULL, use the file provided.
	tempList = newAtomList;//pparse->getAtomList();
	for (i = 0;i<tempList->size;i++){
		if (tempList->atomList[i]->assocAA != NULL){
			delete[](tempList->atomList[i]->assocAA);
			tempList->atomList[i]->assocAA=NULL;
		}
	}

	tempList->highlightAllAAs();				///all atoms are to be used
	tempCaList = new caList(tempList);
	env->setTarget(tempCaList);
	delete(tempList);

	////Generate the Priority Stack for processing
	PriorityStack * stack = new PriorityStack(env);
	///find the seed matches
	findSeedMatches(stack, env->targetList, env);
	env->hasAssocArrays = true;
	///restrict triangles
	env->filterMatchArrays(NUMBESTTRIANGLES);
	if (env->numArrays == 0){
		delete(pparse);
		delete(stack);
		delete(env);
		return -1.00;
	}

	///Start Match Augmentation
	controller = new HierarchialController(env);
	controller->beginTesting();
	controller->term->sortEntries();

	///Output Results
	endTime = clock();
	controller->startTime = startTime;
	controller->endTime = endTime;
	  double finalResult;
	  if(controller->term->nodes[0]->size == env->motif->numResidues){
	    finalResult = controller->term->nodes[0]->RMSD;
	      }
	  else{
	    finalResult = -1;
	  }
	////clean up:
	env->resetEnvironment();
	delete(pparse);
	delete(stack);
	delete(env);
	delete(controller);

	return finalResult;

}
*/

