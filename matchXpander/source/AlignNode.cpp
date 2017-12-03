// AlignNode.cpp: implementation of the AlignNode class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

AlignNode::AlignNode(int ident, int s,  GATenv * e)
{
	size = s;
	isdead = false;
	env = e;
	id = ident;

//	printf("%i\n", id);	////debug

	sourceIds = new int[s];
	targetIds = new int[s];
	currentTier = 0;
	startingIndex = 0;

	int i = 0;
	tempMatrix = new double[16];
	for(i = 0; i<15; i++){
		tempMatrix[i] = 0.0;
	}
	tempMatrix[0] = 1.0;
	tempMatrix[5] = 1.0;
	tempMatrix[10] = 1.0;
	tempMatrix[15] = 1.0;

	accumulationMatrix = NULL;

//	if(id  == 117){
//		printf("OMFG!\n");
//	}


}

AlignNode::~AlignNode()
{
//	printf("\t%i\n", id);
	delete[](sourceIds);
	delete[](targetIds);
	delete[](tempMatrix);

	if(accumulationMatrix != NULL){
		delete[](accumulationMatrix);
	}

}


///only call this after both trans and rot are set.
void AlignNode::initializeAccumulationMatrix(double * trans, double * rot){	
	
	int i = 0;

	///debug
	if(accumulationMatrix != NULL){
		delete[](accumulationMatrix);
	}
	///debug

	accumulationMatrix = new double[16];
	for(i = 0; i<16; i++){
		accumulationMatrix[i] = 0.0;	
	}
	accumulationMatrix[0] = rot[0];
	accumulationMatrix[1] = rot[1];
	accumulationMatrix[2] = rot[2];
	accumulationMatrix[4] = rot[3];
	accumulationMatrix[5] = rot[4];
	accumulationMatrix[6] = rot[5];
	accumulationMatrix[8] = rot[6];
	accumulationMatrix[9] = rot[7];
	accumulationMatrix[10] = rot[8];
	accumulationMatrix[12] = trans[0];
	accumulationMatrix[13] = trans[1];
	accumulationMatrix[14] = trans[2];
	accumulationMatrix[15] = 1.0;
	
}


///only call this after both trans and rot are set.
void AlignNode::initializeAccumulationMatrix(double * accum)
{
	int i = 0;
	///debug
	if(accumulationMatrix != NULL){
		delete[](accumulationMatrix);
	}
	///debug

	accumulationMatrix = new double[16];
	for(i = 0; i<16; i++){
		accumulationMatrix[i] = accum[i];
		tempMatrix[i] = 0.0;
	}
	tempMatrix[0] = 1.0;
	tempMatrix[5] = 1.0;
	tempMatrix[10] = 1.0;
	tempMatrix[15] = 1.0;

}



////Tkaes a match array from the GAT file and inserts into the data struct
////CAll this only after initializeAccumulationMatrix has been called.
void AlignNode::setGeometricPts(int * matchArray)
{
	int i = 0;
	for(i = 0; i<size; i++){
		sourceIds[i] = matchArray[2*i+0];
		targetIds[i] = matchArray[2*i+1];
	}

}



////If you make a copy, you need to assign it an ID
AlignNode * AlignNode::copy(int id)
{
	AlignNode * result = new AlignNode(id, size, env);

	result->initializeAccumulationMatrix(accumulationMatrix);

	int i = 0;
	for(i = 0; i<size; i++){
		result->sourceIds[i] = sourceIds[i];
		result->targetIds[i] = targetIds[i];	
	}

	result->startingIndex = startingIndex;
	result->currentTier = currentTier;

	for(i = 0; i<16; i++){
		result->tempMatrix[i] = tempMatrix[i];
	}
	
	result->RMSD = RMSD;

	return result;
}


//////////////////////////////////////////////////////////
////Specify the index of the source you want to test, and the index of the target you want to test.
////Returns the RMSD if all points are within EXTMATCHTHRESHOLD.  Returns -1 if not.
double AlignNode::testPt(int source, int target)
{
	double result = -1.0;
	int i = 0;

	/*/////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////DEBUG  TESTING do not delete
	printf("now testing:");
	for(i = 0; i<size; i++){
		printf("(%i, %i) ", env->sourceList->acids[ sourceIds[i] ]->Calpha->residue_number, env->targetList->acids[ targetIds[i] ]->Calpha->residue_number);
	}
	printf("(%i, %i)\n", env->sourceList->acids[ source ]->Calpha->residue_number, env->targetList->acids[ target ]->Calpha->residue_number);
	///////////////////////////////////////////////DEBUG  TESTING do not delete
	*//////////////////////////////////////////////////////////////////////////

	for(i = 0; i<size; i++){
		if(source == sourceIds[i] || target == targetIds[i]){  return -1.0;  }
	}

	int * newSources = new int[size+1];
	int * newTargets = new int[size+1];

	for(i = 0; i<size; i++){
		newSources[i] = sourceIds[i];
		newTargets[i] = targetIds[i];
	}
	newSources[size] = source;
	newTargets[size] = target;

	///calculate RMSD
	double * distances = generateSideChainRMSD(newSources, newTargets, size+1, env->motif, env->targetList, NULL, tempMatrix);

	///return the RMSD if all goes well
	result = distances[0];

	/*/////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////DEBUG  TESTING do not delete
	printf("RMSD: %lf\n", result);
	///////////////////////////////////////////////DEBUG  TESTING do not delete
	*//////////////////////////////////////////////////////////////////////////

	////if a pair of points is too far away, then return -1.
	for(i = 0; i<size+1; i++){
		if(distances[i+1] > EXTMATCHTHRESHOLD){
			result = -1.0;
			break;
		}
	}
	
	delete[](distances);
	delete[](newSources);
	delete[](newTargets);
	return result;
}


//////////////////////////////////
///this doesnt multiply in the tempMatrix while testing, whereas the sanityCheck does.
bool AlignNode::rmsdCheck(void)
{
	double * temp = generateSideChainRMSD(sourceIds, targetIds, size, env->motif, env->targetList, &RMSD, tempMatrix);
	delete[](temp);

	return true;
}


//////////////////////////////////
/////Adds the index of a pt and a hole, incremets size, 
void AlignNode::addGeometricPt(int peg, int hole)
{
	int i = 0;
	int * newArray1 = new int[size+1];
	int * newArray2 = new int[size+1];

	for(i = 0; i<size; i++){
		if(sourceIds[i] != peg){ newArray1[i] = sourceIds[i]; }
		else{
			delete[](newArray1);
			delete[](newArray2);
			return;
		}
		if(targetIds[i] != hole){ newArray2[i] = targetIds[i]; }
		else{
			delete[](newArray1);
			delete[](newArray2);
			return;
		}
	}

	newArray1[size] = peg;
	newArray2[size] = hole;

	delete[](sourceIds);
	delete[](targetIds);

	sourceIds = newArray1;
	targetIds = newArray2;

	for(i = 0; i<16; i++){
		accumulationMatrix[i] = tempMatrix[i];
		tempMatrix[i] = 0;
	}
	tempMatrix[0] = 1.0;
	tempMatrix[5] = 1.0;
	tempMatrix[10] = 1.0;
	tempMatrix[15] = 1.0;

	size++;
}


//////////////////////////////////
///This is true if val is in source IDs.  False if not.
bool AlignNode::element(int val)
{
	int i = 0;

	for(i = 0; i<size; i++){
		if(sourceIds[i] == val){
			return true;
		}
	}

	return false;
}


//////////////////////////////////
////printf's results of this alignment node.
void AlignNode::printResults(void)
{
	int i = 0;

	getBestFitRMSD();

	printf("ALIGN NODE: %i\n", id);
	printf("SIZE: %i\n", size);
	printf("RMSD: %lf\n", RMSD);
	printf("Source: ");
	for(i = 0; i<size; i++){
		printf("%i ", sourceIds[i]);
	}
	printf("\n");

	printf("Target: ");
	for(i = 0; i<size; i++){
		printf("%i ", targetIds[i]);
	}
	printf("\n");

	printf("Matrix: \n");
	printMatrix(accumulationMatrix, 16);

	printf("\n");
}


///outputs accumulated translation/rotation info
double * AlignNode::getTotalRotation(void)
{
	double * result = new double[9];
	result[0] = accumulationMatrix[0];
	result[1] = accumulationMatrix[1];
	result[2] = accumulationMatrix[2];
	result[3] = accumulationMatrix[4];
	result[4] = accumulationMatrix[5];
	result[5] = accumulationMatrix[6];
	result[6] = accumulationMatrix[8];
	result[7] = accumulationMatrix[9];
	result[8] = accumulationMatrix[10];

	return result;
}


///outputs accumulated translation/rotation info
double * AlignNode::getTotalTranslation(void)
{
	double * result = new double[3];
	result[0] = accumulationMatrix[12];
	result[1] = accumulationMatrix[13];
	result[2] = accumulationMatrix[14];

	return result;
}


//////////////////////////////////
///////////PROBABLY UNNECESSARY
double AlignNode::getBestFitRMSD(void)
{
	double result;	///where the actual RMSD will go
	double * results = generateSideChainRMSD(sourceIds, targetIds, size, env->motif, env->targetList, &RMSD, accumulationMatrix);

	///set the results
	result = results[0];
	delete[](results);

	RMSD = result;

	return (double) result;
}


void AlignNode::sortMatchIds(void)
{
	smallSort(sourceIds, targetIds, size);

}


/*
////this is packed so that from zero to node->size are thecoords, and the queryPt coords
//// are stored in node->size+1
double * * AlignNode::genCoords(int queryPt)
{
	int i = 0;
	double * * result = new double*[size+1];
	
	for(i = 0; i<size; i++){
		result[i] = transformVector3x4(accumulationMatrix, env->sourceList->acids[sourceIds[i]]->Calpha->coords);
	}
	
	result[size] = transformVector3x4(accumulationMatrix, env->sourceList->acids[queryPt]->Calpha->coords);
	return result;
}
*/
