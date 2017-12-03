// PriorityStack.cpp: implementation of the PriorityStack class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

PriorityStack::PriorityStack(GATenv * e)
{
	env = e;
	currentLevel = 0;


	int i = 0;
	int j = 0;

	///first parse the evolutionary data to get the teir info
	///Get this info from the first source pattern;
	caList * tempList = env->sourceList;
	int num = tempList->numAcids;
	int * temp = new int[num];
	int * tempIndices = new int[num];
	int tally = 1;
	int tally2 = 0;
	int tempVal;
	int * *indexBuffer = new int*[num];
	int * tempCrap = new int[num];
	sizes = new int[num];

	for(i = 0; i<num; i++){
		temp[i] = tempList->acids[i]->significance;		///thsi is the significance
		tempIndices[i] = i;										///this gives you the caList index;
	}

	//now sort them
	smallSort(temp, tempIndices, num);

	///count the number of teirs;
	tempVal = temp[0];
	for(i = 0; i<num; i++){
		if( tempVal < temp[i]){
			tempVal = temp[i];
			///copy tempCrap into an adequately sized array, since we've just finished off a tier;
			indexBuffer[tally-1] = new int[tally2];
			for(j = 0; j<tally2; j++){
				indexBuffer[tally-1][j] = tempCrap[j];
			}
			///remember the size
			sizes[tally-1] = tally2;
			tally++;
			tally2 = 0;

			///Store this tempIndex too
			tempCrap[tally2] = tempIndices[i];
			tally2++;
		}
		else{
			///store the index
			tempCrap[tally2] = tempIndices[i];
			tally2++;
		}
	}
	///copy tempCrap into an adequately sized array;  (this gets the last set, since we've just finished off a tier)
	indexBuffer[tally-1] = new int[tally2];
	for(j = 0; j<tally2; j++){
		indexBuffer[tally-1][j] = tempCrap[j];
		sizes[tally-1] = tally2;
	}

	////Set the number of tiers;
	height = tally;

	////Now add these pointers into the PriorityStack
	flags = new int*[tally];
	for(i = 0; i<height; i++){
		flags[i] = indexBuffer[i];
	}

	delete[](indexBuffer);
	delete[](temp);
	delete[](tempIndices);
	delete[](tempCrap);

}

PriorityStack::~PriorityStack()
{
	delete[](sizes);
	int i = 0;
	for(i = 0; i<height; i++){
		delete[](flags[i]);
	}
	delete[](flags);

	/////do not delete env - it deletes you.
	////GATenv * env;
}



void PriorityStack::nextLevel()
{
	currentLevel++;
}


int * PriorityStack::currentFlags()
{
	return flags[currentLevel];


}

int PriorityStack::currentSize()
{
	return sizes[currentLevel];

}


void PriorityStack::resetStack()
{
	int i = 0;
	int j = 0;
	for(i = 0; i<height; i++){
		for(j = 0; j<sizes[i]; j++){
			flags[i][j] = abs(flags[i][j]);
		}
	}

}


bool PriorityStack::changeLevel(int i)
{
	if(i < height && i >= 0){
		currentLevel = i;
		return true;
	}
	return false;

}

void PriorityStack::printStatus()
{
	int i = 0;
	printf("sizes: ");
	for(i = 0; i<height; i++){
		printf("%i ", sizes[i]);
	}
	printf("\n");

}



