// ActiveList.cpp: implementation of the ActiveList class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ActiveList::ActiveList(GATenv * e)
{
	int i = 0;
	int counter = 0;

	env = e;

	for(i = 0; i<env->numArrays; i++){
		if(env->matchArraySizes[i] > 0){
			counter++;
		}
	}

	size = counter;			////imitially we will fill it with each match array
	capacity = size * 2;			////Capacity doubles when sizes approaches 90% of capacity

	list = new AlignNode*[capacity];

	count = 1;

	///fill the list
	for (i = 0; i<capacity; i++){
		if(i < size){
			list[i] = new AlignNode(count, env->matchArraySizes[i], env);
			////The reason this _can_ be simply "i" is because we do not consider any sources otehr than the 0th.
			list[i]->initializeAccumulationMatrix(env->transArray[i], env->rotArray[i]);
			list[i]->setGeometricPts(env->matchArray[i]);
			list[i]->currentTier = 0;
			list[i]->startingIndex = 0;
			list[i]->RMSD = env->RMSDarray[i];
		
			count++; ////this gives unique identifiers out to new AlignNodes
		}
		else{
			list[i] = NULL;
		}
	}

}


/////This makes ActiveLists preceeding the original one
ActiveList::ActiveList(GATenv * e, int c)
{
	env = e;
	count = c;

	size = 0;							////imitially we will fill it with each match array
	capacity = size * 2 +10;			////Capacity doubles when sizes approaches 90% of capacity, with a constant initial size

	list = new AlignNode*[capacity];

	int i = 0;
	for(i = 0; i<capacity; i++){
		list[i] = NULL;
	}

}


ActiveList::~ActiveList()
{
	///do not delete the env; it deletes you.
///	GATenv * env;
	int i = 0;
	for(i = 0; i<size; i++){
		delete(list[i]);	
	}
	delete[](list);
}

int ActiveList::getNewId()
{
	int result = count;
	count++;

	return result;
}


void ActiveList::expand()
{
	int i = 0;
	capacity = 2*capacity;
	AlignNode * *newNodes = new AlignNode*[capacity];

	for(i = 0; i<capacity; i++){
		if(i<size){
			newNodes[i] = list[i];
		}
		else{
			newNodes[i] = NULL;
		}
	}

	delete[](list);
	list = newNodes;

}


AlignNode * ActiveList::remove(int ind)
{
	int i = 0;
	AlignNode * tempNode;
	AlignNode * result = NULL;
	bool foundit = false;
	int targetIndex = 0;

	for(i = 0; i<size; i++){
		tempNode = list[i];
		if(tempNode == NULL){
			printf("THIS SHOULD NEVER HAPPEN! NULL FOUND DURING REMOVE!\n");
		}
		if(tempNode->id == ind){
			result = tempNode;
			foundit = true;
			targetIndex = i;
			list[i] = NULL;
			break;		
		}
	}

	if(foundit){
		for(i = targetIndex; i<size-1; i++){
			list[i] = list[i+1];
		}
		list[size-1] = NULL;
		size--;
	}

	return result;

}


//This just gives you "any" one, so that you can work with it
AlignNode * ActiveList::remove()
{
	AlignNode * result = list[size-1];
	list[size-1] = NULL;
	size--;
	return result;
}



void ActiveList::addNode(AlignNode * n)
{
	if ((float) size >= 0.7*((float)capacity) ){
		expand();
		list[size] = n;
		size++;
	}
	else{
		list[size] = n;
		size++;
	}


}



