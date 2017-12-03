// TermList.cpp: implementation of the TermList class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TermList::TermList(int c, GATenv * e)
{
	env = e;
	size = 0;
	capacity = c;
	nodes = new AlignNode*[c];

	int i = 0;
	for(i = 0; i<c; i++){
		nodes[i] = NULL;
	}

	idListSize = 0;
	idListCap = 20;
	idList = new int[idListCap];

}


TermList::~TermList()
{
	////Do not delete the env.  it deletes you.
	////GATenv * e

	int i = 0;
	for(i = 0; i<size; i++){
		delete(nodes[i]);
	}
	delete[](nodes);
	delete[](idList);
}

void TermList::printIdList()
{
	int i = 0;
	for(i = 0; i<idListSize; i++){
		printf("      %i\n", idList[i]);
	}
	
}

void TermList::addNode(AlignNode * n)
{
	int i = 0;
	int j = 0;
	AlignNode * tempNode;

	int * tempList = NULL;
	if(idListSize < idListCap){
		idList[idListSize] = n->id;
		idListSize++;
	}
	else{
		tempList = new int[idListCap*2];
		for(i = 0; i<idListSize; i++){
			tempList[i] = idList[i];
		}
		tempList[idListSize] = n->id;
		delete[](idList);
		idList = tempList;
		idListSize++;
		idListCap = idListCap*2;
	}

	////////////////////////////////////////////////////////
	///check for duplicates, and eliminate://///////////////
	#ifdef ELIMINATE_DUPLICATES							////
														////
	for(i = 0; i<size; i++){							////
		if( checkDuplicate(n, nodes[i]) ){				////
			///leak the node for now					////
			delete(n);									////
			n = NULL;									////
			return;										////
		}												////
	}													////
														////
	#endif												////
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////
	#ifdef MINIMIZE_BEST_RMSD							////
														////
//	n->getBestFitRMSD();								////
														////
	#endif												////
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////

	bool inserted = false;
	/////If there's space, keep it
	if(size == 0){
		nodes[0] = n;
		size++;
	}
	else if(size < capacity){
		for(i = 0; i<size+1; i++){
			if( nodes[i] == NULL ){
				nodes[i] = n;
				break;
			}
			if( (nodes[i]->size < n->size) || ( (nodes[i]->size == n->size) && (nodes[i]->RMSD > n->RMSD)) ){
				for(j = size; j>i; j--){
					nodes[j] = nodes[j-1];
				}
				nodes[i] = n;
				break;			
			}
		}
		size++;	
	}
	else{///If not, toss the worst one  (size == capacity)
		for(i = 0; i<capacity; i++){
			tempNode = nodes[i];
			if( (n->size > tempNode->size) || (n->RMSD < tempNode->RMSD && n->size == tempNode->size) ){
				delete(nodes[size-1]);				////bug is the zero node
				
				////move everything else down
				for(j = size-1; j>i; j--){
					nodes[j] = nodes[j-1];
				}
				
				////Insert the new one
				nodes[i] = n;
				inserted = true;
				break;
			}
		}
		///if it wasnt inserted then it is worse than all the others.
		if(!inserted){
			delete(n);
		}
	}
}

bool TermList::checkDuplicate(AlignNode * node1, AlignNode * node2)
{
	if(node1->size != node2->size){
		return false;
	}

	int size = node1->size;
	int i = 0;
	int j = 0;
	int tempSource;
	int tempTarget;
	bool foundSame;

	for(i = 0; i<size; i++){
		foundSame = false;
		tempSource = node1->sourceIds[i];
		tempTarget = node1->targetIds[i];
		for(j = 0; j<size; j++){
			if( (node2->sourceIds[j] == tempSource) && (node2->targetIds[j] == tempTarget) ){
				foundSame = true;
				break;
			}
		}
		if(foundSame){
			continue;
		}
		return false;
	}

	return true;

}

void TermList::sortEntries(void)
{
	int i = 0;

	for(i = 0; i<size; i++){
		nodes[i]->sortMatchIds();
	}
}








































