// HierarchialController.cpp: implementation of the HierarchialController class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

HierarchialController::HierarchialController(GATenv * e)
{
	env = e;
	int i = 0;

	if(!e->hasAssocArrays){

		printf("MatchXpander doesnt have the necessary information to do match augmentation\n");
		printf("MatchXpander will now quit.\n");
		exit(0);
	}

	stack = new PriorityStack(env);
	active = new ActiveList(env);
	newList = new ActiveList(env, active->count);			///this guarantees that all id's are unique
	tierNodes = new AlignNode*[stack->height];

	////Clear the TierNodes
	for(i = 0; i<stack->height; i++){
		tierNodes[i] = NULL;
	}

	term = new TermList(env->numArrays, env);
	startTime = 0;
	endTime = 0;

}

HierarchialController::~HierarchialController()
{
	int i = 0;
	for(i = 0; i<stack->height; i++){
		delete(tierNodes[i]);
	}
	delete[](tierNodes);
	delete(stack);
	delete(active);
	delete(newList);
	delete(term);
}


///tests a particular Index target, with a particular node
bool HierarchialController::testIndex(int ind, AlignNode * node)
{
	bool result = false;
	AlignNode * nodeCopy;
	AlignNode * newNode;
	set_t vals;
	int i = 0;
	double info;

	for(i = 0; i<node->size; i++){
		if(ind == node->sourceIds[i]){
			return true;		///nothing got newed up to here, so its ok.
		}
	}

	///Get similar points from the target
	vals = env->query(ind, node->accumulationMatrix);

	///Find out if it can be pegged If so, add it to the list
	if(size_set(vals) == 0){
		free_set(vals);
		return false;
	}

	////Add the first one onto yourself, and continue, and duplicate off the rest
	nodeCopy = node->copy(newList->getNewId());

	///test out the original align Node, and the result will be the result for that one.
	info = node->testPt(ind, vals[0]);
	if(info >= 0){
		node->addGeometricPt(ind, vals[0]);
		node->RMSD = info;
		result = true;
	}
	else{
		///this is not deleted here because it is handled in a higher level if it is false.
		//printf("negative\n");	///debug
	}

	///add other matches to the newList
	for(i = 1; i<size_set(vals); i++){

//		printf("size_set: %i\n", size_set(vals));

		newNode = nodeCopy->copy(newList->getNewId());	
		info = newNode->testPt(ind, vals[i]);
		if(info > 0){
			////mult in the new matrix
			newNode->addGeometricPt(ind, vals[i]);
			newNode->RMSD = info;

			if(newNode->size == env->sourceList->numAcids){			////new
				newNode->isdead = true;
				term->addNode(newNode); 
				continue;
			}
			newList->addNode(newNode);
			result = true;
		}
		else{
			delete(newNode);
//			printf("negative\n");	///debug
		}
	}

	delete(nodeCopy);
	free_set(vals);

	return result;
}


/////Test a while list of values with testIndex
bool HierarchialController::testArray(int * array, int length, AlignNode * tempNode)
{
	int i = 0;																						////
//	int j = 0;
	bool termStatus = false;  ///true for terminate															////
	int tempStartingIndex = tempNode->startingIndex;														////
	bool verdict;																							////
	int * tempHash = new int[length];
	for(i = 0; i<length; i++){
		tempHash[i] = abs(array[i]);
	}
																											////
	///Test each one of the points																			////
	//positives first																						////
	for(i = tempNode->startingIndex; i<length; i++){														////
		tempNode->startingIndex = i;																		////
		if(array[i] >= 0){			
			////we mark these source points so we dont check them below										////
			tempHash[i] = -tempHash[i];
			if ( !tempNode->element( array[i] ) ){															////
				verdict = testIndex(array[i], tempNode);													////
				if( verdict ){																				////
					termStatus = true;																		////
					array[i] = -array[i];																	////
				}																							////
			}																								////
			else{																							////
				termStatus = true;																			////
			}																								////
		}																									////
	}																										////
																											////
	tempNode->startingIndex = tempStartingIndex;															////
	///Then test negatives - the ones that have already been matched										////
	for(i = tempNode->startingIndex; i<length; i++){														////
		tempNode->startingIndex = tempStartingIndex;														////
		if(array[i] < 0 && tempHash[i] >= 0){																					////
			if( !tempNode->element( abs(array[i]) ) ){														////
				verdict = testIndex(-array[i], tempNode);													////
				if( verdict ){																				////
					termStatus = true;																		////
				}																							////
			}																								////
			else{																							////
				termStatus = true;																			////
			}																								////
		}																									////
	}																										////
																											////
	delete[](tempHash);																						////
	///If termStatus is still false, send it to the terminated list											////
	if(!termStatus){																						////
//		if( active->remove(tempNode->id) == NULL){															////
//			newList->remove(tempNode->id);																	////
//		}																									////
		tempNode->isdead = true;																			////
		term->addNode(tempNode);																			////
	}																										////
																											////
	return termStatus;																						////
}



void HierarchialController::beginTesting()
{
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
//	int m = 0;
	int tempVal;
	AlignNode * node = NULL;
	AlignNode * tempNode = NULL;
	bool verdict;

	////test each seed
	i = 0;

	while(active->size>0){
			////reset the flags, reset the TeirNodes;
		stack->resetStack();
		for(l = 0; l<stack->height; l++){
			if(tierNodes[l] != NULL){
				delete(tierNodes[l]);
				tierNodes[l] = NULL;
			}
		}

		////test each array level - this generates the tierNodes
		node = active->remove(active->list[0]->id);
		multiLevelTest(node, true);
		node = NULL;

		////Then test out any remaining on the newList
		while(newList->size > 0){
			node = newList->remove();
			multiLevelTest(node, false);	////multilevelTest deals with node, so I dont need to delete it
		}
		node = NULL;

									//////BLOCKED OUT FOR DEBUG ONLY								////
		////Then attempt mop-up alternatives with those stored on the tierNodes						////
		for(j = 0; j<stack->height; j++){															////
			for(k = 0; k<stack->sizes[j]; k++){														////
				if(tierNodes[j] == NULL){															////
					break;																			////
				}																					////
				if(stack->flags[j][k] < 0){															////
					continue;																		////
				}																					////
				tempNode = tierNodes[j]->copy(newList->getNewId());									/////
				tempVal = stack->flags[j][k];														////
																									////
				////we only want to mop alternatives that have not been tried						////
				if(tempVal >= 0){																	////
					verdict = testIndex(tempVal, tempNode);											////
					if(verdict){																	////
						tempNode->currentTier = j;													////
						tempNode->startingIndex = k;												////
						multiLevelTest(tempNode, false);		///multilevelTest Doesnt deal with this node properly, bc its in neither list
						while(newList->size > 0){
							node = newList->remove();
							multiLevelTest(node, false);
						}
						tempNode = NULL;
					}
					else{
						delete(tempNode);
						tempNode = NULL;
					}
				}
				else{
				//deletion - tried already
					delete(tempNode);
					tempNode = NULL;
				}
			}
		}
		i++;
									//////BLOCKED OUT FOR DEBUG ONLY
	}
//	///Set all RMSDs
//	for(i = 0; i<term->size; i++){
//		term->nodes[i]->rmsdCheck();
//	}

	////Debug
#ifndef NO_OUTPUT
	printResults();
#endif

}

///tests all levels for a single node.  Prefer unMatched target indices over matched ones.
void HierarchialController::multiLevelTest(AlignNode * node, bool initial)
{
	int i = 0;																				////
//	int j = 0;																				////
	int * tempArray;																		////
	int arraySize;																			////
	AlignNode * tierNode;																	////
	bool infoVal = false;																	////
																							////
	if(node->isdead){																		////
		printf("ERROR: dead AlignNode being tested\n");										////
	}																						////
																							////
	////test each array level																////
	for(i = node->currentTier; i<stack->height; i++){										////
		////if testArray returns false, then this seed is done - terminate					////
		if(!stack->changeLevel(i)){															////
			printf("ERROR: Attempted to change to non-existant Level!\n");					////
			exit(0);																		////
		}																					////
		node->startingIndex = 0;															////
		tempArray = stack->currentFlags();													////
		arraySize = stack->currentSize();													////
																							////
		///first copy the node into the tierNodes											////
		if(initial){																		////
			tierNode = node->copy(newList->getNewId());										////
			tierNodes[i] = tierNode;														////
		}																					////
																
		///Now test the array
		infoVal = testArray(tempArray, arraySize, node);

		///If no matching was successful, testArray will terminate it
		if(!infoVal){
			node = NULL;
			break;
		}
	}


	///if the node still has not been eliminated by the last level, then terminiate it
	if(infoVal){
//		if( active->remove(node->id) == NULL){
//			if(	newList->remove(node->id) == NULL){
//			}
//		}
		node->isdead = true;
		term->addNode(node);
	}

}




void HierarchialController::printResults()
{
	int i = 0;
	for(i = 0; i<term->size; i++){
		term->nodes[i]->printResults();
	}

}













