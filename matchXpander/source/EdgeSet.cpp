// EdgeSet.cpp: implementation of the EdgeSet class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

EdgeSet::EdgeSet()
{
	size = 0;
	edgeRoot = alloc_set(SP_MAP);

}

EdgeSet::~EdgeSet()
{

	int i;

	int size = size_set(edgeRoot);
	set_t tempSet;

	for(i = 0; i<size; i++){
		tempSet = (set_t) ith_map_set(edgeRoot, i);

		free_set(tempSet);
	}
	
	free_set(edgeRoot);
}

int EdgeSet::addEdge(int a, int b)
{
	if (a == b) {
		return -1;
	}

	///insert both ways  first:
	set_t localSet = (set_t) mapsto_set(edgeRoot, a);
	if(localSet == NIL){
		localSet = alloc_set(0);
		edgeRoot = put_set(edgeRoot, a);
		localSet = put_set(localSet, b);
		edgeRoot = associate_set(edgeRoot, a, (ptr_t) localSet);
	}
	else{
		if( contains_set(localSet, b) ){
			return -1;
		}
		localSet = put_set(localSet, b);
		edgeRoot = associate_set(edgeRoot, a, (ptr_t) localSet);
	}

	///insert both ways  second:
	localSet = (set_t) mapsto_set(edgeRoot, b);
	if(localSet == NIL){
		localSet = alloc_set(0);
		edgeRoot = put_set(edgeRoot, b);
		localSet = put_set(localSet, a);
		edgeRoot = associate_set(edgeRoot, b, (ptr_t) localSet);
	}
	else{
		if( contains_set(localSet, a) ){
			return -1;
		}
		localSet = put_set(localSet, a);
		edgeRoot = associate_set(edgeRoot, b, (ptr_t) localSet);
	}
	return 1;
}


bool EdgeSet::contains(int a, int b)
{
	bool result = false;


	set_t localSet = (set_t) mapsto_set(edgeRoot, a);
	if(localSet != NIL){
		if( contains_set(localSet, b) ){
			result = true;
		}
	}
	localSet = (set_t) mapsto_set(edgeRoot, b);
	if(localSet != NIL){
		if( contains_set(localSet, a) ){
			result = true;
		}
	}

	return result;
}


set_t EdgeSet::getNeighbors(int a)
{
	int localSize = 0;
	int i = 0;

	set_t result = alloc_set(0);
	set_t localSet = (set_t) mapsto_set(edgeRoot, a);

	if(localSet != NIL){
		localSize = size_set(localSet);
		for(i = 0; i<localSize; i++){
			result = put_set(result, localSet[i]);
		}
	}
	else{
		free_set(result);
		result = NULL;
	}

	return result;
}



