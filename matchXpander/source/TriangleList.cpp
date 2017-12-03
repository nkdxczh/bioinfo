// TriangleList.cpp: implementation of the TriangleList class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

TriangleList::TriangleList()
{
	size = 0;
	capacity = 40;
	list = new int*[capacity];
}

TriangleList::~TriangleList()
{
	int i = 0; 
	for(i = 0; i<size; i++){
		delete[](list[i]);
	}
	delete[](list);

}

void TriangleList::expand()
{
	capacity = capacity * 2;
	int * *tempVals = new int*[capacity];

	int i =0;
	for(i = 0; i<size; i++){
		tempVals[i] = list[i];
	}

	delete[](list);
	list = tempVals;

}

void TriangleList::addTriangle(int s1, int t1, int s2, int t2, int s3, int t3)
{
///////REDUNDANT edges are eliminated at the EdgeHash level
//	if( t1 == t2 || t1 == t3 || t2 == t3 ){
//		printf("OMFG\n");
//	}
///////REDUNDANT edges are eliminated at the EdgeHash level


	if( ((double)(size+1)) > (double) (.75 * (double) capacity) ){
		expand();
	}
	int * sources = new int[3];
	int * targets = new int[3];
	sources[0] = s1;
	sources[1] = s2;
	sources[2] = s3;
	targets[0] = t1;
	targets[1] = t2;
	targets[2] = t3;

	smallSort(sources, targets, 3);

	int * newArray = new int[6];
	newArray[0] = sources[0];
	newArray[1] = targets[0];
	newArray[2] = sources[1];
	newArray[3] = targets[1];
	newArray[4] = sources[2];
	newArray[5] = targets[2];

	list[size] = newArray;
	size++;

	delete[](sources);
	delete[](targets);

}

void TriangleList::print()
{
	int i = 0;
		printf("Size: %i\n", size);
	for(i = 0; i<size; i++){
		printf("%i:  (%i, %i), (%i, %i), (%i, %i)\n", i, list[i][0], list[i][1], list[i][2], list[i][3], list[i][4], list[i][5]);
	}



}



