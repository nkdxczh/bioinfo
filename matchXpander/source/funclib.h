///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///Processing Functions for Geometric Hashing Code, by Brian Chen//////
///////////////////////////////////////////////////////////////////////
////////////////////Header File: funclib.h/////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
/*********************************************************************/

#if !defined(_FUNCLIB_H_)
#define _FUNCLIB_H_

int translateAtomCode(char * code);
char * encodeAtomCode(int input);

char translateAAcode(char * code);

void smallSort(int * objs, int * refs, int size);
void smallSort(double * objs, int * refs, int size);
void smallSortDescending(int * objs, int * refs, int size);

void mergeSort(double * objs, int * refs, int size);

bool compareLabels(char * label1, char * label2);

//char * stringCopy(char* str, int start, int end);

//int * getMatchingPoints(GeoHashMotif * motif, int * indices, int * targetPts, );

bool subChar(char * set, char val);


void insertionSortAscending(int * list, int size);


#endif



