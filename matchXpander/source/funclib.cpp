///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///Processing Functions for Geometric Hashing Code, by Brian Chen//////
///////////////////////////////////////////////////////////////////////
//////////////////////Implementation: funclib.cpp//////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
/*********************************************************************/

#include "StdAfx.h"


/*
/////////////////////////////////////////////////////////////////////////////////////////////////
/////This function is for getting values out of the file when the position is fixed//////////////
/////Unlike strcpy which copies everything; I want to just copy specifics////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
char * stringCopy(char* str, int start, int end)
{
	int size = end - start + 1;
	char * result = new char[size + 1];

	for(int i = 0; i<size; i++){
		result[i] = str[start+i];	
	}
	result[size] = NULL;

	return result;

}
*/

int translateAtomCode(char * code)
{
	int result;
	if      (strcmp(code, "C")   == 0) result = 1; 	
	else if (strcmp(code, "CA")  == 0) result = 2;	
	else if (strcmp(code, "CB")  == 0) result = 3;	
	else if (strcmp(code, "CD")  == 0) result = 4;	
	else if (strcmp(code, "CD1") == 0) result = 5;	
	else if (strcmp(code, "CD2") == 0) result = 6;	
	else if (strcmp(code, "CE")  == 0) result = 7;	
	else if (strcmp(code, "CE1") == 0) result = 8;	
	else if (strcmp(code, "CE2") == 0) result = 9;	
	else if (strcmp(code, "CE3") == 0) result = 10;	
	else if (strcmp(code, "CG")  == 0) result = 11;	
	else if (strcmp(code, "CG1") == 0) result = 12;	
	else if (strcmp(code, "CG2") == 0) result = 13;	
	else if (strcmp(code, "CH2") == 0) result = 14;	
	else if (strcmp(code, "CZ")  == 0) result = 15;	
	else if (strcmp(code, "CZ2") == 0) result = 16;	
	else if (strcmp(code, "CZ3") == 0) result = 17;	
	else if (strcmp(code, "N")   == 0) result = 18;	
	else if (strcmp(code, "ND1") == 0) result = 19;	
	else if (strcmp(code, "ND2") == 0) result = 20;	
	else if (strcmp(code, "NE")  == 0) result = 21;	
	else if (strcmp(code, "NE1") == 0) result = 22;	
	else if (strcmp(code, "NE2") == 0) result = 23;	
	else if (strcmp(code, "NH1") == 0) result = 24;	
	else if (strcmp(code, "NH2") == 0) result = 25;	
	else if (strcmp(code, "NZ")  == 0) result = 26;	
	else if (strcmp(code, "O")   == 0) result = 27;	
	else if (strcmp(code, "OD1") == 0) result = 28;	
	else if (strcmp(code, "OD2") == 0) result = 29;	
	else if (strcmp(code, "OE1") == 0) result = 30;	
	else if (strcmp(code, "OE2") == 0) result = 31;	
	else if (strcmp(code, "OG")  == 0) result = 32;	
	else if (strcmp(code, "OH")  == 0) result = 33;	
	else if (strcmp(code, "SD")  == 0) result = 34;	
	else if (strcmp(code, "OG1")  == 0) result = 35;	
	else if (strcmp(code, "H")  == 0) result = 36;	
	else if (strcmp(code, "HG")  == 0) result = 37;	
	else if (strcmp(code, "HD2")  == 0) result = 38;	
	else if (strcmp(code, "HE")  == 0) result = 39;	
	else if (strcmp(code, "HH1")  == 0) result = 40;	
	else if (strcmp(code, "HH2")  == 0) result = 41;	
	else if (strcmp(code, "HZ")  == 0) result = 42;	
	else if (strcmp(code, "HE2")  == 0) result = 43;	
	else if (strcmp(code, "HA")  == 0) result = 44;	
	else if (strcmp(code, "HB")  == 0) result = 45;	
	else if (strcmp(code, "HD")  == 0) result = 46;	
	else if (strcmp(code, "HD1")  == 0) result = 47;	
	else if (strcmp(code, "HE1")  == 0) result = 48;	
	else if (strcmp(code, "HH")  == 0) result = 49;	
	else if (strcmp(code, "OXT")  == 0) result = 50;	
	else if (strcmp(code, "HG1")  == 0) result = 51;	
	else if (strcmp(code, "HG2")  == 0) result = 52;	
	else if (strcmp(code, "SG")  == 0) result = 53;	
	else if (strcmp(code, "HE3")  == 0) result = 54;	
	else if (strcmp(code, "HZ2")  == 0) result = 55;	
	else if (strcmp(code, "HZ3")  == 0) result = 56;	
	else {
		printf("ERROR: Unidentified Atom Tag Specified: %s\n", code);
		exit(1);
	}
	return result;
}

char * encodeAtomCode(int input)
{
	char * result = new char[4];
	if      (input == 1) 	{strcpy(result, "C")  ;}
	else if (input == 2)	{strcpy(result, "CA") ;}
	else if (input == 3)	{strcpy(result, "CB") ;}
	else if (input == 4)	{strcpy(result, "CD") ;}
	else if (input == 5)	{strcpy(result, "CD1");}
	else if (input == 6)	{strcpy(result, "CD2");}
	else if (input == 7)	{strcpy(result, "CE") ;}
	else if (input == 8)	{strcpy(result, "CE1");}
	else if (input == 9)	{strcpy(result, "CE2");}
	else if (input == 10)	{strcpy(result, "CE3");}
	else if (input == 11)	{strcpy(result, "CG") ;}
	else if (input == 12)	{strcpy(result, "CG1");}
	else if (input == 13)	{strcpy(result, "CG2");}
	else if (input == 14)	{strcpy(result, "CH2");}
	else if (input == 15)	{strcpy(result, "CZ") ;}
	else if (input == 16)	{strcpy(result, "CZ2");}
	else if (input == 17)	{strcpy(result, "CZ3");}
	else if (input == 18)	{strcpy(result, "N")  ;}
	else if (input == 19)	{strcpy(result, "ND1");}
	else if (input == 20)	{strcpy(result, "ND2");}
	else if (input == 21)	{strcpy(result, "NE") ;}
	else if (input == 22)	{strcpy(result, "NE1");}
	else if (input == 23)	{strcpy(result, "NE2");}
	else if (input == 24)	{strcpy(result, "NH1");}
	else if (input == 25)	{strcpy(result, "NH2");}
	else if (input == 26)	{strcpy(result, "NZ") ;}
	else if (input == 27)	{strcpy(result, "O")  ;}
	else if (input == 28)	{strcpy(result, "OD1");}
	else if (input == 29)	{strcpy(result, "OD2");}
	else if (input == 30)	{strcpy(result, "OE1");}
	else if (input == 31)	{strcpy(result, "OE2");}
	else if (input == 32)	{strcpy(result, "OG") ;}
	else if (input == 33)	{strcpy(result, "OH") ;}
	else if (input == 34)	{strcpy(result, "SD") ;}
	else if (input == 35)   {strcpy(result, "OG1") ;}  
	else if (input == 36)	{strcpy(result, "H") ;}  
	else if (input == 37)	{strcpy(result, "HG") ;}  
	else if (input == 38)	{strcpy(result, "HD2") ;}  
	else if (input == 39) 	{strcpy(result, "HE") ;}
	else if (input == 40) 	{strcpy(result, "HH1") ;}
	else if (input == 41) 	{strcpy(result, "HH2") ;}
	else if (input == 42) 	{strcpy(result, "HZ") ;}
	else if (input == 43) 	{strcpy(result, "HE2") ;}
	else if (input == 44) 	{strcpy(result, "HA") ;}
	else if (input == 45) 	{strcpy(result, "HB") ;}
	else if (input == 46) 	{strcpy(result, "HD") ;}
	else if (input == 47) 	{strcpy(result, "HD1") ;}
	else if (input == 48) 	{strcpy(result, "HE1") ;}
	else if (input == 49) 	{strcpy(result, "HH") ;}
	else if (input == 50) 	{strcpy(result, "OXT") ;}
	else if (input == 51) 	{strcpy(result, "HG1") ;}
	else if (input == 52) 	{strcpy(result, "HG2") ;}
	else if (input == 53) 	{strcpy(result, "SG") ;}
	else if (input == 54) 	{strcpy(result, "HE3") ;}
	else if (input == 55) 	{strcpy(result, "HZ2") ;}
	else if (input == 56) 	{strcpy(result, "HZ3") ;}
	else {
		printf("ERROR: Unidentified Atom number Specified: %i\n", input);
		exit(1);
	}
	return result;
}



char translateAAcode(char * code)
{
	char result;
	
	if		(strcmp(code, "GLY") == 0) result = 'G';		//	G - Glycine (Gly) 
	else if (strcmp(code, "PRO") == 0) result = 'P'; 	//	P - Proline (Pro) 
	else if (strcmp(code, "ALA") == 0) result = 'A';		//	A - Alanine (Ala) 
	else if (strcmp(code, "VAL") == 0) result = 'V';		//	V - Valine (Val) 
	else if (strcmp(code, "LEU") == 0) result = 'L';		//	L - Leucine (Leu) 
	else if (strcmp(code, "ILE") == 0) result = 'I';		//	I - Isoleucine (Ile) 
	else if (strcmp(code, "MET") == 0) result = 'M';		//	M - Methionine (Met) 
	else if (strcmp(code, "CYS") == 0) result = 'C';		//	C - Cysteine (Cys) 
	else if (strcmp(code, "PHE") == 0) result = 'F';		//	F - Phenylalanine (Phe) 
	else if (strcmp(code, "TYR") == 0) result = 'Y';		//	Y - Tyrosine (Tyr) 
	else if (strcmp(code, "TRP") == 0) result = 'W';		//	W - Tryptophan (Trp) 
	else if (strcmp(code, "HIS") == 0) result = 'H';		//	H - Histidine (His) 
	else if (strcmp(code, "LYS") == 0) result = 'K';		//	K - Lysine (Lys) 
	else if (strcmp(code, "ARG") == 0) result = 'R';		//	R - Arginine (Arg) 
	else if (strcmp(code, "GLN") == 0) result = 'Q';		//	Q - Glutamine (Gln) 
	else if (strcmp(code, "ASN") == 0) result = 'N';		//	N - Asparagine (Asn) 
	else if (strcmp(code, "GLU") == 0) result = 'E';		//	E - Glutamic Acid (Glu) 
	else if (strcmp(code, "ASP") == 0) result = 'D';		//	D - Aspartic Acid (Asp) 
	else if (strcmp(code, "SER") == 0) result = 'S';		//	S - Serine (Ser) 
	else if (strcmp(code, "THR") == 0) result = 'T';		//	T - Threonine (Thr) 
	else {
		printf("ERROR: Unidentified Amino Acid Tag Specified\n");
		result = 'X';
	}

	return result;

}


/////////////////////////////
///We sort objs, and re-arrange refs in the same way we re-arranged objs.
void smallSort(int * objs, int * refs, int size){
	int min = INT_MAX;
	int tempRef = 0;
	int tempVal;
	int i = 0;
	int j = 0;

	for(i = 0; i<size; i++){
		//Reset the min val, to keep finding the min
		min = INT_MAX;
		//find the smallest value between the val after i to the end
		for(j = i; j<size; j++){
			if(objs[j] <= min){ 
				tempRef = j; 
				min = objs[j];
			}
		}
		//reArrange the objs array;
		tempVal = objs[i];
		objs[i] = objs[tempRef];
		objs[tempRef] = tempVal;
		//reArrange the refs array;
		tempVal = refs[i];
		refs[i] = refs[tempRef];
		refs[tempRef] = tempVal;
	}
}


/////////////////////////////
///We sort objs, and re-arrange refs in the same way we re-arranged objs.
void smallSort(double * objs, int * refs, int size){
	double min = HUGE_VAL;
	int tempRef = 0;
	double tempVal1;
	int tempVal2;
	int i = 0;
	int j = 0;

	for(i = 0; i<size; i++){
		//Reset the min val, to keep finding the min
		min = HUGE_VAL;
		//find the smallest value between the val after i to the end
		for(j = i; j<size; j++){
			if(objs[j] <= min){ 
				tempRef = j; 
				min = objs[j];
			}
		}
		//reArrange the objs array;
		tempVal1 = objs[i];
		objs[i] = objs[tempRef];
		objs[tempRef] = tempVal1;
		//reArrange the refs array;
		tempVal2 = refs[i];
		refs[i] = refs[tempRef];
		refs[tempRef] = tempVal2;
	}
}


/////////////////////////////
///We sort objs, and re-arrange refs in the same way we re-arranged objs.
void smallSortDescending(int * objs, int * refs, int size){
	int max = 0;
	int tempRef = 0;
	int tempVal;
	int i = 0;
	int j = 0;

	for(i = 0; i<size; i++){
		//Reset the min val, to keep finding the min
		max = 0;
		//find the smallest value between the val after i to the end
		for(j = i; j<size; j++){
			if(objs[j] >= max){ 
				tempRef = j; 
				max = objs[j];
			}
		}
		//reArrange the objs array;
		tempVal = objs[i];
		objs[i] = objs[tempRef];
		objs[tempRef] = tempVal;
		//reArrange the refs array;
		tempVal = refs[i];
		refs[i] = refs[tempRef];
		refs[tempRef] = tempVal;
	}
}


bool compareLabels(char * label1, char * label2)
{
	bool result = false;

	int * lab1 = new int[4];
	int * lab2 = new int[4];
	int * order1 = new int[4];
	int * order2 = new int[4];

	order1[0] = 0;		///dont care about the order, as long as it gets sotred
	order1[1] = 1;
	order1[2] = 2;
	order1[3] = 3;

	order2[0] = 0;
	order2[1] = 1;
	order2[2] = 2;
	order2[3] = 3;

	lab1[0] = (int) label1[0];
	lab1[1] = (int) label1[1];
	lab1[2] = (int) label1[2];
	lab1[3] = (int) label1[3];

	lab2[0] = (int) label2[0];
	lab2[1] = (int) label2[1];
	lab2[2] = (int) label2[2];
	lab2[3] = (int) label2[3];

	smallSort(lab1, order1, 4);		///dont care about the order, as long as it gets sotred
	smallSort(lab2, order2, 4);

	if ( lab1[0] == lab2[0] &&
		lab1[1] == lab2[1] && 
		lab1[2] == lab2[2] && 
		lab1[3] == lab2[3] ){
		result = true;
	}

	delete(lab1);
	delete(lab2);
	delete(order1);
	delete(order2);

	return result;

}



////////////////////////////////////
//If this char is in the set, return true;
//else false;
bool subChar(char * set, char val)
{
	int len = strlen(set);
	int i = 0;
	for (i = 0; i<len; i++){
		if (val == set[i]){
			return true;
		}
	}
	return false;
}

///simple insertion sort.  n^2 runtime - no good for large application
void insertionSortAscending(int * list, int size)
{
	int i = 0;
	int j = 0;
	int tempVal = 0;
	int tempVal2 = 0;
	int pointer = 0;

	for (i = 0; i<size; i++){
		tempVal = list[i];
		tempVal2 = list[i];
		pointer = i;

		for(j = i; j<size; j++){
			if(list[j] < tempVal2){
				tempVal2 = list[j];
				pointer = j;
			}
		}

		list[i] = tempVal2;
		list[pointer] = tempVal;
	}
}


void mergeSort(double * objs, int * refs, int size)
{
	///if its big, split.
	if(size > 10){
		int size1 = ((int) ( ((double)size)/2.0) );
		int size2 = size - size1;

		double * objs1 = new double[size1];
		int * refs1 = new int[size1];
		double * objs2 = new double[size2];
		int * refs2 = new int[size2];

		int i = 0;

		////split into two problems.
		for(i = 0; i<size1; i++){
			objs1[i] = objs[i];
			refs1[i] = refs[i];
		}
		for(i = 0; i<size2; i++){
			objs2[i] = objs[size1+i];
			refs2[i] = refs[size1+i];
		}

		///sort them
		mergeSort(objs1, refs1, size1);
		mergeSort(objs2, refs2, size2);

		///merge the two solutions
		int counter1 = 0;
		int counter2 = 0;
		int caseFlag;
		for(i = 0; i<size; i++){
			caseFlag = 0;
			if(counter1 < size1 && counter2 < size2){
				caseFlag = 0;
			}else if(counter1 < size1 && counter2 >= size2){
				caseFlag = 1;
			}else if(counter1 >= size1 && counter2 < size2){
				caseFlag = 2;
			}else {
				caseFlag = 3;	///error
			}
			switch(caseFlag){
				case 0 : 
					if( (objs1[counter1] <= objs2[counter2]) ){
						objs[i] = objs1[counter1];
						refs[i] = refs1[counter1];
						counter1++;
					}else{
						objs[i] = objs2[counter2];
						refs[i] = refs2[counter2];
						counter2++;
					}break;
				case 1 :
					objs[i] = objs1[counter1];
					refs[i] = refs1[counter1];
					counter1++;
					break;
				case 2 :
					objs[i] = objs2[counter2];
					refs[i] = refs2[counter2];
					counter2++;
					break;
				case 3 :
					printf("ERROR: mergeSort is buggy\n");
					exit(0);
			}
		}

		delete[](objs1);
		delete[](refs1);
		delete[](objs2);
		delete[](refs2);

	}
	else{////otherwise insertion sort it.
		smallSort(objs, refs, size);		
	}




}








