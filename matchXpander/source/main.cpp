////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
//Implementation for the Main methods for the Geometric Hashing Database search
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

#include "StdAfx.h"

/// SETTINGS 
//
//-batchList -outputStats 16pk-1php-input.FAT blist
//-batchList -outputStats 1a3k-1qmj-input.FAT blist
//-FAT 16pk-1php-input.FAT pdb1a81L.pdb 16pk-pdb1a81L
//-FAT 1a3k-1qmj-input.FAT pdb1a81L.pdb 1a3k-pdb1a81L 
//-batchList -outputStats 1aky-1ak2-input.FAT blist
//


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
///Prints out current thresholds
void threshOutput(void)
{
	double temp;

	printf("\n");
	printf("\n");
	printf("Match Augmentation for Protein Structure Pattern Matching Software\n");
	printf("by Brian Chen, Rice University Computer Science, 2003\n");
	printf("Version: 7\n");
	printf("\n");
	printf("Please report any bugs to brianyc@cs.rice.edu\n");
	printf("\n");
	printf("\n");
	printf("\n");
	printf("Current Threshold Settings:\n");
	printf("==========================================================\n");
	temp = EXTMATCHTHRESHOLD;
	printf("ExtMatchThreshold:                               %lf\n", temp);
	temp = TARGETGRIDSIZE;
	printf("TargetGridSize:                                  %lf\n", temp);
	temp = NUMBESTTRIANGLES;
	printf("NumBestTriangles:                                %i\n", (int) temp);
	temp = NUMRIGIDTRANSFORMS;
	printf("NumRigidTransforms:                              %i\n", (int) temp);
	temp = RANGESEARCHSIZE;
	printf("RangeSearchSize:                                 %lf\n", temp);
#ifdef MINIMIZE_BEST_RMSD
	printf("MINIMIZE_BEST_RMSD:                              true\n");
#else 
	printf("MINIMIZE_BEST_RMSD:                              false\n");
#endif/////////////////////////////////////////////////////////
#ifdef MATCH_BY_LEAST_RMSD
	printf("Match_By_Least_RMSD:                             true\n");
#else 
	printf("Match_By_Least_RMSD:                             false\n");
#endif/////////////////////////////////////////////////////////
#ifdef ALIGN_BEST_RMSD
	printf("ALIGN_BEST_RMSD:                                 true\n");
#else 
	printf("ALIGN_BEST_RMSD:                                 false\n");
#endif/////////////////////////////////////////////////////////
#ifdef ELIMINATE_DUPLICATES
	printf("ELIMINATE_DUPLICATES:                            true\n");
#else 
	printf("ELIMINATE_DUPLICATES:                            false\n");
#endif/////////////////////////////////////////////////////////
#ifdef NO_OUTPUT
	printf("No_Output:                                       true\n");
#else 
	printf("No_Output:                                       false\n");
#endif/////////////////////////////////////////////////////////
	printf("==========================================================\n");
	printf("\n");
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

void defaultOutput()
{
	printf("\n");
	printf("\n");
	printf("Match Augmentation for Protein Structure Pattern Matching Software\n");
	printf("by Brian Chen, Rice University Computer Science, 2007\n");
	printf("Version: 3.0\n");
	printf("\n");
	printf("Usage:\n");
	printf("=====================================================================\n");
	printf("Threshold Information\n");
	printf("PROMPT> matchXpander -?\n");
	printf("\n");
	printf("NOTE: IN ALL CASES BELOW, FAT files and MASH files can be used interchangably.\n");
	printf("\n");
	printf("Batch Mode with Batch List: (motif from a .FAT with a list of pdb targets)\n");
	printf("-individuals for XGAT/XSTAT files, -outputStats for .META files, -checkPoint to quit if META exists.\n");
	printf("PROMPT> matchXpander -batchList <InputFile.FAT> <BatchList> (-individuals) (-outputStats) (-checkPoint)\n");
	printf("\n");
	printf("Batch Mode with Binary Batch\n");
	printf("PROMPT> matchXpander -binBatch <InputFile.FAT> <BINfile> <TOCfile> <outputFileName>\n");
	printf("\n");
	printf("Batch Mode with Fat and Pdb Input\n");
	printf("PROMPT> matchXpander -FAT <InputFile.FAT> <TargetFile.pdb> <OutputFileHeader>\n");
	printf("=====================================================================\n");
	printf("\n");
	printf("Please report any bugs to brianyc@cs.rice.edu\n");
	printf("\n");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

void debugStuff()
{


}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                  //                      
//  ###   ###    ###    ######  ##   ##      ###   ###  #######  ######  ##   ##   #####   #####    //                                                                                
//  #### ####   ## ##     ##    ###  ##      #### ####  ##         ##    ##   ##  ##   ##  ##  ##   //                                                                                                            
//  ## ### ##  ##   ##    ##    #### ##      ## ### ##  ######     ##    #######  ##   ##  ##   ##  //                                                                                                              
//  ##  #  ##  #######    ##    ## ####      ##  #  ##  ##         ##    ##   ##  ##   ##  ##   ##  //                                                                  
//  ##     ##  ##   ##    ##    ##  ###      ##     ##  ##         ##    ##   ##  ##   ##  ##  ##   //                                                            
//  ##     ##  ##   ##  ######  ##   ##      ##     ##  #######    ##    ##   ##   #####   #####    //                                                                         
//                                                                                                  //                      
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
	int i = 0;

//	debugStuff();

	///eliminate erroneous cases, process threshold requests
	if( (argc < 2) || (argc > 7) || (argc == 3) ){ defaultOutput(); exit(1); }
	if(strcmp(argv[1], "-?") == 0){ threshOutput(); exit(1); }

	////process switches
	bool individualsOn = false;	bool outputStatsOn = false;	bool checkPointOn = false;
	for(i = 0; i<argc; i++){
		if( strcmp("-individuals", argv[i])==0 ){ individualsOn = true; }
		if( strcmp("-outputStats", argv[i])==0 ){ outputStatsOn = true; }
		if( strcmp("-checkPoint",  argv[i])==0 ){ checkPointOn  = true; }
	}
	
	///check for no input on batchList requests:
	if( ( strcmp("-batchList", argv[1]) == 0 ) && (!individualsOn & !outputStatsOn) ){
		printf("Please select either -individuals or -outputStats for batch processing\n");
		exit(1);
	}
	
	///parse FAT or MASH file, and pass data in.
	AtomList * motifAlist;
	bool isMASHformat = testMASHformat(argv[2]);
	if(isMASHformat){ motifAlist = parseMASHfilePoints(argv[2]); }
	else{
		FATparser * fparse = new FATparser(argv[2]);
		fparse->readFile();
		motifAlist = fparse->sources[0];
	}
	
	///run the batch run
	if( strcmp("-batchList", argv[1]) == 0 ){
		initBatchMXexperiment(motifAlist, argv[3], individualsOn, outputStatsOn, checkPointOn);
	}

	///Run the binary batch run
	if (strcmp("-binBatch", argv[1]) == 0){
		initBinBatchExperiment(motifAlist, argv[3], argv[4], argv[5]);
	}

	////Run the MX experiment with the given target in the FAT file.  
	////(send NULL to the MXexperiment for the Target.)
	if (strcmp("-FAT", argv[1]) == 0){
		initMXExperiment(motifAlist, argv[4], argv[3]);
	}


	return 1;

}
