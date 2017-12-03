
////////////////////////////////////////////////////////
//Execution_H.  This file defines functions which operate MatchXpander
//				Which are called from the main method.
////////////////////////////////////////////////////////
#if !defined(_EXECUTION_H_)
#define _EXECUTION_H_

void initBatchMXexperiment(AtomList * aMotif, char * batchListFile, bool ind, bool out, bool checkPoint);
//void initBatchMXexperiment(char * inputFile, char * batchListFile, bool ind, bool out, bool checkPoint);

void initMXExperiment(AtomList * aMotif, char * output, char * targ);
//void initMXExperiment(char * input, char * output, char * targ);

void initBinBatchExperiment(AtomList * aMotif, char * binFile, char * binToc, char * outputFileName);
//void initBinBatchExperiment(char * fatFile, char * binFile, char * binToc, char * outputFileName);

///deprecated
//double initMXExperimentMPI(char * target, AtomList * motif);

#endif




