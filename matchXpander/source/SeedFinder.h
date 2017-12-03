//
//Declarations of SeedFinder.c

#if !defined(_SEEDFINDER_H_)
#define _SEEDFINDER_H_

void findSeedMatches(PriorityStack * stack, caList * target, GATenv * env);
double * generateSideChainRMSD(int * sourceArray, int * targArray, int numRes, Motif * motif, caList * target, double * RMSD, double * transform);
void generateSideChainRMSD(int * matchArray, int numRes, Motif * motif, caList * target, double * RMSD, double * transform);
double * alignThreePoints(double ** sourceCoords, double ** targCoords, int step);
double alignThreePoints(double ** sourceCoords, double ** targCoords, double * transform);
int * generateAtomMatchListings( int numSourcePts, PdbAtom * *source, aminoAcid * target);

#endif



