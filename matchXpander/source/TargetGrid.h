// TargetGrid.h: interface for the TargetGrid class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_GEOHASH_TARGETGRID_)
#define _GEOHASH_TARGETGRID_

class TargetGrid  
{
public:
	TargetGrid( AtomList * a, int * ids );		///takes an atomList and an array of ints of the residue INDEX (not number) associated for each atom.
	virtual ~TargetGrid();

	void setupTable(int * ids);					///sets up the table with the atomList atoms, and this set of integers which specifies
												///the set of residue INDICES each atom comes from

	set_t querySphere(double x, double y, double z, set_t set, double radius, set_t output);
	int * queryShell(double x, double y, double z, set_t set, double radius);
																			////Tests similarity within a shell of thickness EXTMATCHTHRESHOLD
																			////and radius radius.
	void currentBoxPlus();
	void currentBoxMinus();
//	void moveTestProbe(double x, double y, double z);
	void renderMe();


//////////////////////////////////////DEPRECATED////////////////////////////////////////
//////////////////////////////////////DEPRECATED////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
//	int query(double x, double y, double z, char * assocAA);				////Tests pointwize similarity of two aligned pts, quickly
//	int query(double x, double y, double z, char * assocAA, char atom);		////Tests pointwize similarity of two aligned ATOMS
//	IntStack * query(double x, double y, double z, char * assocAA, char atom, double thresh);	////Tests pointwize similarity of two aligned ATOMS
//	IntStack * testQuery(double ** coords, int size, char * assocAA, char atom, double thresh);	////Tests pointwize similarity of two aligned ATOMS
////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////DEPRECATED////////////////////////////////////////
//////////////////////////////////////DEPRECATED////////////////////////////////////////


//	GeoHashMotif * motif;
	AtomList * list;
	int * *grid;

	int xneg;						///boundary values
	int xpos;
	int yneg;
	int ypos;
	int zneg;
	int zpos;	

	int xdim;						///size values
	int ydim;
	int zdim;

	int num;						///number of points in Constellation

	int currentBox;					///DEBUGGING index
	double testProbeCoordx;			///test probe moves around the space, and activates anyting it gets near (nearest).
	double testProbeCoordy;			///
	double testProbeCoordz;			///
	char * testProbeAssocAA;		///
	int currentProbeFind;			///The current AA found by the probe.

	double xcentroid;				///for rendering (debug)
	double ycentroid;
	double zcentroid;

};

#endif


