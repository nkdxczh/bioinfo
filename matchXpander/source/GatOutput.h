// GatOutput.h: interface for the GatOutput class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_GEOHASH_GATOUTPUT_)
#define _GEOHASH_GATOUTPUT_

class GatOutput  
{
public:
	GatOutput(GATenv * e, HierarchialController * c);
	virtual ~GatOutput();

	FILE * currentOutputFile;
	GATenv * env;
	HierarchialController * controller;

	//////MetaCall for GAT file output
	void outputGatFiles(char * idHeader);
	void outputPartialGatFiles(char * idHeader);
	void outputXGatFiles(char * idHeader);

	//////HeaderPrinter for the GAT format
	void printGlobals();
	//////sourceProtein printer (loop to get them all)
	void printSourceProtein();
	//////TargetMolecule printer
	void printTargetMol();
	//////Match Array Printer
	void printMatchArray();
	//////RMSD Array Printer
	void printINFOArray();

	//////XGAT Match Array Printer
	void printXGatMatchArray();
	//////XGAT RMSD Array Printer
	void printXGatINFOArray();

};

#endif



