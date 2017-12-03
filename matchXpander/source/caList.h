// caList.h: interface for the caList class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_CALIST_H_)
#define _CALIST_H_

class caList  
{
public:
	caList();
	caList(AtomList * aList);
	virtual ~caList();
	void addAcid(aminoAcid * a);
	void setName(char * n);
	void calcMaxSeq();
	void renderMe(bool background, bool colored, bool atoms, bool colormode, bool searchArea, int et, bool matches, bool nohydro);
	void setSourceHighlightBin(int num, int * vals);
	void setTargetHighlightBin(int num, int * vals);
	void clearHighlightBin();
	void setHighlightsForTarget();
	void toString();
	caList * subSet(int size, int * vals);

	aminoAcid * *acids;
	atomBag * bag;
	int numAcids;

//	bool atomsOrAcids;
//	bool backgroundOn;
//	bool solidBackground;
//	int etLev;
//	bool searchArea;

	char * name;
	int maxSeq;
	bool * highlights;	///Instantiate the Highlight Bin on CalcMaxSeq - after all is done

};



#endif



