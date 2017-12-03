// AtomList.h: interface for the AtomList class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_GEOHASH_ATOMLIST_)
#define _GEOHASH_ATOMLIST_

class AtomList  
{
public:
	AtomList(int s);
	AtomList(int s, PdbAtom * *atoms);
	virtual ~AtomList();
	PdbAtom * atomNumber(int i);

	AtomList * copy(void);

	void toString(void);							////output coords
	void highlightAA(int index);
	void highlightUnderThresh(int thresh);
	void highlightAllAAs(void);
	void highlightNoAAs(void);
	void highlightAssociatedAAs(bool info);
	void highlightAssociatedAAs(int * info);
	int * getHighlightedAtoms();
	int * getHighlightedAtomNumbers();

	int * getMXmapping(bool AAorAtomMode);		///returns an array with a mapping, first char is the #.
	int setMXthresh(bool AAorAtomMode);			///finds and highlights everything at appropriate MX level, by dehighlighting inappropriate

	char * name;
	int size;										////num of PdbAtoms stored in AtomList
	PdbAtom * *atomList;							////array of PdbAtoms
	int maxSeq;										////
	bool * highlights;								////

	///Selection Box Controls
	void moveSelectionBox(float dx, float dy, float dz);
	void resizeSelectBox(float dr);
	void highlightInsideBox(bool yesNo);			///either selects or deselects inside box

	///Selection Box 
	float selectx;
	float selecty;
	float selectz;
	float radius;
	bool renderSelectBox;

	////////////////////////////////////////////////////////
	///////////////Rendering Methods////////////////////////
	////////////////////////////////////////////////////////
	#ifdef INTERACTIVEMODE

	///If atomsOrAAs is true, we pick atoms. if false, we pick AAs.
	void processHits(GLint hits, GLuint buffer[], bool atomsOrAAs, float * rot, float * pos);

	////OpenGl calls to actually draw the molecule
	void renderMolecule(int etDisplay, int pickMode);

	#endif
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////

};

#endif //!defined _GEOHASH_ATOMLIST_



