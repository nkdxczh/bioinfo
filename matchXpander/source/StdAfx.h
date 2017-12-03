// stdafx.h : include file for standard system include files,
//  or project specific include files that are used frequently, but
//      are changed infrequently
//
#if !defined(_STDAFX_H_)
#define _STDAFX_H_

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////        INCLUDES        //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

/*--------------------- Libraries ---------------------*/
#include <stdio.h>
#include <stdlib.h>
//#include <GL/glut.h>
//#include <GL/gl.h>
//#include <GL/glui.h>			////Rendering is disabled as a part of this version.  Do not compile Renderer.cpp (old version)
#include <string.h>
#include <assert.h>
#include <math.h>
//#include <glui.h>				////Rendering is disabled as a part of this version.  Do not compile Renderer.cpp (old version)
#include <limits.h>
//#include <fstream.h>
#include <time.h>

/*----------------- Math & Helpers --------------------*/
#include "_lalgebra.h"
#include "mathlib.h"
#include "funclib.h"

/*------------------- Set Library ---------------------*/
#include "set.h"
#include "prime.h"
#include "conf.h"
#include "defs.h"

/*-------------- Data Containing Classes --------------*/
#include "PdbAtom.h"
#include "atomBag.h"
#include "aminoAcid.h"
#include "AtomList.h"
#include "caList.h"
#include "Motif.h"
#include "TargetGrid.h"
#include "LatticeHash.h"
#include "GATenv.h"
#include "AlignNode.h"
#include "EdgeSet.h"
#include "TriangleList.h"

/*----------------- Parsing Headers -------------------*/
///#include "GATparser.h"			///DEPRECATED
#include "FATparser.h"
#include "pdbParser.h"
#include "MASHparser.h"

/*--------------- BitHandling Headers -----------------*/
#include "BitFileProcessing.h"

/*----------------- Functional Classes ----------------*/
#include "ActiveList.h"
#include "TermList.h"
#include "PriorityStack.h"

/*------------------ Control Classes ------------------*/
#include "HierarchialController.h"
#include "SeedFinder.h"
//#include "pegRobot.h"			////this is no longer used

/*----------------- Rendering Classes -----------------*/
//#include "Renderer.h"			////Rendering is disabled as a part of this version.  Do not compile Renderer.cpp (old version)

/*------------------ Output Classes -------------------*/
#include "GatOutput.h"
#include "StatOutput.h"
#include "BatchList.h"
//#include "Renderer.h"

/*------------------ Execution Funcs-------------------*/
#include "Execution.h"

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////        DEFINITIONS        /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

/*-------------Algorithmic #defs -----------*/

#define EXTMATCHTHRESHOLD 7				////error permissible for C-Alpha differences in non-subpattern match
#define NUMRIGIDTRANSFORMS 21			///The number of trans/rots to expect from the pegRobot
#define RANGESEARCHSIZE 10.0				///tbe size of the accessibility sphere
#define NUMBESTTRIANGLES 30				///the number of triangles to be considered in 
#define TARGETGRIDSIZE 5				////size of the grid on which the target is gridded

#define MINIMIZE_BEST_RMSD				////sorts by least fit RMSD instead of alignment RMSD
#define ALIGN_BEST_RMSD					////Minimize alignments with least RMSD after new alignments have been expanded.  
										///////NOTE: Code not using this is still functional, but no longer supported.
#define ELIMINATE_DUPLICATES			////self explanatory

//#define ASSOCIATE_SIDECHAIN_ATOMS		////If this is defined, automatically match sidechain atoms if other members
										////of the same amino acid are matched.
				
#define RMSD_CLAMP	.0000000001			////Identical structures frustrate a numerical issue in _lalgebra.cpp which
										////causes zero-RMSD situations to turn into epsilon-RMSD situations.  The result
										////is transformation matrices which are inaccurate, propagating error through MX
										////This #define finds cases with RMSD less than .0000000001, and sets them to zero
										////And remaps the transformation matrix to the identity.
										/////////The caveat is that if in fact these are structures with actual RMSD this
										/////////low, and yet are not identical, then the identity transformation will 
										/////////be incorrect. However, we do not anticipate for this to occur in real life.
/*------------Text Parsing #defs ----------*/

#define NUMCOLUMNS 300					////parsing param for max # of columns

/*------------- Rendering #defs -----------*/
#define ATOMRADIUS 1.7					////atom radius in angstroms
#define ATOMVERTSUBS 6					////# of vertical subdivisions
#define ATOMHORISUBS 5					////# of horizontal subdivisions

#define HIGHLIGHTEDATOMRADIUS 2.0		////atom radius in angstroms of highlighted Atom
#define HIGHLIGHTEDATOMVERTSUBS 5		////# of vertical subdivisions of highlighted Atom
#define HIGHLIGHTEDATOMHORISUBS 3		////# of horizontal subdivisions of highlighted Atom

#define CLIPPINGVOL 75000.0f			///the depth of the back clipping plane of the view frustrum
#define POINTSIZE 1.0					///size to render points
#define LINEWIDTH 2.0					///size to render lines
#define ZDISPLACEMENT -200.0			///standard Z displacement for easy rendering

/*------------- Picking #defs ------------*/

#define BUFSIZE 512						///buffer size for picking

/*------------- Output #defs -------------*/

#define NO_OUTPUT						///turn this on to turn off verbose output while operating.  For debug purposes.

/*----------------------------------------*/






#endif
