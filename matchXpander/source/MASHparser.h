// MASHparser.h : 
//  Parses MASH files instead of FAT files
//

#ifndef MASH_PARSER_H
#define MASH_PARSER_H

#include "StdAfx.h"

///this determines if the file is in the MASH format by testing the first line.
bool testMASHformat(char * fileName);

///This parses the MASH file for the points, and returns an atomList
AtomList * parseMASHfilePoints(char * fileName);









#endif
