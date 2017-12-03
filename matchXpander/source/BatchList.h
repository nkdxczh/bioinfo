// BatchList.h: interface for the BatchList class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BATCHLIST_H__8B1BE518_399A_486E_A7A1_4C6FE2E9DD59__INCLUDED_)
#define AFX_BATCHLIST_H__8B1BE518_399A_486E_A7A1_4C6FE2E9DD59__INCLUDED_

class BatchList  
{
public:
	BatchList(char * fname, Motif * m);
	BatchList(Motif * m, int thisSize);
	virtual ~BatchList();

	int size;
	int currentFile;
	char ** fileList;
	FILE * file;

	int years;
	int days;
	int hours;
	int minutes;
	int seconds;
	int ms;

	long startTime;
	long endTime;

	///Output info
	char ** expList;
	char ** matchSizes;
	char ** matchRMSDs;
	char ** matchListings;

	///Global Statistics
	int sourceSize;
	int numMatchedPerfect;
	double avgPointsMatched;
	double avgMatchRMSD;

	///input Data
	Motif * motif;

	///Methods
	char * getNextFile();
	bool finished();
	void parse();
	void printOutMetaStats(char * fname);
	void beginTime(long val);
	void updateTime(long val);
	void printTime();


};

#endif


