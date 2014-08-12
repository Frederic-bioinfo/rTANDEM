/*
loadmzparser: object for loading data using mzParser.

Purpose: Default mzXML and mzML handlers in X!Tandem are poorly implemented and cannot
read these files if, among other things, they are compressed. Instead, intercept loading
with mzParser.

Pros: smaller, faster, and easier to use than pWiz. This code stands independent of
X!Tandem code, so easier merging when X!Tandem updates. Much smaller code footprint than
previous implementations of this functionality. No need to drag boost into X!Tandem. Easier
to compile.

Cons: does not support all the formats of pWiz; However, X!Tandem supports many of these
other, less useful formats anyway.

Author: You can blame Mike Hoopmann at ISB for this.

License: The same artistic license on X!Tandem. Rather than reprint it here (AGAIN), just
look at any other source file in X!Tandem.
*/

#ifndef LOADMZPARSER_H
#define LOADMZPARSER_H

#include "loadmspectrum.h"
#include "mzParser.h"

class loadmzparser : public loadmspectrum
{
 public:
	loadmzparser(vector<mspectrum>& _vS, mspectrumcondition& _sC, mscore& _m);
	virtual ~loadmzparser();

	virtual bool get();
	int guessCharge();
	virtual bool open(string &_s);
	virtual bool open_force(string &_s);

private:

	void processSpectrum(int charge);

	vector<mspectrum>* vSpec; 
	mspectrumcondition* SpecCond;
	mscore* m_pScore; 
	long m_lLoaded;

	BasicSpectrum s;
	MzParser* mzp;

};


#endif
