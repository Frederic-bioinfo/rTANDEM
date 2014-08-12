/* See loadmzparser.h for license and information */

#include "stdafx.h"
#include <algorithm>
#include "mspectrum.h"
#include "mspectrumcondition.h"
#include "loadmzparser.h"

loadmzparser::loadmzparser(vector<mspectrum>& _vS, mspectrumcondition& _sC, mscore& _m){
	vSpec=&_vS;
	SpecCond=&_sC;
	m_pScore=&_m;
	mzp = new MzParser(&s);
	m_lLoaded=0;
}
loadmzparser::~loadmzparser(){
	vSpec=NULL;
	SpecCond=NULL;
	m_pScore=NULL;
	delete mzp;
}

bool loadmzparser::get(){

	while(mzp->readSpectrum()){
		if(s.getMSLevel()!=2) continue;

		if(s.getPrecursorCharge()!=0) {
			processSpectrum(s.getPrecursorCharge());
		} else {
			// Guess charge state
			if( guessCharge() == 1)	{
				processSpectrum(1);
			}	else { // Multiple charge, most likely 2 or 3
				processSpectrum(2);
				s.setScanNum(s.getScanNum() + 100000000);
				processSpectrum(3);
				s.setScanNum(s.getScanNum() + 100000000);	
				processSpectrum(4);
				s.setScanNum(s.getScanNum() - 200000000);	
			}
		}
	}
	return true;

}

int loadmzparser::guessCharge(){
	// All this small routine does is trying to guess the charge state of the precursor
	// from the ratio of the integrals of the intensities below and above m_precursorMz
	float intBelow = 0;
	float intTotal = 0;

	for(int i=0; i<s.getPeaksCount(); i++)	{
		intTotal += (float)s[i].intensity;
		if(s[i].mz < s.getPrecursorMZ()) intBelow += (float)s[i].intensity;   
	}

	// There is no particular reason for the 0.95. It's there just
	// to compensate for the noise.... 
	if(intTotal == 0.0 || intBelow/intTotal > 0.95) return 1;
	else return 2;
}

bool loadmzparser::open(string &_s){
	return mzp->load(&_s[0]);
}
bool loadmzparser::open_force(string &_s){
	return false;
}

void loadmzparser::processSpectrum(int charge){
	mspectrum m;
	m.m_tId = s.getScanNum();
	char rtStr[16];
	sprintf(rtStr,"PT%.2lfS",s.getRTime(false));
	m.m_strRt = rtStr; //wants RT as a string...
	m.m_uiType = 0;
	if(s.getActivation()==CID) m.m_uiType = I_Y|I_B;
	else if(s.getActivation()==ETD || s.getActivation()==ETDSA)	m.m_uiType = I_C|I_Z;
	m.m_dMH = (s.getPrecursorMZ() - 1.007276)*(float)charge + 1.007276;
	m.m_fZ = (float) charge;

	mi miCurrent;
	m.clear_intensity_values();
	for(int i=0; i<s.getPeaksCount(); i++)	{
		miCurrent.m_fM = (float)s[i].mz;
		miCurrent.m_fI = (float)s[i].intensity;
		// Only push if both mass and intensity are non-zero.
		if (miCurrent.m_fM && miCurrent.m_fI ) m.m_vMI.push_back(miCurrent);
	}

	if(SpecCond->condition(m, *m_pScore) ){ // || true to force
		vSpec->push_back(m);
		m_lLoaded++;
		if(m_lLoaded == 2000)	{
			cout << ".";
			cout.flush();
			m_lLoaded = 0;
		}
	}
}

