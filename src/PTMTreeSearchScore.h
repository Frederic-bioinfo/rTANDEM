/* PTMTreeSearch plugin
For furter details see or if PTMTreeSearch is used in a study please cite the following article:
PTMTreeSearch: a novel two-stage tree-search algorithm with pruning rules for the identification of post-translational modification of proteins in MS/MS spectra
Attila Kertesz-Farkas, Beata Reiz, Roberto Vera, Michael P Myers, Sandor Pongor
Bioinformatics 30 (2), 234-241, 2014 

FILE version: 2013-Febr-11  // Modified for rTANDEM on January 2014
*/

#include "mscore.h"
#include "mscore_tandem.h"
#include "mprocess.h"
//#include "MinMaxHeap.hpp"

typedef vector<float> vPeaks;
typedef vector<int> vPeakIdx;

#define MAX_NODE_NUM 3000000
#define MAX_PEPT_LEN 256
#define MAX_ION_TYPE 50
#define MAX_BATCH	1


class ionScoreType{
public:
	double m_dScore;
	unsigned int m_uiMatchNum;
	float		 m_fSumInt;	
	unsigned int m_uiLastPeakIdx;
};


class PTMState{
public:
	double*			m_dKey;			//as a key for the node ordering
	unsigned int	m_iLevel;
	ionScoreType*	m_pIonScores;
	unsigned int	m_iPTMCnt;
	double			m_dPTMMass;
	unsigned int*	m_uipPTMMap;
	char*			m_bPeakMatchHist;
	unsigned int	m_uiStateCnt;
	double			m_dLastPTM;
	double*			m_dRatio;				//can be deleted
	unsigned int*	m_uiMatchNum;
	unsigned int*	m_uiActiveSpectra;
	string			m_peptide;
	unsigned int	m_uiSpectraNum;
	

	PTMState(unsigned int uiIonTypeNum, unsigned long lPeptLen, unsigned int uiSpectraNum){
		m_iLevel = 0;
		m_iPTMCnt = 0;
		m_dPTMMass = 0.0;
		m_uiStateCnt = 0;
		m_dLastPTM = 0.0;
		m_uiSpectraNum	 = uiSpectraNum;
 
		m_bPeakMatchHist= new char[lPeptLen*uiSpectraNum];
		memset(m_bPeakMatchHist, 0, sizeof(char)*lPeptLen*uiSpectraNum);

		m_uipPTMMap = new unsigned int[lPeptLen];
		memset(m_uipPTMMap, 0, sizeof(unsigned int)*lPeptLen);

		m_pIonScores = new ionScoreType[uiIonTypeNum*uiSpectraNum];
		memset(m_pIonScores, 0, sizeof(ionScoreType)*uiIonTypeNum*uiSpectraNum);
		
		m_dKey = new double[uiSpectraNum];
		memset(m_dKey, 0, sizeof(double)*uiSpectraNum);

		m_uiActiveSpectra = new unsigned int[uiSpectraNum];
		memset(m_uiActiveSpectra, 0, sizeof(unsigned int)*uiSpectraNum);

		m_dRatio = new double[uiSpectraNum];
		memset(m_dRatio, 0, sizeof(double)*uiSpectraNum);

		m_uiMatchNum = new unsigned int[uiSpectraNum];
		memset(m_uiMatchNum, 0, sizeof(unsigned int)*uiSpectraNum);
	}
	PTMState(){
	}
	~PTMState(){
	}
	void release(){
		delete m_bPeakMatchHist;
		delete m_uipPTMMap;
		delete m_pIonScores;
		delete m_dKey;
		delete m_uiActiveSpectra;
		delete m_dRatio;
		delete m_uiMatchNum;
	}
	bool operator <(const PTMState &tn) const{
		return m_dKey[0] < tn.m_dKey[0];
	}	
};
//for Gready approach
//bool grt(const PTMState &a, const PTMState b) {
//	return a.m_dKey[0] > b.m_dKey[0];
//}

typedef vector<float> spectrumInt;

class PTMTreeSearchScore : public mscore_tandem {
public:
	PTMTreeSearchScore();
	~PTMTreeSearchScore();
	enum {
		match_O	= 0x00,
		match_B = 0x01,
		match_Y = 0x02,
	} MatchType;
	unsigned long mconvert(double _m, const double _z);
	unsigned long mconvert(double _m, const long _c);
	unsigned int	m_uiIonTypeNum;
	unsigned int	m_uiPtmBound;
	double**		m_pdMods;
	unsigned int*	m_uipModNum;
	bool			m_bRefine;
	PTMState*		m_pBestGoal;
	PTMState**		m_BestGoal;
	PTMState**		m_DecoyBestGoal;
	merrors			m_errValues;
	double*			m_pdTheoPeaks[MAX_ION_TYPE];
	float*			m_pfTheoIntens[MAX_ION_TYPE];
	int				m_iChargeState[MAX_ION_TYPE];
	double*			m_pdPeaks;
	float*			m_pfIntens;
	double			m_dMatchMassTol;
	bool			m_bHitCheck;
	double*			m_dFixResidues;
	double 			dMs[MAX_PEPT_LEN];
	double 			dKs[MAX_PEPT_LEN]; 
	bool			m_bTrypticDigestion;
	unsigned long 	m_lRound;
	double 			m_dHitRatioTh;
	double			m_dPTMEntropyTh;
	
	unsigned int 	m_uiSpectraNum;
	long*			m_lChargeStateLimit;
	unsigned int*	m_uiActiveSpectra;
	double*			m_dSpectraPM;
	unsigned int	m_CurrentSpectra;
	double*			m_dDPM;
	unsigned int*	m_uiSpectrumID;
	vector<mspectrum>* m_vSpectra;
	vector<spectrumInt> m_vfSpectraInt;
	spectrumInt spectraSumInt;	
	vector<mi>**	m_vSpectraEqualS;
	unsigned int	m_uiMaxSpectraNum;
	double*			m_dAllMatch;
	double*			m_dBestScore;
	unsigned int	m_uiBatchCnt;
	long lChargeLimit;	

	
	bool add_dA(const unsigned long _t,const long _c);
	bool add_dB(const unsigned long _t,const long _c);
	bool add_dC(const unsigned long _t,const long _c);
	bool add_dY(const unsigned long _t,const long _c);
	bool add_dX(const unsigned long _t,const long _c);
	bool add_dZ(const unsigned long _t,const long _c);

	void set_parent_tolerance(double _a, double _b);
	bool get_aa(vector<maa> &_m,const size_t _a,double &_d);
	double seq_mh()	{
		return m_dSeqMH + m_pBestGoal->m_dPTMMass;
	}
	bool permute();	
	double m_dHyperValue;
	double m_dScoreValue;
	unsigned int i;
	inline double scoring(ionScoreType* _ionScore){
/*		m_dHyperValue = 1.0;
		m_dScoreValue = 0.0;
		for (i = 0; i < m_uiIonTypeNum; i++){
			m_dHyperValue *= hfactor(_ionScore[i].m_uiMatchNum);
			m_dScoreValue += _ionScore[i].m_dScore;
		}
		return m_dHyperValue*m_dScoreValue;
*/		
		// count only the number of matches.
		m_dHyperValue = 0.0;
		for (i = 0; i < m_uiIonTypeNum; i++){
			m_dHyperValue += _ionScore[i].m_uiMatchNum;
		}
		return m_dHyperValue;
	}

	double calculateMEM(char* bPeakMatchHist, unsigned int* PTMMap);
	
	bool AnchorFixedModification();
	bool ResetBestNodes();
protected:
	float score(const size_t _i);
};

/*
 * mscorefactory_ptmsearch implements a factory for creating ptmsearchscore instances.
 */
class mscorefactory_PTMTreeSearch : public mpluginfactory
{
public:
    mscorefactory_PTMTreeSearch();

    virtual mplugin* create_plugin();
};
