/* PTMTreeSearch plugin
For furter details see or if PTMTreeSearch is used in a study please cite the following article:
PTMTreeSearch: a novel two-stage tree-search algorithm with pruning rules for the identification of post-translational modification of proteins in MS/MS spectra
Attila Kertesz-Farkas, Beata Reiz, Roberto Vera, Michael P Myers, Sandor Pongor
Bioinformatics 30 (2), 234-241, 2014 

FILE version: 2013-Febr-11  // Modified for rTANDEM on January 2014
*/

#include <cstring> 
#include "stdafx.h"
#include "msequence.h"
#include "mspectrum.h"
#include "msequtilities.h"
#include "xmlparameter.h"
#include "PTMTreeSearchScore.h"
#include <float.h>
#include <algorithm>
#include "mscore.h"
#include "stack_ptmtreesearch.h"
//#include "priority_deque"

// Factory instance, registers itself with the mscoremanager.
static mscorefactory_PTMTreeSearch factory;
	
mscorefactory_PTMTreeSearch::mscorefactory_PTMTreeSearch()
{
    mscoremanager::register_factory("ptmtreesearch-score", this);
}

mplugin* mscorefactory_PTMTreeSearch::create_plugin()
{
    return new PTMTreeSearchScore();
}

/*
 * mconvert converts from mass and charge to integer ion value
 * for mi vector.
 */
unsigned long PTMTreeSearchScore::mconvert(double _m, const double _z)
{
  return (unsigned long)((m_pSeqUtilFrag->m_dProton + _m/_z)*m_dWE);
}
unsigned long mscore::mconvert(double _m, const long _c)
{
/*
 * calculate the conversion factor between an m/z value and its integer value
 * as referenced in m_vsmapMI
 */
	const double dZ = (double)_c;
	return (unsigned long)((m_pSeqUtilFrag->m_dProton + _m/dZ)*m_dWidth/m_dErr);
}


PTMTreeSearchScore::PTMTreeSearchScore(){
//	mscore_tandem::mscore_tandem();
	m_bRefine = false;
	m_uiPtmBound = 0;
	m_dMatchMassTol = 0.4;

	m_uiIonTypeNum = 0;
	unsigned int i;
	for (i = 0; i < MAX_ION_TYPE; i++){
		m_pdTheoPeaks[i] = new double[MAX_PEPT_LEN];
		m_pfTheoIntens[i] = new float[MAX_PEPT_LEN];
	}
	m_dFixResidues = new double[MAX_PEPT_LEN];
	m_bTrypticDigestion = true;

	m_CurrentSpectra = 0;
	m_BestGoal = NULL;
	m_uiSpectraNum = 1;
	m_uiMaxSpectraNum = 0;
	ResetBestNodes();
	m_uiSpectraNum = 0;
	m_uiBatchCnt = 0;
	lChargeLimit = 0;
}  

PTMTreeSearchScore::~PTMTreeSearchScore(){
	unsigned int i;
	for (i = 0; i < MAX_ION_TYPE; i++){
		delete m_pdTheoPeaks[i];
		delete m_pfTheoIntens[i];
	}
	for (i = 0; i < m_uiMaxSpectraNum; i++){
		m_BestGoal[i]->release();
		delete m_BestGoal[i];
	}
	delete m_BestGoal;
	delete m_dFixResidues;	
}
bool PTMTreeSearchScore::permute(){
	if (m_bRefine == 0)
		return mscore_tandem::permute();
	return false;
}

bool PTMTreeSearchScore::ResetBestNodes(){
	unsigned int s;
	if (m_uiSpectraNum > m_uiMaxSpectraNum){
		if (m_BestGoal != NULL){
			for (s = 0; s < m_uiMaxSpectraNum; s++){
				m_BestGoal[s]->release();
				delete m_BestGoal[s];
			}
			delete m_BestGoal;
		}

		m_uiMaxSpectraNum = m_uiSpectraNum;
			
		m_BestGoal = new PTMState*[m_uiMaxSpectraNum];
		for (s = 0; s < m_uiMaxSpectraNum; s++){
			m_BestGoal[s] = new PTMState();
			m_BestGoal[s]->m_bPeakMatchHist= new char[MAX_PEPT_LEN];
			m_BestGoal[s]->m_uipPTMMap = new unsigned int[MAX_PEPT_LEN];
			m_BestGoal[s]->m_pIonScores = new ionScoreType[MAX_ION_TYPE];			
			m_BestGoal[s]->m_dKey = new double[1];
//			m_BestGoal[s]->m_bHighestPeakMatch = new bool[1];
			m_BestGoal[s]->m_uiActiveSpectra = new unsigned int[1];
			m_BestGoal[s]->m_dRatio = new double[1];
			m_BestGoal[s]->m_uiMatchNum = new unsigned int[1];
		}
	}
	//reset values to zero
	for (s = 0; s < m_uiSpectraNum; s++){
		m_BestGoal[s]->m_iLevel = 0;
		m_BestGoal[s]->m_iPTMCnt = 0;
		m_BestGoal[s]->m_dPTMMass = 0.0;
		m_BestGoal[s]->m_uiStateCnt = 0;
		m_BestGoal[s]->m_dLastPTM = 0.0;
		m_BestGoal[s]->m_dKey[0] = 0.0;
		m_BestGoal[s]->m_uiActiveSpectra[0] = 0;
		m_BestGoal[s]->m_dRatio[0] = 0.0;
		m_BestGoal[s]->m_uiMatchNum[0] = 0;
		memset(m_BestGoal[s]->m_bPeakMatchHist, 0, sizeof(char)*MAX_PEPT_LEN);
		memset(m_BestGoal[s]->m_uipPTMMap, 0, sizeof(unsigned int)*MAX_PEPT_LEN);
		memset(m_BestGoal[s]->m_pIonScores, 0, sizeof(ionScoreType)*MAX_ION_TYPE);
	}
	m_pBestGoal = m_BestGoal[0];
	return true;
}

float PTMTreeSearchScore::score(const size_t _i){
	if (m_bRefine == 0)
		return mscore_tandem::score(_i);
	m_fScore = -1.0;
	m_fHyper = -1.0;
	double dFactor = 1.0;
/*
 * return -1000.0 if there is no sequence available
 */
	if(m_pSeq == NULL)
		return -1000.0;

	if (m_lSeqLength < 2){
		return -10000.0;
	}
	prescore(_i);	
/*
 * initialize values for the protein modeling session
 */
	double dScore = (float)0.0;
	double dValue = (float)0.0;
	unsigned long lType;
	//theoretical peak generation from the database peptide sequence
	long a = 0;
	unsigned int s;

	if (m_CurrentSpectra == 0){
		//delete bestnodes;
		m_uiSpectraNum 		= m_State.m_lEqualsS;
		if (m_uiSpectraNum > MAX_BATCH){
			m_uiSpectraNum = MAX_BATCH;
		}
		if (m_uiSpectraNum > m_State.m_lEqualsS - m_uiBatchCnt*MAX_BATCH){
			m_uiSpectraNum = m_State.m_lEqualsS - m_uiBatchCnt*MAX_BATCH;
		}
		
		ResetBestNodes();
		m_dMatchMassTol = m_dErr;
		vector<mi>::iterator it;
		vector<mi>::iterator itEnd;
		vector<mi>::iterator itBegin;
		vector<mi>::reverse_iterator rit;
		vector<mi>::reverse_iterator ritEnd; //= spectra->rend();
		vector<mi>::reverse_iterator ritBegin; //= spectra->rbegin();
		
		float fMaxIntens = 0;

		bool active;
		unsigned int uiActiveSpectraCnt;
		unsigned int uiActiveSpectraNum;

		m_lChargeStateLimit = new long[m_uiSpectraNum];
		m_uiActiveSpectra 	= new unsigned int[m_uiSpectraNum];
		m_dSpectraPM 		= new double[m_uiSpectraNum];
		m_dDPM				= new double[m_uiSpectraNum];
		m_vSpectraEqualS	= new vector<mi>*[m_uiSpectraNum];
		m_dAllMatch			= new double[m_uiSpectraNum];
		m_dBestScore		= new double[m_uiSpectraNum];
		m_uiSpectrumID		= new unsigned int[m_uiSpectraNum];		
		
		uiActiveSpectraCnt 	= 0;
		lChargeLimit = 0;
		unsigned int uiSpectrumID;

		for (s = 0; s < m_uiSpectraNum; s++){
			uiSpectrumID = m_State.m_plEqualsS[s + m_uiBatchCnt*MAX_BATCH];
			m_uiSpectrumID[s] = uiSpectrumID;
						
			m_uiActiveSpectra[uiActiveSpectraCnt] = s;
			
			m_dSpectraPM[s] = m_vSpec[uiSpectrumID].m_fMH;
			m_dBestScore[s] = 0.0;
			if (m_vSpectra->at(uiSpectrumID).m_fScore > 0.0){
				m_dBestScore[s] = m_vSpectra->at(uiSpectrumID).m_fHyper/m_vSpectra->at(uiSpectrumID).m_fScore;
			}
			
			m_lChargeStateLimit[s] = (long)m_vSpec[uiSpectrumID].m_fZ;
//			long lChargeLimit = (long)m_vSpec[uiSpectrumID].m_fZ;
			// 2006.12.01: to reduce problems with high charge states, the hyperscore is
			// calculated on the 2 best fragment ion charge states. in the
			// previous versions the hyperscore used all available charge states, giving
			// a considerable advantage to highly charged parent ions
			if(m_lChargeStateLimit[s] == 1)	{
				m_lChargeStateLimit[s] = 2;
			}
			if((m_lType & T_C) || (m_lType & T_Z))	{
				if(m_lChargeStateLimit[s] > 2)	{
					m_lChargeStateLimit[s]--;
				}
			}		
			if (m_errValues.check(m_dSpectraPM[s],m_dSeqMH)){
				continue;
			}
			if (lChargeLimit < m_lChargeStateLimit[s]){
				lChargeLimit = m_lChargeStateLimit[s];
			}	
//			lChargeLimit = 3;
//			if (lChargeLimit < m_lChargeStateLimit[s]){
//				m_lChargeStateLimit[s] = lChargeLimit;
//			}	
			m_dDPM[s] = m_dSpectraPM[s] - m_dSeqMH;	
			m_vSpectraEqualS[s] = &m_vSpectra->at(uiSpectrumID).m_vMI;
			m_dAllMatch[s] = (double)(spectraSumInt[uiSpectrumID]);
			
			++uiActiveSpectraCnt;
		}
		uiActiveSpectraNum = uiActiveSpectraCnt;

		//generate theoretical peaks
 		unsigned int cnt = 0;
		memset(m_dFixResidues, 0, sizeof(double)*(m_lSeqLength+1));
		AnchorFixedModification();
		for (lType = T_Y; lType < m_lType+1; lType *= 2)	{
			if(lType & m_lType)	{
				for (a = 1; a < lChargeLimit; a++){
					m_iChargeState[cnt] = a;
					memset(&m_pdTheoPeaks[cnt][0],0,sizeof(double)*(m_lSeqLength+1));
					memset(&m_pfTheoIntens[cnt][0],0,sizeof(float)*(m_lSeqLength+1));
					//b-ion guide
//					m_pdPeaks = m_pdTheoPeaks[cnt];
//					m_pfIntens =m_pfTheoIntens[cnt]; 
					//y-ion guide
					m_pdTheoPeaks[cnt][0] = 0.0;
					m_pfTheoIntens[cnt][0] = 0.0;
					m_pdPeaks = m_pdTheoPeaks[cnt]+1;
					m_pfIntens =m_pfTheoIntens[cnt]+1; 
					if (T_Y & lType){
						add_dY(lType,a);
						m_iChargeState[cnt] *= -1;
					} else if (T_B & lType){
						add_dB(lType,a);
					} else if (T_X & lType){
						add_dX(lType,a);
						m_iChargeState[cnt] *= -1;
					} else if (T_A & lType){
						add_dA(lType,a);
					} else if (T_Z & lType){
						add_dZ(lType,a);
						m_iChargeState[cnt] *= -1;
					} else if (T_C & lType){
						add_dC(lType,a);
					}
					cnt++;
				}
			}
		}
		m_uiIonTypeNum = cnt;

		//PTMSearch local variables
		unsigned int state_cnt=1;
		unsigned int uiLevel;
		unsigned int uiModIt;
		double dPTM;
		double dPTMMass;
		double dPTMMiss;
		unsigned int uiPTMCnt;
		double dMod;
		double dPTMEntropy;
		unsigned int uiPtmBound;
		unsigned int uiMatchNumPoss;
		unsigned int MatchOffset;
		unsigned int ScoreOffset;
		unsigned int currentSpectra;
		double dBestBestKey = 0.0;	
		unsigned int newstate_size;
		char cBYMatch;
		unsigned int uiLengthLevel;
		unsigned int uiActiveLengthOffset;
		
		if (m_uiPtmBound < (m_lSeqLength)/3){
			uiPtmBound = m_uiPtmBound;
		} else {
			uiPtmBound = (m_lSeqLength)/3;
		}
		
		PTMState* state;
		PTMState* newstate;
		ionScoreType* currentScores;
		ionScoreType* stateScores;
		
		//ptmsearch init. create root state and initialize the best_goal	
		stack <PTMState, vector<PTMState> > PTMStateDeque;
//		int m_iQueueLimit = m_lSeqLength;
//		priority_deque <PTMState, deque<PTMState> > PTMStateDeque;

		PTMState root(m_uiIonTypeNum, m_lSeqLength, uiActiveSpectraNum);
		
		memcpy(root.m_uiActiveSpectra, m_uiActiveSpectra, sizeof(unsigned int)*uiActiveSpectraNum);
		//y-ion guide
		root.m_iLevel = m_lSeqLength-1;
		
		PTMStateDeque.push(root);
		
		newstate = new PTMState(m_uiIonTypeNum, m_lSeqLength, uiActiveSpectraNum);
		newstate_size = uiActiveSpectraNum;
		
		float fMatchedInt;
		while (!PTMStateDeque.empty()){
		
			state = PTMStateDeque.access_top();

			uiLevel = state->m_iLevel;
			uiModIt = (state->m_uipPTMMap[uiLevel])++;
			uiPTMCnt = state->m_iPTMCnt; 

			if (uiModIt > 0){
				if (uiPTMCnt == uiPtmBound){
					state->release();
					PTMStateDeque.pop();
					continue;
				}
				uiPTMCnt++;
			}
			if (m_dFixResidues[uiLevel] != 0.0 && uiModIt > 0){
				state->release();
				PTMStateDeque.pop();
				continue;
			}
			if (m_dFixResidues[uiLevel] != 0.0 && uiModIt == 0){
				state->m_dLastPTM = m_dFixResidues[uiLevel];
			}
			dPTM = m_pdMods[m_pSeq[uiLevel]][uiModIt];

			//check whether the the current modification mass would neutralize the previous modification
			if ( uiPTMCnt > 0 && fabs((float)(state->m_dLastPTM + dPTM)) < m_dMatchMassTol){
				if (uiModIt == m_uipModNum[m_pSeq[uiLevel]]){ //this was the last modification, so delete this state
					state->release();
					PTMStateDeque.pop();
				}
				continue;
			}

			//calculate the matches and scores
			uiActiveSpectraNum = state->m_uiSpectraNum;
			if (newstate_size < uiActiveSpectraNum){
				newstate->release();
				delete newstate;
				newstate = new PTMState(m_uiIonTypeNum, m_lSeqLength, uiActiveSpectraNum);
				newstate_size = uiActiveSpectraNum;
			}
			memcpy(newstate->m_pIonScores, state->m_pIonScores, sizeof(ionScoreType)*m_uiIonTypeNum*uiActiveSpectraNum);
			
			currentScores = newstate->m_pIonScores - m_uiIonTypeNum;
			stateScores = state->m_pIonScores;
			
			uiActiveSpectraCnt = 0;
			dPTMMass  = state->m_dPTMMass + dPTM;
			//b-ion guide
//			uiLengthLevel = 2*(m_lSeqLength - uiLevel - 2);
			//y-ion guide
			uiLengthLevel = 2*(uiLevel-1);
			
			uiActiveLengthOffset = 0;
			for (s = 0; s < uiActiveSpectraNum; s++){
								
				currentSpectra = state->m_uiActiveSpectra[s];
				currentScores += m_uiIonTypeNum;

				dPTMMiss = m_dDPM[currentSpectra] - dPTMMass;

				itBegin = m_vSpectraEqualS[currentSpectra]->begin();
				itEnd = m_vSpectraEqualS[currentSpectra]->end();
				ritBegin = m_vSpectraEqualS[currentSpectra]->rbegin();
				ritEnd = m_vSpectraEqualS[currentSpectra]->rend();
				
				cBYMatch = match_O;
				fMatchedInt = 0.0;
				
				for (i = 0; i < m_uiIonTypeNum; i++){
					fMatchedInt += currentScores[i].m_fSumInt;									
					if (m_pdTheoPeaks[i][uiLevel] == 0.0) continue;
					if (abs(m_iChargeState[i]) >= m_lChargeStateLimit[currentSpectra]) continue;
					if (m_iChargeState[i] > 0){
						//b-ion guide
//						dMod = dPTMMass/m_iChargeState[i];
						//y-ion guide
						dMod = dPTMMiss/m_iChargeState[i];
						for (rit = ritBegin + currentScores[i].m_uiLastPeakIdx; rit != ritEnd; ++rit){
							if (fabs((float)(rit->m_fM - m_pdTheoPeaks[i][uiLevel] - dMod)) <  m_dMatchMassTol){
								++currentScores[i].m_uiMatchNum;
								currentScores[i].m_dScore += rit->m_fI*m_pfTheoIntens[i][uiLevel];
								currentScores[i].m_fSumInt += rit->m_fI;
								fMatchedInt += rit->m_fI;								
								cBYMatch |= match_B;
								break;
							}
							if (m_pdTheoPeaks[i][uiLevel] + dMod > rit->m_fM)
								break;
							++currentScores[i].m_uiLastPeakIdx;
						}
					} else {
						//b-ion guide
//						dMod = dPTMMiss/(abs(m_iChargeState[i]));
						//y-ion guide
						dMod = dPTMMass/(abs(m_iChargeState[i]));
						//Notice that this is reverse iterator!
						for (it = itBegin + currentScores[i].m_uiLastPeakIdx; it < itEnd; ++it){
							if (fabs((double)(it->m_fM - m_pdTheoPeaks[i][uiLevel]-dMod)) < m_dMatchMassTol){
								++currentScores[i].m_uiMatchNum;
								currentScores[i].m_dScore += it->m_fI*m_pfTheoIntens[i][uiLevel];
								currentScores[i].m_fSumInt += it->m_fI;								
								fMatchedInt += it->m_fI;
								cBYMatch |= match_Y;
								break;
							}
							if (m_pdTheoPeaks[i][uiLevel] + dMod < it->m_fM)
								break;
							++currentScores[i].m_uiLastPeakIdx;
						}
					}
					if (uiModIt > 0){
						stateScores[i].m_uiLastPeakIdx = currentScores[i].m_uiLastPeakIdx;
					}								
				}
				stateScores+=m_uiIonTypeNum;
				// A modified peak has to have a match
				if (m_bHitCheck == true && uiModIt >0 ){
					if (uiLevel == 0){
//						if (!cBYMatch){ //b-ion guide
						if (!(state->m_bPeakMatchHist[s*m_lSeqLength + 1] )){ //y-ion guide
							continue;
						}
					} else if (uiLevel == m_lSeqLength-1){
//						if ( !(state->m_bPeakMatchHist[s*m_lSeqLength + uiLevel -1] ) ){ //b-ion guide
						if ( !cBYMatch ){ //y-ion guide
							continue;
						}					
					} else				
//						if ( !(cBYMatch & match_B) && (!(cBYMatch & match_Y) && !(state->m_bPeakMatchHist[s*m_lSeqLength + uiLevel -1] & match_Y)) ){ //b-ion guide
						if ( !(cBYMatch & match_B) && (!(cBYMatch & match_Y) && !(state->m_bPeakMatchHist[s*m_lSeqLength + uiLevel +1] & match_B)) ){ //y-ion guide
							continue;
						}
				}
				//optimistic estimation about the HitRatio at the end.
				newstate->m_uiMatchNum[uiActiveSpectraCnt] = state->m_uiMatchNum[s];
				if (cBYMatch & match_B){
					++newstate->m_uiMatchNum[uiActiveSpectraCnt];
				}
				if (cBYMatch & match_Y){
					++newstate->m_uiMatchNum[uiActiveSpectraCnt];
				}

				uiMatchNumPoss = newstate->m_uiMatchNum[uiActiveSpectraCnt] + uiLengthLevel;
				if (  uiMatchNumPoss < 3){
					continue;  
				}
				if (hfactor(uiMatchNumPoss) < m_dBestScore[currentSpectra]){
					continue;
				}
				if (dBestBestKey*0.6 > uiMatchNumPoss){
					continue;
				}
				if (m_lRound > 0 && uiMatchNumPoss < m_BestGoal[currentSpectra]->m_dKey[0]){
					continue;
				}
//				double dMatchNumPoss= (double)((fMatchedInt + m_vfSpectraInt[m_uiSpectrumID[currentSpectra]][2*m_lSeqLength]-m_vfSpectraInt[m_uiSpectrumID[currentSpectra]][2*(m_lSeqLength-uiLevel)]));
				unsigned int uiMaxPeakLen = 2*m_lSeqLength-2;
				if (uiMaxPeakLen >= m_vfSpectraInt[m_uiSpectrumID[currentSpectra]].size()){
					uiMaxPeakLen = m_vfSpectraInt[m_uiSpectrumID[currentSpectra]].size()-1;
				}
				double dMatchNumPoss= (double)(fMatchedInt + m_vfSpectraInt[m_uiSpectrumID[currentSpectra]][uiMaxPeakLen]*2*uiLevel/(uiMaxPeakLen));
				newstate->m_dRatio[uiActiveSpectraCnt] = (double)dMatchNumPoss/m_vfSpectraInt[m_uiSpectrumID[currentSpectra]][uiMaxPeakLen];

				if (  newstate->m_dRatio[uiActiveSpectraCnt] < m_dHitRatioTh){
					continue;
				}
				// Peak with highest intensity always has to mutch to a theoretical peak.
				newstate->m_uiActiveSpectra[uiActiveSpectraCnt] = currentSpectra;//state->m_uiActiveSpectra[s];
				newstate->m_dKey[uiActiveSpectraCnt] = newstate->m_uiMatchNum[uiActiveSpectraCnt];//scoring(currentScores);
//				memcpy(&newstate->m_bPeakMatchHist[uiActiveLengthOffset], &state->m_bPeakMatchHist[s*m_lSeqLength],sizeof(char)*(uiLevel)); 
				memcpy(&newstate->m_bPeakMatchHist[uiActiveLengthOffset], &state->m_bPeakMatchHist[s*m_lSeqLength],sizeof(char)*(m_lSeqLength)); 
				newstate->m_bPeakMatchHist[uiActiveLengthOffset + uiLevel] = cBYMatch;
				if (uiActiveSpectraCnt < s){
					memcpy(newstate->m_pIonScores+(uiActiveSpectraCnt*m_uiIonTypeNum), currentScores, sizeof(ionScoreType)*m_uiIonTypeNum);
				}				
				++uiActiveSpectraCnt;
				uiActiveLengthOffset += m_lSeqLength;
			}
			
			if (uiActiveSpectraCnt == 0){
				if (uiModIt == m_uipModNum[m_pSeq[uiLevel]]){
					state->release();
					PTMStateDeque.pop();
				}				
				continue;	
			}
//			if (uiLevel == m_lSeqLength-1){ //last level, end of recursion for b-ion guide
			if (uiLevel == 0){ //last level, end of recursion for y-ion guide
				for (s = 0; s < uiActiveSpectraCnt; s++){

					currentSpectra = newstate->m_uiActiveSpectra[s];
					
					if (m_errValues.check(m_dSpectraPM[currentSpectra],m_dSeqMH+dPTMMass) && uiPTMCnt > 0){		//mass check and at least there is one PTM.
						if (newstate->m_dKey[s] > m_BestGoal[currentSpectra]->m_dKey[0] || 										// better score
							(newstate->m_dKey[s] == m_BestGoal[currentSpectra]->m_dKey[0] && uiPTMCnt < m_BestGoal[currentSpectra]->m_iPTMCnt ) || 	// equal score but less PTMs
							(newstate->m_dKey[s] == m_BestGoal[currentSpectra]->m_dKey[0] && uiPTMCnt == m_BestGoal[currentSpectra]->m_iPTMCnt && newstate->m_dRatio[s] > m_BestGoal[currentSpectra]->m_dRatio[0])  // equal score and eqaul PTMs but better hit ratio
							){ 
							// calculate statistics for the modifications 
							MatchOffset = s*m_lSeqLength;
							
							dPTMEntropy = calculateMEM(&newstate->m_bPeakMatchHist[MatchOffset], state->m_uipPTMMap);
							if ( dPTMEntropy > m_dPTMEntropyTh || m_lRound > 0){
								m_BestGoal[currentSpectra]->m_dKey[0] = newstate->m_dKey[s];
								memcpy(m_BestGoal[currentSpectra]->m_bPeakMatchHist, &newstate->m_bPeakMatchHist[MatchOffset], sizeof(char)*m_lSeqLength);
								memcpy(m_BestGoal[currentSpectra]->m_uipPTMMap, state->m_uipPTMMap, sizeof(unsigned int)*m_lSeqLength);
								memcpy(m_BestGoal[currentSpectra]->m_pIonScores,  newstate->m_pIonScores+(s*m_uiIonTypeNum), sizeof(ionScoreType)*m_uiIonTypeNum);
								m_BestGoal[currentSpectra]->m_dPTMMass = dPTMMass;
								m_BestGoal[currentSpectra]->m_iPTMCnt  = uiPTMCnt;
								m_BestGoal[currentSpectra]->m_dRatio[0]   = newstate->m_dRatio[s];
								m_BestGoal[currentSpectra]->m_uiMatchNum[0]   = newstate->m_uiMatchNum[s];
								m_BestGoal[currentSpectra]->m_dLastPTM = dPTMEntropy + 1.0;
								m_BestGoal[currentSpectra]->m_peptide = m_pSeq;

								if (dBestBestKey < newstate->m_dKey[s]){
									dBestBestKey = newstate->m_dKey[s];
								}	
							}
						}
					}
				}
				if (uiModIt == m_uipModNum[m_pSeq[uiLevel]]){
					state->release();
					PTMStateDeque.pop();
				}
			} else {  //generate the next child

				newstate->m_iLevel = uiLevel-1;
				newstate->m_iPTMCnt = uiPTMCnt;
				newstate->m_dPTMMass = dPTMMass;
				newstate->m_uiStateCnt = state_cnt++;
				newstate->m_uiSpectraNum = uiActiveSpectraCnt;
				memcpy(newstate->m_uipPTMMap, state->m_uipPTMMap, sizeof(unsigned int)*m_lSeqLength);

				if (uiModIt>0){
					newstate->m_dLastPTM = dPTM;
				} else {
					newstate->m_dLastPTM = state->m_dLastPTM;
				}

				if (uiModIt == m_uipModNum[m_pSeq[uiLevel]]){
					state->release();
					PTMStateDeque.pop();
				}
				PTMStateDeque.push(*newstate);
				delete newstate;
				newstate = new PTMState(m_uiIonTypeNum, m_lSeqLength, uiActiveSpectraCnt);
				newstate_size = uiActiveSpectraCnt;
			}
/*			//using with greedy approach
			if (PTMStateDeque.size() > m_iQueueLimit && PTMStateDeque.size() > 2){
				make_heap(PTMStateDeque.begin(), PTMStateDeque.end(), grt);
				PTMStateDeque.access_top()->release();
				PTMStateDeque.pop();
				make_heap(PTMStateDeque.begin(), PTMStateDeque.end());
			}
*/			if (state_cnt ==  MAX_NODE_NUM){  //hard limit
				while (!PTMStateDeque.empty()){
					state = PTMStateDeque.access_top();
					state->release();
					PTMStateDeque.pop();
				}
				break;
			}
		}
		newstate->release();
		delete newstate;
		delete m_uiActiveSpectra;
		delete m_dSpectraPM;
		delete m_dDPM;
		delete m_vSpectraEqualS;
		delete m_dAllMatch;
		delete m_dBestScore;
	} //END 	if (m_CurrentSpectra == 0){

	m_pBestGoal = m_BestGoal[m_CurrentSpectra];
	unsigned long lValue = 0;
	unsigned long lValueTotal = 0;
	unsigned long lS = S_Y;	
	dScore = 0;
	i = 0;
	for (lType = T_Y; lType < m_lType+1; lType *= 2){
		lValueTotal = 0;
		dValue = 0.0;
		if(lType & m_lType)	{
			for (a = 1; a < lChargeLimit; a++){
				if (a < m_lChargeStateLimit[m_CurrentSpectra]) {
					dScore += m_pBestGoal->m_pIonScores[i].m_dScore;
					dValue += m_pBestGoal->m_pIonScores[i].m_dScore;
					lValueTotal += m_pBestGoal->m_pIonScores[i].m_uiMatchNum;
					dFactor *= hfactor(m_pBestGoal->m_pIonScores[i].m_uiMatchNum);
				}
				i++;
			}
		}
		m_pfScore[lS] = (float) dValue;
		m_plCount[lS] = lValueTotal;
		lS++;
	}

	dScore *= sfactor();
	m_fScore = (float)dScore;
	dFactor *= dScore;
	if(dFactor > FLT_MAX)	{
		m_fHyper = FLT_MAX;
	}
	else	{
		m_fHyper = (float)dFactor;
	}

/*
 * returning 1.0 for a zero score makes the logic in mprocess easier. see mprocess:create_score
 * to see why.
 */
	if(dScore == 0.0)	{
		dScore = 1.0;
	}

	char *pS = strchr(m_pSeq,'s');
	char *pT = strchr(m_pSeq,'t');

	if(m_bPhosphoBias && m_fHyper < FLT_MAX && (pS || pT))	{
		int iST = 0;
		char *pV = strstr(m_pSeq,"sP");
		while(pV)	{
			iST++;
			pV++;
			pV = strstr(pV,"sP");
		}
		pV = strstr(m_pSeq,"tP");
		while(pV)	{
			iST++;
			pV++;
			pV = strstr(pV,"tP");
		}
		double dV = (double)m_fHyper*(1.0 + 0.001*iST);
		if(dV < FLT_MAX)	{
			 m_fHyper = (float)dV;
		}
		double dNeutral = 0.0;
		unsigned long lNeutral = 0;
		m_dWE = (double)(m_dWidth/m_dErr);
		if((pS && m_pSeqUtilFrag->m_bPhosphoSerine) || (pT && m_pSeqUtilFrag->m_bPhosphoThreonine))	{
			nMap::iterator itMap = m_seqUtil.m_mapNeutralLoss.find(m_lId);
			if(itMap != m_seqUtil.m_mapNeutralLoss.end())	{
				double dHyper = (double)m_fHyper * (double)itMap->second;
				if(dHyper < FLT_MAX)	{
					m_fHyper = (float)dHyper;
				}
			}
			else	{
				dNeutral = m_dSeqMH - (79.966331 + m_seqUtil.m_dWater) - m_seqUtil.m_dProton;
#ifdef PLUGGABLE_SCORING
				lNeutral = mconvert(dNeutral,(long)m_vSpec[m_lId].m_fZ);
#else
				lNeutral = mconvert(dNeutral,(double)m_vSpec[m_lId].m_fZ);
#endif
				float fV = ion_check(lNeutral,m_lId);
				if(fV >= 20.0)	{
					double dHyper = (double)m_fHyper * (double)1.001;
					if(dHyper < FLT_MAX)	{
						m_fHyper = (float)dHyper;
					}
					m_seqUtil.m_mapNeutralLoss.insert(nMap::value_type(m_lId,(float)1.001));
				}
				else	{
					m_seqUtil.m_mapNeutralLoss.insert(nMap::value_type(m_lId,(float)1.0));
				}
			}
		}
	}
	
	m_CurrentSpectra++;
	if (m_CurrentSpectra == m_uiSpectraNum){
		m_CurrentSpectra = 0;
		delete m_lChargeStateLimit;
		if (m_uiSpectraNum + m_uiBatchCnt*MAX_BATCH < m_State.m_lEqualsS){
			++m_uiBatchCnt;
		} else {
			m_uiBatchCnt = 0;
		}
	}
	return (float) dScore;
}
double PTMTreeSearchScore::calculateMEM(char* bPeakMatchHist, unsigned int* PTMMap){
	
	double entropy = 0.0;
	unsigned int i;
	unsigned int cnt = 0;
	int hit, Bhit, Yhit;
	int len;
	int modNum;
	double coeff = 0.0;
	
	//b peaks
	Bhit = 0;
	Yhit = 0;	
	hit = 0;
	len = 0;
	for (i = 0; i < m_lSeqLength-1; ++i){
		if (PTMMap[i] > 1){
			dKs[cnt] = 0.0;
			dMs[cnt] = 0.0;
			if (hit > 0){
				dKs[cnt] = (double)hit;
				dMs[cnt] = (double)hit/len;
				coeff += dMs[cnt];
			}
			++cnt;	
			hit = 0;
			len = 0;			
		} 
		++len;
		if (bPeakMatchHist[i] & match_B ){
			++hit;
			++Bhit;
		}
		if (bPeakMatchHist[i] & match_Y){
			++Yhit;
		}		
	}
	dKs[cnt] = 0.0;
	dMs[cnt] = 0.0;
	if (hit>0){
		dKs[cnt] = (double)hit;
		dMs[cnt] = (double)hit/len;
		coeff += dMs[cnt];
	}
	if (coeff > 0.0){
		if (Yhit >0){
			coeff *= 2;
		}
		for (i = 0; i <= cnt; i++){
			dMs[i] /= coeff;
		}
	}

	//y peaks
	++cnt;
	modNum = cnt;
	hit = 0;
	len = 0;
	coeff = 0.0;
	for (i = 0; i < m_lSeqLength-1; ++i){
		if (PTMMap[i] > 1){
			dKs[cnt] = 0.0;
			dMs[cnt] = 0.0;		
			if (hit > 0){
				dKs[cnt] = (double)hit;
				dMs[cnt] = (double)hit/len;
				coeff += dMs[cnt];
			}
			++cnt;	
			hit = 0;
			len = 0;			
		} 
		++len;
		if (bPeakMatchHist[i] & match_Y){
			++hit;
		}
	}
	dKs[cnt] = 0.0;
	dMs[cnt] = 0.0;
	if (hit>0){
		dKs[cnt] = (double)hit;
		dMs[cnt] = (double)hit/len;
		coeff += dMs[cnt];
	}
	
	if (coeff > 0.0){
		if (Bhit>0) {
			coeff *= 2;
		}
		for (i = modNum; i <= cnt; i++){
			dMs[i] /= coeff;
		}
	}
	
	dMs[0] += dMs[cnt];
	dKs[0] += dKs[cnt];	
	dKs[cnt] = 0.0;
	dMs[modNum-1] += dMs[modNum];
	dKs[modNum-1] += dKs[modNum];	
	dKs[modNum] = 0.0;
	
#ifdef MSVC
	//entropy;
	for (i = 0; i < cnt; i++){
		if (dKs[i] < 0.1) continue; 
		entropy += dMs[i]*log(dMs[i])/log(2.0);
	}
	return entropy/(log(1.0/((double)(cnt-1)))/log(2.0));
#else
	//entropy;
	for (i = 0; i < cnt; i++){
		if (dKs[i] < 0.1) continue; 
		entropy += dMs[i]*log2(dMs[i]);
	}
	return entropy/log2(1.0/((double)(cnt-1)));
#endif
}

void PTMTreeSearchScore::set_parent_tolerance(double _a, double _b){

	vector<mspectrumdetails>::iterator itDetails;
	vector<mspectrumdetails>::iterator itEnd = m_vDetails.end();

	for (itDetails = m_vDetails.begin(); itDetails != itEnd; itDetails++){
		itDetails->m_dL -= (float)_b;
		itDetails->m_dU += (float)_a;
	}
}

bool PTMTreeSearchScore::get_aa(vector<maa> &_m,const size_t _a,double &_d){

	bool rValue = mscore_tandem::get_aa(_m, _a, _d);
	if (m_pBestGoal->m_iPTMCnt != 0){
		double dDelta = 0.0;
		maa aaValue;
		for (unsigned long i=0; i <m_lSeqLength; i++){
			if (m_pBestGoal->m_uipPTMMap[i] > 1){
				aaValue.m_cRes = m_pSeq[i];
				aaValue.m_lPos = (long)_a + i;
				aaValue.m_dMod = m_pdMods[m_pSeq[i]][m_pBestGoal->m_uipPTMMap[i]-1];
				aaValue.m_strId = m_pBestGoal->m_uipPTMMap[i]-1;
				aaValue.m_fPval = m_pBestGoal->m_dLastPTM;
				aaValue.m_dRatio = m_pBestGoal->m_dRatio[0];
				aaValue.m_dMatch = m_pBestGoal->m_uiMatchNum[0];
				aaValue.m_dPos = i;
				aaValue.m_strId = m_pBestGoal->m_peptide;
				dDelta += aaValue.m_dMod;
				_m.push_back(aaValue);
			}
		}
		if (_d != 1000000.0){
			_d += dDelta;
		}
	}
	return rValue;
}

bool PTMTreeSearchScore::AnchorFixedModification(){
	if(m_bIsN)	{
		m_dFixResidues[0] += m_pSeqUtilFrag->m_fNT;
	}
	m_dFixResidues[0] += (m_pSeqUtilFrag->m_dCleaveN - m_pSeqUtilFrag->m_dCleaveNdefault);	
	if(m_Term.m_lN)	{
		m_dFixResidues[0] += m_pSeqUtilFrag->m_pdAaMod['['];
	}
	m_dFixResidues[0] += m_pSeqUtilFrag->m_pdAaFullMod['['];

	unsigned int a = 0;
	char tC;
	while(a < m_lSeqLength)	{
		tC = m_pSeq[a];
		m_dFixResidues[a] += m_pSeqUtilFrag->m_pdAaMod[tC] + m_pSeqUtilFrag->m_pdAaFullMod[tC];
		if(m_pSeqUtilFrag->m_bPrompt)	{
			m_dFixResidues[a] += m_pSeqUtilFrag->m_pdAaPrompt[tC];
		}
		if (m_pSeqUtilFrag->m_bSequenceMods)	{
			SMap::iterator itSeq = m_pSeqUtilFrag->m_mapMods.find(m_tSeqPos+a);
			if(itSeq != m_pSeqUtilFrag->m_mapMods.end()){
				m_dFixResidues[a] += itSeq->second;
			}
		}
		a++;		
	}
/*
 * deal with non-hydrolytic cleavage and deal with protein C-teminus
 */
	a = m_lSeqLength-1;
	m_dFixResidues[a] += (m_pSeqUtilFrag->m_dCleaveC - m_pSeqUtilFrag->m_dCleaveCdefault);
	if(m_Term.m_lC)	{
		m_dFixResidues[a] += m_pSeqUtilFrag->m_pdAaMod[']'];
	}
	m_dFixResidues[a] += m_pSeqUtilFrag->m_pdAaFullMod[']'];

	if(m_bIsC)	{
		m_dFixResidues[a] += m_pSeqUtilFrag->m_fCT;
	}	
	return true;
}
/*
 * create list of non-zero predicted intensity values for a-ions and their
 * integer converted m/z values
 */
bool PTMTreeSearchScore::add_dA(const unsigned long _t,const long _c)
{
	unsigned long a = 0;
/*
 * get the conversion factor between a straight sequence mass and an a-ion
 */
	double dValue = m_pSeqUtilFrag->m_dA;
/*
 * deal with protein N-terminus
 */
	if(m_bIsN)	{
		dValue += m_pSeqUtilFrag->m_fNT;		
	}
/*
 * deal with non-hydrolytic cleavage
 */
	dValue += (m_pSeqUtilFrag->m_dCleaveN - m_pSeqUtilFrag->m_dCleaveNdefault);
	if(m_Term.m_lN)	{
		dValue += m_pSeqUtilFrag->m_pdAaMod['['];
	}
	dValue += m_pSeqUtilFrag->m_pdAaFullMod['['];
//	unsigned long lValue = 0;
/*
 * calculate the conversion factor between an m/z value and its integer value
 * as referenced in m_vsmapMI
 */
	size_t tC = 0;
	float *pfScore = m_pSeqUtilFrag->m_pfAScore;
	unsigned long lCount = 0;
/*
 * from N- to C-terminus, calcuate fragment ion m/z values and store the results
 * look up appropriate scores from m_pSeqUtilFrag->m_pfAScore
 */
	const unsigned long tPos = (unsigned long) m_tSeqPos;
	m_dWE = m_dWidth/m_dErr;
	const double dZ = (double)_c;
	while(a < m_lSeqLength)	{
		tC = m_pSeq[a];
//#ifdef PLUGGABLE_SCORING
//		dValue += m_pSeqUtilFrag->getAaMass((char)tC, tPos+a);
//		lValue = mconvert(dValue, _c);
//#else
		dValue += m_pSeqUtilFrag->m_pdAaMass[tC];
		dValue += m_pSeqUtilFrag->m_pdAaMod[tC];
		dValue += m_pSeqUtilFrag->m_pdAaFullMod[tC];
		m_dFixResidues[a] = m_pSeqUtilFrag->m_pdAaMod[tC] + m_pSeqUtilFrag->m_pdAaFullMod[tC];
		if(m_pSeqUtilFrag->m_bPrompt)	{
			dValue += m_pSeqUtilFrag->m_pdAaPrompt[tC];
		}
		if (m_pSeqUtilFrag->m_bSequenceMods)	{
			SMap::iterator itSeq = m_pSeqUtilFrag->m_mapMods.find(tPos+a);
			if(itSeq != m_pSeqUtilFrag->m_mapMods.end()){
				dValue += itSeq->second;
			}
		}
//		lValue = mconvert(dValue, dZ);
//#endif
		//modified by AKF
		m_pdPeaks[a] = m_pSeqUtilFrag->m_dProton + dValue/dZ;   
		m_pfIntens[a] = pfScore[tC];							

//		m_plSeq[lCount] = lValue;
//		m_pfSeq[lCount] = pfScore[tC];
		lCount++;
		a++;
	}
/*
 * set the next integer mass value to 0: this marks the end of the array 
 */
//	m_pdPeaks[a] = 0.0;
//	m_plSeq[lCount] = 0;
	return true;
}
/*
 * create list of non-zero predicted intensity values for b-ions and their
 * integer converted m/z values
 */
bool PTMTreeSearchScore::add_dB(const unsigned long _t,const long _c)
{
	unsigned long a = 0;
/*
 * get the conversion factor between a straight sequence mass and a b-ion
 */
	double dValue = m_pSeqUtilFrag->m_dB;
/*
 * deal with protein N-terminus
 */
	if(m_bIsN)	{
		dValue += m_pSeqUtilFrag->m_fNT;		
	}
/*
 * deal with non-hydrolytic cleavage
 */
	dValue += (m_pSeqUtilFrag->m_dCleaveN - m_pSeqUtilFrag->m_dCleaveNdefault);
	if(m_Term.m_lN)	{
		dValue += m_pSeqUtilFrag->m_pdAaMod['['];
	}
	dValue += m_pSeqUtilFrag->m_pdAaFullMod['['];
/*
 * calculate the conversion factor between an m/z value and its integer value
 * as referenced in m_vsmapMI
 */
	long lCount = 0;
	float *pfScore = m_pSeqUtilFrag->m_pfBScore;
	float *pfScorePlus = m_pSeqUtilFrag->m_pfYScore;
/*
 * from N- to C-terminus, calcuate fragment ion m/z values and store the results
 * look up appropriate scores from m_pSeqUtilFrag->m_pfBScore
 */
	const unsigned long tPos = (unsigned long) m_tSeqPos;
	size_t tC = 0;
	m_dWE = m_dWidth/m_dErr;
	const double dZ = (double)_c;
	while(a < m_lSeqLength-1)	{
		tC = m_pSeq[a];
//#ifdef PLUGGABLE_SCORING
//		dValue += m_pSeqUtilFrag->getAaMass((char)tC, tPos+a);
//		lValue = mconvert(dValue, _c);
//#else
		dValue +=m_pSeqUtilFrag-> m_pdAaMass[tC];
		dValue += m_pSeqUtilFrag->m_pdAaMod[tC];
		dValue += m_pSeqUtilFrag->m_pdAaFullMod[tC];
		if(m_pSeqUtilFrag->m_bPrompt)	{
			dValue += m_pSeqUtilFrag->m_pdAaPrompt[tC];
		}
		if (m_pSeqUtilFrag->m_bSequenceMods)	{
			SMap::iterator itSeq = m_pSeqUtilFrag->m_mapMods.find(tPos+a);
			if(itSeq != m_pSeqUtilFrag->m_mapMods.end())
				dValue += itSeq->second;
		}
//		lValue = mconvert(dValue, dZ);
//#endif
		//modified by AKF
		m_pdPeaks[a] = m_pSeqUtilFrag->m_dProton + dValue/dZ;
		m_pfIntens[a] = pfScore[tC]*pfScorePlus[m_pSeq[a+1]];
		if(a == 1)	{
			if(m_pSeq[1] == 'P')	{
				m_pfIntens[a] *= 10;
			}
			else	{
				m_pfIntens[a] *= 3;
			}
		}

/*		m_plSeq[lCount] = lValue;
//		m_pfSeq[lCount] = pfScore[tC]*pfScorePlus[m_pSeq[a+1]];
		if(a == 1)	{
			if(m_pSeq[1] == 'P')	{
				m_pfSeq[lCount] *= 10;
			}
			else	{
				m_pfSeq[lCount] *= 3;
			}
		}
*/		lCount++;
		a++;
	}
//	m_plSeq[lCount] = 0;
	return true;
}
/*
 * create list of non-zero predicted intensity values for c-ions and their
 * integer converted m/z values
 */
bool PTMTreeSearchScore::add_dC(const unsigned long _t,const long _c)
{
	unsigned long a = 0;
/*
 * get the conversion factor between a straight sequence mass and a b-ion
 */
	double dValue = m_pSeqUtilFrag->m_dC;
/*
 * deal with protein N-terminus
 */
	if(m_bIsN)	{
		dValue += m_pSeqUtilFrag->m_fNT;
	}
/*
 * deal with non-hydrolytic cleavage
 */
	dValue += (m_pSeqUtilFrag->m_dCleaveN - m_pSeqUtilFrag->m_dCleaveNdefault);
	if(m_Term.m_lN)	{
		dValue += m_pSeqUtilFrag->m_pdAaMod['['];
	}
	dValue += m_pSeqUtilFrag->m_pdAaFullMod['['];

//	unsigned long lValue = 0;
/*
 * calculate the conversion factor between an m/z value and its integer value
 * as referenced in m_vsmapMI
 */
	size_t tC = 0;
	long lCount = 0;
	float *pfScore = m_pSeqUtilFrag->m_pfBScore;
	float *pfScorePlus = m_pSeqUtilFrag->m_pfYScore;
/*
 * from N- to C-terminus, calcuate fragment ion m/z values and store the results
 * look up appropriate scores from m_pSeqUtilFrag->m_pfBScore
 */
	m_dWE = m_dWidth/m_dErr;
	double dZ = (double)_c;
	const unsigned long tPos = (unsigned long) m_tSeqPos;
	while(a < m_lSeqLength-2)	{
		tC = m_pSeq[a];
//#ifdef PLUGGABLE_SCORING
//		dValue += m_pSeqUtilFrag->getAaMass((char)tC, tPos+a);
//		lValue = mconvert(dValue, _c);
//#else
		dValue +=m_pSeqUtilFrag-> m_pdAaMass[tC];
		dValue += m_pSeqUtilFrag->m_pdAaMod[tC];
		dValue += m_pSeqUtilFrag->m_pdAaFullMod[tC];
		if(m_pSeqUtilFrag->m_bPrompt)	{
			dValue += m_pSeqUtilFrag->m_pdAaPrompt[tC];
		}
		if (m_pSeqUtilFrag->m_bSequenceMods)	{
			SMap::iterator itSeq = m_pSeqUtilFrag->m_mapMods.find(tPos+a);
			if(itSeq != m_pSeqUtilFrag->m_mapMods.end())
				dValue += itSeq->second;
		}
//		lValue = mconvert(dValue, dZ);
//#endif
		//modified by AKF
		m_pdPeaks[a] = m_pSeqUtilFrag->m_dProton + dValue/dZ;
		m_pfIntens[a] = pfScore[tC]*pfScorePlus[m_pSeq[a+1]];

//		m_plSeq[lCount] = lValue;
//		m_pfSeq[lCount] = pfScore[tC]*pfScorePlus[m_pSeq[a+1]];
		lCount++;
		a++;
	}
	m_plSeq[lCount] = 0;
	return true;
}
/*
 * create list of non-zero predicted intensity values for x-ions and their
 * integer converted m/z values
 */
bool PTMTreeSearchScore::add_dX(const unsigned long _t,const long _c)
{
	long a = m_lSeqLength - 1;
/*
 * get the conversion factor between a straight sequence mass and an x-ion
 */
	double dValue = m_pSeqUtilFrag->m_dX;
/*
 * deal with non-hydrolytic cleavage
 */
	dValue += (m_pSeqUtilFrag->m_dCleaveC - m_pSeqUtilFrag->m_dCleaveCdefault);
	if(m_Term.m_lC)	{
		dValue += m_pSeqUtilFrag->m_pdAaMod[']'];
	}
	dValue += m_pSeqUtilFrag->m_pdAaFullMod[']'];
/*
 * deal with protein C-teminus
 */
	if(m_bIsC)	{
		dValue += m_pSeqUtilFrag->m_fCT;	
	}
//	unsigned long lValue = 0;
/*
 * calculate the conversion factor between an m/z value and its integer value
 * as referenced in m_vsmapMI
 */
	size_t tC = 0;
	unsigned long lCount = 0;
	float fSub = 0.0;
	float *pfScore = m_pSeqUtilFrag->m_pfXScore;
/*
 * from C- to N-terminus, calcuate fragment ion m/z values and store the results
 * look up appropriate scores from m_pSeqUtilFrag->m_pfAScore
 */
	m_dWE = m_dWidth/m_dErr;
	double dZ = (double)_c;
	const unsigned long tPos = (unsigned long) m_tSeqPos;
	while(a > 0)	{
		tC = m_pSeq[a];
//#ifdef PLUGGABLE_SCORING
//		dValue += m_pSeqUtilFrag->getAaMass((char)tC, tPos+a);
//		lValue = mconvert(dValue, _c);
//#else
		dValue +=m_pSeqUtilFrag-> m_pdAaMass[tC];
		dValue += m_pSeqUtilFrag->m_pdAaMod[tC];
		dValue += m_pSeqUtilFrag->m_pdAaFullMod[tC];
		if(m_pSeqUtilFrag->m_bPrompt)	{
			dValue += m_pSeqUtilFrag->m_pdAaPrompt[tC];
		}
		if (m_pSeqUtilFrag->m_bSequenceMods)	{
			SMap::iterator itSeq = m_pSeqUtilFrag->m_mapMods.find(tPos+a);
			if(itSeq != m_pSeqUtilFrag->m_mapMods.end())
				dValue += itSeq->second;
		}
//		lValue = mconvert(dValue, dZ);
//#endif
		//modified by AKF
		m_pdPeaks[a-1] = m_pSeqUtilFrag->m_dProton + dValue/dZ;
		m_pfIntens[a-1] = pfScore[tC];

//		m_plSeq[lCount] = lValue;
//		m_pfSeq[lCount] = pfScore[tC];
		lCount++;
		a--;
	}
/*
 * set the next integer mass value to 0: this marks the end of the array 
 */
//	m_plSeq[lCount] = 0;
	return true;
}
/*
 * create list of non-zero predicted intensity values for y-ions and their
 * integer converted m/z values
 */
bool PTMTreeSearchScore::add_dY(const unsigned long _t,const long _c)
{
	long a = m_lSeqLength - 1;
/*
 * get the conversion factor between a straight sequence mass and a y-ion
 */
	double dValue = m_pSeqUtilFrag->m_dY;
//	unsigned long lValue = 0;
/*
 * deal with non-hydrolytic cleavage
 */
	dValue += (m_pSeqUtilFrag->m_dCleaveC - m_pSeqUtilFrag->m_dCleaveCdefault);
	if(m_Term.m_lC)	{
		dValue += m_pSeqUtilFrag->m_pdAaMod[']'];
	}
	dValue += m_pSeqUtilFrag->m_pdAaFullMod[']'];

/*
/*
 * deal with protein C-teminus
 */
	if(m_bIsC)	{
		dValue +=  m_pSeqUtilFrag->m_fCT;	
	}
	unsigned long lCount = 0;
	float fSub = 0.0;
	float *pfScore = m_pSeqUtilFrag->m_pfYScore;
	float *pfScoreMinus = m_pSeqUtilFrag->m_pfBScore;
/*
 * from C- to N-terminus, calcuate fragment ion m/z values and store the results
 * look up appropriate scores from m_pSeqUtilFrag->m_pfAScore
 */
	long tPos = (unsigned long) m_tSeqPos;
	size_t tC = 0;
	m_dWE = m_dWidth/m_dErr;
	double dZ = (double)_c;
	while(a > 0)	{
		tC = m_pSeq[a];
//#ifdef PLUGGABLE_SCORING
//		dValue += m_pSeqUtilFrag->getAaMass((char)tC, tPos+a);
//		lValue = mconvert(dValue, _c);
//#else
		dValue +=m_pSeqUtilFrag-> m_pdAaMass[tC];
		dValue += m_pSeqUtilFrag->m_pdAaMod[tC];
		dValue += m_pSeqUtilFrag->m_pdAaFullMod[tC];
		if(m_pSeqUtilFrag->m_bPrompt)	{
			dValue += m_pSeqUtilFrag->m_pdAaPrompt[tC];
		}
		if (m_pSeqUtilFrag->m_bSequenceMods)	{
			SMap::iterator itSeq = m_pSeqUtilFrag->m_mapMods.find(tPos+a);
			if(itSeq != m_pSeqUtilFrag->m_mapMods.end())
				dValue += itSeq->second;
		}
//		lValue = mconvert(dValue, dZ);
//#endif
		if(_t == 0)	{
			if(a < 5)	{
				//modified by AKF
				m_pdPeaks[a-1] = m_pSeqUtilFrag->m_dProton + dValue/dZ;
				m_pfIntens[a-1] = pfScore[tC]*pfScoreMinus[m_pSeq[a-1]];

//				m_plSeq[lCount] = lValue;
//				m_pfSeq[lCount] = pfScore[tC]*pfScoreMinus[m_pSeq[a-1]];
				lCount++;
			}
		}
		else	{
			m_pdPeaks[a-1] = m_pSeqUtilFrag->m_dProton + dValue/dZ;
			m_pfIntens[a-1] = pfScore[tC]*pfScoreMinus[m_pSeq[a-1]];
			if(a == 2)	{
				if(m_pSeq[1] == 'P')	{
					m_pfIntens[a-1] *= 10;
				}
				else	{
					m_pfIntens[a-1] *= 3;
				}
			}

/*			m_plSeq[lCount] = lValue;
//			m_pfSeq[lCount] = pfScore[tC]*pfScoreMinus[m_pSeq[a-1]];
			if(a == 2)	{
				if(m_pSeq[1] == 'P')	{
					m_pfSeq[lCount] *= 10;
				}
				else	{
					m_pfSeq[lCount] *= 3;
				}
			}
*/			lCount++;
		}
		a--;
	}
/*
 * set the next integer mass value to 0: this marks the end of the array 
 */
//	m_plSeq[lCount] = 0;
	return true;
}

/*
 * create list of non-zero predicted intensity values for y-ions and their
 * integer converted m/z values
 */
bool PTMTreeSearchScore::add_dZ(const unsigned long _t,const long _c)
{
	long a = m_lSeqLength - 1;
/*
 * get the conversion factor between a straight sequence mass and a y-ion
 */
	double dValue = m_pSeqUtilFrag->m_dZ;
//	unsigned long lValue = 0;
/*
 * deal with non-hydrolytic cleavage
 */
	dValue += (m_pSeqUtilFrag->m_dCleaveC - m_pSeqUtilFrag->m_dCleaveCdefault);
	if(m_Term.m_lC)	{
		dValue += m_pSeqUtilFrag->m_pdAaMod[']'];
	}
	dValue += m_pSeqUtilFrag->m_pdAaFullMod[']'];
/*
/*
 * deal with protein C-teminus
 */
	if(m_bIsC)	{
		dValue +=  m_pSeqUtilFrag->m_fCT;	
	}
	size_t tC = 0;
	unsigned long lCount = 0;
	float fSub = 0.0;
	float *pfScore = m_pSeqUtilFrag->m_pfYScore;
	float *pfScoreMinus = m_pSeqUtilFrag->m_pfBScore;
/*
 * from C- to N-terminus, calcuate fragment ion m/z values and store the results
 * look up appropriate scores from m_pSeqUtilFrag->m_pfAScore
 */
	m_dWE = m_dWidth/m_dErr;
	double dZ = (double)_c;
	const unsigned long tPos = (unsigned long) m_tSeqPos;
	while(a > 0)	{
		tC = m_pSeq[a];
//#ifdef PLUGGABLE_SCORING
//		dValue += m_pSeqUtilFrag->getAaMass((char)tC, tPos+a);
//		lValue = mconvert(dValue, _c);
//#else
		dValue +=m_pSeqUtilFrag-> m_pdAaMass[tC];
		dValue += m_pSeqUtilFrag->m_pdAaMod[tC];
		dValue += m_pSeqUtilFrag->m_pdAaFullMod[tC];
		if(m_pSeqUtilFrag->m_bPrompt)	{
			dValue += m_pSeqUtilFrag->m_pdAaPrompt[tC];
		}
		if (m_pSeqUtilFrag->m_bSequenceMods)	{
			SMap::iterator itSeq = m_pSeqUtilFrag->m_mapMods.find(tPos+a);
			if(itSeq != m_pSeqUtilFrag->m_mapMods.end())
				dValue += itSeq->second;
		}
//		lValue = mconvert(dValue, dZ);
//#endif
		//modified by AKF
		m_pdPeaks[a-1] = m_pSeqUtilFrag->m_dProton + dValue/dZ;
		m_pfIntens[a-1] = pfScore[tC]*pfScoreMinus[m_pSeq[a-1]];
//		m_plSeq[lCount] = lValue;
//		m_pfSeq[lCount] = pfScore[tC]*pfScoreMinus[m_pSeq[a-1]];
		lCount++;
		a--;
	}
/*
 * set the next integer mass value to 0: this marks the end of the array 
 */
//	m_plSeq[lCount] = 0;
	return true;
}
