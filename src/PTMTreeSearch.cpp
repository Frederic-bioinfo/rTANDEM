/* PTMTreeSearch plugin
For furter details see or if PTMTreeSearch is used in a study please cite the following article:
PTMTreeSearch: a novel two-stage tree-search algorithm with pruning rules for the identification of post-translational modification of proteins in MS/MS spectra
Attila Kertesz-Farkas, Beata Reiz, Roberto Vera, Michael P Myers, Sandor Pongor
Bioinformatics 30 (2), 234-241, 2014 

FILE version: 2013-Febr-11  // Modified for rTANDEM on January 2014
*/

#include <cstring> 
#include "stdafx.h"
#include <sys/timeb.h>
#include <ctime>
#include "PTMTreeSearch.h"
#include "mprocess.h"
#include "mscore.h"
#include "PTMTreeSearchScore.h"
#include "msequenceserver.h"
#include "msequencecollection.h"

//#ifdef PARALELL_PTMTREESEARCH
//#define USE_MPI
#ifdef USE_MPI
#include "ownercompute.h"
#include "mpi.h"
#endif


bool reverseOrderFNC(float a,float b){
  return a>b;
}

PTMTreeSearch::PTMTreeSearch(){
  m_pProcess=NULL;
}

PTMTreeSearch::~PTMTreeSearch(){
}

bool PTMTreeSearch::set_mprocess(mprocess *_p){
  if(m_pProcess != NULL)
    delete m_pProcess;
  m_pProcess = _p;
  return true;
}

/*
 * process PTMTreeSearch
 */
bool PTMTreeSearch::refine()
{
  string strKey;
  string strValue;
  size_t tActiveNow = 0;
  size_t a = 0;
  size_t tPips = 0;
  
  //PTMTreeSearch parameters
  unsigned int uiPtmBound;
  double dMassLowerBound;
  double dMassUpperBound;
  bool bHitCheck;
  bool bKeepFixedMods;
  double dHitRatioTH;
  double dPTMEntropyTH;
  double dPTMMinEval;
  
  if(m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF)	{
    Rprintf("\tPTMTreeSearch\n");
    //cout << "\tPTMTreeSearch "<<endl;
  }
  
  //Load parameters
  strKey = "refine, tic percent";
  m_pProcess->m_xmlValues.get(strKey,strValue);
  double dTicPercent = atof(strValue.c_str());
  if(dTicPercent == 0)	{
    dTicPercent = 20.0;
  }
  size_t tTicMax = (size_t)((double)m_pProcess->m_vseqBest.size()*dTicPercent/100.0);
  if(tTicMax < 1)	{
    tTicMax = 1;
  }
  strKey = "scoring, maximum missed cleavage sites";
  m_pProcess->m_xmlValues.get(strKey,strValue);
  m_pProcess->m_tMissedCleaves = atoi(strValue.c_str());
  //	if(m_pProcess->m_tMissedCleaves < 2)	{
  //		m_pProcess->m_tMissedCleaves = 2;
  //	}
  //	m_pProcess->m_tMissedCleaves = 1;
  
  strKey = "protein, cleavage site";
  m_pProcess->m_xmlValues.get(strKey,strValue);
  m_pProcess->m_Cleave.load(strValue);
  strKey = "protein, cleavage semi";
  m_pProcess->m_xmlValues.get(strKey,strValue);
  if(strValue == "yes")	{
    m_pProcess->m_semiState.activate(true);
  }
  m_pProcess->m_semiState.activate(false);
  
  /*	m_pProcess->m_pScore->set_pam(false);
	strKey = "refine, cleavage semi";
	m_pProcess->m_xmlValues.get(strKey,strValue);
	if(strValue == "yes")	{
	m_pProcess->m_semiState.activate(true);
	}	else	{
	m_pProcess->m_semiState.activate(false);
	}
  */	m_pProcess->m_semiState.activate(false);
  
  //Load PTMTreeSearch Paramteres
  strKey = "refine, PTMTreeSearch 1r ptm bound";
  if (!m_pProcess->m_xmlValues.get(strKey,strValue)){
    strValue = "1";
  }
  uiPtmBound = atoi(strValue.c_str());
  
  strKey = "refine, PTMTreeSearch mass lower bound";
  if (!m_pProcess->m_xmlValues.get(strKey,strValue)){
    strValue = "50";
  }
  dMassLowerBound = atof(strValue.c_str());
  
  strKey = "refine, PTMTreeSearch mass upper bound";
  if (!m_pProcess->m_xmlValues.get(strKey,strValue)){
    strValue = "100";
  }
  dMassUpperBound = atof(strValue.c_str());
  
  bKeepFixedMods = false;
  strKey = "refine, PTMTreeSearch keep fixed modifications";
  if (!m_pProcess->m_xmlValues.get(strKey,strValue)){
    strValue = "no";
  }
  if (strValue == "yes"){
    bKeepFixedMods = true;
  }
  
  bHitCheck = false;
  strKey = "refine, PTMTreeSearch 1r NPM";
  if (!m_pProcess->m_xmlValues.get(strKey,strValue)){
    strValue = "yes";
  }
  if (strValue == "yes"){
    bHitCheck = true;
  }
  dHitRatioTH = 0.20;  //Score ratio SR
  strKey = "refine, PTMTreeSearch 1r SR";
  if (m_pProcess->m_xmlValues.get(strKey,strValue)){
    dHitRatioTH = atof(strValue.c_str());
  }
  
  dPTMEntropyTH = 0.8;
  strKey = "refine, PTMTreeSearch 1r MPE";
  if (m_pProcess->m_xmlValues.get(strKey,strValue)){
    dPTMEntropyTH = atof(strValue.c_str());
  }
  
  dPTMMinEval = -3.0;
  strKey = "refine, PTMTreeSearch 1r min eval";
  if (m_pProcess->m_xmlValues.get(strKey,strValue)){
    dPTMMinEval = atof(strValue.c_str());
  }
  
  //load the PTM DB
  Modification mod;
  mod.m_strName = "no modif";
  for (int i = 0; i <255; i++){
    m_vmodModificationsDB[i].push_back(mod);
  }
  bool inputDB = false;
  
  //Load modifications from the specified XML file from OMSSA project
  strKey = "refine, PTMTreeSearch modif file";
  if (m_pProcess->m_xmlValues.get(strKey,strValue)){
    SAXModsHandler mods_input(strValue,&m_vmodModificationsDB[0]);
    if(!mods_input.load()){
      if (m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF){
	///cout << "\t\tThe input parameter file \"" << strValue << "\" could not be located.\nCheck the file path name and try again.\nSkipping PTMTreeSearch.\n";
	Rprintf("\t\tThe input parameter file \"%s", strValue.c_str());
	Rprintf("\" could not be located.\nCheck the file path name and try again.\nSkipping PTMTreeSearch.\n");
      }
    } else {
      inputDB = true;
    }
  }
  
  //Load modifications from RESID
  strKey = "refine, PTMTreeSearch resid modif file";
  if (m_pProcess->m_xmlValues.get(strKey,strValue)){
    SAXResidHandler residmods_input(strValue,&m_vmodModificationsDB[0]);
    if(!residmods_input.load()){
      if (m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF){
	//cout << "\t\tThe input parameter file \"" << strValue << "\" could not be located.\nCheck the file path name and try again.\n";
	Rprintf("\t\tThe input parameter file \"%s",strValue.c_str());
	Rprintf("\" could not be located.\nCheck the file path name and try again.\n");
      }
    } else {
      inputDB = true;
    }
  }
  //Load modifications from unimod from MASCOT xml file format.
  strKey = "refine, PTMTreeSearch unimod modif file";
  if (m_pProcess->m_xmlValues.get(strKey,strValue)){
    SAXUnimodHandler unimods_input(strValue,&m_vmodModificationsDB[0]);
    if(!unimods_input.load()){
      if (m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF){
	
	//cout << "\t\tThe input parameter file \"" << strValue << "\" could not be located.\nCheck the file path name and try again.\n");
	Rprintf("\t\tThe input parameter file \"%s", strValue.c_str());
	Rprintf("\" could not be located.\nCheck the file path name and try again.\n");
      }
    } else {
      inputDB = true;
    }
  }
  
  //Load modifications from uniprot knowledge
  strKey = "refine, PTMTreeSearch uniprot modif file";
  if (m_pProcess->m_xmlValues.get(strKey,strValue)){
    if (!openUniprotPTMs(strValue.c_str(), &m_vmodModificationsDB[0])){
      if (m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF){
	//cout << "\t\tFailure with uniprot file.\n";
	Rprintf("\t\tFailure with uniprot file.\n");
      }
    } else {
      inputDB = true;
    }
  }
  if (!inputDB && (m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF)){
    //cout << "\t\tNo modification DB has been loaded successfully.\nSkipping PTMTreeSearch.\n";
    Rprintf("\t\tNo modification DB has been loaded successfully.\nSkipping PTMTreeSearch.\n");
    return true;
  }
  
  if (bKeepFixedMods == true){
    for  (int i = 0; i < 255; i++){
      if (m_pProcess->m_pScore->m_seqUtil.m_pdAaFullMod[i] != 0.0){
	m_vmodModificationsDB[i].clear();
	mod.m_strName = "fixed modif";
	mod.m_dAveMass = m_pProcess->m_pScore->m_seqUtil.m_pdAaFullMod[i];
	mod.m_dMonoMass = m_pProcess->m_pScore->m_seqUtilAvg.m_pdAaFullMod[i];
	mod.m_strResidues = (char)i;
	m_vmodModificationsDB[i].push_back(mod);
      }
    }
  }
  
  //remove the double modifications from the PTM DB
  unsigned int modNum[255];
  double* mods[255];
  unsigned int modCnt = 0;
  for (unsigned int i = 0; i < 255; i++){
    for (unsigned int j = 0; j < m_vmodModificationsDB[i].size(); j++){
      m_vmodModificationsDB[i][j].m_bType = m_pProcess->m_pScore->m_pSeqUtilFrag->m_calc.getMassType() == masscalc::monoisotopic;
      for (unsigned int k = 0; k < j; k++){
	if (floor(m_vmodModificationsDB[i][k].m_dMonoMass*10000) == floor(m_vmodModificationsDB[i][j].m_dMonoMass*10000)){
	  m_vmodModificationsDB[i].erase(m_vmodModificationsDB[i].begin()+j);
	  j--;
	  break;
	}
      }
    }
    modNum[i] = m_vmodModificationsDB[i].size()-1;
    mods[i] = new double[modNum[i]+1];
    modCnt += modNum[i];
    //		vector<Modification>::iterator modIt;
    sort(m_vmodModificationsDB[i].begin() + 1,  m_vmodModificationsDB[i].end());
    for (unsigned int j = 0; j < m_vmodModificationsDB[i].size(); j++){		
      if (m_pProcess->m_pScore->m_pSeqUtilFrag->m_calc.getMassType() == masscalc::monoisotopic){
	mods[i][j] = m_vmodModificationsDB[i][j].m_dMonoMass;
      } else {
	mods[i][j] = m_vmodModificationsDB[i][j].m_dAveMass;
      }
    }		
  }
  
  if (m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF){
    //cout << "\t\tNumber of unique modifications that have been loaded is:" << modCnt << endl;
    //cout << "\t\t";
    Rprintf("\t\tNumber of unique modifications that have been loaded is: %lu \n", modCnt);
    Rprintf("\t\t");
  }
  
  //Turn off other modifications in X!tandem
  strValue = "";
  m_pProcess->m_pScore->m_seqUtil.modify_all(strValue);
  m_pProcess->m_pScore->m_seqUtilAvg.modify_all(strValue);
  m_pProcess->m_pScore->m_seqUtil.modify_maybe(strValue);
  m_pProcess->m_pScore->m_seqUtilAvg.modify_maybe(strValue);
  m_pProcess->m_pScore->m_seqUtil.modify_motif(strValue);
  m_pProcess->m_pScore->m_seqUtilAvg.modify_motif(strValue);
  m_pProcess->m_strLastMods.clear();
  
  vector <string> vstrModifications = m_pProcess->m_vstrModifications;
  vector <string> vstrMods = m_pProcess->m_vstrMods;
  m_pProcess->m_vstrModifications.clear();
  m_pProcess->m_vstrMods.clear();
  
  spectrumInt spectraIntTMP;
  int i;
  for (a = 0; a < m_pProcess->m_vSpectra.size(); ++a){
    for (i = 0; i < m_pProcess->m_vSpectra[a].m_vMI.size(); ++i){
      spectraIntTMP.push_back(m_pProcess->m_vSpectra[a].m_vMI[i].m_fI);
    }
    sort(spectraIntTMP.begin(),spectraIntTMP.end(),reverseOrderFNC);
    for (i = 1; i < spectraIntTMP.size(); ++i){
      spectraIntTMP[i] += spectraIntTMP[i-1];
    }
    ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->spectraSumInt.push_back(spectraIntTMP[i-1]);
    ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->m_vfSpectraInt.push_back(spectraIntTMP);
    spectraIntTMP.clear();
  }
  /*
   * score the spectra against each sequence
   */
  ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->set_parent_tolerance(dMassLowerBound,dMassUpperBound);
  ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->sort_details();
  ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->m_bRefine = true;
  ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->m_uiPtmBound = uiPtmBound;
  ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->m_errValues = m_pProcess->m_errValues;
  ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->m_vSpectra = &m_pProcess->m_vSpectra;
  ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->m_pdMods = mods;
  ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->m_uipModNum = modNum;
  ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->m_bHitCheck = bHitCheck;
  ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->m_dHitRatioTh = dHitRatioTH;
  ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->m_dPTMEntropyTh = dPTMEntropyTH;
  ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->m_lRound = 0;
  m_pProcess->m_bCrcCheck = true;
  /*	if (m_pProcess->m_Cleave.m_lType & 0x02){ //tryptic digestion
	((PTMTreeSearchScore*)(m_pProcess->m_pScore))->m_bTrypticDigestion = true;
	} else {
	((PTMTreeSearchScore*)(m_pProcess->m_pScore))->m_bTrypticDigestion = false;
	}	
  */
  vector<mspectrum> vSpectra;
  vSpectra = m_pProcess->m_vSpectra;
  
  for (a = 0;  a < m_pProcess->m_vseqBest.size(); a++){
    m_pProcess->score(m_pProcess->m_vseqBest[a]);
    
    tPips++;
    if(tPips == tTicMax)	{
      if(m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF)	{
	Rprintf("."); //cout << ".";
	m_pProcess->m_prcLog.log(".");
      }
      tPips = 0;
    }
  }
  
  if(m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF)	{
    Rprintf("done.\n"); //cout << "done.\n";
  }
  /*
    strKey = "refine, PTMTreeSearch 2r";
    if (!m_pProcess->m_xmlValues.get(strKey,strValue)){
    strValue = "yes";
    }
    if (strValue != "yes"){
    for (int i = 0; i < 255; i++){
    delete [] mods[i];
    }
    m_pProcess->m_vstrModifications = vstrModifications;
    m_pProcess->m_vstrMods = vstrMods;
    ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->set_parent_tolerance(-dMassLowerBound,-dMassUpperBound);
    ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->sort_details();		
    return true;
    }
  */
  //2nd round of PTMTreeSearch, 	
  if(m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF)	{
    //cout << "\tPTMTreeSearch 2nd round" << endl;
    Rprintf("\tPTMTreeSearch 2nd round\n");
  }
  //calculate evalues.
  strKey = "refine, maximum valid expectation value";
  m_pProcess->m_xmlValues.get(strKey,strValue);
  double dMaxExpect = 0.1;
  if(strValue.size() > 0)	{
    double dMaxExpect = atof(strValue.c_str());
  }
  
  vector<mspectrum>::iterator itS = m_pProcess->m_vSpectra.begin();
  vector<maa>::iterator itAa;
  vector<maa>::iterator itAaEnd;
  int iModSize;
  double expect;
  
  while (itS != m_pProcess->m_vSpectra.end())	{
    itS->m_hHyper.model();
    itS->m_hHyper.set_protein_factor(1.0);
    if(!itS->m_vseqBest.empty() && !itS->m_vseqBest[0].m_vDomains.empty() && !itS->m_vseqBest[0].m_vDomains[0].m_vAa.empty())	{
      
      expect = log10(itS->m_hHyper.expect(m_pProcess->m_pScore->hconvert(itS->m_vseqBest[0].m_vDomains[0].m_fHyper)));
      itAa =		itS->m_vseqBest[0].m_vDomains[0].m_vAa.begin();
      itAaEnd =	itS->m_vseqBest[0].m_vDomains[0].m_vAa.end();
      if (expect < 3){
	for (; itAa != itAaEnd; itAa++){
	  iModSize = floor(itAa->m_dMod*10000);
	  for (unsigned int k = 1; k < m_vmodModificationsDB[itAa->m_cRes].size(); k++){
	    if (floor((mods[itAa->m_cRes][k])*10000) == iModSize){
	      if ((float)itAa->m_fPval>0.00001){
		m_vmodModificationsDB[itAa->m_cRes][k].m_uiCount++;
		if (m_vmodModificationsDB[itAa->m_cRes][k].m_dTarget < itAa->m_fPval-1.0)
		  m_vmodModificationsDB[itAa->m_cRes][k].m_dTarget = itAa->m_fPval-1.0;
		if (m_vmodModificationsDB[itAa->m_cRes][k].m_dTargetEval > expect)
		  m_vmodModificationsDB[itAa->m_cRes][k].m_dTargetEval = expect;
	      }
	      break;
	    }
	  }
	}
      }
    }
    itS++;
  }
#ifdef USE_MPI
  
  //merge modifications from different processes.
  MPI_Status	status;
  double	 	*sendbuffer, *recvbuffer;
  int package_size = 7;
  sendbuffer  = new double[modCnt*package_size];
  recvbuffer	= new double[modCnt*package_size];
  int			buffsize;
  int			ret;
  
  //slave processes send the ptm counts to the master process.
  if (MPI::COMM_WORLD.Get_rank() == 0){
    for (int j = 1; j < MPI::COMM_WORLD.Get_size(); j++){
      ret = MPI_Recv(recvbuffer, modCnt*package_size, MPI_DOUBLE, j, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      if (ret != MPI_SUCCESS){	return false;	}
      for (int k = 0; k < modCnt*package_size; k+=package_size){
	m_vmodModificationsDB[(int)floor(recvbuffer[k])][(int)floor(recvbuffer[k+1])].m_uiCount += recvbuffer[k+2];
	if (m_vmodModificationsDB[(int)floor(recvbuffer[k])][(int)floor(recvbuffer[k+1])].m_dTarget < recvbuffer[k+3])
	  m_vmodModificationsDB[(int)floor(recvbuffer[k])][(int)floor(recvbuffer[k+1])].m_dTarget = recvbuffer[k+3];
	if (m_vmodModificationsDB[(int)floor(recvbuffer[k])][(int)floor(recvbuffer[k+1])].m_dTargetEval > recvbuffer[k+4])
	  m_vmodModificationsDB[(int)floor(recvbuffer[k])][(int)floor(recvbuffer[k+1])].m_dTargetEval = recvbuffer[k+4];				
      }
    }
    vector <double> vdTargetEval;
    vdTargetEval.resize(modCnt);
    for (int i = 0; i < 255; i++){
      for (unsigned int j = 0; j < m_vmodModificationsDB[i].size(); j++){
	if ( m_vmodModificationsDB[i][j].m_uiCount > 0 ){
	  vdTargetEval.push_back(m_vmodModificationsDB[i][j].m_dTargetEval);
	}
      }
    }
    double dPTMth = dPTMMinEval;
    
    for (int i = 0; i < 255; i++){
      for (unsigned int j = 0; j < m_vmodModificationsDB[i].size(); j++){
	if ( m_vmodModificationsDB[i][j].m_dTargetEval < dPTMth ){
	  m_vmodModificationsDB[i][j].m_uiCount = 2;
	} else {
	  m_vmodModificationsDB[i][j].m_uiCount = 0;
	}
      }
    }
  } else {
    unsigned int cnt = 0;
    for (int j = 0; j < 255; j++){
      for (int k = 1; k < m_vmodModificationsDB[j].size(); k++){
	sendbuffer[cnt++] = j;
	sendbuffer[cnt++] = k;
	sendbuffer[cnt++] = m_vmodModificationsDB[j][k].m_uiCount;
	sendbuffer[cnt++] = m_vmodModificationsDB[j][k].m_dTarget;
	sendbuffer[cnt++] = m_vmodModificationsDB[j][k].m_dTargetEval;
	cnt++;
	cnt++;
      }
    }
    ret = MPI_Send(sendbuffer, modCnt*package_size, MPI_DOUBLE, 0,0,MPI_COMM_WORLD);
    if (ret != MPI_SUCCESS){	return false;	}
  }	
  
  ret = MPI_Barrier(MPI_COMM_WORLD);
  if (ret != MPI_SUCCESS){	return false;	}
  
  //master process sends the ptmdetails of the modifications to the slave processes.
  if (MPI::COMM_WORLD.Get_rank() == 0){
    unsigned cnt = 0;
    for (int j = 0; j < 255; j++){
      for (int k = 1; k < m_vmodModificationsDB[j].size(); k++){
	sendbuffer[cnt++] = j;
	sendbuffer[cnt++] = k;
	sendbuffer[cnt++] = m_vmodModificationsDB[j][k].m_uiCount;
	cnt++;
	cnt++;
	cnt++;
	cnt++;
      }
    }
    for (int j = 1; j < MPI::COMM_WORLD.Get_size(); j++){
      ret = MPI_Send(sendbuffer, modCnt*package_size, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
      if (ret != MPI_SUCCESS){	return false;	}
    }
  } else {
    ret = MPI_Recv(recvbuffer, modCnt*package_size, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    
    for (int k = 0; k < modCnt*package_size; k+=package_size){
      m_vmodModificationsDB[(int)floor(recvbuffer[k])][(int)floor(recvbuffer[k+1])].m_uiCount = recvbuffer[k+2];
      
    }
  }	
  
  delete [] sendbuffer;
  delete [] recvbuffer;
#endif
  
  for (unsigned int i = 0; i < 255; i++){
    for (unsigned int k = 1; k < m_vmodModificationsDB[i].size(); k++){
      if (m_vmodModificationsDB[i][k].m_uiCount < 1){
	m_vmodModificationsDB[i].erase(m_vmodModificationsDB[i].begin()+k);
	k--;
      }
    }
  }
  
  //all processes have the same modifications with the same frequency.
  //filter the nonsignificant modifications.
  for (int i = 0; i < 255; i++){
    delete [] mods[i];
  }
  
  modCnt = 0;
  for (int i = 0; i < 255; i++){
    modNum[i] = m_vmodModificationsDB[i].size()-1;
    mods[i] = new double[modNum[i]+1];
    modCnt += modNum[i];
    for (unsigned int j = 0; j < m_vmodModificationsDB[i].size(); j++){
      if (m_pProcess->m_pScore->m_pSeqUtilFrag->m_calc.getMassType() == masscalc::monoisotopic){
	mods[i][j] = m_vmodModificationsDB[i][j].m_dMonoMass;
      } else {
	mods[i][j] = m_vmodModificationsDB[i][j].m_dAveMass;
      }
    }
  }
  
  if(m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF)	{
    //    cout << "\t\tNumber of unique modifications have been found is:" << modCnt << endl << "\t\t";
    Rprintf("\t\tNumber of unique modifications have been found is: %lu \n", modCnt);
    Rprintf("\t\t");
  }
  
  strKey = "refine, PTMTreeSearch 2r";
  if (!m_pProcess->m_xmlValues.get(strKey,strValue)){
    strValue = "yes";
  }
  if (strValue != "yes"){
    for (int i = 0; i < 255; i++){
      delete [] mods[i];
    }
    m_pProcess->m_vstrModifications = vstrModifications;
    m_pProcess->m_vstrMods = vstrMods;
    ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->set_parent_tolerance(-dMassLowerBound,-dMassUpperBound);
    ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->sort_details();		
    return true;
  }
  
  //2nd round of PTMTreeSearch.
  m_pProcess->m_vSpectra = vSpectra;
  ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->m_lRound = 1;
  
  strKey = "refine, PTMTreeSearch 2r ptm bound";
  if (m_pProcess->m_xmlValues.get(strKey,strValue)){
    uiPtmBound = atoi(strValue.c_str());
  }
  if(m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF)	{
    //		cout << "\t\tNumber of unique modifications that have been loaded is:" << modCnt << endl << "\t\t"; /// commented out before rTANDEM
    //cout << "\t\tPTMBound = " << uiPtmBound << endl << "\t\t";
    Rprintf("\t\tPTMBound = %lu \n",(unsigned long)uiPtmBound);
    Rprintf("\t\t");
  }
  
  bHitCheck = false;
  int noptmhitlimit = 0;
  strKey = "refine, PTMTreeSearch 2r no NPM limit";
  if (m_pProcess->m_xmlValues.get(strKey,strValue)){
    noptmhitlimit = atoi(strValue.c_str());
  }
  if (modCnt > noptmhitlimit){
    bHitCheck = true;
  }	
  
  dHitRatioTH = 0.05;
  strKey = "refine, PTMTreeSearch 2r SR";
  if (m_pProcess->m_xmlValues.get(strKey,strValue)){
    dHitRatioTH = atof(strValue.c_str());
  }
  
  dPTMEntropyTH = 0.0;
  strKey = "refine, PTMTreeSearch 2r MPE";
  if (m_pProcess->m_xmlValues.get(strKey,strValue)){
    dPTMEntropyTH = atof(strValue.c_str());
  }
  ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->m_uiPtmBound = uiPtmBound;
  ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->m_bHitCheck = bHitCheck;
  ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->m_dHitRatioTh = dHitRatioTH;
  ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->m_dPTMEntropyTh = dPTMEntropyTH;
  m_pProcess->m_pScore->set_pam(false);
  strKey = "refine, cleavage semi";
  m_pProcess->m_xmlValues.get(strKey,strValue);
  if(strValue == "yes")	{
    m_pProcess->m_semiState.activate(true);
  }	else	{
    m_pProcess->m_semiState.activate(false);
  }
  m_pProcess->m_semiState.activate(false);
  
  
  for (a = 0;  a < m_pProcess->m_vseqBest.size(); a++){
    
    m_pProcess->score(m_pProcess->m_vseqBest[a]);
    tPips++;
    if(tPips == tTicMax)	{
      if(m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF)	{
	Rprintf("."); //cout << ".";
	m_pProcess->m_prcLog.log(".");
      }
      tPips = 0;
    }
  }
  
  //clear modifications
  for (int i = 0; i < 255; i++){
    delete [] mods[i];
  }
  m_pProcess->m_vstrModifications = vstrModifications;
  m_pProcess->m_vstrMods = vstrMods;
  ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->set_parent_tolerance(-dMassLowerBound,-dMassUpperBound);
  ((PTMTreeSearchScore*)(m_pProcess->m_pScore))->sort_details();	
  
  if(m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF)	{
    Rprintf(" done.\n"); //cout << " done.\n";
  }

  //Send the results via email
  /*
  char emailcommand[2048];
  string strEmailAddress;
  strKey = "output, path";
  if(m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF)	{
    if (m_pProcess->m_xmlValues.get(strKey,strValue)){
      strKey = "refine, PTMTreeSearch emailaddress";

      if (m_pProcess->m_xmlValues.get(strKey,strEmailAddress)){
	
	sprintf(emailcommand, "echo \"PTMTreeSearch results can be seen at:\n\n http://net.icgeb.org/thegpm-cgi/plist.pl?npep=0&path=/thegpm-cgi/%s&proex=-1&ltype= \" | mailx -s \"Peptide identification results with PTMTreeSearch\" -r noreply@icgeb.org -S smtp=smtp://smtp.icgeb.trieste.it:25 %s\n", strValue.c_str(),strEmailAddress.c_str() );
	system(emailcommand);
      }
    }
  }
  */
  return true;
}


bool PTMTreeSearch::openUniprotPTMs(const char* filename, vector<Modification>* _mods){
  
  FILE *f = fopen(filename, "r");
  if (f == NULL) {
    return false;
  }
  char row[2048]; 
  
  unsigned int m_uiModType;	//type of the modification
  double m_dMonoMass;			//monoisotopic mass
  double m_dAveMass;			//average mas
  string m_strResidues;		//the residues that can be modified
  string m_strUnimod;			//unimod id
  string m_strName;			//name of the modification
  unsigned int m_uiId;		//unique id of the modification
  
  Modification mod;
  char resid;
  unsigned int len;
  m_uiModType = 0;
  
  while (true){
    if (fgets(row, 2048, f) == NULL) break;
    if (row[0] == '-') continue;
    if (row[0] == ' ') continue;
    row[2] = '\0';
    if (strcmp(row, "AC")==0){
      m_strName = row+5;
      m_strName.append(";");
    }
    if (strcmp(row, "DR")==0){
      m_strName.append(row+5);
      m_strName.append(";");
    }
    if (strcmp(row, "MM")==0){
      m_dMonoMass = atof(row+5);
    }
    if (strcmp(row, "MA")==0){
      m_dAveMass = atof(row+5);
    }
    if (strcmp(row, "PA")==0){
      m_uiModType++;
      if (strcmp(row+5, "Amino acid side chain.\n")==0){
	m_uiModType += 10;
      }
    }
    if (strcmp(row, "PP")==0){
      m_uiModType++;
      if (strcmp(row+5, "Anywhere.\n")==0){
	m_uiModType += 10;
      }
    }
    if (strcmp(row, "TG")==0){
      char* tok = strtok(row+5,"-.");
      while (tok != NULL){
	if (strcmp(tok, "Isoleucine") == 0){	m_strResidues.append("I");	} else
	  if (strcmp(tok, "Leucine") == 0){		m_strResidues.append("L");	} else 
	    if (strcmp(tok, "Valine") == 0){		m_strResidues.append("V");	} else
	      if (strcmp(tok, "Phenylalanine") == 0){	m_strResidues.append("F");	} else
		if (strcmp(tok, "Methionine") == 0){	m_strResidues.append("M");	} else
		  if (strcmp(tok, "Cysteine") == 0){		m_strResidues.append("C");	} else
		    if (strcmp(tok, "Alanine") == 0){		m_strResidues.append("A");	} else
		      if (strcmp(tok, "Glycine") == 0){		m_strResidues.append("G");	} else
			if (strcmp(tok, "Proline") == 0){		m_strResidues.append("P");	} else
			  if (strcmp(tok, "Threonine") == 0){		m_strResidues.append("T");	} else
			    if (strcmp(tok, "Serine") == 0){		m_strResidues.append("S");	} else
			      if (strcmp(tok, "Tyrosine") == 0){		m_strResidues.append("Y");	} else
				if (strcmp(tok, "Tryptophan") == 0){	m_strResidues.append("W");	} else
				  if (strcmp(tok, "Glutamine") == 0){		m_strResidues.append("Q");	} else
				    if (strcmp(tok, "Asparagine") == 0){	m_strResidues.append("N");	} else
				      if (strcmp(tok, "Histidine") == 0){		m_strResidues.append("H");	} else
					if (strcmp(tok, "Glutamate") == 0){		m_strResidues.append("E");	} else
					  if (strcmp(tok, "Aspartate") == 0){		m_strResidues.append("D");	} else
					    if (strcmp(tok, "Lysine") == 0){		m_strResidues.append("K");	} else
					      if (strcmp(tok, "Arginine") == 0){		m_strResidues.append("R");	}
	tok = strtok(NULL, "-.");
      }
    }
    if (strcmp(row, "//")==0){
      if ((m_uiModType > 20 || m_uiModType == 0) && m_dAveMass > 0 && m_dMonoMass > 0){
	len = m_strResidues.length();
	for (size_t i = 0; i < len; i++){
	  mod.m_dAveMass = m_dAveMass;
	  mod.m_dMonoMass = m_dMonoMass;
	  mod.m_strName = m_strName;
	  mod.m_strUnimod = m_strUnimod;
	  mod.m_uiModType = 0;
	  resid = m_strResidues[i];
	  mod.m_strResidues = resid;
	  m_vmodModificationsDB[resid].push_back(mod);
	}
      }
      //restore variables to default;
      m_uiId			= 0;
      m_uiModType		= 0;
      m_strName		= "";
      m_dMonoMass		= 0.0;
      m_dAveMass		= 0.0;
      m_strResidues	= "";
      m_strUnimod		= "";		 
    }
  }
  return true;
}

/*
 * PTMTreeSearchmanager contains static short-cuts for dealing with PTMTreeSearch
 * plug-ins.
 */
const char* PTMTreeSearchmanager::TYPE = "refinement,PTMTreeSearch algorithm";

/*
 * create_PTMTreeSearch creates the correct PTMTreeSearch object for the given set of XML
 * parameters.
 */
PTMTreeSearch* PTMTreeSearchmanager::create_PTMTreeSearch(XmlParameter &_x)
{
  string strValue;
  string strKey = TYPE;
  if (!_x.get(strKey,strValue))
    strValue = "tandem";
  return (PTMTreeSearch*) mpluginmanager::get().create_plugin(TYPE, strValue.data());
}

/*
 * register_factory registers a factory with the mpluginmanager for
 * creating mpmods derived objects.
 */
void PTMTreeSearchmanager::register_factory(const char* _spec, mpluginfactory* _f)
{
  mpluginmanager::get().register_factory(TYPE, _spec, _f);
}

// Factory instance, registers itself with the PTMTreeSearchmanager.
static PTMTreeSearchfactory_tandem factory;

PTMTreeSearchfactory_tandem::PTMTreeSearchfactory_tandem()
{
  PTMTreeSearchmanager::register_factory("tandem", this);
}

mplugin* PTMTreeSearchfactory_tandem::create_plugin()
{
  return new PTMTreeSearch();
}


SAXModsHandler::SAXModsHandler(const string& _p, vector<Modification> *_mods)
{
  m_strXmlPath	= _p;
  m_vmodModificationsDB = _mods;
  m_uiInputType	= 0;
  m_uiId			= 0;
  m_uiModType		= 0;
  m_dMonoMass		= 0.0;
  m_dAveMass		= 0.0;
  m_strName		= "";
  m_strResidues	= "";
  m_strUnimod		= "";		 
}

SAXModsHandler::~SAXModsHandler()
{
}

void SAXModsHandler::startElement(const XML_Char *el, const XML_Char **attr)
{
  if(isElement("MSMod", el)){
    m_uiInputType = 1; //id node.
    return;
  }
  if(isElement("MSModType", el)){
    m_uiInputType = 2; //modification type
    return;
  }
  if(isElement("MSModSpec_name", el)){
    m_uiInputType = 3; //modification name
    return;
  }
  if(isElement("MSModSpec_monomass", el)){
    m_uiInputType = 4; //modification monomass
    return;
  }
  if(isElement("MSModSpec_averagemass", el)){
    m_uiInputType = 5; //modification average mass
    return;
  }
  if(isElement("MSModSpec_residues_E", el)){
    m_uiInputType = 6; //modification residues
    return;
  }
  if(isElement("MSModSpec_unimod", el)){
    m_uiInputType = 7; //modification unimod id
    return;
  }
}

void SAXModsHandler::endElement(const XML_Char *el)
{
  size_t len;
  Modification mod;
  char resid;
  if(isElement("MSModSpec", el)){
    if (m_uiModType == 0){
      len = m_strResidues.length();
      for (size_t i = 0; i < len; i++){
	mod.m_dAveMass = m_dAveMass;
	mod.m_dMonoMass = m_dMonoMass;
	mod.m_strName = m_strName;
	mod.m_strUnimod = m_strUnimod;
	mod.m_uiModType = m_uiModType;
	resid = m_strResidues[i];
	mod.m_strResidues = resid;
	m_vmodModificationsDB[resid].push_back(mod);
      }
    }
    
    //restore variables to default;
    m_uiId			= 0;
    m_uiModType		= 0;
    m_strName		= "";
    m_dMonoMass		= 0.0;
    m_dAveMass		= 0.0;
    m_strResidues	= "";
    m_strUnimod		= "";		 
  }
  if (m_uiInputType > 0)
    m_uiInputType = 0;
  
}
// changed in 2006.09.15 to improve handling of characters
// suggested by Brendan Maclean
void SAXModsHandler::characters(const XML_Char *s, int len)
{
  if (m_uiInputType == 0)
    return;
  
  if (len > 1024) 
    len = 1024;
  
  XML_Char ss[1024];
  memcpy(ss,s, len);
  ss[len] = '\0';
  
  if (m_uiInputType == 1){
    m_uiId = atoi(s);
  } else if (m_uiInputType == 2){
    m_uiModType = atoi(ss);
  }else if (m_uiInputType  == 3){
    m_strName = ss;
  } else if (m_uiInputType == 4){
    m_dMonoMass = atof(ss);
  } else if (m_uiInputType == 5){
    m_dAveMass = atof(ss);
  } else if (m_uiInputType == 6){
    m_strResidues.append(s,len);
  } else if (m_uiInputType == 7){
    m_strUnimod = ss;
  }
}

bool SAXModsHandler::load(){
  ifstream ifXml;
  /*
   * open an input stream, return false if it can't be opened
   */
  ifXml.open(m_strXmlPath.c_str());
  if(ifXml.fail()){
    //cout << "\nFailed to open: \"" << m_strXmlPath.c_str() << "\"\n"; // commente out before rTANDEM
    return false;
  }
  setFileName( m_strXmlPath.data() );
  
  parse();
  return true;
}


SAXResidHandler::SAXResidHandler(const string& _p, vector<Modification> *_mods)
{
  m_strXmlPath	= _p;
  m_vmodModificationsDB = _mods;
  m_uiInputType	= 0;
  m_uiId			= 0;
  m_strModType	= "";
  m_dMonoMass		= 0.0;
  m_dAveMass		= 0.0;
  m_strName		= "";
  m_strResidues	= "";
  m_strUnimod		= "";		 
}

SAXResidHandler::~SAXResidHandler()
{
}

void SAXResidHandler::startElement(const XML_Char *el, const XML_Char **attr)
{
  if(isElement("Entry", el)){
    m_uiInputType = 1; //id node.
    m_uiId = atoi(attr[1]+2);
    return;
  }
  if(isElement("Code", el)){
    m_uiInputType = 2; //modification name
    return;
  }
  if(isElement("SequenceSpec", el)){
    m_uiInputType = 3; //modification residues
    return;
  }
  if(isElement("Condition", el)){
    m_uiInputType = 4; //modification unimod id
    return;
  }
  if(isElement("CorrectionBlock", el)){
    m_uiInputType = 5; //modification type
    return;
  }
  if (m_uiInputType == 5)
    if(isElement("Weight", el)){
      if (isAttr("type", attr[0]) && isAttr("chemical", attr[1]) ){
	m_uiInputType = 6; //modification monomass
	return;
      }
      if (isAttr("type", attr[0]) && isAttr("physical", attr[1]) ){
	m_uiInputType = 7; //modification monomass
	return;
      }
    }
  
}

void SAXResidHandler::endElement(const XML_Char *el)
{
  size_t pos = 0;
  size_t len;
  Modification mod;
  char resid;
  if(isElement("Entry", el)){
    if (m_strModType.empty() == true){
      while ((pos = m_strResidues.find(',')>m_strResidues.npos)){
	m_strResidues.erase(pos,1);
      }
      while ((pos = m_strResidues.find(' ')>m_strResidues.npos)){
	m_strResidues.erase(pos,1);
      }
      len = m_strResidues.length();
      for (size_t i = 0; i < len; i++){
	mod.m_dAveMass = m_dAveMass;
	mod.m_dMonoMass = m_dMonoMass;
	mod.m_strName = m_strName;
	mod.m_strUnimod = m_strUnimod;
	mod.m_uiModType = 0;
	resid = m_strResidues[i];
	mod.m_strResidues = resid;
	m_vmodModificationsDB[resid].push_back(mod);
      }
    }
    
    //restore variables to default;
    m_uiId			= 0;
    m_strModType	= "";
    m_strName		= "";
    m_dMonoMass		= 0.0;
    m_dAveMass		= 0.0;
    m_strResidues	= "";
    m_strUnimod		= "";		 
  }
  if(!isElement("CorrectionBlock", el) && m_uiInputType>=5){
    return;
  }
  if (m_uiInputType > 0)
    m_uiInputType = 0;
  
}
// changed in 2006.09.15 to improve handling of characters
// suggested by Brendan Maclean
void SAXResidHandler::characters(const XML_Char *s, int len)
{
  if (m_uiInputType == 0)
    return;
  
  if (len > 1024) 
    len = 1024;
  
  XML_Char ss[1024];
  memcpy(ss,s, len);
  ss[len] = '\0';
  
  if (m_uiInputType == 1){
  }else if (m_uiInputType  == 2){
    m_strName = ss;
  } else if (m_uiInputType == 3){
    m_strResidues.append(s,len);
  } else if (m_uiInputType == 4){
    m_strModType = ss;
  } else if (m_uiInputType == 6){
    double mass = atof(ss);
    if ((int)floor(mass*10000) != 0)
      m_dMonoMass = mass;
    m_uiInputType = 5;
  } else if (m_uiInputType == 7){
    double mass = atof(ss);
    if ((int)floor(mass*10000) != 0)
      m_dAveMass = mass;
    m_uiInputType = 5;
  }
}

bool SAXResidHandler::load(){
  ifstream ifXml;
  /*
   * open an input stream, return false if it can't be opened
   */
  ifXml.open(m_strXmlPath.c_str());
  if(ifXml.fail()){
    //		cout << "\nFailed to open: \"" << m_strXmlPath.c_str() << "\"\n";
    return false;
  }
  setFileName( m_strXmlPath.data() );
  
  parse();
  return true;
}

SAXUnimodHandler::SAXUnimodHandler(const string& _p, vector<Modification> *_mods)
{
  m_strXmlPath	= _p;
  m_vmodModificationsDB = _mods;
  m_uiInputType	= 0;
  m_uiId			= 0;
  m_strModType	= "";
  m_dMonoMass		= 0.0;
  m_dAveMass		= 0.0;
  m_strName		= "";
  m_strResidues	= "";
  m_strUnimod		= "";		 
}

SAXUnimodHandler::~SAXUnimodHandler()
{
}

void SAXUnimodHandler::startElement(const XML_Char *el, const XML_Char **attr)
{
  if(isElement("umod:mod", el)){
    for (int i = 0; attr[i]; i+=2){
      if (isAttr("record_id",attr[i])){
	m_uiId = atoi(attr[i+1]);
	m_strName.append("UnimodRecordID_");
	m_strName.append(attr[i+1]);
      }
    }
    return;
  }
  if(isElement("umod:specificity", el)){
    bool position = false;
    bool classification = false;
    char site = '\0';
    for (int i = 0; attr[i]; i+=2){
      if (isAttr("position", attr[i])){
	if (strcmp(attr[i+1], "Anywhere") == 0){
	  position = true;
	}
      }
      if (isAttr("classification", attr[i])){
	if (strcmp(attr[i+1], "Post-translational") == 0){
	  classification = true;
	}
      }
      if (isAttr("site", attr[i])){
	site = attr[i+1][0];
      }
    }
    if (classification && position)
      m_strResidues += site;
    return;
  }
  if(isElement("umod:delta", el)){
    for (int i = 0; attr[i]; i+=2){
      if (isAttr("mono_mass", attr[i])){
	m_dMonoMass = atof(attr[i+1]);
      }
      if (isAttr("avge_mass", attr[i])){
	m_dAveMass = atof(attr[i+1]);
      }
    }
    return;
  }
}

void SAXUnimodHandler::endElement(const XML_Char *el)
{
  size_t pos = 0;
  size_t len;
  Modification mod;
  char resid;
  if(isElement("umod:mod", el)){
    len = m_strResidues.length();
    for (size_t i = 0; i < len; i++){
      mod.m_dAveMass = m_dAveMass;
      mod.m_dMonoMass = m_dMonoMass;
      mod.m_strName = m_strName;
      mod.m_strUnimod = m_strUnimod;
      mod.m_uiModType = 0;
      resid = m_strResidues[i];
      mod.m_strResidues = resid;
      m_vmodModificationsDB[resid].push_back(mod);
    }
    //restore variables to default;
    m_uiId			= 0;
    m_strModType	= "";
    m_strName		= "";
    m_dMonoMass		= 0.0;
    m_dAveMass		= 0.0;
    m_strResidues	= "";
    m_strUnimod		= "";		 
  }
}
// changed in 2006.09.15 to improve handling of characters
// suggested by Brendan Maclean
void SAXUnimodHandler::characters(const XML_Char *s, int len)
{
}

bool SAXUnimodHandler::load(){
  ifstream ifXml;
  /*
   * open an input stream, return false if it can't be opened
   */
  ifXml.open(m_strXmlPath.c_str());
  if(ifXml.fail()){
    //		cout << "\nFailed to open: \"" << m_strXmlPath.c_str() << "\"\n";
    return false;
  }
  setFileName( m_strXmlPath.data() );
  
  parse();
  return true;
}
