/* PTMTreeSearch plugin
For furter details see or if PTMTreeSearch is used in a study please cite the following article:
PTMTreeSearch: a novel two-stage tree-search algorithm with pruning rules for the identification of post-translational modification of proteins in MS/MS spectra
Attila Kertesz-Farkas, Beata Reiz, Roberto Vera, Michael P Myers, Sandor Pongor
Bioinformatics 30 (2), 234-241, 2014 

FILE version: 2013-Febr-11  // Modified for rTANDEM on January 2014
*/

#ifndef PTMTREESEARCH_H
#define PTMTREESEARCH_H

//#define PARALELL_PTMTREESEARCH

#include "mrefine.h"
#include "saxhandler.h"

class mpeptide;
class msequence;
class mprocess;
class Modification{
public:
	double m_dTarget;
	double m_dTargetEval;
	unsigned int m_uiCount;
	unsigned int m_uiModType;
	double m_dMonoMass;
	double m_dAveMass;
	string m_strResidues;
	string m_strName;
	string m_strUnimod;
	bool m_bType;		//true = monoisotpic, false = average isotopic.
	Modification(){
		m_dMonoMass = 0.0;
		m_dAveMass  = 0.0;
		m_uiModType = 0;
		m_uiCount	= 0;
		m_dTarget = 0.0;
		m_dTargetEval = 30.0;
		m_bType = true;
	}
	bool operator <(const Modification &m) const{
		
		return m_bType ? m_dMonoMass < m.m_dMonoMass : m_dAveMass < m.m_dAveMass;
	}
};

class PTMTreeSearch : public mrefine
{
public:
	PTMTreeSearch();
	~PTMTreeSearch(void);
	bool refine();			//entry point for potential modifications processing
	bool set_mprocess(mprocess *_p); //set the m_pProcess object
	mprocess *m_pProcess;	//pointer to the parent mprocess class
	vector<Modification> m_vmodModificationsDB[255];  //used to store all information about a modification
	vector<double> m_vdModifications[255];		//used only the modification mass for ptm identification in order to carry out calculations faster
	bool openUniprotPTMs(const char* filename, vector<Modification>* _mods);
};

/*
*	ptmsearchmanager is used to organize and facilitate the creation of PTMTreeSearch objects
*/
class PTMTreeSearchmanager
{
public:
	static const char* TYPE;
	static PTMTreeSearch* create_PTMTreeSearch(XmlParameter &_x);
	static void register_factory(const char* _spec, mpluginfactory* _f);
};

/*
 * mrefinefactory_tandem implements a factory for creating PTMTreeSearch
 * instances.
 */
class PTMTreeSearchfactory_tandem : public mpluginfactory
{
public:
	PTMTreeSearchfactory_tandem();
	virtual mplugin* create_plugin();
};


/**
* Uses eXpat SAX parser to parse Modifications Input data.
*/
class SAXModsHandler : public SAXHandler
{
public:
	SAXModsHandler(const string& _p, vector<Modification>* _map);
	~SAXModsHandler();
	bool load();
protected:
	// -----------------------------------------------------------------------
	//  Overrides of SAXHandler functions
	// -----------------------------------------------------------------------
	/**
	* Receive notification of the start of an element.
	*
	* <p>By default, do nothing.  Application writers may override this
	* method in a subclass to take specific actions at the start of
	* each element (such as allocating a new tree node or writing
	* output to a file).</p>
	*/
	void startElement(const XML_Char *el, const XML_Char **attr);

	/**
	* Receive notification of the end of an element.
	*
	* <p>By default, do nothing.  Application writers may override this
	* method in a subclass to take specific actions at the end of
	* each element (such as finalising a tree node or writing
	* output to a file).</p>
	*/
	void endElement(const XML_Char *el);

	/**
	* Receive notification of character data inside an element.
	*
	* <p>By default, do nothing.  Application writers may override this
	* method to take specific actions for each chunk of character data
	* (such as adding the data to a node or buffer, or printing it to
	* a file).</p>
	*/
	void characters(const XML_Char *s, int len);
	string m_strXmlPath;
	string m_strKey;
	vector<Modification>* m_vmodModificationsDB;
	unsigned int m_uiInputType;			//type of the XML node. 
	unsigned int m_uiModType;	//type of the modification
	double m_dMonoMass;			//monoisotopic mass
	double m_dAveMass;			//average mas
	string m_strResidues;		//the residues that can be modified
	string m_strUnimod;			//unimod id
	string m_strName;			//name of the modification
	unsigned int m_uiId;		//unique id of the modification
};
/**
* Uses eXpat SAX parser to parse Modifications Input data.
*/
class SAXResidHandler : public SAXHandler
{
public:
	SAXResidHandler(const string& _p, vector<Modification>* _map);
	~SAXResidHandler();
	bool load();
protected:
	// -----------------------------------------------------------------------
	//  Overrides of SAXHandler functions
	// -----------------------------------------------------------------------
	/**
	* Receive notification of the start of an element.
	*
	* <p>By default, do nothing.  Application writers may override this
	* method in a subclass to take specific actions at the start of
	* each element (such as allocating a new tree node or writing
	* output to a file).</p>
	*/
	void startElement(const XML_Char *el, const XML_Char **attr);

	/**
	* Receive notification of the end of an element.
	*
	* <p>By default, do nothing.  Application writers may override this
	* method in a subclass to take specific actions at the end of
	* each element (such as finalising a tree node or writing
	* output to a file).</p>
	*/
	void endElement(const XML_Char *el);

	/**
	* Receive notification of character data inside an element.
	*
	* <p>By default, do nothing.  Application writers may override this
	* method to take specific actions for each chunk of character data
	* (such as adding the data to a node or buffer, or printing it to
	* a file).</p>
	*/
	void characters(const XML_Char *s, int len);
	string m_strXmlPath;
	string m_strKey;
	vector<Modification>* m_vmodModificationsDB;
	unsigned int m_uiInputType;			//type of the XML node. 
	string m_strModType;		//type of the modification
	double m_dMonoMass;			//monoisotopic mass
	double m_dAveMass;			//average mas
	string m_strResidues;		//the residues that can be modified
	string m_strUnimod;			//unimod id
	string m_strName;			//name of the modification
	unsigned int m_uiId;		//unique id of the modification
};

class SAXUnimodHandler : public SAXHandler
{
public:
	SAXUnimodHandler(const string& _p, vector<Modification>* _map);
	~SAXUnimodHandler();
	bool load();
protected:
	// -----------------------------------------------------------------------
	//  Overrides of SAXHandler functions
	// -----------------------------------------------------------------------
	/**
	* Receive notification of the start of an element.
	*
	* <p>By default, do nothing.  Application writers may override this
	* method in a subclass to take specific actions at the start of
	* each element (such as allocating a new tree node or writing
	* output to a file).</p>
	*/
	void startElement(const XML_Char *el, const XML_Char **attr);

	/**
	* Receive notification of the end of an element.
	*
	* <p>By default, do nothing.  Application writers may override this
	* method in a subclass to take specific actions at the end of
	* each element (such as finalising a tree node or writing
	* output to a file).</p>
	*/
	void endElement(const XML_Char *el);

	/**
	* Receive notification of character data inside an element.
	*
	* <p>By default, do nothing.  Application writers may override this
	* method to take specific actions for each chunk of character data
	* (such as adding the data to a node or buffer, or printing it to
	* a file).</p>
	*/
	void characters(const XML_Char *s, int len);
	string m_strXmlPath;
	string m_strKey;
	vector<Modification>* m_vmodModificationsDB;
	unsigned int m_uiInputType;			//type of the XML node. 
	string m_strModType;		//type of the modification
	double m_dMonoMass;			//monoisotopic mass
	double m_dAveMass;			//average mas
	string m_strResidues;		//the residues that can be modified
	string m_strUnimod;			//unimod id
	string m_strName;			//name of the modification
	unsigned int m_uiId;		//unique id of the modification
};

#endif

