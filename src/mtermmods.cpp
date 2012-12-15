/*
 Copyright (C) 2005 Ronald C Beavis, all rights reserved
 X! tandem 
 This software is a component of the X! proteomics software
 development project

Use of this software governed by the Artistic license, as reproduced here:

The Artistic License for all X! software, binaries and documentation

Preamble
The intent of this document is to state the conditions under which a
Package may be copied, such that the Copyright Holder maintains some 
semblance of artistic control over the development of the package, 
while giving the users of the package the right to use and distribute 
the Package in a more-or-less customary fashion, plus the right to 
make reasonable modifications. 

Definitions
"Package" refers to the collection of files distributed by the Copyright 
	Holder, and derivatives of that collection of files created through 
	textual modification. 

"Standard Version" refers to such a Package if it has not been modified, 
	or has been modified in accordance with the wishes of the Copyright 
	Holder as specified below. 

"Copyright Holder" is whoever is named in the copyright or copyrights 
	for the package. 

"You" is you, if you're thinking about copying or distributing this Package. 

"Reasonable copying fee" is whatever you can justify on the basis of 
	media cost, duplication charges, time of people involved, and so on. 
	(You will not be required to justify it to the Copyright Holder, but 
	only to the computing community at large as a market that must bear 
	the fee.) 

"Freely Available" means that no fee is charged for the item itself, 
	though there may be fees involved in handling the item. It also means 
	that recipients of the item may redistribute it under the same
	conditions they received it. 

1. You may make and give away verbatim copies of the source form of the 
Standard Version of this Package without restriction, provided that 
you duplicate all of the original copyright notices and associated 
disclaimers. 

2. You may apply bug fixes, portability fixes and other modifications 
derived from the Public Domain or from the Copyright Holder. A 
Package modified in such a way shall still be considered the Standard 
Version. 

3. You may otherwise modify your copy of this Package in any way, provided 
that you insert a prominent notice in each changed file stating how and 
when you changed that file, and provided that you do at least ONE of the 
following: 

a.	place your modifications in the Public Domain or otherwise make them 
	Freely Available, such as by posting said modifications to Usenet 
	or an equivalent medium, or placing the modifications on a major 
	archive site such as uunet.uu.net, or by allowing the Copyright Holder 
	to include your modifications in the Standard Version of the Package. 
b.	use the modified Package only within your corporation or organization. 
c.	rename any non-standard executables so the names do not conflict 
	with standard executables, which must also be provided, and provide 
	a separate manual page for each non-standard executable that clearly 
	documents how it differs from the Standard Version. 
d.	make other distribution arrangements with the Copyright Holder. 

4. You may distribute the programs of this Package in object code or 
executable form, provided that you do at least ONE of the following: 

a.	distribute a Standard Version of the executables and library files, 
	together with instructions (in the manual page or equivalent) on 
	where to get the Standard Version. 
b.	accompany the distribution with the machine-readable source of the 
	Package with your modifications. 
c.	give non-standard executables non-standard names, and clearly 
	document the differences in manual pages (or equivalent), together 
	with instructions on where to get the Standard Version. 
d.	make other distribution arrangements with the Copyright Holder. 

5. You may charge a reasonable copying fee for any distribution of 
this Package. You may charge any fee you choose for support of 
this Package. You may not charge a fee for this Package itself. 
However, you may distribute this Package in aggregate with other 
(possibly commercial) programs as part of a larger (possibly 
commercial) software distribution provided that you do not a
dvertise this Package as a product of your own. You may embed this 
Package's interpreter within an executable of yours (by linking); 
this shall be construed as a mere form of aggregation, provided that 
the complete Standard Version of the interpreter is so embedded. 

6. The scripts and library files supplied as input to or produced as 
output from the programs of this Package do not automatically fall 
under the copyright of this Package, but belong to whomever generated 
them, and may be sold commercially, and may be aggregated with this 
Package. If such scripts or library files are aggregated with this 
Package via the so-called "undump" or "unexec" methods of producing 
a binary executable image, then distribution of such an image shall 
neither be construed as a distribution of this Package nor shall it 
fall under the restrictions of Paragraphs 3 and 4, provided that you 
do not represent such an executable image as a Standard Version of 
this Package. 

7. C subroutines (or comparably compiled subroutines in other languages) 
supplied by you and linked into this Package in order to emulate 
subroutines and variables of the language defined by this Package 
shall not be considered part of this Package, but are the equivalent 
of input as in Paragraph 6, provided these subroutines do not change 
the language in any way that would cause it to fail the regression 
tests for the language. 

8. Aggregation of this Package with a commercial distribution is always 
permitted provided that the use of this Package is embedded; that is, 
when no overt attempt is made to make this Package's interfaces visible 
to the end user of the commercial distribution. Such use shall not be 
construed as a distribution of this Package. 

9. The name of the Copyright Holder may not be used to endorse or promote 
products derived from this software without specific prior written permission. 

10. THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED 
WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF 
MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE. 

The End 
*/

// File version: 2005-09-01

#include "stdafx.h"
#include <sys/timeb.h>
#include <ctime>
#include "mtermmods.h"
#include "mprocess.h"
#include "mscore.h"

mtermmods::mtermmods(){
	m_pProcess=NULL;
}
mtermmods::~mtermmods(){
}

bool mtermmods::set_mprocess(mprocess *_p){
	if(m_pProcess != NULL)
		delete m_pProcess;
	m_pProcess = _p;
	return true;
}

/*
 * process N or C terminal modifications using cleavage at any site and 50 missed cleavages
 */
bool mtermmods::refine(){
	vector<mspectrum> vspRollback;
	string strKey;
	string strValue;
	string strTemp;
	size_t a = 0;
	size_t tActiveNow = 0;

/*
 *  set the maximum expectation value
 */
	strKey = "refine, maximum valid expectation value";
	m_pProcess->m_xmlValues.get(strKey,strValue);
	if(strValue.size() > 0)	{
		m_dMaxExpect = atof(strValue.c_str());
	}
	if(m_pProcess->m_bRefineCterm)
		strKey = "refine, potential C-terminus modifications";
	else
		strKey = "refine, potential N-terminus modifications";

	m_pProcess->m_xmlValues.get(strKey,strValue);
	if(m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF)	{
		if(m_pProcess->m_bRefineCterm)
//			cout << "\tmodified C-terminus ";
			Rprintf("\tmodified C-terminus ");
		else
//			cout << "\tmodified N-terminus ";
			Rprintf("\tmodified N-terminus ");
		//cout.flush();
		m_pProcess->m_prcLog.log("modified terminus");
	}
	m_pProcess->create_rollback(vspRollback);
	strTemp = "[X]|[X]";
	m_pProcess->m_Cleave.load(strTemp);
	m_pProcess->m_tMissedCleaves = 50;
/*
 * score the spectra against each sequence
 */
	m_pProcess->score_terminus(strValue);
/*
 * update active state of spectra
 */
	m_pProcess->load_best_vector();
	while(a < m_pProcess->m_vSpectra.size())	{
		if(!m_pProcess->m_vSpectra[a].m_bActive)
			tActiveNow++;
		a++;
	}
	if(tActiveNow >= m_pProcess->m_tActive)	{
		if(m_pProcess->m_bRefineCterm)
			m_pProcess->m_tRefineCterminal = tActiveNow - m_pProcess->m_tActive;
		else
			m_pProcess->m_tRefineNterminal = tActiveNow - m_pProcess->m_tActive;
	}
	m_pProcess->rollback(vspRollback,m_dMaxExpect,0.1);
	m_pProcess->m_tActive = tActiveNow;
	if(m_pProcess->m_lThread == 0 || m_pProcess->m_lThread == 0xFFFFFFFF)	{
//		cout << " done.\n";
		Rprintf(" done.\n");
		m_pProcess->m_prcLog.log("done");
	}
	//cout.flush();

	return true;
}

/*
 * mpmodsmanager contains static short-cuts for dealing with mpmods
 * plug-ins.
 */
const char* mtermmodsmanager::TYPE = "refinement, terminal modifications algorithm";

/*
 * create_mtermmods creates the correct mtermmods object for the given set of XML
 * parameters.
 */
mtermmods* mtermmodsmanager::create_mtermmods(XmlParameter &_x)
{
	string strValue;
	string strKey = TYPE;
	if (!_x.get(strKey,strValue))
		strValue = "tandem";
	return (mtermmods*) mpluginmanager::get().create_plugin(TYPE, strValue.data());
}

/*
 * register_factory registers a factory with the mpluginmanager for
 * creating mtermmods derived objects.
 */
void mtermmodsmanager::register_factory(const char* _spec, mpluginfactory* _f)
{
	mpluginmanager::get().register_factory(TYPE, _spec, _f);
}


// Factory instance, registers itself with the mtermmodsmanager.
static mtermmodsfactory_tandem factory;
	
mtermmodsfactory_tandem::mtermmodsfactory_tandem()
{
	mtermmodsmanager::register_factory("tandem", this);
}

mplugin* mtermmodsfactory_tandem::create_plugin()
{
	return new mtermmods();
}
