/*
 * Portions are Copyright (c) 2003-2006 Fred Hutchinson Cancer Research Center
 * Additional code Copyright (c) 2010-2011 Institute for Systems Biology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * Adapted from TPP project (tools.proteomecenter.org/wiki/) on 2013-11-26 
 */

#ifndef MSCORE_HRK_H
#define MSCORE_HRK_H
#include "mscore.h" //added rTANDEM
#include "mscore_k.h" //Needed for miLookup class
class mscore_hrk : public mscore {
 protected:
  friend class mscorefactory_hrk;
  mscore_hrk(void); // Should only be created through mscorefactory_tandem
 public:
  ~mscore_hrk(void);	//Destructor
  bool add_mi(mspectrum &_s);
  bool add_mi_hr(mspectrum &_s);
  bool clear();
  bool load_param(XmlParameter &_x);	// allows score object to issue warnings, or set variables based on xml
  unsigned long mconvert(double _m, const long _c); // convert mass to integer ion m/z for mi vector
  bool precondition(mspectrum &_s);	// called before spectrum conditioning
  void prescore(const size_t _i);	// called before scoring
  void report_score(char* _buff, float _h);	// format hyper score for output
  double sfactor();	// factor applied to final convolution score
 protected:
  bool add_A(const unsigned long _t,const long _c);
  bool add_B(const unsigned long _t,const long _c);
  bool add_C(const unsigned long _t,const long _c);
  bool add_Y(const unsigned long _t,const long _c);
  bool add_X(const unsigned long _t,const long _c);
  bool add_Z(const unsigned long _t,const long _c);
  double dot(unsigned long *_v);	// this is where the real scoring happens
  double dot_hr(unsigned long *_v);
  float ion_check(const unsigned long _v,const size_t _d);
 protected:
  unsigned long imass(double _m){
    return (unsigned long)((_m/m_dIsotopeCorrection) + 0.5);
  }
  //Data Members
 protected:
  double	m_dIsotopeCorrection;
  int	m_maxEnd;
  miLookup	m_miUsed;
  vector<vmiType> m_vmiType;
  vector<int>	m_vmiUsed;
};
/*
 * mscorefactory_hrk implements a factory for creating mscore_hrk instances.
 */
class mscorefactory_hrk : public mpluginfactory {
 public:
  mscorefactory_hrk();
  mplugin* create_plugin();
};
#endif 
