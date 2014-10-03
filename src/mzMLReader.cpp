#include "mzParser.h"

using namespace std;

int main(int argc, char* argv[]){

	if(argc!=2) {
	  Rprintf("USAGE: mzMLReader <mzML file>\n");
	  //exit(0);
	  return(0);
	}

	MSDataFile* msd;
	ChromatogramListPtr sl;
	ChromatogramPtr s2;
	string st=argv[1];
	vector<TimeIntensityPair> pairs;

	msd = new MSDataFile(argv[1]);
	if(!msd->run.chromatogramListPtr->get()) Rprintf("WTF\n");

	sl = msd->run.chromatogramListPtr;
	for (int j=1; j<(int)sl->size(); j++) {
		s2 = sl->chromatogram(j, true);
		Rprintf("%i\t", j);
		Rprintf("%s\n", s2->id.c_str());
		// cout << j << "\t" << s2->id << endl;
		s2->getTimeIntensityPairs(pairs);
		for(int k=0;k<(int)pairs.size();k++) {
		  Rprintf("%lu ", pairs[k].time);
		  Rprintf("%lu\n", pairs[k].intensity);
		  //cout << pairs[k].time << " " << pairs[k].intensity << endl;
		}
	}

	//exit(1);
	return(1);

	BasicSpectrum s;
	BasicChromatogram chromat;
	MzParser sax(&s,&chromat);
	sax.load(argv[1]);

	bool bLastWasSpectrum=true;
	char c='a';
	char str[256];
	int num;
	while(c!='x'){

		if(bLastWasSpectrum){
		  Rprintf("\nCurrent specturm:\n");
		  Rprintf("\tScan number: %i\n", s.getScanNum());
		  Rprintf("\tRetention Time: %ul\n", s.getRTime());
		  Rprintf("\tMS Level: %i\n", s.getMSLevel());
		  Rprintf("\tNumber of Peaks: %i\n", s.size());
		  //	cout << "\nCurrent spectrum:" << endl;
		  //	cout << "\tScan number: " << s.getScanNum() << endl;
		  //	cout << "\tRetention Time: " << s.getRTime() << endl;
		  //	cout << "\tMS Level: " << s.getMSLevel() << endl;
		  //	cout << "\tNumber of Peaks: " << s.size() << endl;
		} else {
			chromat.getIDString(str);
			Rprintf("\nCurrent chromatogram:\n");
			Rprintf("\tID: %s\n", str);
			Rprintf("\tNumber of Peaks: %i\n", chromat.size());
			// cout << "\nCurrent chromatogram:" << endl;
			// cout << "\tID: " << str << endl;
			//  cout << "\tNumber of Peaks: " << chromat.size() << endl;
		}
		Rprintf("\nMenu:\n\t'c' to grab a new chromatogram\n\t's' to grab a new spectrum\n\t'p' to show peaks\n\t'x' to exit\n");
		Rprintf("Please enter your choice: ");
		// cout << "\nMenu:\n\t'c' to grab a new chromatogram\n\t's' to grab a new spectrum\n\t'p' to show peaks\n\t'x' to exit" << endl;
		// cout << "Please enter your choice: ";
		cin >> c;

		switch(c){
			case 'c':
				if(sax.highChromat()==0){
				  Rprintf("No chromatograms in the file.\n");
					//cout << "No chromatograms in the file." << endl;
				} else {
				  Rprintf("Enter a number from 0 to %i: ", (int) sax.highChromat()-1);
				  //	cout << "Enter a number from 0 to " << sax.highChromat()-1 << ": ";
					cin >> str;
					num=(int)atoi(str);
					if(num<0 || num>sax.highChromat()-1) {
					  Rprintf("Bad number! BOOOOO!\n");
					  //cout << "Bad number! BOOOOO!" << endl;
					} else {
					  if(!sax.readChromatogram(num)) Rprintf("Chromatogram number not in file.\n");
						else bLastWasSpectrum=false;
					}
				}
				break;
			case 'p':
				if(bLastWasSpectrum){
					for(unsigned int i=0;i<s.size();i++) printf("%.6lf\t%.1lf\n",s[i].mz,s[i].intensity);
				} else {
					for(unsigned int i=0;i<chromat.size();i++) printf("%.6lf\t%.1lf\n",chromat[i].time,chromat[i].intensity);
				}
				break;
			case 's':
			  Rprintf("Enter a number from %i", (int) sax.lowScan());
			  Rprintf("to %i: ", (int)sax.highScan());
//			cout << "Enter a number from " << sax.lowScan() << " to " << sax.highScan() << ": ";
				cin >> str;
				num=(int)atoi(str);
				if(num<sax.lowScan() || num>sax.highScan()) {
				  Rprintf("Bad number! BOOOOO!\n");
				  //cout << "Bad number! BOOOOO!" << endl;
				} else {
				  if(!sax.readSpectrum(num)) Rprintf("pectrum number not in file.\n"); 
					  //cout << "Spectrum number not in file." << endl;
					else bLastWasSpectrum=true;
				}
				break;
			case 'x':
				break;
			default:
			  Rprintf("\nInvalid command!\n"); 
			  //	cout << "\nInvalid command!" << endl;
				break;
		}
	}


	return 0;
}
