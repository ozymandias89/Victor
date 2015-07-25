/*
 * MobiSaver.h
 *
 *  Created on: 22 lug 2015
 *      Author: riccardo
 */

#ifndef MOBI_SOURCES_MOBISAVER_H_
#define MOBI_SOURCES_MOBISAVER_H_

#include <Saver.h>
#include <ProteinModels.h>
#include <vector>
#include <iterator>
#include <iostream>
#include <fstream>
#include <unistd.h>

using namespace Victor::Biopool;

namespace Victor {
namespace Mobi {

class MobiSaver : public Saver{
public:
	MobiSaver(ProteinModels __protein,string __output, bool __verbose = false);
	virtual ~MobiSaver();


	void mob_eveScalD(vector <double> __everageDistance);
	void mob_eveScalD_filtSecS(vector <double> __everageDistance, vector<char> __mobSecS);
	void mob_stanD(vector <double> __Scale_distance);
	void mob_stanD_withMask(vector <double> __Scale_distance);
	void mob_aPHI(vector <double> __angle_PHI);
	void mob_aPSI(vector <double> __angle_PSI);
	void mob_SecS(vector <char> __Mob_SecStructure);
	void allMobility(vector<double> __everageDistance,
			vector<double> __Scale_distance, vector<double> __angle_PHI,
			vector<double> __angle_PSI, vector<char> __Mob_SecStructure);

private:
	bool verbose;
	string out;
};
}
}

#endif /* MOBI_SOURCES_MOBISAVER_H_ */
