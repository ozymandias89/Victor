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


	void save_mob_everageScalD(vector <double> __everageDistance);
	void save_mob_standardDeviation(vector <double> __Scale_distance);
	void save_mob_angle_PHI(vector <double> __angle_PHI);
	void save_mob_angle_PSI(vector <double> __angle_PSI);
	void save_mob_SecondaryStructure(vector <char> __Mob_SecStructure);
	void save_allMobility(vector<double> __everageDistance,
			vector<double> __Scale_distance, vector<double> __angle_PHI,
			vector<double> __angle_PSI, vector<char> __Mob_SecStructure);

private:
	bool verbose;
	string out;
};
}
}

#endif /* MOBI_SOURCES_MOBISAVER_H_ */
