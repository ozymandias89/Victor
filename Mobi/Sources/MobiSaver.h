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


	virtual void save_mob_everageScalD();

private:
	bool verbose;
	string out;
};
}
}

#endif /* MOBI_SOURCES_MOBISAVER_H_ */
