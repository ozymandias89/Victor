/* This file is part of Victor
Victor is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
Victor is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with Victor. If not, see <http://www.gnu.org/licenses/>.
*/
/**
* @file PreteinModel.h
* @author Riccardo Zanella
* @date Lug 2015
* @version 0.1
*/
#ifndef MOBI_SOURCES_PROTEINMODEL_
#define MOBI_SOURCES_PROTEINMODEL_

#include<Protein.h>
#include<PdbLoader.h>
#include<AtomCode.h>


using namespace Victor::Biopool;

namespace Victor { namespace Mobi {
/**
* @brief Extends Protein class with functionalities related to manipulation of NMR model and models comparations.
* Scale Distance metric is provided.
*/
class ProteinModels : public Protein{


public:

	 // CONSTRUCTORS/DESTRUCTOR:
	ProteinModels()	;

	ProteinModels(const Protein& _orig);

	virtual ~ProteinModels();



	// PREDICATES:

	//void load(PdbLoader& pl, char chain, vector<unsigned int> models);

	virtual void load(PdbLoader& pl);

	unsigned int selectModels(PdbLoader& pl);

	void loadSameModels(PdbLoader& pl);

	//Spacer& getModel(unsigned int _model);

	// MODIFIERS:

	void setVerbose() {
		verbose = true;
	}

	void save(ProteinModels& prot, string outputFile);

	vector <Spacer> getModels();

	// MODIFIERS:
private:
	bool verbose;
	vector <Spacer> models;

};
}
}
#endif /* MOBI_SOURCES_PROTEINMODEL_ */
