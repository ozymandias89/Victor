/*
 * ScaleDistance.cc
 *
 *  Created on: 11 lug 2015
 *      Author: riccardo
 */
#include <ScaleDistance.h>

using namespace Victor;
using namespace Victor::Mobi;
using namespace Victor::Biopool;

ScaleDistance::ScaleDistance(const ProteinModels& modelli, bool verbose) :
	verbose(verbose) {

	//vector <int>* ScaleD = new vector<int>(5);

	Spacer primo_modello;
	Spacer secondo_modello;

	unsigned int u = 0;
	while (u <((modelli.size())/2)) {
		if (verbose)
		cout << "#modello numero: " << u << endl;
		primo_modello = modelli.models[u];
		this-> getCaAtom(primo_modello, true);
		u++;
		cout << "#modello numero: " << u << endl;
		secondo_modello = modelli.models[u];
		this-> getCaAtom(secondo_modello, false);
		u++;
	}

	//cout << "##### Atomo:" << this->CaVector1[0].getCode() << endl;




}



void ScaleDistance :: getCaAtom (Spacer& s, bool flag){

	for (unsigned int u=0; u < s.sizeAmino(); u++) {
		AminoAcid a = s.getAmino(u);
		if (verbose)
			cout << "#Aminoacido numero: " << u << ":" << endl;

		unsigned int i = 1;

		if (flag)
			this->CaVector1.push_back(a.getAtom(i));
		else
			this->CaVector2.push_back(a.getAtom(i));


		if (a.getAtom(i).getCode()!= CA)
			ERROR("L'atomo trovato non è un carbonio alfa.", exception);

	}

}

ScaleDistance::~ScaleDistance() {
}

