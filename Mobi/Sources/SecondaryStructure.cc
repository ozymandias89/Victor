/*
 * SecondaryStructure.cc
 *
 *  Created on: 21 lug 2015
 *      Author: riccardo
 */


#include "SecondaryStructure.h"

using namespace Victor;
using namespace Victor::Mobi;
using namespace Victor::Biopool;

SecondaryStructure::SecondaryStructure(const ProteinModels& modelli ,bool verbose) :
	 verbose(verbose), models(modelli.original_models) {}

SecondaryStructure::~SecondaryStructure() {}



vector<char> SecondaryStructure::getMobilitySecondaryStructure() {

	mobility.clear();
	unsigned int num_amino = 0;
	if (models.size() != 0)
		num_amino = models[0].sizeAmino();
	else
		ERROR(
				"Nessun modello presente nella proteina, si ricorda di caricare dei modelli.",
				exception);

	//create a support vector to get secondary structure
	vector<set<char> > vettore_di_supporto(num_amino, set<char>());

	//create a vector for manipulated all amiocid's chars
	vector<vector<char> > sec_Structures(num_amino,
			vector<char>(models.size()));

	if (verbose) {
		cout << "\nStart procedure select secondary structure... " << endl;
		cout << "Model size: " << models.size() << endl;
	}
	//per ogni modello
	for (unsigned int i = 0; i < models.size(); i++) {

		if (verbose)
			cout << "Models number: " << i <<endl;
		models[i].setDSSP(false);
		vettore_di_supporto = models[i].getDSSP();
		//cout << "support vector size: " << vettore_di_supporto.size() << endl;

		//scorro il vettore di set
		for (unsigned int j = 0; j < num_amino; j++) {
			//mi tiro fuori il primo valore char e lo metto in un vettore di char
			sec_Structures[j][i] = (*(vettore_di_supporto[j].begin()));
			if (verbose)
				cout << sec_Structures[j][i];
		}
		if (verbose)
		cout  << endl;
	}






		char test;
		bool flag;

		for (unsigned int i = 0; i < num_amino; i++) {

			test = sec_Structures[i][0];
			flag = true;

			//se la prima lettera è vuota oppure è una S

			if ((test != 'G') && (test != 'H') && (test != 'I') && (test != 'T')
					&& (test != 'E') && (test != 'B')) {
				//verifico che siano tutte S e/o vuote

				for (unsigned int j = 0; j < models.size(); j++) {

					if ((sec_Structures[i][j] != 'G')
							&& (sec_Structures[i][j] != 'H')
							&& (sec_Structures[i][j] != 'I')
							&& (sec_Structures[i][j] != 'T')
							&& (sec_Structures[i][j] != 'E')
							&& (sec_Structures[i][j] != 'B')) {
					} else {

						mobility.push_back('M');
						flag = false;
						break;

					}

				}
				if (flag) {
					mobility.push_back('c');
				}
				//la prima lettera è una lettera di struttura
			} else {

				//verifico che siano tutte uguali a test
				for (unsigned int j = 0; j < models.size(); j++) {

					if ((sec_Structures[i][j] == test)) {
					} else {

						mobility.push_back('M');
						flag = false;
						break;

					}

				}
				if (flag) {
					mobility.push_back('.');
				}

			}

		}



	return mobility;
}

