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

ScaleDistance::ScaleDistance(const ProteinModels& modelli, bool standardDeviation ,bool verbose) :
		standardDev(standardDeviation), verbose(verbose), models(modelli.models) {

	Spacer* primo_modello = new Spacer;
	Spacer* secondo_modello = new Spacer;

	unsigned int num_atomi = 0;
	double ScalD, distance= -1.0;

	if (models.size() != 0)
		num_atomi = models[0].sizeAmino();
	else
		ERROR("Nessun modello presente nella proteina", exception);




	vector<vector<double> > dist_from_Ca_atoms(num_atomi,
			vector<double>(0.0));

	cout << "################# " << dist_from_Ca_atoms.size() << "   " << dist_from_Ca_atoms[0].size() << endl;


	if (verbose){
		cout << "Numero di modelli presenti: " << models.size() << endl;
		cout << "Ogni modello ha un numero di atomi pari a: " << num_atomi << endl;
	}
	//scorro i modelli e per ognuno creo il vettore di atomi
	unsigned int u = 0;

	while (u < (models.size())) {
		//if (verbose)
		//cout << "#modello numero: " << u << endl;
		if (verbose)
				cout << "#Distanza tra modello: " << u;

		primo_modello = &models[u];
		getCaAtom(primo_modello, true);
		u++;
		//if (verbose)
		//cout << "#modello numero: " << u << endl;

		if (verbose)
						cout << " e modello: " << u << endl;

		secondo_modello = &models[u];
		getCaAtom(secondo_modello, false);
		u++;


		//cout << "GRandezza vettore " << dist_from_Ca_atoms[0].size() << endl;

		for (unsigned int i = 0; i < num_atomi; i++) {
			distance = CaVector1[i].distance(CaVector2[i]);
			ScalD = 1.0 / (pow (( 1.0 + (distance/4.0) ), 2.0));
			dist_from_Ca_atoms[i].push_back(ScalD);
			if (verbose)
			cout << dist_from_Ca_atoms[i].back() << endl;
		}

		//cout << "GRAndezza vettore " << dist_from_Ca_atoms[0].size() << endl;

		CaVector1.clear();
		CaVector2.clear();

	}



		int count;
		double sum;
		for (unsigned int i = 0; i < num_atomi; i++) {
			count=0;
			sum=0;
			cout << "Grandezza vettore " << dist_from_Ca_atoms[0].size() << endl;
			vector<double>::iterator atom = dist_from_Ca_atoms[i].begin();
			while (atom != dist_from_Ca_atoms[i].end()){

				sum = sum + (*atom);

				atom++;
				count++;
			}
			//cout << "SOMMA " << sum << " COUNT " << count << endl;
			ScD.push_back((sum / count));
			cout << "media per il " << i << " atomo = " << ScD.back() << endl;


		}

		if (standardDev){


			double standDev=-1;

			cout << "#############################" << endl;

			vector<double>::iterator everage = ScD.begin();

			for (unsigned int i = 0; i < num_atomi; i++) {
				sum=0;
				count=0;
				vector<double>::iterator atom = dist_from_Ca_atoms[i].begin();


				while (atom != dist_from_Ca_atoms[i].end()){
					cout << "Stampo distanza fra atomi " << *atom << endl;
					sum += (pow ( ((*atom) - (*everage) ) , 2.0 ));
							atom++;
							count++;

				}

				standDev = sqrt(sum / count);

				cout << "sum : " << sum << " COUNT " << count << endl;
				cout << "media : " << *everage << endl;
				cout << "DEVSTADARD per il " << i << " atomo = " << standDev << endl;
				*everage = standDev;

				everage++;


			}


		}






}

ScaleDistance::~ScaleDistance() {
}


void ScaleDistance::getCaAtom(Spacer* s, bool flag) {

//per un modello fissato creo il vettore di atomi
	AminoAcid* a = new AminoAcid;

	for (unsigned int u = 0; u < s->sizeAmino(); u++) {
		a = &(s->getAmino(u));


		unsigned int i = 1;

		if (flag)
			this->CaVector1.push_back(a->getAtom(i));
		else
			this->CaVector2.push_back(a->getAtom(i));

		if (a->getAtom(i).getCode() != CA)
			ERROR("L'atomo trovato non è un carbonio alfa.", exception);

	}

}



vector <double> ScaleDistance::get_ScalDist(){
	return ScD;
}
