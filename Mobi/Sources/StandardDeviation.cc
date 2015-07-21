/*
 * ScaleDistance.cc
 *
 *  Created on: 11 lug 2015
 *      Author: riccardo
 */
#include "StandardDeviation.h"

using namespace Victor;
using namespace Victor::Mobi;
using namespace Victor::Biopool;

StandardDeviation::StandardDeviation(const ProteinModels& modelli ,bool verbose) :
	 verbose(verbose), original_models(modelli.original_models) ,models(modelli.models) {}

StandardDeviation::~StandardDeviation() {
}


void StandardDeviation::getCaAtom(Spacer* s, bool flag) {

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

vector <double> StandardDeviation::get_everage_distance(){

	dist_everage.clear();

	Spacer* primo_modello = new Spacer;
	Spacer* secondo_modello = new Spacer;

			unsigned int num_atomi = 0;
			double ScalD, distance= -1.0;

			if (models.size() != 0)
				num_atomi = models[0].sizeAmino();
			else
				ERROR("Nessun modello presente nella proteina, si ricorda di eseguire TmImpose.", exception);




			vector<vector<double> > dist_from_Ca_atoms(num_atomi,
					vector<double>(0.0));



			if (verbose){
				cout << "Numero di modelli presenti: " << models.size() << endl;
				cout << "Ogni modello ha un numero di atomi pari a: " << num_atomi << endl;
			}
			//scorro i modelli e per ognuno creo il vettore di atomi
			unsigned int u = 0;

			while (u < (models.size())) {

				primo_modello = &models[u];
				getCaAtom(primo_modello, true);
				u++;

				secondo_modello = &models[u];
				getCaAtom(secondo_modello, false);
				u++;


				for (unsigned int i = 0; i < num_atomi; i++) {
					distance = CaVector1[i].distance(CaVector2[i]);
					ScalD = 1.0 / (1.0 + ( pow ((distance/4.0), 2.0)));
					dist_from_Ca_atoms[i].push_back(ScalD);
					}


				CaVector1.clear();
				CaVector2.clear();

			}



				int count;
				double sum;
				for (unsigned int i = 0; i < num_atomi; i++) {
					count=0;
					sum=0;
					vector<double>::iterator atom = dist_from_Ca_atoms[i].begin();
					while (atom != dist_from_Ca_atoms[i].end()){

						sum = sum + (*atom);

						atom++;
						count++;
					}
					dist_everage.push_back((sum / count));


				}

return dist_everage;

}

vector <double> StandardDeviation::get_standard_deviation(){

	dist_everage.clear();
	ScD.clear();

	Spacer* primo_modello = new Spacer;
		Spacer* secondo_modello = new Spacer;

		unsigned int num_atomi = 0;
		double ScalD, distance= -1.0;

		if (models.size() != 0)
			num_atomi = models[0].sizeAmino();
		else
			ERROR("Nessun modello presente nella proteina, si ricorda di eseguire TmImpose.", exception);




		vector<vector<double> > dist_from_Ca_atoms(num_atomi,
				vector<double>(0.0));



		if (verbose){
			cout << "Numero di modelli presenti: " << models.size() << endl;
			cout << "Ogni modello ha un numero di atomi pari a: " << num_atomi << endl;
		}
		//scorro i modelli e per ognuno creo il vettore di atomi
		unsigned int u = 0;

		while (u < (models.size())) {

			primo_modello = &models[u];
			getCaAtom(primo_modello, true);
			u++;

			secondo_modello = &models[u];
			getCaAtom(secondo_modello, false);
			u++;


			for (unsigned int i = 0; i < num_atomi; i++) {
				distance = CaVector1[i].distance(CaVector2[i]);
				ScalD = 1.0 / (1.0 + ( pow ((distance/4.0), 2.0)));
				dist_from_Ca_atoms[i].push_back(ScalD);
				}


			CaVector1.clear();
			CaVector2.clear();

		}



			int count;
			double sum;
			for (unsigned int i = 0; i < num_atomi; i++) {
				count=0;
				sum=0;
				vector<double>::iterator atom = dist_from_Ca_atoms[i].begin();
				while (atom != dist_from_Ca_atoms[i].end()){

					sum = sum + (*atom);

					atom++;
					count++;
				}
				dist_everage.push_back((sum / count));


			}



				double standDev=-1;


				vector<double>::iterator everage = dist_everage.begin();

				for (unsigned int i = 0; i < num_atomi; i++) {
					sum=0;
					count=0;
					vector<double>::iterator atom = dist_from_Ca_atoms[i].begin();


					while (atom != dist_from_Ca_atoms[i].end()){
						sum += (pow ( ((*atom) - (*everage) ) , 2.0 ));
								atom++;
								count++;

					}

					standDev = sqrt(sum / count);

					ScD.push_back(standDev);

					everage++;


				}


	return ScD;
}


vector<double> StandardDeviation::get_StandarDev_angle_PHI() {

	double num_amino = 0;
	double sum = 0;

	//cout << "modelli"  << original_models.size()  << "aminoacidi" << original_models[0].sizeAmino() << endl;

	if (original_models.size() < 2)
		ERROR("Modelli insufficienti", exception)
	else
		num_amino = original_models[0].sizeAmino();

	angle_PHI.clear();
	vector<vector<double> > misure_PHI(num_amino,
			vector<double>(original_models.size()));
	vector<double> everage(num_amino);

//calcolo angoli
	for (unsigned int i = 0; i < original_models.size(); i++)
		for (unsigned int j = 0; j < num_amino; j++)
			misure_PHI[j][i] = original_models[i].getAmino(j).getPhi();

//calcolo media angoli
	for (unsigned int j = 0; j < num_amino; j++) {
		sum = 0;

		for (unsigned int i = 0; i < original_models.size(); i++) {
			sum = sum + misure_PHI[j][i];
		}

		everage[j] = (sum / (original_models.size()));

	}

	//calcolo stadDev
	double standDev = -1;

	for (unsigned int j = 0; j < num_amino; j++) {
		sum = 0;
		for (unsigned int i = 0; i < original_models.size(); i++) {

			sum += (pow((misure_PHI[j][i] - everage[j]), 2.0));
		}

		standDev = sqrt(sum / original_models.size());
		angle_PHI.push_back(standDev);
	}

	return angle_PHI;

}


vector<double> StandardDeviation::get_StandarDev_angle_PSI() {

	double num_amino = 0;
	double sum = 0;

	//cout << "modelli"  << original_models.size()  << "aminoacidi" << original_models[0].sizeAmino() << endl;

	if (original_models.size() < 2)
		ERROR("Modelli insufficienti", exception)
	else
		num_amino = original_models[0].sizeAmino();

	angle_PSI.clear();
	vector<vector<double> > misure_PSI(num_amino,
			vector<double>(original_models.size()));
	vector<double> everage(num_amino);

//calcolo angoli
	for (unsigned int i = 0; i < original_models.size(); i++)
		for (unsigned int j = 0; j < num_amino; j++)
			misure_PSI[j][i] = original_models[i].getAmino(j).getPsi();

//calcolo media angoli
	for (unsigned int j = 0; j < num_amino; j++) {
		sum = 0;

		for (unsigned int i = 0; i < original_models.size(); i++) {
			sum = sum + misure_PSI[j][i];
		}

		everage[j] = (sum / (original_models.size()));

	}

	//calcolo stadDev
	double standDev = -1;

	for (unsigned int j = 0; j < num_amino; j++) {
		sum = 0;
		for (unsigned int i = 0; i < original_models.size(); i++) {

			sum += (pow((misure_PSI[j][i] - everage[j]), 2.0));
		}

		standDev = sqrt(sum / original_models.size());
		angle_PSI.push_back(standDev);
	}

	return angle_PSI;

}




