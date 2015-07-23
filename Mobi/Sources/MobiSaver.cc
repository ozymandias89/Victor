/*
 * MobiSaver.cc
 *
 *  Created on: 22 lug 2015
 *      Author: riccardo
 */

#include "MobiSaver.h"



using namespace Victor;
using namespace Victor::Biopool;
using namespace Victor::Mobi;

const string OUT = "ResultMobi.fasta";

MobiSaver::MobiSaver(ProteinModels prot, string output, bool verbose): Saver(), verbose(verbose) {
	// TODO Auto-generated constructor stub

	//Se il file non esiste lo creo inserendo la sequenza aminoacidica in testa altrimenti non creo nulla
	out = output + OUT;

	if (access(out.c_str(), F_OK) != 0){

		if (verbose)
				cout<<"\n ### Create fasta file in output ###"<<endl;

		ofstream fout;			//stream in output
		vector<char>::iterator walk = prot.sequence.begin();

		fout.open(out.c_str());
				if (!fout) {
							ERROR("Error to create file!", error);
						} else {
							fout << "> sequence" << endl;
							while (walk != prot.sequence.end()) {
								fout << *(walk);
								walk++;
							}


							fout.close();

							if (verbose)
								cout << "file success created!" << endl;
						}

	}


}

MobiSaver::~MobiSaver() {
	// TODO Auto-generated destructor stub
}



void MobiSaver::save_mob_everageScalD(vector <double> everageDistance){

	if (everageDistance.size()==0)
		ERROR ("Vector everage distance empty!", exception);

	ofstream fout (out.c_str(), ofstream::app);

	if (!fout)
		ERROR("Could not open file for writing.", error);

	if (verbose)
		cout<<"\n ### Salvataggio sul file della mobilità data la distanza media ###"<<endl;

	fout << endl;
	fout << "> Average scale distance 0" << endl;

	for (vector<double>::iterator walk = everageDistance.begin();
			walk != everageDistance.end(); walk++) {
		if (*(walk) < 0.85)
			fout << "M";
		else
			fout << ".";
	}

	fout.close();

	if (verbose)
		cout << "Mobility from average scaled distance saved!" << endl;



}


void MobiSaver::save_mob_standardDeviation(vector<double> Scale_distance) {

	if (Scale_distance.size() == 0)
		ERROR("Vector standard deviation empty!", exception);

	ofstream fout(out.c_str(), ofstream::app);

	if (!fout)
		ERROR("Could not open file for writing.", error);

	if (verbose)
		cout
				<< "\n ### Salvataggio sul file della mobilità data dalla standard deviation ###"
				<< endl;

	fout << endl;
	fout << "> Standard deviation distance 1" << endl;

	for (vector<double>::iterator walk = Scale_distance.begin();
			walk != Scale_distance.end(); walk++) {
		if (*(walk) > 0.09)
			fout << "M";
		else
			fout << ".";
	}

	fout.close();

	if (verbose)
		cout << "Mobility from standard deviation saved!" << endl;

}


void MobiSaver::save_mob_angle_PHI(vector<double> angle_PHI){

	if (angle_PHI.size() == 0)
			ERROR("Vector angle PHI empty!", exception);

		ofstream fout(out.c_str(), ofstream::app);

		if (!fout)
			ERROR("Could not open file for writing.", error);

		if (verbose)
			cout
					<< "\n ### Salvataggio sul file della mobilità data dall'angle_PHI ###"
					<< endl;

		fout << endl;
		fout << "> Phi angle 2" << endl;

		for (vector<double>::iterator walk = angle_PHI.begin();
				walk != angle_PHI.end(); walk++) {
			if (*(walk) > 20)
				fout << "M";
			else
				fout << ".";
		}

		fout.close();

		if (verbose)
			cout << "Mobility from angle_PHI saved!" << endl;

}

void MobiSaver::save_mob_angle_PSI(vector<double> angle_PSI){

	if (angle_PSI.size() == 0)
			ERROR("Vector angle PSI empty!", exception);

		ofstream fout(out.c_str(), ofstream::app);

		if (!fout)
			ERROR("Could not open file for writing.", error);

		if (verbose)
			cout
					<< "\n ### Salvataggio sul file della mobilità data dall'angle_PSI ###"
					<< endl;

		fout << endl;
		fout << "> Psi angle 3" << endl;

		for (vector<double>::iterator walk = angle_PSI.begin();
				walk != angle_PSI.end(); walk++) {
			if (*(walk) > 20)
				fout << "M";
			else
				fout << ".";
		}

		fout.close();

		if (verbose)
			cout << "Mobility from angle_PSI saved!" << endl;

}

void MobiSaver::save_mob_SecondaryStructure(vector<char> Mob_SecStructure){

	if (Mob_SecStructure.size() == 0)
			ERROR("Vector Mob_SecStructure empty!", exception);

		ofstream fout(out.c_str(), ofstream::app);

		if (!fout)
			ERROR("Could not open file for writing.", error);

		if (verbose)
			cout
					<< "\n ### Salvataggio sul file della mobilità data la struttura secondaria ###"
					<< endl;

		fout << endl;
		fout << "> Secondary structure 4" << endl;

		for (vector<char>::iterator walk = Mob_SecStructure.begin();
				walk != Mob_SecStructure.end(); walk++)
				fout << *(walk);

		fout.close();

		if (verbose)
			cout << "Mobility from sacondary structure saved!" << endl;

}

void MobiSaver::save_allMobility(vector<double> everageDistance,
		vector<double> Scale_distance, vector<double> angle_PHI,
		vector<double> angle_PSI, vector<char> Mob_SecStructure) {

	save_mob_everageScalD(everageDistance);
	save_mob_standardDeviation(Scale_distance);
	save_mob_angle_PHI(angle_PHI);
	save_mob_angle_PSI(angle_PSI);
	save_mob_SecondaryStructure(Mob_SecStructure);

}
