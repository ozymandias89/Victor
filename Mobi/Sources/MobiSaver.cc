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



void MobiSaver::save_mob_everageScalD (){

	ofstream fout;			//stream in output


	if (verbose)
		cout<<"\n ###Salvataggio sul file della mobilità data la distanza media ###"<<endl;






}
