#include <iostream>
#include <ProteinModels.h>
#include <String2Number.h>
#include <PdbSaver.h>

using namespace Victor;
using namespace Victor::Mobi;
using namespace Victor::Biopool;

const string PDB = ".pdb";

ProteinModels::ProteinModels() : Protein(), verbose(false) {}

ProteinModels::ProteinModels(const Protein& _orig) : Protein(_orig), verbose(false) {}

ProteinModels::~ProteinModels() {}

// PREDICATES:

/**
 *   Return a pointer to the Spacer with the requested chainID
 * @param c (char), the chainID
 * @return A pointer to the Spacer
 */
void ProteinModels::load(PdbLoader& pl) {
	if (!verbose)
		pl.setNoVerbose();
//Load Models
	for (unsigned int i = 1; i <= pl.getMaxModels(); i++) {
		if (verbose)
			cout << ">>>model#" << i << endl;
		pl.setModel(i);
		pl.checkModel();
		this->Protein::load(pl);
	}
	if (verbose)
	 		 cout << "\n" << "Modelli caricati nella proteina" << endl;

}


/**
 *   Return a pointer to the Spacer with the requested chainID
 * @param c (char), the chainID
 * @return A pointer to the Spacer
 */
unsigned int ProteinModels::selectModels(PdbLoader& pl){

	 unsigned int modelNum;
	 string input;


	 cout << "Quanti ne vuoi caricare? Inserisci un valore compreso tra 2 e ";
	 cout << pl.getMaxModels() << endl;


	 getline(cin, input);


	 modelNum = stouiDEF(input);

	 while (modelNum > pl.getMaxModels() || modelNum < 2)
		{
		 cerr << "Number of models selected out of bounds! Reinput please..." << endl;
		 getline(cin, input);
		 modelNum = stouiDEF(input);

		}

	 return modelNum;

 }


/**
 *   Return a pointer to the Spacer with the requested chainID
 * @param c (char), the chainID
 * @return A pointer to the Spacer
 */
void ProteinModels::loadSameModels(PdbLoader& pl){

		 unsigned int modelNum = selectModels(pl);

 		 for (unsigned int i=1; i<= modelNum; i++)
 		     {
 			 	 if (verbose)
 			 		 cout << "\t>>>model#" << i << endl;
 		    	 pl.setModel(i);		  //dal pdb loader scelgo il modello da caricare nella proteina
 		    	 pl.checkModel();
 		    	 this->Protein::load(pl);          // creates the Protein object
 		     }
 		 if (verbose)
 		 cout << "\n" << "Modelli caricati nella proteina" << endl;

 }

/*procedura che stampa su file il vettore di spacer*/
void ProteinModels::printModels(string outputFile){

	ofstream fout;			//stream in output
	string outputFile_1;
	PdbSaver ps(fout);
	const string VECTOR = "Vector";

	unsigned int i = 0;

	vector<Spacer>::iterator walk = models.begin();
	while (walk != models.end()) {

		outputFile_1 = outputFile + VECTOR + (itosDEF(i)) + PDB;
		if (verbose)
			cout << outputFile_1 << endl;

		fout.open(outputFile_1.c_str());
		if (!fout) {
			ERROR("Could not open file for writing.", exception);
		} else {

			ps.saveSpacer(*walk);
			walk++;
			i++;
			ps.endFile();
			fout.close();

			if (verbose)
				cout << "file chiuso" << endl;
		}
	}

}

void ProteinModels::save(string outputFile){

	ofstream fout;			//stream in output
	string outputFile_1;
	PdbSaver ps(fout);

	if (verbose)
				cout << "Salvataggio dei modelli su file" << endl;
	// Open the proper output stream (file or stdout)
	for (unsigned int i = 0; i < this->size(); i++) {


		outputFile_1 = outputFile + (itosDEF(i)) + PDB;

		if (verbose)
			cout << outputFile_1 << endl;

		fout.open(outputFile_1.c_str());

		if (!fout) {
			ERROR("Could not open file for writing.", exception);
		} else
			ps.saveSpacer(*(this->getSpacer(i)));

		ps.endFile();
		fout.close();

		if (verbose)
			cout << "file chiuso" << endl;
	}

}

void ProteinModels::addModels(Spacer& sp){
	this->models.push_back(sp);

}


vector<Spacer> ProteinModels::getModels() {
	return this->models;
}
