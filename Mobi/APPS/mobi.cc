/*  This file is part of Victor.

 Victor is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Victor is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with Victor.  If not, see <http://www.gnu.org/licenses/>.
 */

// --*- C++ -*------x-----------------------------------------------------------
//
//
// Description:     Mobi è un software che riceve in input
//					un file .pdb creato con tecnologia NMR.
//	 				Questo file contiene molti modelli della stessa proteina,
//					li confronta e restituisce in un file di output
//	 				le parti mobili della proteina.
//
//
// -----------------x-----------------------------------------------------------
#include <ProteinModels.h>
#include <MobiSaver.h>
#include <SecondaryStructure.h>
#include <TmScore.h>
#include <string>
#include <GetArg.h>
#include <PdbSaver.h>
#include <StandardDeviation.h>

using namespace Victor;
using namespace Victor::Biopool;
using namespace Victor::Mobi;

void sShowHelp() {

	cout << "NAME" << endl;
	cout << "\t mobi" << endl;

	cout << "SYNOPSIS" << endl;
	cout << "\t mobi [FILE]... [OPTION]..." << endl;

	cout << "DESCRIPTION \n" << endl;
	cout
			<< "\t Mobi è un software che riceve in input un file .pdb creato con tecnologia NMR."
			<< endl;
	cout
			<< "\t Questo file contiene molti modelli della stessa proteina, li confronta e restituisce in un file di output "
			<< endl;
	cout << "\t le parti mobili della proteina." << endl;

	cout << "\n" << endl;

	cout << "\t -v  verbose output. \n" << endl;
	cout << "\t -o [FILE OUTPUT]  Output to file (default stdout)\n" << endl;
	cout << "\t -h  help" << endl;

	exit(EXIT_SUCCESS);

}

int main(int argc, char* argv[]) {

	bool v;
	string inputFile, outputFile, outputFile_1;
	string input;

	//guide with -h option
	if (getArg("h", argc, argv)) {
		sShowHelp();
		return 1;
	}
	if (argc == 1) {
		sShowHelp();
		return 1;
	}

	cout << "Welcome to Mobi!" << endl;

	//loan pdb file from input
	inputFile = argv[1];

	ifstream inFile(inputFile.c_str());

	if (!inFile)
		ERROR("Error opening input .pdb file.", exception);

	getArg("o", outputFile, argc, argv, "!");
	v = getArg("v", argc, argv);

	PdbLoader pl(inFile);    // creates the PdbLoader object

	ProteinModels prot;


	// Protein prot;

	if (v) {
		pl.setVerbose();
		prot.setVerbose();
	} else
		pl.setNoVerbose();

	//method that ask how much models load in protein

	do {
		cout << "Questo file pdb contiene ";
		cout << pl.getMaxModels() << endl;

		cout << "Vuoi caricare tutti i modelli? [y/n]" << endl;
		getline(cin, input);
		//		cin >> input;

	} while ((strcmp(input.c_str(), "y")) != 0
			&& (strcmp(input.c_str(), "n")) != 0); //(!cin.fail() && input != 'y' && input != 'n');

	if ((strcmp(input.c_str(), "y")) == 0) {
		prot.load(pl);
	} else {
		prot.loadSameModels(pl);
	}

	if (outputFile != "!") {
		outputFile = "./Mobi/data/" + outputFile;
	} else
		outputFile = "./Mobi/data/stdout";

	// Mi salvo i vari modelli
	prot.save(outputFile);

	/////////////////////////////////////////////////////////////////////

//	unsigned int d = 0;
//
//	TmScore tm("./Mobi/data/TMscore", outputFile, v);
//
//	Protein* traslata = new Protein();
//
//	for (unsigned int i = 0; i < prot.size()-1; i++)
//		for (unsigned int j = i+1; j < prot.size(); j++)
//			{
//				traslata = tm.TmImpose(outputFile + (itosDEF(i)),
//						outputFile + (itosDEF(j)));
//				if (traslata != NULL) {
//
//					prot.addModels(*(traslata->getSpacer(d)));
//				} else
//					ERROR("Errore nella creazione della proteina traslata",
//							exeption);
//
//				prot.addModels(*(prot.getSpacer(j)));
//
//
//			}
//
//  //prot.printModels(outputFile);
//
//	StandardDeviation std(prot, v);
//	vector<double> SD;
//
//	vector<double> EVERAGE;
//
//	EVERAGE = std.get_everage_distance();
//
//	cout << "EVERAGE" << endl;
//		for (vector<double>::iterator everage = EVERAGE.begin(); everage != EVERAGE.end();
//				everage++) {
//			cout << *everage << endl;
//
//		}
//
//
//	SD = std.get_standard_deviation();
//	cout << "STANDARD DEVIATION" << endl;
//	for (vector<double>::iterator walk = SD.begin(); walk != SD.end();
//			walk++) {
//		cout << *walk << endl;
//	}


///////////////////////////////////////////////////




//	StandardDeviation std(prot, v);
//
//	vector<double> ANGLE_PHI;
//	ANGLE_PHI=std.get_StandarDev_angle_PHI();
//
//
//	cout << "STANDARD DEVIATION" << endl;
//	for (vector<double>::iterator walk = ANGLE_PHI.begin(); walk != ANGLE_PHI.end();
//			walk++) {
//		cout << *walk << endl;
//
//	}
//
//
//	vector<double> ANGLE_PSI;
//		ANGLE_PSI=std.get_StandarDev_angle_PSI();
//
//
//		cout << "STANDARD DEVIATION" << endl;
//		for (vector<double>::iterator walk = ANGLE_PSI.begin(); walk != ANGLE_PSI.end();
//				walk++) {
//			cout << *walk << endl;
//
//		}


///////////////////////////////////////////////////////////////////////////////////

//	SecondaryStructure sstr (prot, v);
//
//	vector<char> MOB;
//	MOB = sstr.getMobilitySecondaryStructure();
//
//	cout << "\nSECONDARY STRUCTURE" << endl;
//	for (vector<char>::iterator walk = MOB.begin(); walk != MOB.end(); walk++) {
//		cout << *walk << " ";
//
//	}
//
//	cout << endl;


/////////////////////////////////////////////////////////////////////////////////////

    MobiSaver* saver = new MobiSaver(prot, outputFile, v);
    saver->save_mob_everageScalD();
    delete saver;

	prot.remove(outputFile);

	return 0;
}
