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

/*
 @file    Mobi.cc
 @author  Riccardo Zanella, riccardozanella89@gmail.com
 @version 1.0
 */

/* --*- C++ -*------x-----------------------------------------------------------
 *
 *
 * Description:     Mobi is a software that take in input .pdb file created
 *					Whit NMR technology.
 *	 				this file contained more models of the same protein,
 *					he take it, then he imposed all models whit all models
 *					and return mobility part of the protein in a .pdb file
 *
 *
 * -----------------x-----------------------------------------------------------
 */

// Includes:
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
			<< "\t Mobi is a software that take in input a .pdb file created by NMR tecnology."
			<< endl;
	cout
			<< "\t this file contained more models of the same protein, then Mobi imposed models and return "
			<< endl;
	cout << "\t  aminoacid's mobility." << endl;

	cout << "\n" << endl;

	cout << "\t -v  verbose output. \n" << endl;
	cout << "\t -o [FILE OUTPUT]  Output to file (default stdout)\n" << endl;
	cout << "\t -h  help" << endl;

	exit(EXIT_SUCCESS);

}

int main(int argc, char* argv[]) {

	bool v;
	string inputFile, outputFile, input;

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

	if (v) {
		pl.setVerbose();
		prot.setVerbose();
	} else
		pl.setNoVerbose();

	// --------------------------------------------------
	// 1. ask how much models load in protein
	// --------------------------------------------------

	do {
		cout << "Questo file pdb contiene ";
		cout << pl.getMaxModels() << endl;

		cout << "Vuoi caricare tutti i modelli? [y/n]" << endl;
		getline(cin, input);
		//		cin >> input;

	} while ((strcmp(input.c_str(), "y")) != 0
			&& (strcmp(input.c_str(), "n")) != 0); //(!cin.fail() && input != 'y' && input != 'n');

	if (!v)
		cout << "\nLOAD..." << endl;

	if ((strcmp(input.c_str(), "y")) == 0) {
		prot.load(pl);
	} else {
		prot.loadSameModels(pl);
	}

	// --------------------------------------------------
	// 2. set path to output file
	// --------------------------------------------------

	if (outputFile != "!") {
		outputFile = "./Mobi/data/" + outputFile;
	} else
		outputFile = "./Mobi/data/stdout";

	// --------------------------------------------------
	// 3. method that save in separated file each models
	// --------------------------------------------------

	prot.save(outputFile);

	/////////////////////////////////////////////////////////////////////

	//SUPER IMPOSITION
	unsigned int d = 0;

	// ------------------------------------------------------
	// 4. TmScore object make super imposition of each models
	// ------------------------------------------------------

	TmScore tm("./Mobi/data/TMscore", outputFile, v);

	Protein* traslata = new Protein();

	for (unsigned int i = 0; i < prot.size() - 1; i++)
		for (unsigned int j = i + 1; j < prot.size(); j++) {
			traslata = tm.TmImpose(outputFile + (itosDEF(i)),
					outputFile + (itosDEF(j)));
			if (traslata != NULL) {

				prot.addModels(*(traslata->getSpacer(d)));
			} else
				ERROR("Errore nella creazione della proteina traslata",
						exeption);

			prot.addModels(*(prot.getSpacer(j)));

		}

	// ----------------------------------------------------------------------------------------
	// 5. Standard Deviation object calculate metrics (everage ScalD, StanDevScalD ,phi, psi)
	// ----------------------------------------------------------------------------------------

	StandardDeviation std(prot, v);

	vector<double> ever;

	ever = std.get_everage_distance();

	//print vector
//	cout << "EVERAGE" << endl;
//		for (vector<double>::iterator everage = ever.begin(); everage != ever.end();
//				everage++) {
//			cout << *everage << endl;
//
//		}

	vector<double> SD;
	SD = std.get_standard_deviation();

	//print vector
//	cout << "STANDARD DEVIATION" << endl;
//	for (vector<double>::iterator walk = SD.begin(); walk != SD.end();
//			walk++) {
//		cout << *walk << endl;
//	}

///////////////////////////////////////////////////

	//ANGLE

	vector<double> ANGLE_PHI;
	ANGLE_PHI = std.get_StandarDev_angle_PHI();

//	cout << "ANGLE PHI" << endl;
//	for (vector<double>::iterator walk = ANGLE_PHI.begin(); walk != ANGLE_PHI.end();
//			walk++) {
//		cout << *walk << endl;
//
//	}

	vector<double> ANGLE_PSI;
	ANGLE_PSI = std.get_StandarDev_angle_PSI();

//		cout << "ANGLE PSI" << endl;
//		for (vector<double>::iterator walk = ANGLE_PSI.begin(); walk != ANGLE_PSI.end();
//				walk++) {
//			cout << *walk << endl;
//
//		}

///////////////////////////////////////////////////////////////////////////////////

	//SECONDARY STRUCTURE

	SecondaryStructure sstr(prot, v);

	vector<char> MOB;
	MOB = sstr.getMobilitySecondaryStructure();

	//print vector
//	cout << "\nSECONDARY STRUCTURE" << endl;
//	for (vector<char>::iterator walk = MOB.begin(); walk != MOB.end(); walk++) {
//		cout << *walk << " ";
//
//	}
//
//	cout << endl;

/////////////////////////////////////////////////////////////////////////////////////

	// --------------------------------------------------
	// 6. MobiSaver object save mobility in output file
	// --------------------------------------------------

	MobiSaver* saver = new MobiSaver(prot, outputFile, v);
	saver->allMobility(ever, SD, ANGLE_PHI, ANGLE_PSI, MOB);

	delete saver;

	prot.remove(outputFile);

	return 0;
}
