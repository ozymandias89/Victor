/*
 @file    TestMobi.cc
 @author  Riccardo Zanella, riccardozanella89@gmail.com
 @version 1.0
 */

/* This file is part of Victor.
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

#include <iostream>
#include <cppunit/TestSuite.h>
#include <cppunit/ui/text/TestRunner.h>

#include <TestProteinModels.h>

using namespace Victor::Mobi;
using namespace Victor::Biopool;
using namespace std;

int main() {
	CppUnit::TextUi::TestRunner runner;

	cout << "Creating Test Suites:" << endl;
	runner.addTest(TestProteinModels::suite());

//	runner.addTest(TestMobiProtein::suite());
//	runner.addTest(TestTM::suite());
//	runner.addTest(TestMobiMethods::suite());

	cout << "Running the unit tests. " << endl;
	runner.run();
	return 0;
}




/**
 test object protein models. if i == 0 (load and save)
 if i == 1 (load same models and save)
 if i == 2 (load same models and save and impose and save vector of models)
 @param (string output) , name of output file
 @param int i , select type of test
 @return void
 */
/*
void TestMobi::TestProteinModels(string output, int i) {

	// --------------------------------------------------
	// 1. Initialization
	// --------------------------------------------------
	ifstream inputFile(output.c_str());
	ProteinModels* protein_models = new ProteinModels();
	protein_models->setVerbose();
	PdbLoader pl(inputFile);
	pl.setVerbose();

	if (i == 0) {
		// --------------------------------------------------
		// 2. Load
		// --------------------------------------------------
		protein_models->load(pl);

		// --------------------------------------------------
		// 3. Save
		// --------------------------------------------------
		protein_models->save(output);
	}
	if (i == 1) {
		// --------------------------------------------------
		// 2. Load
		// --------------------------------------------------
		protein_models->loadSameModels(pl);

		// --------------------------------------------------
		// 3. Save
		// --------------------------------------------------
		protein_models->save(output);

	}


	if (i == 2) {
		// --------------------------------------------------
		// 2. Load
		// --------------------------------------------------
		protein_models->loadSameModels(pl);

		// --------------------------------------------------
		// 3. Save
		// --------------------------------------------------
		protein_models->save(output);


		//SUPER IMPOSITION
			unsigned int d = 0;

			// ------------------------------------------------------
			// 4. TmScore object make super imposition of each models
			// ------------------------------------------------------

			TmScore tm("./Mobi/data/TMscore", output, true);

			Protein* traslata = new Protein();

			for (unsigned int i = 0; i < protein_models->size() - 1; i++)
				for (unsigned int j = i + 1; j < protein_models->size(); j++) {
					traslata = tm.TmImpose(output + (itosDEF(i)),
							output + (itosDEF(j)));
					if (traslata != NULL) {

						protein_models->addModels(*(traslata->getSpacer(d)));
					} else
						ERROR("Error in the creation of shift models.",
								exeption);

					protein_models->addModels(*(protein_models->getSpacer(j)));

				}

			// ------------------------------------------------------
			// 5. save model vector that contain couple of models,
			//    the first rototraslate, the second is the model
			// ------------------------------------------------------

			protein_models->printModels(output);

	}

}
*/
