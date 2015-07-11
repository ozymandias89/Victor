/* This file is part of Victor
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
/**
 * @file TmScore.h
 * @author Riccardo Zanella
 * @date Lug 2015
 * @version 0.1
 */

#ifndef MOBI_SOURCES_TMSCOREBIN_H_
#define MOBI_SOURCES_TMSCOREBIN_H_

#include <iostream>
#include <string>
#include <ProteinModels.h>


using namespace Victor::Biopool;
using namespace Victor::Mobi;
using namespace std;


extern const string TMTMP_OUT;


namespace Victor {

namespace Mobi {

/**
 * @brief TMScore functionalities through external binary.
 */

class TmScore {

public:
	/**
	 * @brief Setup.
	 * @param _binary (string) full path to binary TMScore file, must have execution permission
	 * @param _tmp (string) full path to temp dir, must have write permission
	 */
	TmScore(string _binary = "TMScore", string output = ".", bool _verbose = false);
	/**
	 * @brief Given two pdb files containing each a model of the same protein,
	 * call TMScore binary to superimpose the first over the second. The superimposed (rotated/traslated)
	 * model is then loaded in a ProteinModels using the double pointer provided.
	 * @param modelFile (string) full path to model#1 file
	 * @param nativeFile (string) full path to model#2 file
	 * @param imposedModel(ProteinModel**) double pointer of type ProteinModel, as output
	 */

	Protein* TmImpose(string modelFile, string nativeFile);//, ProteinModels** imposedModel);
	/**
	 * @brief Given a ProteinModel call TMScore binary to superimpose two models contained in it.
	 * The superimposed (rotated/traslated) model is then loaded in a ProteinModel using the double pointer provided.
	 * @param prot(ProteinModel&) reference to ProteinModel object
	 * @param model (unsigned int) model#1 name in ProteinModel object
	 * @param native (unsigned int) model#2 name in ProteinModel object
	 * @param imposedModel(ProteinModel**) double pointer of type ProteinModel, as output
	 */
//	virtual void TMImpose(ProteinModel& prot, unsigned int model,
//			unsigned int native, ProteinModel** imposedModel);
//	/**
//	 * @brief the TMScore output is not pdb conformant. This static method read the output and fix it in a memory buffer.
//	 * Then loads a ProteinModel with the superimposed model only.
//	 * @param pdbFile (string) full path to TMScore output
//	 * @param imposedModel (ProteinModel**) double pointer of type PRoteinModel, as output
//	 */
	//Spacer* spacerFromTMOutput(string pdbFile);
	/**
//	 * Set verbosity.
//	 * Verbosity is applied on PdbLoader instances.
//	 * @param v (bool) new verbosity option
//	 */
	void setVerbose(bool v) {
		verbose = v;
	}
	/**
	 * Default deconstructor
	 */
	~TmScore();

private:
	string binary;
	const string outputTmScore;
	bool verbose;
};
}
}
#endif /* MOBI_SOURCES_TMSCOREBIN_H_ */
