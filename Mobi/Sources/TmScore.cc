/* 	This file is part of Victor

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

/**
 * @file TmScore.cc
 * @author Riccardo Zanella
 * @date Lug 2015
 * @version 0.1
 */
#include <TmScore.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <Debug.h>
#include <stdio.h>
#include <errno.h>

using namespace Victor::Mobi;
using namespace Victor::Biopool;
using namespace std;

const string TMTMP_OUT = "TMScore.pdb";


TmScore::TmScore(string _binary, string output, bool verbose) :
		binary(_binary), outputTmScore(output + TMTMP_OUT), verbose(verbose) {}

//TmScore::TmScore(string _binary, string output, bool verbose) :
//		binary(_binary), tmp(output.substr(output.length() - 1, 1) == "/" ?
//						output : output + "/"), verbose(verbose) {
//}

TmScore::~TmScore(){}



Protein* TmScore::TmImpose(string modelFile, string nativeFile) {

	Protein* prot = NULL;
	//se i file sono accessibili in lettura
	if (access(modelFile.c_str(), R_OK) == 0
			&& access(nativeFile.c_str(), R_OK) == 0) {

		//se l'eseguibile ha i permessi di esecuzione
		if (access(binary.c_str(), X_OK) == 0) {
			pid_t pid;
			if (this->verbose)
				cout << "Accesso consentito." << endl;
			//Pipeing TMS output
			int pipefd[2];
			pipe(pipefd);
			//Forking child for TMS
			pid = fork();
			if (pid < 0) {
				ERROR("Unable to fork child process", error);
			} else {
				if (pid == 0) {
					close(pipefd[0]);						//close input pipe
					if (dup2(pipefd[1], 1) < 0)
						ERROR("Unable to redirect output TMScore", error);//stdout
					if (dup2(pipefd[1], 2) < 0)
						ERROR("Unable to redirect error output TMScore", error);//stderr
					close(pipefd[1]);						//close output pipe
					if (this->verbose)
						cout << "Eseguo TMScore" << endl;
					if (execl(binary.c_str(), binary.c_str(), modelFile.c_str(),
							nativeFile.c_str(), "-o", outputTmScore.c_str(),
							NULL))
						ERROR("Unable to exec", error);
				} else {
					//Read output from child pipe

					//////////////////////////////////////////////////////////////////////////
					int status;

//					ifstream inFile(nativeFile.c_str());
//					ofstream out ((outputTmScore+"11111").c_str());
//					stringstream buffer;
//					string line;
//
//					if (verbose) {
//						cout << "Conversione file native file per renderlo leggibile" << endl;
//					}
//					bool flag=false;
//
//					while (inFile) {
//						line = readLine(inFile);
//						if (line.substr(0, 6) == "ATOM  " && !flag) {
//							buffer << "MODEL        1" << endl;
//							out << "MODEL        1" << endl;
//							buffer << line << endl;
//							out << line << endl;
//							flag=true;
//							} else if (line.substr(0, 6) == "ATOM  " && flag) {
//							buffer << line << endl;
//							out << line << endl;
//						} else if (line.substr(0, 3) == "TER") {
//							buffer << line << endl;
//							buffer << "ENDMDL" << endl;
//							out << line << endl;
//							out << "ENDMDL" << endl;
//
//							break;
//						}
//					}
//					if (verbose) {
//						cout << "Fine lettura file nativo" << endl;
//					}
//
//					buffer.clear();
//					out.close();
//
//					cout << "AAAAAAAAAAAAAAAAAAAAAAAA" << endl;
//
//					PdbLoader loader (buffer);
//
//				    if (!verbose)
//				    	loader.setNoVerbose();
//
//				    loader.setModel(1);
//				    cout << "BBBBBBBBBBBBBBBBBBBBBBBBBB" << endl;
//				    loader.checkModel();
//				    cout << "CCCCCCCCCCCCCCCCCCCCCCCCCCCCc" << endl;
//					Protein temp;
//					temp.load(loader);
//					unsigned int d=0;
//					sp[0]=*(temp.getSpacer(d));

					wait(&status);

					remove(outputTmScore.c_str());

					if (access((outputTmScore + "_atm").c_str(), R_OK) != 0)
						ERROR("Cannot read pdb file to fix", exception);

					ifstream inFile((outputTmScore + "_atm").c_str());
					stringstream buffer;
					string line;

					if (verbose) {
						cout
								<< "Conversione file rototraslato per renderlo leggibile"
								<< endl;
					}
					while (inFile) {
						line = readLine(inFile);
						if (line.substr(0, 16) == "REMARK  TM-score") {
							buffer << line << endl;
							buffer << "MODEL        1" << endl;

						} else if (line.substr(0, 6) == "ATOM  ") {
							buffer << line << endl;
						} else if (line.substr(0, 3) == "TER") {
							buffer << line << endl;
							buffer << "ENDMDL" << endl;
							break;
						}
					}
					if (verbose) {
						cout << "Fine lettura file rototraslato" << endl;
					}

					buffer.clear();

					//Load from memory buffer
					PdbLoader pl(buffer);

					if (!verbose)
						pl.setNoVerbose();

					pl.setModel(1);
					pl.checkModel();

					prot = new Protein();

					if (verbose)
						cout << "### Caricamento proteina rototraslata ###" <<endl;

					pl.loadProtein(*prot);




				}
			}
		} else
			ERROR("No access to " + binary + "  binary!", exception);
	} else
		ERROR("No access to pdb files " + modelFile + " or " + nativeFile,
				exception);
	return prot;
}

//void TmScore::TMImpose(ProteinModels& prot, unsigned int model, unsigned int native, ProteinModels** imposedModel) {
//
//	stringstream sstm;
//	sstm << "TMScore between models " << model << " and " << native;
//	DEBUG_MSG(sstm.str());
//
//	//Save spacers in pdb files
//	std::ofstream fout;
//	PdbSaver ps(fout);
//	fout.open((tmp + TMTMP_IN1).c_str());
//	ps.saveSpacer(prot.getModel(model));
//	ps.endFile();
//	fout.close();
//	fout.open((tmp + TMTMP_IN2).c_str());
//	ps.saveSpacer(prot.getModel(native));
//	ps.endFile();
//	fout.close();
//
//	//Call TMScore binary
//	return TMImpose((tmp + TMTMP_IN1), (tmp + TMTMP_IN2), imposedModel);
//}
//
//Spacer* TmScore::spacerFromTMOutput(string pdbFile) {
//	if (access(pdbFile.c_str(), R_OK) != 0)
//		ERROR("Cannot read pdb file to fix", exception);
//
//	ifstream inFile(pdbFile.c_str());
//	stringstream buffer;
//	string line;
//	while (inFile) {
//		line = readLine(inFile);
//		if (line.substr(0, 16) == "REMARK  TM-score") {
//			buffer << line << endl;
//			if (line.substr(6, 8) == "TM-score")
//				buffer << "MODEL        1" << endl;
//		} else if (line.substr(0, 6) == "ATOM  ") {
//			buffer << line << endl;
//		} else if (line.substr(0, 3) == "TER") {
//			buffer << line << endl;
//			buffer << "ENDMDL" << endl;
//			break;
//		}
//	}
//	buffer.clear();
//	buffer.seekg(0);
//	//Load from memory buffer
//	PdbLoader pl(buffer);
//	if (!verbose)
//		pl.setNoVerbose();
//	//Load and return protein object
//	//*imposedModel = new ProteinModels();
//	pl.setModel(1);
//	pl.checkModel();
//	Protein prot;
//	pl.loadProtein(prot);
//	unsigned int d=0;
//	return prot.getSpacer(d);
//}

//char buffer[2048];
//					bool flag=true;
//					//la read implementa un sync al suo interno, questo mi permette di aspettare l'output del figlio senza pensieri
//					while(flag){
//					int a=read(pipefd[0], buffer, sizeof(buffer));
//					if ( a> 0) {
//						if (this->verbose) cout << "Accesso alla read." <<endl;
//						char * tok;
//						tok = strtok(buffer, "\n\r");			//spezzetto la stringa in token
//						while (tok != NULL) {
//							//cerco il valore del tm score
//							if (string(tok).find("TM-score") == 0)
//								if (string(tok).find("=") > 0){
//									//stringa valore TM-score
//									score_tmp=string(tok).substr(string(tok).find("=") + 2,6);
//									//converto in double il valore del TM-score
//									score = std::strtod(score_tmp.c_str(), NULL);
//									flag=false;
//								}
//							//termino ricerca
//							tok = strtok(NULL, "\n\r");
//						}
//					}
//					if (a == 0) {
//						ERROR ("Connessione chiusa", exeption);
//					}
//					if (a == -1) {
//						cout << "Errore nella read" << endl;
//						strerror(errno);
//					}
////					if (score < 0)
////						ERROR("TM-score failed or non recognised output",exception)
////						else if (this->verbose) cout << "valore dello score: " << score <<endl;
//				}
//					if (this->verbose) cout << "valore dello score: " << score <<endl;
