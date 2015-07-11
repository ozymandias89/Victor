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
 #include <TmScore.h>
 #include <string>
 #include <GetArg.h>
 #include <PdbSaver.h>

 

  using namespace Victor;
  using namespace Victor::Biopool;
  using namespace Victor::Mobi;



 void sShowHelp(){

	 cout << "NAME" << endl;
	 cout << "\t mobi" << endl;

	 cout << "SYNOPSIS" << endl;
	 cout << "\t mobi [FILE]... [OPTION]..." << endl;

	 cout << "DESCRIPTION \n" << endl;
	 cout << "\t Mobi è un software che riceve in input un file .pdb creato con tecnologia NMR." << endl;
	 cout << "\t Questo file contiene molti modelli della stessa proteina, li confronta e restituisce in un file di output " << endl;
	 cout << "\t le parti mobili della proteina." << endl;

	 cout << "\n" <<endl;

	 cout << "\t -v  verbose output. \n" <<endl;
	 cout << "\t -o [FILE OUTPUT]  Output to file (default stdout)\n" <<endl;
	 cout << "\t -h  help" << endl;

	 exit (EXIT_SUCCESS);


 }




  int main( int argc, char* argv[] ) {

	  bool v;
	  string inputFile, outputFile, outputFile_1;
	  string input;



	  //guide with -h option
	  if (getArg("h", argc, argv)) {
	          sShowHelp();
	          return 1;
	      }
	 if (argc == 1){
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

     if (v){
    	 pl.setVerbose();
    	 prot.setVerbose();
     }else
    	 pl.setNoVerbose();

     //method that ask how much models load in protein

     do {
		cout << "Questo file pdb contiene ";
		cout << pl.getMaxModels() << endl;

		cout << "Vuoi caricare tutti i modelli? [y/n]" << endl;
		getline(cin, input);
		//		cin >> input;

	} while ((strcmp(input.c_str(), "y"))!=0 && (strcmp(input.c_str(), "n"))!=0);//(!cin.fail() && input != 'y' && input != 'n');


     if ((strcmp(input.c_str(), "y"))==0){
    	// cin.ignore();
//    	 cout << "1" << endl;
    	 prot.load(pl);
//    	 cout << "2" << endl;
     }
     else{
    	 //cin.ignore();
//    	 cout << "1111" << endl;
    	 prot.loadSameModels(pl);
//    	 cout << "22" << endl;
     }

     if (outputFile != "!"){
    	outputFile="./Mobi/data/" + outputFile;
  	  }	else outputFile="./Mobi/data/stdout";


     // Mi salvo i vari modelli
       prot.save(prot, outputFile);



       /////////////////////////////////////////////////////////////////////

       TmScore tm("./Mobi/data/TMscore", outputFile, v);


       Protein* sp = new Protein();

//for (unsigned int i=0; i<prot.size(); i++)
//for (unsigned int j=0; j<prot.size(); j++)
//if (i!=j) {
//	cout << outputFile + (itosDEF(i)) << outputFile + (itosDEF(j)) << endl;

       sp=tm.TmImpose(outputFile + (itosDEF(0)), outputFile + (itosDEF(1)));



	if (sp != NULL) {

		ofstream pippo;			//stream in output

		unsigned int d = 0;
		PdbSaver ps(pippo);
		pippo.open((outputFile + "111").c_str());
		ps.saveSpacer(*(sp->getSpacer(d)));
		ps.endFile();
		pippo.close();
	}
    			   //    			   tm.TmImpose(outputFile + (itosDEF(0)), outputFile + (itosDEF(1)));
//    			   sp=*(tm.spacerFromTMOutput(tm.outputTmScore + "_atm"));
//    			   prot.models.push_back(sp);
//    			   cout << prot.models.size()<< endl;
//    		//}

       return 0;
  }
