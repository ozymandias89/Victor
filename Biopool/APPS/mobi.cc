 #include <PdbLoader.h>
 #include <Protein.h>
 #include <iostream>

 

  using namespace Victor::Biopool;
  using namespace Victor; 

 

  int main( int argc, char* argv[] ) {

 

    // string inputFile = "/home/riccardo/mobi/1AB2_input.pdb";
    string inputFile = "/home/riccardo/Victor/Biopool/Tests/data/3DFR.pdb";


     ifstream inFile( inputFile.c_str() );

     PdbLoader pl(inFile);    // creates the PdbLoader object

 

     Protein prot;            

     prot.load( pl );         // creates the Protein object

	cout << "!!!!!!!!!!!!!!!!!" << endl;
        
	cout << prot.getChainLetter (0) << endl;
	

  }
