
/*
 @file    test.cc
 @author  Riccardo Zanella, riccardozanella89@gmail.com
 @version 1.0
 */


// Includes:
#include <TestProteinModels.h>

using namespace Victor;
using namespace Victor::Mobi;

int main(int argc, char* argv[]) {

	if (argc!=3)
		ERROR("Give me me output path! and the integer to test. (see signature test)", error);


	string output = argv[1];
	int i =atoi(argv[2]);

	//TestMobi* test = new TestMobi;

	//test->TestProteinModels(output,i);



	return 0;
}
