/*
 @file    TestProteinModels.h
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

#ifndef MOBI_TESTS_TESTPROTEINMODELS_H_
#define MOBI_TESTS_TESTPROTEINMODELS_H_

#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include "ProteinModels.h"

using namespace std;
using namespace Victor::Mobi;


class TestProteinModels : public CppUnit::TestFixture {
public:
	TestProteinModels(){}
	virtual ~TestProteinModels(){}


	 static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite(
				"TestProteinModels");

		suiteOfTests->addTest(
				new CppUnit::TestCaller<TestProteinModels>(
						"Test1 - Populate the ProteinModels",
						&TestProteinModels::testPopulation));

		return suiteOfTests;
	}

	 /** @brief Setup method. */
	void setUp() {
	}
	/** @brief Teardown method. */
	void tearDown() {
	}

protected:
	 /** @brief Test for parse/write. */
	void testPopulation() {

		cout << "Prima prova" << endl;
		TestProteinModels* test = new TestProteinModels();
		delete test;
	}

	//void TestProteinModels(string __output, int __i);

};

#endif /* MOBI_TESTS_TESTPROTEINMODELS_H_ */
