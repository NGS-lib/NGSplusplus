#include <gtest/gtest.h>
#include <string>
#include <sstream>
//TODO : Split tests by type of file??
/**< To test private member directly, changed private for public */
//#define private public
#include "IO/Parser/uParser.h" 
//#include "uParser.h"

using namespace std;
using namespace NGS;

/**< Fixtures */
class CustomConstructorTests_ValidList: public ::testing::Test {
public:
	CustomConstructorTests_ValidList() {
		m_fieldsList.push_back("CHR");
		m_fieldsList.push_back("STRAND");
		m_fieldsList.push_back("START_POS");
		m_fieldsList.push_back("END_POS");
		m_fieldsList.push_back("NAME");
		m_fieldsList.push_back("JUNK");
	}
	vector<string> m_fieldsList;
};

class CustomConstructorTests_EmptyList: public ::testing::Test {
public:
	CustomConstructorTests_EmptyList() { }
	vector<string> m_fieldsList;
};

class CustomConstructorTests_InvalidListMandatoryCHR: public ::testing::Test {
public:
	CustomConstructorTests_InvalidListMandatoryCHR() {
		m_fieldsList.push_back("STRAND");
		m_fieldsList.push_back("START_POS");
		m_fieldsList.push_back("END_POS");
		m_fieldsList.push_back("NAME");
		m_fieldsList.push_back("JUNK");
	}
	vector<string> m_fieldsList;
};

class CustomConstructorTests_InvalidListMandatorySTART_POS: public ::testing::Test {
public:
	CustomConstructorTests_InvalidListMandatorySTART_POS() {
		m_fieldsList.push_back("CHR");
		m_fieldsList.push_back("STRAND");
		m_fieldsList.push_back("END_POS");
		m_fieldsList.push_back("NAME");
		m_fieldsList.push_back("JUNK");
	}
	vector<string> m_fieldsList;
};

class CustomConstructorTests_InvalidListEND_POS: public ::testing::Test {
public:
	CustomConstructorTests_InvalidListEND_POS() {
		m_fieldsList.push_back("CHR");
		m_fieldsList.push_back("STRAND");
		m_fieldsList.push_back("START_POS");
		m_fieldsList.push_back("NAME");
		m_fieldsList.push_back("JUNK");
	}
	vector<string> m_fieldsList;
};
