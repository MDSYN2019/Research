#pragma once

#ifndef _StringTools_
#define _StringTools_

#include <string>  
#include <iostream>  
#include <cwctype>  
#include <vector>
#include <sstream>

using namespace std;

namespace StringTools {

	string trimString(string& s);
	void splitString(std::string& s, std::string delim, std::vector< std::string >* ret);
}

#endif
