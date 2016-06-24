#pragma once
#include <string>
#include<vector>
using namespace std;
using std::string;

#ifndef _MESSAGES_INCLUDED_
#define _MESSAGES_INCLUDED_

class Messages
{
public:
	void NotImplementedElement(int elemType, string whereHapp);
	~Messages();
	string logMessage(string);


};
#endif


//http ://www.tutorialspoint.com/cplusplus/cpp_date_time.htm