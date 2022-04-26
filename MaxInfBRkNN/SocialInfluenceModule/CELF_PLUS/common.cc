#include "common.h"

float getCurrentMemoryUsage() {

	string pid = intToStr(unsigned(getpid()));
	string outfile = "temp/tmp_" + pid + ".txt";
	string command = "pmap " + pid + " | grep -i Total | awk '{print $2}' > " + outfile;
	//	cout << "Command to be executed: " << command << endl;
	system(command.c_str());

	string mem_str;
	ifstream ifs(outfile.c_str());
	std::getline(ifs, mem_str);
	ifs.close();

	mem_str = mem_str.substr(0, mem_str.size()-1);
	float mem = (float)strToInt(mem_str);

//	cout << "Memory used at this point : " << mem/1024 << " MB" << endl;
	return mem/1024;
	
	return 0;
}

/**** 
  Function for deleting pointers in map.
	example to delete a map : 
	for_each( StringMap.begin(),
              StringMap.end(),
              DeleteMap<int,CString*>());
  ****/
template<class A, class B>
struct DeleteMap
{
    // Overloaded () operator.
    // This will be called by for_each() function.
    bool operator()(pair<A,B> x) const
    {
        // Assuming the second item of map is to be
        // deleted. Change as you wish.
        delete x.second;
        return true;
    }
};


