#pragma once

//#include <windows.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include "serialize.h"
#include "resultInfo.h"
using namespace std;

#define DEFAULT_CACHESIZE 1024
#define DEFAULT_BLOCKLENGTH	4096

//-----------
// LRU buffer
//-----------
enum AccessMode {isData,isIndex};

int SEQ_PAGE_ACCESS=0;
int RAN_PAGE_ACCESS=0;
int INDEX_PAGE_ACCESS=0;
int DATA_PAGE_ACCESS=0;
int CACHE_ACCESSED=0;

int blocklength=4096*4;//DEFAULT_BLOCKLENGTH;
int cachesize;	 	// max. number of cache blocks

// old: int *cache_cont; store transBlockID
int *cache_block,*cache_user; 	// block Id, user (file ID) separate !

int *lastTime;		// >=0 means used
int nextTime;
char **cache;		// cache cont

//int LastBlockIdStored;	// for determining SEQ or RAN
int LastBlockId,LastUserId;


typedef struct {
	int *cache_block_self,*cache_user_self; 	// block Id, user (file ID) separate !
	int *lastTime_self;		// >=0 means used
	int nextTime_self;
	char **cache_self;
	int LastBlockId_self;
	int LastUserId_self;
}PriveCacheData;

int getBlockLength() {
	return blocklength;
}

void InitCache(int csize) {
	assert(blocklength>0);	assert(csize>0);
	LastBlockId=-2;	LastUserId=-2;
	nextTime=0;
	cachesize=csize;
	cache_block = new int[cachesize];	cache_user = new int[cachesize];

	cache = new char*[cachesize];
	lastTime=new int[cachesize];
	for (int i=0;i<cachesize;i++) {
		lastTime[i]=-1;
		cache_block[i]=-1;  cache_user[i]=-1;
		cache[i] = new char[blocklength];
	}
}

void RefreshStat() {
	SEQ_PAGE_ACCESS=0;
	RAN_PAGE_ACCESS=0;
	INDEX_PAGE_ACCESS=0;
	DATA_PAGE_ACCESS=0;
	CACHE_ACCESSED=0;
}

void RefreshCache() {	// call before query execution
	LastBlockId=-2;	LastUserId=-2;
	nextTime=0;
	for (int i=0;i<cachesize;i++) {
		lastTime[i]=-1;
		cache_block[i]=-1;  cache_user[i]=-1;
	}
}

void DestroyCache() {
	// no need to flush the cache (since the query alg. are read-only )
	delete[] cache_block;	delete[] cache_user;
	delete[] lastTime;
	for (int i=0;i<cachesize;i++) delete[] cache[i];
	delete[] cache;
}


void DestroyPrivateCache(int *cache_block_self, int *cache_user_self, int *lastTime_self, int nextTime_self,char **cache_self) {
	// no need to flush the cache (since the query alg. are read-only )
	delete[] cache_block_self;	delete[] cache_user_self;
	delete[] lastTime_self;
	for (int i=0;i<cachesize;i++) delete[] cache_self[i];
	delete[] cache_self;
}

// light-weight get function
char* getCacheBlockRef(int UserId,int BlockId) {
	CACHE_ACCESSED++;
	for (int i=0;i<cachesize;i++)
		if ((cache_block[i]==BlockId)&&(cache_user[i]==UserId)&&
			(lastTime[i]>=0)) {
			lastTime[i]=nextTime++;
			return cache[i];
		}
	return NULL;	// NULL if non-existent
}

bool getCacheBlock(char* buffer,int UserId,int BlockId) {
	CACHE_ACCESSED++;
	int _t = cachesize;
	for (int i=0;i<cachesize;i++)
		if ((cache_block[i]==BlockId)&&(cache_user[i]==UserId)&&
			(lastTime[i]>=0)) {
			memcpy(buffer,cache[i],blocklength);
			lastTime[i]=nextTime++;
			//cout<<"block id="<<BlockId;
			//cout<<", hit Cache!"<<endl;
			return true;
		}
	return false;
}

bool getPrivateCacheBlock(char* buffer,int UserId,int BlockId,PriveCacheData& pcd) {
	CACHE_ACCESSED++;
	for (int i=0;i<cachesize;i++)
		if ((pcd.cache_block_self[i]==BlockId)&&(pcd.cache_user_self[i]==UserId)&&
			(lastTime[i]>=0)) {
			memcpy(buffer,pcd.cache_self[i],blocklength);
			lastTime[i]=nextTime++;
			return true;
		}
	return false;
}


// user's responsibility
// the place for counting number of page accesses
void storeCacheBlock(char* buffer,int UserId,int BlockId,AccessMode mode) {
	int index=-1;
	for (int i=0;i<cachesize;i++)	// search for empty block
		if (lastTime[i]<0) {index=i;	break;}

	if (index<0) {
		index=0;	// full, evict LRU block
		for (int i=0;i<cachesize;i++)
			if (lastTime[i]<lastTime[index]) index=i;
	}
	//int i = 0；
	//memcpy(&i,buffer, sizeof(int));

	memcpy(cache[index],buffer,blocklength);     //-------------------------问题出在这里
	cache_block[index]=BlockId;	cache_user[index]=UserId;
	lastTime[index]=nextTime++;

	if (mode==isData)
		DATA_PAGE_ACCESS++;
	else if (mode==isIndex)
		INDEX_PAGE_ACCESS++;

	if ((LastUserId==UserId)&&(abs(LastBlockId-BlockId)<=1))
		SEQ_PAGE_ACCESS++;
	else
		RAN_PAGE_ACCESS++;
	LastBlockId=BlockId;	LastUserId=UserId;
}





void printPageAccess() {
	printf("Cache requests:\t%d\n",CACHE_ACCESSED);
	printf("Index page accesses:\t%d\n",INDEX_PAGE_ACCESS);
	printf("Data page accesses:\t%d\n",DATA_PAGE_ACCESS);
	//printf("Seq.:\t%d\n",SEQ_PAGE_ACCESS);
	//printf("Ran.:\t%d\n",RAN_PAGE_ACCESS);
	printf("Total page accesses:\t%d\n",INDEX_PAGE_ACCESS+DATA_PAGE_ACCESS);
}

//----------
// FreqCache
//----------

struct FreqCache {
	char* buffer;
	int UserId,BlockId;
};
void InitFreqCache(FreqCache& fc) {
	fc.buffer=new char[blocklength];	// assume BlkLen inited
	fc.UserId=-1;	fc.BlockId=-1;
}

// assume buf of size BlkLen
void storeFreqCache(FreqCache& fc,char* buf,int uid,int bid) {
	memcpy(fc.buffer,buf,blocklength);
	fc.UserId=uid;	fc.BlockId=bid;
}

bool inFreqCache(FreqCache& fc,int uid,int bid) {
	//cout<<"inFreqCache"<<endl;
	return (fc.UserId==uid&&fc.BlockId==bid);
}

void DestroyFreqCache(FreqCache& fc) {
	delete[] fc.buffer;
}




//sigmod 2018 work

class IOcontroller
{
public:
	static void mkdir_absence(const char* outFolder)
	{
#if defined(_WIN32)
		CreateDirectoryA(outFolder, nullptr); // can be used on Windows
#else
		mkdir(outFolder, 0733); // can be used on non-Windows
#endif
	}

	/// Save a serialized file
	template <class T>
	static void save_file(const std::string filename, const T& output)
	{
		std::ofstream outfile(filename, std::ios::binary);
		if (!outfile.eof() && !outfile.fail())
		{
			StreamType res;
			serialize(output, res);
			outfile.write(reinterpret_cast<char*>(&res[0]), res.size());
			outfile.close();
			res.clear();
			std::cout << "Save file successfully: " << filename << '\n';
		}
		else
		{
			std::cout << "Save file failed: " + filename << '\n';
			exit(1);
		}
	}

	/// Load a serialized file
	template <class T>
	static void load_file(const std::string filename, T& input)
	{
		std::ifstream infile(filename, std::ios::binary);
		if (!infile.eof() && !infile.fail())
		{
			infile.seekg(0, std::ios_base::end);
			const std::streampos fileSize = infile.tellg();
			infile.seekg(0, std::ios_base::beg);
			std::vector<uint8_t> res(fileSize);
			infile.read(reinterpret_cast<char*>(&res[0]), fileSize);
			infile.close();
			input.clear();
			auto it = res.cbegin();
			input = deserialize<T>(it, res.cend());
			res.clear();
		}
		else
		{
			std::cout << "Cannot open file: " + filename << '\n';
			exit(1);
		}
	}

	/// Save graph structure to a file
	static void save_graph_struct(const std::string graphName, const Graph& vecGraph, const bool isReverse)
	{
		std::string postfix = ".vec.graph";
		if (isReverse) postfix = ".vec.rvs.graph";
		const std::string filename = graphName + postfix;
		save_file(filename, vecGraph);
	}

	/// Load graph structure from a file
	static void load_graph_struct(const std::string graphName, Graph& vecGraph, const bool isReverse)
	{
		std::string postfix = ".vec.graph";
		if (isReverse) postfix = ".vec.rvs.graph";
		const std::string filename = graphName + postfix;
		load_file(filename, vecGraph);
	}

	/// Get out-file name
	static std::string get_out_file_name(const std::string graphName, const std::string algName, const int seedsize,
		const std::string probDist, const float probEdge)
	{
		if (probDist == "UNI")
		{
			return graphName + "_" + algName + "_k" + std::to_string(seedsize) + "_" + probDist + std::
				to_string(probEdge);
		}
		return graphName + "_" + algName + "_k" + std::to_string(seedsize) + "_" + probDist;
	}

	/// Print the results
	static void write_result(const std::string& outFileName, const TResult& resultObj, const std::string& outFolder)
	{
		const auto approx = resultObj.get_approximation();
		const auto runTime = resultObj.get_running_time();
		const auto influence = resultObj.get_influence();
		const auto influenceOriginal = resultObj.get_influence_original();
		const auto seedSize = resultObj.get_seed_size();
		const auto RRsetsSize = resultObj.get_RRsets_size();

		std::cout << "   --------------------" << std::endl;
		std::cout << "  |Approx.: " << approx << std::endl;
		std::cout << "  |Time (sec): " << runTime << std::endl;
		std::cout << "  |Influence: " << influence << std::endl;
		std::cout << "  |Self-estimated influence: " << influenceOriginal << std::endl;
		std::cout << "  |#Seeds: " << seedSize << std::endl;
		std::cout << "  |#RR sets: " << RRsetsSize << std::endl;
		std::cout << "   --------------------" << std::endl;
		mkdir_absence(outFolder.c_str());
		std::ofstream outFileNew(outFolder + "/" + outFileName);
		if (outFileNew.is_open())
		{
			outFileNew << "Approx.: " << approx << std::endl;
			outFileNew << "Time (sec): " << runTime << std::endl;
			outFileNew << "Influence: " << influence << std::endl;
			outFileNew << "Self-estimated influence: " << influenceOriginal << std::endl;
			outFileNew << "#Seeds: " << seedSize << std::endl;
			outFileNew << "#RR sets: " << RRsetsSize << std::endl;
			outFileNew.close();
		}
	}

	/// Print the seeds
	static void write_order_seeds(const std::string& outFileName, const TResult& resultObj, const std::string& outFolder)
	{
		auto vecSeed = resultObj.get_seed_vec();
		mkdir_absence(outFolder.c_str());
		const auto outpath = outFolder + "/seed";
		mkdir_absence(outpath.c_str());
		std::ofstream outFile(outpath + "/seed_" + outFileName);
		for (auto i = 0; i < vecSeed.size(); i++)
		{
			outFile << vecSeed[i] << '\n';
		}
		outFile.close();
	}
};

using TIO = IOcontroller;
using PIO = std::shared_ptr<IOcontroller>;
