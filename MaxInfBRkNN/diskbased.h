#ifndef  DISKBASED_H

#include "btree.h"
#include "utility.h"
#include "netshare.h"
#include "blk_file.cc"
#include "IOcontroller.h"


int PtNum;
int PtFileSize,AdjFileSize,GNFileSize;
int PLFileSize, POIIdxFileSize, PIVFileSize;
int V2PFileSize, NVDGFileSize;
int V2PIDXFileSize, NVDGIDXFileSize;

FILE *PtFile,*AdjFile;
FILE *GtreeNodeFile;  // the file of .gtree-disk
//--------------------------------------0928
FILE *VertexHashFile;
FILE *VertexHashIDXFile;

FILE *NVD_ADJGraph;
FILE *NVD_ADJGraphIDXFile;

/*FILE *POILeafFile;  //根据leaf node为单位分组写入
FILE *POIIdxFile;   //基于以上形式的索引地址
FILE *POIInvFile;   //纯粹的只有兴趣点信息的倒排索引*/


//-------数据-------- 地址---------倒排--------
int O2ULFileSize, O2UIdxFileSize,O2UIVFileSize;
FILE *O2ULeafFile;  //双色体对象源数据文件
FILE *O2UIdxFile;  //双色体对象数据块首地址索引文件
FILE *O2UInvFile; //双色体对象关键词信息的倒排索引
int o_baseAddr;   //object content 在 idx中的基地址
int u_baseAddr;	  //user content 在 idx中的基地址
int term_idx_address;



BTree *PtTree;
int i_capacity;
char *gBuffer;	// temp, global buffer for atomic use
FreqCache FC_A,FC_P; // FC for ensureCache exclusive use !
FreqCache FC_GN; // FC of gtree-disk file!
FreqCache FC_PL;
FreqCache FC_poiIdx; // FC of object invert index
FreqCache FC_PIV; // FC of object invert index


//below is for 双色体下问题
FreqCache FC_O2UL;
FreqCache FC_O2UIdx; // FC of object invert index
FreqCache FC_O2UIV; // FC of object invert index


//below is for NVD
FreqCache FC_v2p;
FreqCache FC_NVD_ADJGraph;
FreqCache FC_v2pIdx;
FreqCache FC_NVD_ADJGraphIdx;


int BlkLen;


#define PtTid 	1000
#define PtFid 	2000
#define AdjTid 	3000
#define AdjFid 	4000
#define GNFid 5000
#define O2ULFid 6000
#define O2ULIdxFid 7000
#define PInvFid 8000
#define O2UInvFid 9000

#define V2PFid 10000
#define NVDFid 11000

#define V2PFIDXid 12000
#define NVDFIDXid 13000


#define Ref(obj) ((void*)&obj)

int block_num =0;
int idx_block_num =0;
int data_block_num=0;

struct NodeStruct {
	char level;
	int num_entries;
	int* key;
	int* son;
};

NodeStruct* createNodeStruct() {
	NodeStruct* node=new NodeStruct();
	node->level=-1;		node->num_entries=-1;
	node->key=new int[i_capacity];
	node->son=new int[i_capacity];
	return node;
}

int getFileSize(FILE* f) {	// side effect, setting f to begin
	fseek(f,0,SEEK_END);
	int filesize=ftell(f);
	//fseek(f,0,SEEK_SET);
	return filesize;
}



void ReadIndexBlock(BTree* bt,int block,NodeStruct* node) {
	int UserId=(bt==PtTree)?PtTid:AdjTid;
	if (!getCacheBlock(gBuffer,UserId,block)) {
		CachedBlockFile* cf=bt->file;
		cf->read_block(gBuffer,block);
		storeCacheBlock(gBuffer,UserId,block,isIndex);
	}

	int j=0;	// read node header
	memcpy(&(node->level), &gBuffer[j], sizeof(char));
	j += sizeof(char);
	memcpy(&(node->num_entries), &gBuffer[j], sizeof(int));
	j += sizeof(int);
	int _size = node->num_entries;

	// read node content
	for (int i=0;i<node->num_entries;i++) {
		memcpy(&(node->key[i]),&gBuffer[j],sizeof(int));
		memcpy(&(node->son[i]),&gBuffer[j+sizeof(int)],sizeof(int));
		j += sizeof(int)+sizeof(int);
	}
}

int inline pointQuery(BTree* bt,int key,int& TreeKey) {
	static NodeStruct* node=createNodeStruct();
	TreeKey=-2;
	if (key<0) return -3;	// assume non-negative keys

	ReadIndexBlock(bt,bt->root,node);
	while (true) {
		int curLevel=(int)(node->level);
		int Son=-1;

		if (curLevel>1) {
			for (int i=(node->num_entries-1);i>=0;i--)
				if (key<=node->key[i]) Son=i;
			if (Son<0) return -4;
			ReadIndexBlock(bt,node->son[Son],node);
		} else {  // curLevel is 1
			//printf("%d %d %d\n",key,node->key[0],node->key[node->num_entries-1]);
			for (int i=(node->num_entries-1);i>=0;i--) {
				if (key>=node->key[i]) { // use a different test than above
					TreeKey=node->key[i];
					return node->son[i];
				}
			}
			return -5;
		}
	}
	return -6;
}
enum ReadType{POIDATA, ADJDATA};



//this for access index file including O2UIdx, Gtree and O2U inverted file
char* getFlatBlock(FILE* myfile,int BlockId) {
	int UserId,FileSize;
	FreqCache* FC;

	if (myfile==PtFile) {
		FC=&FC_P;	UserId=PtFid;	FileSize=PtFileSize;   //有可能是FC_P有问题
	}
	else if(myfile==AdjFile){
		FC=&FC_A;	UserId=AdjFid;	FileSize=AdjFileSize;
	}

	else if(myfile==O2ULeafFile){
		FC=&FC_O2UL;  UserId=O2ULFid;	FileSize = O2ULFileSize;
	}


	//检查是否在频繁cache中命中，若是则返回频繁cache中内容
	if (inFreqCache(*FC,UserId,BlockId)) return (FC->buffer);
	//检查是否在cache数组中命中，若是则把cache内容存储到频繁cache中，否则进行IO操作
	if (!getCacheBlock(FC->buffer,UserId,BlockId)) {    //从文件中读取数据
		int readsize=min(BlkLen,FileSize-BlockId*BlkLen);
		fseek(myfile,BlockId*BlkLen,SEEK_SET);
		fread(FC->buffer,readsize,sizeof(char),myfile);
		storeCacheBlock(FC->buffer,UserId,BlockId,isData);
		if(UserId==6000||UserId==8000)block_num++;   //访问双色体数据读取block次数加1 ！
		//cout<<block_num<<endl;
	}//该if语句结束后，更新的都是FC中的内容,更新的内容要么从cache数据中获得，要么从频繁cache中获得
	FC->UserId=UserId;	FC->BlockId=BlockId;
	return (FC->buffer);
}


//this for access index file including O2UIdx, Gtree and O2U inverted file
char* getFlatBlockForIdx(FILE* myfile,int BlockId) {
	int UserId,FileSize;
	FreqCache* FC;

	//if (myfile==POIIdxFile){
		//FC=&FC_poiIdx; UserId = POILeafIdxid; FileSize=POIIdxFileSize;
	//}
	if (myfile==O2UIdxFile){
		FC=&FC_O2UIdx; UserId = O2ULIdxFid; FileSize=O2UIdxFileSize;
	}
	else if(myfile==GtreeNodeFile){
		FC=&FC_GN;  UserId=GNFid;	FileSize = GNFileSize;
		//cout<<"FC=&FC_GN;  UserId=GNFid;\tFileSize = GNFileSize"<<endl;
	}
	else if (myfile==O2UInvFile){
		FC=&FC_O2UIV; UserId= O2UInvFid;  FileSize = O2UIVFileSize;
		//cout<<"FC=&FC_O2UIV; UserId= O2UInvFid;  FileSize = O2UIVFileSize;"<<endl;
	}
	else if(myfile == VertexHashFile){
		FC =&FC_v2p;  UserId= V2PFid;   FileSize = V2PFileSize;
	}
	else if(myfile == NVD_ADJGraph){
		FC = &FC_NVD_ADJGraph; UserId= NVDFid; FileSize = NVDGFileSize;
	}
	else if(myfile == VertexHashIDXFile){
		FC =&FC_v2pIdx;  UserId= V2PFIDXid;   FileSize = V2PIDXFileSize;
	}
	else if(myfile == NVD_ADJGraphIDXFile){
		FC =&FC_NVD_ADJGraphIdx;  UserId= NVDFIDXid;   FileSize = NVDGIDXFileSize;
	}

	//检查是否在频繁cache中命中，若是则返回频繁cache中内容
	if (inFreqCache(*FC,UserId,BlockId)) {
		//cout<<"block id="<<BlockId;
		//cout<<", hit FreqCache"<<endl;
		return (FC->buffer);  //error
	}
	//检查是否在cache数组中命中，若是则把cache内容存储到频繁cache中，否则进行IO操作
	if (!getCacheBlock(FC->buffer,UserId,BlockId)) {    //从文件中读取数据
		int readsize=min(BlkLen,FileSize-BlockId*BlkLen);
		fseek(myfile,BlockId*BlkLen,SEEK_SET);
		fread(FC->buffer,readsize,sizeof(char),myfile);
		storeCacheBlock(FC->buffer,UserId,BlockId,isIndex);
#ifdef IO_TRACK
		cout<<"UserId="<<UserId<<",block id="<<BlockId;
		cout<<", access disk and store cache"<<endl;
#endif
		if(UserId==5000||UserId==9000)block_num++;   //访问GIMTree文件 或 倒排索引表读取block次数加1 ！  //gc
/*#ifdef IO_TRACK
		cout<<block_num<<endl;
#endif*/
	}//该if语句结束后，更新的都是FC中的内容,更新的内容要么从cache数据中获得，要么从频繁cache中获得
	FC->UserId=UserId;	FC->BlockId=BlockId;
	return (FC->buffer);
}

//access the page from index file
void getIdxBlockByjins(FILE* f,int addr, int readSize, void* readData){
	//int readSize = sizeof(int);
	int _address_in_block = addr%BlkLen;
	int end_inblock_addr = addr%BlkLen+readSize-1;
	int first_blockID = addr/BlkLen;
#ifdef IO_TRACK
	cout<<"first_blockID="<<first_blockID<<endl;
#endif

	char* BlockAddr; char full_data[readSize]; char* required_data;
	//数据跨block
	if(end_inblock_addr < BlkLen){ //数据未跨block

#ifdef IO_TRACK
		cout<<"数据未跨block"<<endl;
		cout<<"block id="<<(addr/BlkLen)<<", inaddr="<<(addr%BlkLen)<<endl;
#endif

		BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
		memcpy(full_data,BlockAddr+(addr%BlkLen), readSize);
		addr += readSize;
		//required_data = full_data;

	}

	else{   //数据跨block

#ifdef IO_TRACK
		cout<<"数据跨block"<<endl;
#endif
		int read_size_in_fistBlock = BlkLen-(_address_in_block);
		int read_size_in_secondBlock = readSize - read_size_in_fistBlock;
		int second_blockID = first_blockID + 1;
		char data_part1[read_size_in_fistBlock];char data_part2[read_size_in_secondBlock];
		//读第一部分
		BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
		memcpy(data_part1,BlockAddr+(addr%BlkLen), read_size_in_fistBlock);

		//读第二部分
		BlockAddr = getFlatBlockForIdx(f, addr/BlkLen + 1);
		memcpy(data_part2,BlockAddr, read_size_in_secondBlock);

		addr += readSize;

		for(int i=0;i<read_size_in_fistBlock;i++)
			full_data[i] = data_part1[i];
		for(int i=0;i<read_size_in_secondBlock;i++)
			full_data[i+read_size_in_fistBlock]= data_part2[i];


	}

	memcpy(readData,full_data, readSize);

}

//access the page from data file
void getDataBlockByjins(FILE* f,int addr, int readSize, void* readData){
	//int readSize = sizeof(int);
	int _address_in_block = addr%BlkLen;
	int end_inblock_addr = addr%BlkLen+readSize-1;
	int first_blockID = addr/BlkLen;

	char* BlockAddr; char full_data[readSize]; char* required_data;
	//数据跨block
	if(end_inblock_addr > BlkLen){
		int read_size_in_fistBlock = BlkLen-(_address_in_block);
		int read_size_in_secondBlock = readSize - read_size_in_fistBlock;
		int second_blockID = first_blockID + 1;
		char data_part1[read_size_in_fistBlock];char data_part2[read_size_in_secondBlock];
		//读第一部分
		BlockAddr = getFlatBlock(f, addr/BlkLen);
		memcpy(data_part1,BlockAddr+(addr%BlkLen), read_size_in_fistBlock);

		//读第二部分
		BlockAddr = getFlatBlock(f, addr/BlkLen + 1);
		memcpy(data_part2,BlockAddr, read_size_in_secondBlock);

		addr += readSize;

		for(int i=0;i<read_size_in_fistBlock;i++)
			full_data[i] = data_part1[i];
		for(int i=0;i<read_size_in_secondBlock;i++)
			full_data[i+read_size_in_fistBlock]= data_part2[i];

		//required_data = full_data;

	}

	else{
		BlockAddr = getFlatBlock(f, addr/BlkLen);
		memcpy(full_data,BlockAddr+(addr%BlkLen), readSize);
		addr += readSize;
		//required_data = full_data;
	}

	memcpy(readData,full_data, readSize);

}




void getBichromaticIdx_old(){
	FILE* f = O2UIdxFile;
	int address_logic = 0;

	//读出poi_size;
	int pSize =-1;
	char* BlockAddr=getFlatBlockForIdx(f,address_logic/BlkLen);
	memcpy(&pSize,BlockAddr+(address_logic%BlkLen), sizeof(int));
	address_logic += sizeof(int);

	//读出pbase地址;
	int pBase = -1;
	BlockAddr=getFlatBlockForIdx(f,address_logic/BlkLen);
	memcpy(&pBase,BlockAddr+(address_logic%BlkLen),sizeof(int));
	address_logic += sizeof(int);
	o_baseAddr = pBase;

	//读出uSize地址;
	int uSize = -1;
	BlockAddr=getFlatBlockForIdx(f,address_logic/BlkLen);
	memcpy(&uSize,BlockAddr+(address_logic%BlkLen),sizeof(int));
	address_logic += sizeof(int);


	//读出ubase地址;
	int uBase = -1;
	BlockAddr=getFlatBlockForIdx(f,address_logic/BlkLen);
	memcpy(&uBase,BlockAddr+(address_logic%BlkLen),sizeof(int));
	address_logic += sizeof(int);
	u_baseAddr = uBase;


}

void getBichromaticIdx(){
	FILE* f = O2UIdxFile;
	int address_logic = 0;

	//读出poi_size;
	int pSize =-1;
	char* BlockAddr=getFlatBlockForIdx(f,address_logic/BlkLen);
	memcpy(&pSize,BlockAddr+(address_logic%BlkLen), sizeof(int));
	address_logic += sizeof(int);

	//读出pbase地址;
	int pBase = -1;
	BlockAddr=getFlatBlockForIdx(f,address_logic/BlkLen);
	memcpy(&pBase,BlockAddr+(address_logic%BlkLen),sizeof(int));
	address_logic += sizeof(int);
	o_baseAddr = pBase;

	//读出uSize地址;
	int uSize = -1;
	BlockAddr=getFlatBlockForIdx(f,address_logic/BlkLen);
	memcpy(&uSize,BlockAddr+(address_logic%BlkLen),sizeof(int));
	address_logic += sizeof(int);


	//读出ubase地址;
	int uBase = -1;
	BlockAddr=getFlatBlockForIdx(f,address_logic/BlkLen);
	memcpy(&uBase,BlockAddr+(address_logic%BlkLen),sizeof(int));
	address_logic += sizeof(int);
	u_baseAddr = uBase;


}



// modify address for Adj as well ??
int PTGRP_HEADSIZE=3*sizeof(int)+sizeof(float);
int PTGRP_ITEMSIZE=sizeof(float);
int ADJGRP_HEADSIZE=sizeof(int);
int ADJGRP_ITEMSIZE=4*sizeof(int)+sizeof(float);	// changed  (重要！) , 4个int分别为： Nj, POIDataAddress, POIIdxAddress, OkeyEdgeAddr.


enum VarE {ADJNODE_A,DIST_A,PTKEY_A,PIDXKEY_A,PIDX_A,PT_P,  OKEY_A};

//针对邻接表上的数据内容的地址偏移函数
void getVarE(VarE type,void* buf,int BaseAddr,int pos) {
	// note default values !
	FILE* f=AdjFile;
	if(f==NULL) cout <<"faile to open Adjfile.";
	int size=sizeof(int),addr=-1;
	int VarBase=BaseAddr+ADJGRP_HEADSIZE+pos*ADJGRP_ITEMSIZE;  //邻接点数+ 位置*（Nj_id, dist, poiGrpKey, poiIdxKey）

	// for VarE in AdjFile
	if (type==ADJNODE_A) addr=VarBase;   //要往下偏移的
	if (type==DIST_A) {addr=VarBase+sizeof(int);size=sizeof(float);}
	if (type==PTKEY_A) addr=VarBase+sizeof(int)+sizeof(float);
	if (type==PIDX_A) addr=VarBase+sizeof(int)+sizeof(float)+sizeof(int);
	if (type==OKEY_A) addr=VarBase+sizeof(int)+sizeof(float)+sizeof(int)+ sizeof(int);

	// for VarE in PtFile
	if (type==PT_P) {
		addr=BaseAddr+PTGRP_HEADSIZE+pos*PTGRP_ITEMSIZE;
		f=PtFile;	size=sizeof(float);
	}

		//char* BlockAddr=getFlatBlock(f,addr/BlkLen);
		//memcpy(buf,BlockAddr+(addr%BlkLen),size);
	getDataBlockByjins(f,addr,size,buf);
}




enum FixedF {SIZE_A,NI_P,NJ_P,DIST_P,SIZE_P, POINTA,INVA};

void getFixedF(FixedF type,void* buf,int BaseAddr) {
	// note default values !
	//cout<<"getFixedF函数内！";
	FILE* f=PtFile;
	//if(f==NULL) cout <<"faile to open ptfile.";
	int size=sizeof(int),addr=-1;

	// for FixedF in PtFile, 根据不同的数据类型，进行各自地址的计算
	if (type==NI_P) addr=BaseAddr;
	if (type==NJ_P) addr=BaseAddr+sizeof(int);
	if (type==DIST_P) { addr=BaseAddr+2*sizeof(int); size=sizeof(float); }
	if (type==SIZE_P) { addr=BaseAddr+2*sizeof(int)+sizeof(float);}

	if (type==POINTA) { addr=BaseAddr+3*sizeof(int)+sizeof(float);size=sizeof(float);}
	if (type==INVA) { addr=BaseAddr+3*sizeof(int)+sizeof(float);size=sizeof(bool)*REALTOTALKEY;}
	// for FixedF in AdjFile
	if (type==SIZE_A) {
		addr=BaseAddr;
		f=AdjFile;
		//cout<<"getFixedF函数内部获取邻边信息"<<endl;
	}
	//块号= addr/BlkLen， 根据块号取得块地址
	//cout<<"块号(addr/BlkLen)="<<(addr/BlkLen)<<endl;
		//char* BlockAddr=getFlatBlock(f,addr/BlkLen);
	//cout<<"addr="<<addr<<endl;
	//cout<<"(块内地址)addr%BlkLen"<<(addr%BlkLen)<<endl;
		//memcpy(buf,BlockAddr+(addr%BlkLen),size);
	getDataBlockByjins(f,addr, size,buf);
}



//this is for adjacent list
int getAdjListGrpAddr(int NodeID) {	// using AdjFile
	int addr=sizeof(int)+NodeID*sizeof(int); //计算出保存顶点NodeID的邻接表位置的那条记录所在地址
	int GrpAddr;
	getDataBlockByjins(AdjFile, addr, sizeof(int),Ref(GrpAddr));
	//char* BlockAddr=getFlatBlock(AdjFile,addr/BlkLen); //根据块号(addr/BlkLen)获取对应块地址
	//memcpy(Ref(GrpAddr),BlockAddr+(addr%BlkLen),sizeof(int)); //最终地址ַ= 块地址+ 块内地址(addr%BlkLen)
	return GrpAddr;
}



enum FixedPOI {POI_ID,POI_DIS,POI_KEY,POI_CHECK};

int getPoiIdxAddr(int poi_idxList_key, int poi_idx) {	// using AdjFile
	int addr=poi_idxList_key+poi_idx*sizeof(int); //计算出保存顶点NodeID的邻接表位置的那条记录所在地址
	return addr;
}

void getPOIAddr(void* buf,int BaseAddr) {
	// note default values !
	//cout<<"getFixedF函数内！";
	FILE* f=PtFile;
	int addr = BaseAddr;
	//char* BlockAddr=getFlatBlock(f,addr/BlkLen);
	//memcpy(buf,BlockAddr+(addr%BlkLen), sizeof(int));
	getDataBlockByjins(f, addr, sizeof(int),buf);
}



POI getPOIDATA_has_CheckIN(int BaseAddr) {
	// note default values !
	//cout<<"getFixedF函数内！";
	FILE* f=PtFile;
	int addr = BaseAddr;
	char* BlockAddr=getFlatBlock(f,addr/BlkLen);
	//cout<<"addr="<<addr<<endl;
	//cout<<"(块内地址)addr%BlkLen"<<(addr%BlkLen)<<endl;
	POI p;
	memcpy(&p.id,BlockAddr+(addr%BlkLen), sizeof(int));
	memcpy(&p.dis,BlockAddr+(addr%BlkLen)+sizeof(int), sizeof(float));
	char* address_base = BlockAddr+(addr%BlkLen)+sizeof(int)+sizeof(float);
	int key_size = 0;
	memcpy(&key_size,address_base, sizeof(int));
	//cout<<"key_size="<<key_size<<endl;
	for(int i=0;i<key_size;i++){
		int term =  -1;
		memcpy(&term,address_base+sizeof(int)+i*sizeof(int), sizeof(int));
		//cout<<term<<endl;
		p.keywords.push_back(term);
	}
	address_base += (1+key_size)*sizeof(int);
	int check_usr_num = 0;
	memcpy(&check_usr_num,address_base, sizeof(int));
	cout<<"check_usr_num="<<check_usr_num<<endl;
	return p;
}

POI getPOIDATA(int BaseAddr) {
	// note default values !
	//cout<<"getFixedF函数内！";
	FILE* f=PtFile; POI p;
	int addr = BaseAddr;
	//char* BlockAddr=getFlatBlock(f,addr/BlkLen);
	//memcpy(&p.id,BlockAddr+(addr%BlkLen), sizeof(int));
	getDataBlockByjins(f, addr, sizeof(int),Ref(p.id));
	//更新基址位置,获得dis
	addr = BaseAddr+ sizeof(int);  //id的偏移


	//BlockAddr=getFlatBlock(f,addr/BlkLen);
	//memcpy(&p.dis,BlockAddr+(addr%BlkLen), sizeof(float));
	getDataBlockByjins(f, addr, sizeof(float),Ref(p.dis));
	addr = addr + sizeof(float);  //dis的偏移


	//更新基址位置,获得key_size
	int key_size = 0;
	getDataBlockByjins(f, addr, sizeof(int),Ref(key_size));
	//BlockAddr=getFlatBlock(f,addr/BlkLen);
	//memcpy(&key_size,BlockAddr+(addr%BlkLen), sizeof(int));

	//更新基址位置,获得key_size个term
	addr = addr+sizeof(int);
	for(int i=0;i<key_size;i++){
		int term =  -1;
		//BlockAddr=getFlatBlock(f,addr/BlkLen);
		//memcpy(&term,BlockAddr+(addr%BlkLen), sizeof(int));
		getDataBlockByjins(f, addr, sizeof(int),Ref(term));
		//cout<<term<<endl;
		p.keywords.push_back(term);
		p.keywordSet.insert(term);
		addr += sizeof(int);  //key_size的偏移
	}

	//获得check_size
	int usr_size = 0;
	//BlockAddr=getFlatBlock(f,addr/BlkLen);
	//memcpy(&usr_size,BlockAddr+(addr%BlkLen), sizeof(int));
	getDataBlockByjins(f, addr, sizeof(int),Ref(usr_size));
	addr = addr+sizeof(int);
	for(int i=0;i<usr_size;i++){
		int usr_id =  -1;
		//BlockAddr=getFlatBlock(f,addr/BlkLen);
		//memcpy(&usr_id,BlockAddr+(addr%BlkLen), sizeof(int));
		getDataBlockByjins(f, addr, sizeof(int),Ref(usr_id));
		//cout<<term<<endl;
		p.check_ins.push_back(usr_id);
		addr += sizeof(int);  //key_size的偏移
	}

	return p;
}


//this is for adj file


set<int> getOKeyEdge(int BaseAddr) {
	// note default values !
	set<int> keyword_onEdge;
	FILE* f=AdjFile;
	if(f==NULL) cout <<"faile to open Adjfile.";
	int size=sizeof(int),addr=-1;
	// change the address to data block for keywords on edge
	addr= BaseAddr; // the address of KeySizeOnEdge

	int _oKeySize = -1;
	char* BlockAddr=getFlatBlock(f,addr/BlkLen);
	memcpy(&_oKeySize,BlockAddr+(addr%BlkLen),size);
	addr+= sizeof(int);
	if(_oKeySize>0){
		int _key;
		for(int i=0;i<_oKeySize;i++){
			BlockAddr=getFlatBlock(f,addr/BlkLen);
			memcpy(&_key,BlockAddr+(addr%BlkLen),size);
			addr+= sizeof(int);
			keyword_onEdge.insert(_key);
		}
	}
	return keyword_onEdge;

}

set<int> getOKeyEdge(int Ni, int Adj_th) {
	// note default values !
	set<int> keyword_onEdge;
	FILE* f=AdjFile;
	if(f==NULL) cout <<"faile to open Adjfile.";
	int size=sizeof(int),addr=-1;
	//1. 取得顶点邻接表数据块基地址
	int AdjGrpAddr = getAdjListGrpAddr(Ni);
	//2. 进行地址偏移，确定相关内容的最终地址
	int okeyAddr = -1;
	getVarE(OKEY_A, Ref(okeyAddr),AdjGrpAddr,Adj_th);
	if(okeyAddr==-1)  //该条边无兴趣点，无任何关键字
		return keyword_onEdge;

	// 3. 根据步骤2中地址访问磁盘
	//okeyAddr = 12588320;
	addr= okeyAddr; // the address of KeySizeOnEdge
	int _oKeySize = -1;
	//char* BlockAddr=getFlatBlock(f,addr/BlkLen);
	//memcpy(&_oKeySize,BlockAddr+(addr%BlkLen),size);
	getDataBlockByjins(f,addr, sizeof(int),Ref(_oKeySize));
	addr+= sizeof(int);
	if(_oKeySize>0){
		int _key;
		for(int i=0;i<_oKeySize;i++){
			//BlockAddr=getFlatBlock(f,addr/BlkLen);
			//memcpy(&_key,BlockAddr+(addr%BlkLen),size);
			getDataBlockByjins(f,addr, sizeof(int),Ref(_key));
			addr+= sizeof(int);
			keyword_onEdge.insert(_key);
		}
	}
	return keyword_onEdge;

}

//返回边上兴趣点数据块基地址的指针
int getPOIDataAddrOfEdge(int Ni, int Adj_th) {
	// note default values !
	FILE* f=AdjFile;
	if(f==NULL) cout <<"faile to open Adjfile.";
	int size=sizeof(int),addr=-1;
	//1. 取得顶点邻接表数据块基地址
	int AdjGrpAddr = getAdjListGrpAddr(Ni);
	//2. 进行地址偏移，确定相关内容的最终地址
	int poiDataAddr = -1;
	getVarE(PTKEY_A, Ref(poiDataAddr),AdjGrpAddr,Adj_th);

	return poiDataAddr;

}

//返回边上兴趣点对应的地址索引表基地址的指针
int getPOIIdxAddrOfEdge(int Ni, int Adj_th) {
	// note default values !
	FILE* f=AdjFile;
	if(f==NULL) cout <<"faile to open Adjfile.";
	int size=sizeof(int),addr=-1;
	//1. 取得顶点邻接表数据块基地址
	int AdjGrpAddr = getAdjListGrpAddr(Ni);
	//2. 进行地址偏移，确定相关内容的最终地址
	int poiIdxAddr = -1;
	getVarE(PIDX_A, Ref(poiIdxAddr),AdjGrpAddr,Adj_th);

	return poiIdxAddr;

}

//读取该边的上兴趣点的个数
int getPOISizeOfEdge(int Ni, int Adj_th){
	int poiData_Address = -1;
	poiData_Address = getPOIDataAddrOfEdge(Ni,Adj_th);
	int poi_size = -1;
	getFixedF(SIZE_P,Ref(poi_size),poiData_Address); //读取边上poi个数
	return  poi_size;
}

//根据相关地址， 返回该edge 上第 poi_th个POI
POI getPOIDataOnEdge(int idx_baseAddr, int poi_th){
	int poi_idx_address = getPoiIdxAddr(idx_baseAddr,poi_th);  //索引表所在的地址
	int poi_address;
	getPOIAddr(&poi_address, poi_idx_address);  //根据poi索引地址（poi_idx_address），获取poi数据块地址（poi_address）
	//根据poi数据块地址（poi_address）， 获取poi数据内容
	POI p = getPOIDATA(poi_address);
	return p;
}

POI getPOIDataOnEdge(int Ni, int adj_th, int poi_th){
	int idx_baseAddr = getPOIIdxAddrOfEdge(Ni,adj_th);
	int poi_idx_address = getPoiIdxAddr(idx_baseAddr,poi_th);  //索引表所在的地址
	int poi_address;
	getPOIAddr(&poi_address, poi_idx_address);  //根据poi索引地址（poi_idx_address），获取poi数据块地址（poi_address）
	//根据poi数据块地址（poi_address）， 获取poi数据内容
	POI p = getPOIDATA(poi_address);
	return p;
}





//this is for gree tree index


int getGTreeNodeAddr(int treenode){
	int addr_logic= sizeof(int)+treenode* sizeof(int);
	int addr_real; int readSize = sizeof(int);
	//cout<<"addr_logic="<<addr_logic<<endl;
	getIdxBlockByjins(GtreeNodeFile,addr_logic, readSize,Ref(addr_real));  
	//memcpy(Ref(addr_real),data, sizeof(int));
	//cout<<"addr_real="<<addr_real<<endl;
	return addr_real;

}




enum GIMTree_AccessPattern{OnlyU,OnlyO, Both};

TreeNode getGIMTreeNodeData_old(int node_id, int pattern){

	int BaseAddr = getGTreeNodeAddr(node_id);

	FILE* f = GtreeNodeFile; TreeNode tn;
	int addr = BaseAddr;  //基地址
	//读入border 的个数
	char* BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
	int border_size;
	memcpy(Ref(border_size),BlockAddr+(addr%BlkLen), sizeof(int));
	addr += sizeof(int);
	//依次读入各border id
	for(int i=0;i<border_size;i++){
		//依次读入各border id
		BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
		int border_id;
		memcpy(Ref(border_id),BlockAddr+(addr%BlkLen), sizeof(int));
		tn.borders.push_back(border_id);
		addr += sizeof(int);
	}

	//读入child 的个数
	BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
	int child_size;
	memcpy(Ref(child_size),BlockAddr+(addr%BlkLen), sizeof(int));
	addr += sizeof(int);
	//依次读入各child node id
	for(int i=0;i<child_size;i++){
		BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
		int child_node;
		memcpy(Ref(child_node),BlockAddr+(addr%BlkLen), sizeof(int));
		tn.children.push_back(child_node);
		addr += sizeof(int);
	}

	//读入isleaf flag
	BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
	bool isleaf;
	memcpy(Ref(isleaf),BlockAddr+(addr%BlkLen), sizeof(bool));
	tn.isleaf = isleaf;
	//cout<<"tn.isleaf="<<isleaf<<endl;
	addr += sizeof(bool);

	//读入 leaf 个数
	BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
	int leaf_size;
	memcpy(Ref(leaf_size),BlockAddr+(addr%BlkLen), sizeof(int));
	addr += sizeof(int);
	//依次读入各 leaf id
	for(int i=0;i<leaf_size;i++){
		BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
		int leaf_id;
		memcpy(Ref(leaf_id),BlockAddr+(addr%BlkLen), sizeof(int));
		tn.leafnodes.push_back(leaf_id);
		addr += sizeof(int);
	}

	//读入father
	BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
	int father =-2;
	memcpy(Ref(father),BlockAddr+(addr%BlkLen), sizeof(int));
	tn.father = father;
	addr += sizeof(int);

	//读入minSocialCount
	BlockAddr = getFlatBlockForIdx(f,addr/BlkLen);
	int _count = -1;
	memcpy(Ref(_count),BlockAddr+(addr%BlkLen), sizeof(int));
	tn.minSocialCount = _count;
	addr += sizeof(int);


	//读入entry对应keyword的vocabulary
	if(pattern==OnlyO || pattern==Both){
		//读入 poi vocabulary 中的keyword 的个数
		int keyword_size;
		BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
		memcpy(Ref(keyword_size),BlockAddr+(addr%BlkLen), sizeof(int));
		addr += sizeof(int);
		//针对兴趣点： 依次读入各 <keyword, inneraddress, related_object_cnt>
		int innerAddress =-1;
		int keyword = -2;
		//fseek(f,addr,SEEK_SET);
		for(int i=0;i<keyword_size;i++){

			//read okeyword
			BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
			memcpy(Ref(keyword),BlockAddr+(addr%BlkLen), sizeof(int));
			tn.objectUKeySet.insert(keyword);
			cout<<"keyword="<<keyword<<",addr="<<addr<<endl;
			addr += sizeof(int);

			if(keyword==178){
				cout<<"find 178"<<endl;
			}


			BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
			memcpy(Ref(innerAddress),BlockAddr+(addr%BlkLen), sizeof(int));
			tn.term_innerAddr_OMap[keyword].push_back(innerAddress);
			cout<<" innerAddress="<<innerAddress<<",addr="<<addr<<endl;
			addr += sizeof(int);
			//cout<<tn.term_innerAddr_OMap[keyword][0];


			//read related object_size
			int object_size = -1;
			BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
			memcpy(Ref(object_size),BlockAddr+(addr%BlkLen), sizeof(int));
			tn.term_object_countMap[keyword] = object_size;
			cout<<"object_size="<<object_size<<",addr="<<addr<<endl;
			addr += sizeof(int);



			if(true){  //read  keyword==178
				cout<<"keyword="<<keyword<<",innerAddress="<<innerAddress<<" object_size="<<object_size<<endl;
				//getchar();
			}
		}
	}
	else if(pattern==OnlyU){
		//跳过 object 部分的keyword内容，直接进行地址更改即可
		BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
		int Okeyword_size;  //计算要跳过的字节数
		memcpy(Ref(Okeyword_size),BlockAddr+(addr%BlkLen), sizeof(int));
		addr += sizeof(int);
		addr += sizeof(int)*3*Okeyword_size;

		//读入 usr vocabulary 中keyword 的个数
		BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
		int keyword_size2;
		memcpy(Ref(keyword_size2),BlockAddr+(addr%BlkLen), sizeof(int));
		addr += sizeof(int);
		//针对用户：依次读入各 <keyword, inneraddress>， 并加入 term_innerAddr_UMap 中
		for(int i=0;i<keyword_size2;i++){
			//read ukeyword
			int keyword;
			BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
			memcpy(Ref(keyword),BlockAddr+(addr%BlkLen), sizeof(int));
			tn.userUKeySet.insert(keyword);
			addr += sizeof(int);

			//read innerAddr
			int innerAddress;
			BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
			memcpy(Ref(innerAddress),BlockAddr+(addr%BlkLen), sizeof(int));
			tn.term_innerAddr_UMap[keyword].push_back(innerAddress);
			addr += sizeof(int);

			//read related user_size
			int user_size = -1;
			BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
			memcpy(Ref(user_size),BlockAddr+(addr%BlkLen), sizeof(int));
			tn.term_usr_countMap[keyword] = user_size;
			addr += sizeof(int);

		}
	}


	return  tn;

}

TreeNode getGIMTreeNodeData(int node_id, int pattern){

    if(GTreeCache[node_id].objectUKeySet.size()>0&& pattern==OnlyO)
        return GTreeCache[node_id];
    else if(GTreeCache[node_id].userUKeySet.size()>0 && pattern==OnlyU)
        return GTreeCache[node_id];

	//cout<<"获取node"<<node_id<<"地址"<<endl;
	int node_Addr = getGTreeNodeAddr(node_id);

	FILE* f = GtreeNodeFile; TreeNode tn;
	int addr = node_Addr;  //基地址



	//cout<<"读入border 的个数"<<endl;
	int border_size;
	getIdxBlockByjins(f,addr,sizeof(int),Ref(border_size));
	addr += sizeof(int);
	//依次读入各border id
	for(int i=0;i<border_size;i++){
		//依次读入各border id
		int border_id;
		getIdxBlockByjins(f,addr,sizeof(int),Ref(border_id));
		tn.borders.push_back(border_id);
		addr += sizeof(int);
	}

	//cout<<"读入child 的个数"<<endl;
	int child_size;
	getIdxBlockByjins(f,addr,sizeof(int),Ref(child_size));
	addr += sizeof(int);
	//依次读入各child node id
	for(int i=0;i<child_size;i++){
		int child_node;
		getIdxBlockByjins(f,addr,sizeof(int),Ref(child_node));
		tn.children.push_back(child_node);
		addr += sizeof(int);
	}

	//cout<<"读入isleaf flag"<<endl;
	bool isleaf;
	getIdxBlockByjins(f,addr,sizeof(bool),Ref(isleaf));
	tn.isleaf = isleaf;
	addr += sizeof(bool);

	//读入 leaf 个数
	int leaf_size;
	getIdxBlockByjins(f,addr,sizeof(int),Ref(leaf_size));
	addr += sizeof(int);
	//依次读入各 leaf id
	for(int i=0;i<leaf_size;i++){
		int leaf_id;
		getIdxBlockByjins(f, addr, sizeof(int),Ref(leaf_id));
		tn.leafnodes.push_back(leaf_id);
		addr += sizeof(int);
	}

	//读入father
	int father =-2;
	getIdxBlockByjins(f, addr, sizeof(int),Ref(father));
	tn.father = father;
	addr += sizeof(int);

	//读入minSocialCount
	int _count = -1;
	getIdxBlockByjins(f, addr, sizeof(int),Ref(_count));
	tn.minSocialCount = _count;
	addr += sizeof(int);


	//cout<<"读入entry对应keyword的vocabulary"<<endl;
	if(pattern==OnlyO || pattern==Both){
		//读入 poi vocabulary 中的keyword 的个数
		int keyword_size;
		getIdxBlockByjins(f, addr, sizeof(int),Ref(keyword_size));
		addr += sizeof(int);
		//针对兴趣点： 依次读入各 <keyword, inneraddress, related_object_cnt>

		//fseek(f,addr,SEEK_SET);
		for(int i=0;i<keyword_size;i++){

			//read okeyword
			int keyword =-1;
			getIdxBlockByjins(f, addr, sizeof(int),Ref(keyword));
			tn.objectUKeySet.insert(keyword);
			addr += sizeof(int);
			//cout<<"keyword="<<keyword<<",addr="<<addr<<endl;


			//read innerAddress
			int innerAddress=-100;
			int readSize = sizeof(int);
			getIdxBlockByjins(f,addr,readSize,Ref(innerAddress));
			tn.term_innerAddr_OMap[keyword].push_back(innerAddress);
			addr +=  sizeof(int);
			//cout<<tn.term_innerAddr_OMap[keyword][0];


			//read related object_size
			int object_size = -1;
			getIdxBlockByjins(f,addr, sizeof(int),Ref(object_size));
			tn.term_object_countMap[keyword] = object_size;
			addr += sizeof(int);
			//cout<<"object_size="<<object_size<<",addr="<<addr<<endl;


			//cout<<"keyword="<<keyword<<",innerAddress="<<innerAddress<<" object_size="<<object_size<<endl;
				//getchar();
		}
	}
	else if(pattern==OnlyU){
		//跳过 object 部分的keyword内容，直接进行地址更改即可
		int Okeyword_size;
		getIdxBlockByjins(f,addr, sizeof(int),Ref(Okeyword_size));
		//计算要跳过的字节数
		addr += sizeof(int);
		addr += sizeof(int)*3*Okeyword_size;

		//读入 usr vocabulary 中keyword 的个数
		int keyword_size2;
		getIdxBlockByjins(f,addr, sizeof(int),Ref(keyword_size2));
		addr += sizeof(int);

		//针对用户：依次读入各 <keyword, inneraddress>， 并加入 term_innerAddr_UMap 中
		for(int i=0;i<keyword_size2;i++){
			//read ukeyword
			int keyword;
			getIdxBlockByjins(f,addr, sizeof(int),Ref(keyword));
			tn.userUKeySet.insert(keyword);
			addr += sizeof(int);

			//read innerAddr
			int innerAddress;
			getIdxBlockByjins(f,addr, sizeof(int),Ref(innerAddress));
			tn.term_innerAddr_UMap[keyword].push_back(innerAddress);
			addr += sizeof(int);

			//read related user_size
			int user_size = -1;
			getIdxBlockByjins(f,addr, sizeof(int),Ref(user_size));
			tn.term_usr_countMap[keyword] = user_size;
			addr += sizeof(int);

		}
	}

	GTreeCache[node_id] = tn;

	return  tn;

}

TreeNode getGIMTreeNodeData_fast(int node_id, int pattern, int term_id){

	//cout<<"获取node"<<node_id<<"地址"<<endl;
	int node_Addr = getGTreeNodeAddr(node_id);

	FILE* f = GtreeNodeFile; TreeNode tn;
	int addr = node_Addr;  //基地址


	//cout<<"读入border 的个数"<<endl;
	int border_size;
	getIdxBlockByjins(f,addr,sizeof(int),Ref(border_size));
	addr += sizeof(int);
	//依次读入各border id
	for(int i=0;i<border_size;i++){
		//依次读入各border id
		int border_id;
		getIdxBlockByjins(f,addr,sizeof(int),Ref(border_id));
		tn.borders.push_back(border_id);
		addr += sizeof(int);
	}

	//cout<<"读入child 的个数"<<endl;
	int child_size;
	getIdxBlockByjins(f,addr,sizeof(int),Ref(child_size));
	addr += sizeof(int);
	//依次读入各child node id
	for(int i=0;i<child_size;i++){
		int child_node;
		getIdxBlockByjins(f,addr,sizeof(int),Ref(child_node));
		tn.children.push_back(child_node);
		addr += sizeof(int);
	}

	//cout<<"读入isleaf flag"<<endl;
	bool isleaf;
	getIdxBlockByjins(f,addr,sizeof(bool),Ref(isleaf));
	tn.isleaf = isleaf;
	addr += sizeof(bool);

	//读入 leaf 个数
	int leaf_size;
	getIdxBlockByjins(f,addr,sizeof(int),Ref(leaf_size));
	addr += sizeof(int);
	//依次读入各 leaf id
	for(int i=0;i<leaf_size;i++){
		int leaf_id;
		getIdxBlockByjins(f, addr, sizeof(int),Ref(leaf_id));
		tn.leafnodes.push_back(leaf_id);
		addr += sizeof(int);
	}

	//读入father
	int father =-2;
	getIdxBlockByjins(f, addr, sizeof(int),Ref(father));
	tn.father = father;
	addr += sizeof(int);

	//读入minSocialCount
	int _count = -1;
	getIdxBlockByjins(f, addr, sizeof(int),Ref(_count));
	tn.minSocialCount = _count;
	addr += sizeof(int);


	//cout<<"读入entry对应keyword的vocabulary"<<endl;
	if(pattern==OnlyO){
		//读入 poi vocabulary 中的keyword 的个数
		int keyword_size;
		getIdxBlockByjins(f, addr, sizeof(int),Ref(keyword_size));
		addr += sizeof(int);
		//针对兴趣点： 依次读入各 <keyword, inneraddress, related_object_cnt>

		//fseek(f,addr,SEEK_SET);
		for(int i=0;i<keyword_size;i++){

			//read okeyword
			int keyword =-1;

			getIdxBlockByjins(f, addr, sizeof(int),Ref(keyword));
			tn.objectUKeySet.insert(keyword);
			addr += sizeof(int);

			if(keyword>term_id)
				return tn;

			//if(keyword==336){
				//cout<<"find t"<<keyword<<endl;
				//getchar();
			//}


			//read innerAddress
			int innerAddress=-100;
			int readSize = sizeof(int);
			if(true){  //keyword!=336
				getIdxBlockByjins(f,addr,readSize,Ref(innerAddress));
				tn.term_innerAddr_OMap[keyword].push_back(innerAddress);
				addr +=  sizeof(int);
			}/*else{
				int read_addr = addr;
				fseek(f,read_addr,SEEK_SET);
				fread(Ref(innerAddress),1, sizeof(int),f);
				addr += sizeof(int);

			}*/

			//cout<<"innerAddress="<<tn.term_innerAddr_OMap[keyword][0];


			//read related object_size
			int object_size = -1;
			getIdxBlockByjins(f,addr, sizeof(int),Ref(object_size));
			tn.term_object_countMap[keyword] = object_size;
			addr += sizeof(int);
			//cout<<"object_size="<<object_size<<",addr="<<addr<<endl;


			//cout<<"keyword="<<keyword<<",innerAddress="<<innerAddress<<" object_size="<<object_size<<endl;
			//
		}
		//getchar();
	}
	else if(pattern==OnlyU){
		//跳过 object 部分的keyword内容，直接进行地址更改即可
		int Okeyword_size;
		getIdxBlockByjins(f,addr, sizeof(int),Ref(Okeyword_size));
		//计算要跳过的字节数
		addr += sizeof(int);
		addr += sizeof(int)*3*Okeyword_size;

		//读入 usr vocabulary 中keyword 的个数
		int keyword_size2;
		getIdxBlockByjins(f,addr, sizeof(int),Ref(keyword_size2));
		addr += sizeof(int);

		//针对用户：依次读入各 <keyword, inneraddress>， 并加入 term_innerAddr_UMap 中
		for(int i=0;i<keyword_size2;i++){
			//read ukeyword
			int keyword;
			getIdxBlockByjins(f,addr, sizeof(int),Ref(keyword));
			tn.userUKeySet.insert(keyword);
			addr += sizeof(int);
			if(keyword>term_id)
				return tn;


			//read innerAddr
			int innerAddress;
			getIdxBlockByjins(f,addr, sizeof(int),Ref(innerAddress));
			tn.term_innerAddr_UMap[keyword].push_back(innerAddress);
			addr += sizeof(int);

			//read related user_size
			int user_size = -1;
			getIdxBlockByjins(f,addr, sizeof(int),Ref(user_size));
			tn.term_usr_countMap[keyword] = user_size;
			addr += sizeof(int);

		}
	}

	return  tn;

}



//below is for POI data org Leaf



//below is for bichromatic object(poi, user)
int getO2UAddrByIdxFile(int id, int type){
	int addr_logic; int addr_real;
	int current_base = -1;
	if(type==OnlyO){
		current_base = o_baseAddr;  //for object
	}
	else if(type==OnlyU){
		current_base = u_baseAddr;  //for user
	}

	addr_logic= current_base  + id * sizeof(int);

	//char* BlockAddr=getFlatBlockForIdx(O2UIdxFile,addr_logic/BlkLen);
	//memcpy(Ref(addr_real),BlockAddr+(addr_logic%BlkLen), sizeof(int));
	getIdxBlockByjins(O2UIdxFile,addr_logic, sizeof(int),Ref(addr_real));
	return addr_real;

}

int getO2UOrgLeafID(int id, int type){
	int BaseAddr = getO2UAddrByIdxFile(id, type);
	FILE* f = O2ULeafFile;
	int addr = BaseAddr;  //基地址
	char* BlockAddr = getFlatBlock(f, addr/BlkLen);
	int e_id;
	memcpy(Ref(e_id),BlockAddr+(addr%BlkLen), sizeof(int));
	addr += sizeof(int);
	return e_id;
}


POI getPOIFromO2UOrgLeafData(int id){
	if(POIs[id].id >0)
		return  POIs[id];
	POI poi;FILE* f = O2ULeafFile;

	int BaseAddr = getO2UAddrByIdxFile(id, OnlyO);

	int addr = BaseAddr;  //获得POI数据块的基地址

	int poi_id;
	//char* BlockAddr = getFlatBlock(f, addr/BlkLen);
	//memcpy(Ref(poi_id),BlockAddr+(addr%BlkLen), sizeof(int));
	getDataBlockByjins(f,addr, sizeof(int),Ref(poi_id));
	poi.id = poi_id;
	addr += sizeof(int);

	//Ni
	int Ni=0;
	//BlockAddr = getFlatBlock(f, addr/BlkLen);
	//memcpy(Ref(Ni),BlockAddr+(addr%BlkLen), sizeof(int));
	getDataBlockByjins(f,addr, sizeof(int),Ref(Ni));
	poi.Ni = Ni;
	addr += sizeof(int);

	//Nj
	int Nj=0;
	//BlockAddr = getFlatBlock(f, addr/BlkLen);
	//memcpy(Ref(Nj),BlockAddr+(addr%BlkLen), sizeof(int));
	getDataBlockByjins(f,addr, sizeof(int),Ref(Nj));
	poi.Nj = Nj;
	addr += sizeof(int);

	//dist
	float dist = 0.0;
	//BlockAddr = getFlatBlock(f, addr/BlkLen);
	//memcpy(Ref(dist),BlockAddr+(addr%BlkLen), sizeof(float));
	getDataBlockByjins(f,addr, sizeof(float),Ref(dist));
	poi.dist= dist;
	addr += sizeof(float);

	//dis
	float dis = 0.0;
	//BlockAddr = getFlatBlock(f, addr/BlkLen);
	//memcpy(Ref(dis),BlockAddr+(addr%BlkLen), sizeof(float));
	getDataBlockByjins(f,addr, sizeof(float),Ref(dis));
	poi.dis = dis;
	addr += sizeof(float);

	int key_size = 0;
	//BlockAddr = getFlatBlock(f, addr/BlkLen);
	//memcpy(Ref(key_size),BlockAddr+(addr%BlkLen), sizeof(int));
	getDataBlockByjins(f,addr, sizeof(int),Ref(key_size));
	addr += sizeof(int);

	//依次读入各keyword_size
	for(int i=0;i<key_size;i++){
		//依次读入各keywords
		int keyword;
		//BlockAddr = getFlatBlock(f, addr/BlkLen);
		//memcpy(Ref(keyword),BlockAddr+(addr%BlkLen), sizeof(int));
		getDataBlockByjins(f,addr, sizeof(int),Ref(keyword));
		if(keyword>vocabularySize) break;
		poi.keywords.push_back(keyword);
        poi.keywordSet.insert(keyword);
		addr += sizeof(int);
	}

	int checkin_size = 0;
	//BlockAddr = getFlatBlock(f, addr/BlkLen);
	//memcpy(Ref(checkin_size),BlockAddr+(addr%BlkLen), sizeof(int));
	getDataBlockByjins(f,addr, sizeof(int),Ref(checkin_size));	
	addr += sizeof(int);

	for(int i=0; i<checkin_size;i++){
		//依次读入各check-in
		int check_in_usr;
		//BlockAddr = getFlatBlock(f, addr/BlkLen);
		//memcpy(Ref(check_in_usr),BlockAddr+(addr%BlkLen), sizeof(int));
		getDataBlockByjins(f,addr, sizeof(int),Ref(check_in_usr));
		poi.check_ins.push_back(check_in_usr);
		addr += sizeof(int);
	}
	int _size = POIs.size();
	POIs[id] = poi;
	return  poi;

}

User getUserFromO2UOrgLeafData(int id){
    if(Users[id].keywords.size()>0)
        return Users[id];
	User usr;

	int BaseAddr = getO2UAddrByIdxFile(id, OnlyU);
	FILE* f = O2ULeafFile;
	int addr = BaseAddr;  //获得POI数据块的基地址

	int usr_id;
	//char* BlockAddr = getFlatBlock(f, addr/BlkLen);
	//memcpy(Ref(usr_id),BlockAddr+(addr%BlkLen), sizeof(int));
	getDataBlockByjins(f,addr, sizeof(int),Ref(usr_id));
	usr.id = usr_id;
	addr += sizeof(int);

	//Ni
	int Ni=0;
	//BlockAddr = getFlatBlock(f, addr/BlkLen);
	//memcpy(Ref(Ni),BlockAddr+(addr%BlkLen), sizeof(int));
	getDataBlockByjins(f,addr, sizeof(int),Ref(Ni));
	usr.Ni = Ni;
	addr += sizeof(int);

	//Nj
	int Nj=0;
	//BlockAddr = getFlatBlock(f, addr/BlkLen);
	//memcpy(Ref(Nj),BlockAddr+(addr%BlkLen), sizeof(int));
	getDataBlockByjins(f,addr, sizeof(int),Ref(Nj));
	usr.Nj = Nj;
	addr += sizeof(int);

	//dist
	float dist = 0.0;
	//BlockAddr = getFlatBlock(f, addr/BlkLen);
	//memcpy(Ref(dist),BlockAddr+(addr%BlkLen), sizeof(float));
	getDataBlockByjins(f,addr, sizeof(float),Ref(dist));
	usr.dist= dist;
	addr += sizeof(float);

	//dis
	float dis = 0.0;
	//BlockAddr = getFlatBlock(f, addr/BlkLen);
	//memcpy(Ref(dis),BlockAddr+(addr%BlkLen), sizeof(float));
	getDataBlockByjins(f,addr, sizeof(float),Ref(dis));
	usr.dis = dis;
	addr += sizeof(float);

	int key_size = 0;
	//BlockAddr = getFlatBlock(f, addr/BlkLen);
	//memcpy(Ref(key_size),BlockAddr+(addr%BlkLen), sizeof(int));
	getDataBlockByjins(f,addr, sizeof(int),Ref(key_size));
	addr += sizeof(int);

	//依次读入各keyword_size
	for(int i=0;i<key_size;i++){
		//依次读入各keywords
		int keyword;
		//BlockAddr = getFlatBlock(f, addr/BlkLen);
		//memcpy(Ref(keyword),BlockAddr+(addr%BlkLen), sizeof(int));
		getDataBlockByjins(f,addr, sizeof(int),Ref(keyword));
        if(keyword>vocabularySize) break;
		usr.keywords.push_back(keyword);
		addr += sizeof(int);
	}

	int friend_size = 0;
	//BlockAddr = getFlatBlock(f, addr/BlkLen);
	//memcpy(Ref(friend_size),BlockAddr+(addr%BlkLen), sizeof(int));
	getDataBlockByjins(f,addr, sizeof(int),Ref(friend_size));
	addr += sizeof(int);

	for(int i=0; i<friend_size;i++){
		//依次读入各check-in
		int friend_usr;
		//BlockAddr = getFlatBlock(f, addr/BlkLen);
		//memcpy(Ref(friend_usr),BlockAddr+(addr%BlkLen), sizeof(int));
		getDataBlockByjins(f,addr, sizeof(int),Ref(friend_usr));
		usr.friends.push_back(friend_usr);
		addr += sizeof(int);
	}
    Users[id] = usr;
	return  usr;

}



#define   TRACKTFIDF1

//this is for inverted post list file of poi as well as gtree node
float getTermIDFWeight(int term){

#ifdef TRACKTFIDF
		cout<<"获取t"<<term<<"的weight"<<endl;
#endif
	if(termWeight_Map[term]>0){

#ifdef TRACKTFIDF
		cout<<"t"<<term<<"的weight当前已被缓存！"<<endl;
#endif
		return termWeight_Map[term];
	}



	//int addr_logic = 2*sizeof(int)+term*(sizeof(int)+sizeof(int));
	if(!term>0){
		cout<<"term id error!"<<endl;
		exit(-1);
	}
	int term_record_size = 7;  //size: term_base, term_obase, term_ubase, term_oleafaddr, term_uleafaddr, term_object, term_user
	int addr_logic = sizeof(int)+(term-1)*(sizeof(float)+term_record_size*sizeof(int));  //term_weight + address_count(7)

	int addr_base =0 ; int addr_leaf = -1;

	//读取term的IDF weight
	float term_weight = -1.0;
	//char* BlockAddr = getFlatBlockForIdx(O2UInvFile,addr_logic/BlkLen);
	//memcpy(Ref(term_weight),BlockAddr+(addr_logic%BlkLen), sizeof(float));
	getIdxBlockByjins(O2UInvFile,addr_logic, sizeof(float),Ref(term_weight));
	addr_logic += sizeof(float);


    termWeight_Map[term]=term_weight;
	return term_weight;
}

int term_size = -1;
//得到当前数据集文本关键词词汇（vocabulary)的个数
float getTermSize(){

	if(term_size >0)
		return term_size;

	
	int addr_logic = 0;  //term_weight + address_count(7)
			

	//读取term的IDF weight
	int _size = -1.0;
	
	getIdxBlockByjins(O2UInvFile,addr_logic, sizeof(int),Ref(_size));
	addr_logic += sizeof(float);

	term_size = _size;


	return term_size;
}


vector<int> getTermOjectInvList(int term){


	vector<int> posting_list;
	//int addr_logic = 2*sizeof(int)+term*(sizeof(int)+sizeof(int));
	if(!term>0){
		cout<<"term id error!"<<endl;
		exit(-1);
	}
	int term_record_size = 7;  //size: term_base, term_obase, term_ubase, term_oleafaddr, term_uleafaddr, term_object, term_user
	int addr_start = sizeof(int)+(term-1)*(sizeof(float)+term_record_size*sizeof(int));  //term_size, record_size: term_weight + address_count(7)

	int addr_logic = addr_start+ sizeof(float)+5*sizeof(int);

	int addr_base =0 ; int addr_leaf = -1;

	//读取term对应的（object）posting_list 的地址
	int term_object_inv_addr = -1.0;
	getIdxBlockByjins(O2UInvFile,addr_logic, sizeof(int),Ref(term_object_inv_addr));
	addr_logic = term_object_inv_addr;

	//取posting_list size大小

	int object_size = -1;
	//cout<<"从地址"<<addr_logic<<"处数据块，读取object_size..."<<endl;
	getIdxBlockByjins(O2UInvFile,addr_logic, sizeof(int),Ref(object_size));
	addr_logic += sizeof(int);
	//cout<<"------------------------------读到 object_size="<<object_size<<endl;

	for(int i=0;i<object_size;i++){
		int o_id = -1;
		//cout<<"从地址"<<addr_logic<<"处数据块，读取o_id..."<<endl;
		getIdxBlockByjins(O2UInvFile,addr_logic, sizeof(int),Ref(o_id));
		//cout<<"------------------------------读到 o"<<o_id<<endl;
		posting_list.push_back(o_id);
		addr_logic += sizeof(int);
	}
	//getchar();

	return posting_list;
}

int getTermOjectInvListSize(int term){


	vector<int> posting_list;
	//int addr_logic = 2*sizeof(int)+term*(sizeof(int)+sizeof(int));
	if(!term>0){
		cout<<"t"<<term<<",term id error!"<<endl;
		exit(-1);
	}
	int term_record_size = 7;  //size: term_base, term_obase, term_ubase, term_oleafaddr, term_uleafaddr, term_object, term_user
	int addr_start = sizeof(int)+(term-1)*(sizeof(float)+term_record_size*sizeof(int));  //term_size, record_size: term_weight + address_count(7)

	int addr_logic = addr_start+ sizeof(float)+5*sizeof(int);

	int addr_base =0 ; int addr_leaf = -1;

	//读取term对应的（object）posting_list 的地址
	int term_object_inv_addr = -1.0;
	getIdxBlockByjins(O2UInvFile,addr_logic, sizeof(int),Ref(term_object_inv_addr));
	addr_logic = term_object_inv_addr;

	//取posting_list size大小

	int object_size = -1;
	//cout<<"从地址"<<addr_logic<<"处数据块，读取object_size..."<<endl;
	getIdxBlockByjins(O2UInvFile,addr_logic, sizeof(int),Ref(object_size));
	addr_logic += sizeof(int);
	//cout<<"------------------------------读到 object_size="<<object_size<<endl;
	return object_size;
}


vector<int> getTermUserInvList(int term){


	vector<int> posting_list;
	//int addr_logic = 2*sizeof(int)+term*(sizeof(int)+sizeof(int));
	if(!term>0){
		cout<<"term id error!"<<endl;
		exit(-1);
	}
	int term_record_size = 7;  //size: term_base, term_obase, term_ubase, term_oleafaddr, term_uleafaddr, term_object, term_user
	int addr_start = sizeof(int)+(term-1)*(sizeof(float)+term_record_size*sizeof(int));  //term_size, record_size: term_weight + address_count(7)

	int addr_logic = addr_start+ sizeof(float)+6*sizeof(int);

	int addr_base =0 ; int addr_leaf = -1;

	//读取term对应的（object）posting_list 的地址
	int term_user_inv_addr = -1.0;
	getIdxBlockByjins(O2UInvFile,addr_logic, sizeof(int),Ref(term_user_inv_addr));
	addr_logic = term_user_inv_addr;

	//取posting_list size大小

	int user_size = -1;
	//cout<<"从地址"<<addr_logic<<"处数据块，读取object_size..."<<endl;
	getIdxBlockByjins(O2UInvFile,addr_logic, sizeof(int),Ref(user_size));
	addr_logic += sizeof(int);
	//cout<<"------------------------------读到 object_size="<<object_size<<endl;

	for(int i=0;i<user_size;i++){
		int u_id = -1;
		//cout<<"从地址"<<addr_logic<<"处数据块，读取o_id..."<<endl;
		getIdxBlockByjins(O2UInvFile,addr_logic, sizeof(int),Ref(u_id));
		//cout<<"------------------------------读到 o"<<o_id<<endl;
		posting_list.push_back(u_id);
		addr_logic += sizeof(int);
	}
	//getchar();

	return posting_list;
}




//below for textual similarity computation

//TF-IDF similarity


double tfIDF_term(int term){
    double a = 0.0; double b = 0.0;
    float idfScore = -1;
#ifndef DiskAccess
    int poi_size = POIs.size();
    vector<int> poi_list = invListOfPOI[term];
    idfScore = log2(poi_size*1.0/poi_list.size());
#else
    idfScore = getTermIDFWeight(term);
#endif
    a += idfScore;
    return  a;
}


double tfIDF_termSet(set<int> termSet){
	double a = 0.0; double b = 0.0;
	int poi_size = POIs.size();

	for(int term:termSet) {
		double idfScore = 0;
		idfScore = tfIDF_term(term);
		a += idfScore;
	}
	return  a;
}



double tfIDF(set<int> Qkey, set<int> Pkey){
    double a = 0.0; double b = 0.0;
    int poi_size = POIs.size();

    set<int> interKey ;//= obtain_itersection_jins(Qkey,Pkey);
    set_intersection(Qkey.begin(),Qkey.end(),Pkey.begin(),Pkey.end(),inserter(interKey,interKey.begin()));
    if(interKey.size()>0){
        for(int term:interKey){

            double idfScore = tfIDF_term(term);

            a += idfScore;
        }
        return  a;
    }
    else return 0;

}


double textRelevance(vector<int> Ukey, vector<int> Pkey){
    set<int> Ukeyset; set<int>  Pkeyset;
    for (int term :Ukey)
    	Ukeyset.insert(term);
    for (int term :Pkey)
    	Pkeyset.insert(term);
    return  tfIDF(Ukeyset,Pkeyset);
    //take into account the ratio of overlap keywords
}

//TF-IDF similarity

double textRelevanceSet(set<int> Ukeyset, set<int>  Pkeyset){

    return  tfIDF(Ukeyset,Pkeyset);
    //take into account the ratio of overlap keywords
}





//返回term对应的索引表中地址内容<baseAddress,
						  // obaseAddr, oLeafAddr
						  // ubaseAddr, uLeafAddr>

//p2u_invlist文件中索引列表字段的每行内容
typedef struct {
	int base;
	int obase;
	int oleaf;
	int ubase;
	int uleaf;

}IdxBundle;

//不建议用，会出错
IdxBundle getO2UTermInvPostAddr_old(int term){
	//int addr_logic = term_size+ term_size*(每行内容：term_weight, termbase, obase, oleaf, ubase, uleaf);
	int addr_logic = sizeof(int)+(term-1)*(sizeof(float)+5*sizeof(int));


	//读取term的IDF weight
	float term_weight = -1.0;
	char* BlockAddr = getFlatBlockForIdx(O2UInvFile,addr_logic/BlkLen);
	memcpy(Ref(term_weight),BlockAddr+(addr_logic%BlkLen), sizeof(int));
	termWeight_Map[term]=term_weight;  //cache term weight
	addr_logic += sizeof(float);

	double idfScore = log2(POIs.size()*1.0/invListOfPOI[term].size());



	//读取term的倒排列表内容基地址
	int addr_base =0 ;
	BlockAddr = getFlatBlockForIdx(O2UInvFile,addr_logic/BlkLen);
	memcpy(Ref(addr_base),BlockAddr+(addr_logic%BlkLen), sizeof(int));
	addr_logic += sizeof(int);

	//读取object entry的倒排表内容地址
	int addr_obase = -1;
	BlockAddr = getFlatBlockForIdx(O2UInvFile,addr_logic/BlkLen);
	memcpy(Ref(addr_obase),BlockAddr+(addr_logic%BlkLen), sizeof(int));
	addr_logic += sizeof(int);

	//读取object_leaf 表的地址
	int addr_oleaf = -1;
	BlockAddr = getFlatBlockForIdx(O2UInvFile,addr_logic/BlkLen);
	memcpy(Ref(addr_oleaf),BlockAddr+(addr_logic%BlkLen), sizeof(int));
	addr_logic += sizeof(int);


	//读取user entry的倒排表内容地址
	int addr_ubase =-1;
	BlockAddr = getFlatBlockForIdx(O2UInvFile,addr_logic/BlkLen);
	memcpy(Ref(addr_ubase),BlockAddr+(addr_logic%BlkLen), sizeof(int));
	addr_logic += sizeof(int);

	//读取user_leaf 表的地址
	int addr_uleaf = -1;
	BlockAddr = getFlatBlockForIdx(O2UInvFile,addr_logic/BlkLen);
	memcpy(Ref(addr_uleaf),BlockAddr+(addr_logic%BlkLen), sizeof(int));
	addr_logic += sizeof(int);

	//地址赋值
	IdxBundle idxBundle;
	idxBundle.base = addr_base;
	idxBundle.obase = addr_obase;
	idxBundle.oleaf = addr_oleaf;
	idxBundle.ubase = addr_ubase;
	idxBundle.uleaf = addr_uleaf;

	return idxBundle;
}


//获得单词倒排列表索引表的信息
IdxBundle getO2UTermInvPostAddr(int term){
	//int addr_logic = term_size+ term_size*(每行内容：term_weight, termbase, obase, oleaf, ubase, uleaf);
	int addr_logic = sizeof(int)+(term-1)*(sizeof(float)+inv_idx_record_size*sizeof(int));


	//读取term的IDF weight
	float term_weight = -1.0;
	getIdxBlockByjins(O2UInvFile,addr_logic, sizeof(float),Ref(term_weight));
	termWeight_Map[term]=term_weight;  //cache term weight
	addr_logic += sizeof(float);

	//double idfScore = log2(POIs.size()*1.0/invListOfPOI[term].size());


	//读取term的倒排列表内容基地址
	int addr_base =0 ;
	getIdxBlockByjins(O2UInvFile,addr_logic, sizeof(int), Ref(addr_base));
	addr_logic += sizeof(int);

	//读取object entry的倒排表内容地址
	int addr_obase = -1;
	getIdxBlockByjins(O2UInvFile,addr_logic, sizeof(int),Ref(addr_obase));
	addr_logic += sizeof(int);

	//读取object_leaf 表的地址
	int addr_oleaf = -1;
	getIdxBlockByjins(O2UInvFile,addr_logic, sizeof(int),Ref(addr_oleaf));
	addr_logic += sizeof(int);


	//读取user entry的倒排表内容地址
	int addr_ubase =-1;
	getIdxBlockByjins(O2UInvFile,addr_logic, sizeof(int),Ref(addr_ubase));
	addr_logic += sizeof(int);

	//读取user_leaf 表的地址
	int addr_uleaf = -1;
	getIdxBlockByjins(O2UInvFile,addr_logic, sizeof(int), Ref(addr_uleaf));
	addr_logic += sizeof(int);

	//地址赋值
	IdxBundle idxBundle;
	idxBundle.base = addr_base;
	idxBundle.obase = addr_obase;
	idxBundle.oleaf = addr_oleaf;
	idxBundle.ubase = addr_ubase;
	idxBundle.uleaf = addr_uleaf;

	return idxBundle;
}



vector<int> getObjectTermRelatedLeaves(int term){

	int BaseAddr = getO2UTermInvPostAddr(term).oleaf;


	FILE* f = O2UInvFile; vector<int> related_leaves;
	int addr = BaseAddr;  //基地址
	//读入leaf size
	int related_leaf_size=0;
	char* BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
	memcpy(Ref(related_leaf_size),BlockAddr+(addr%BlkLen), sizeof(int));
	addr+= sizeof(int);
	//读入各个相关leaf node的id
	for(int i=0;i<related_leaf_size;i++){
		int leaf_id = -1;
		char* BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
		memcpy(Ref(leaf_id),BlockAddr+(addr%BlkLen), sizeof(int));
		addr+= sizeof(int);
		related_leaves.push_back(leaf_id);
	}
	return related_leaves;
}

vector<int> getUsrTermRelatedLeaves(int term){

	int BaseAddr = getO2UTermInvPostAddr(term).uleaf;


	FILE* f = O2UInvFile; vector<int> related_leaves;
	int addr = BaseAddr;  //基地址
	//读入leaf size
	int related_leaf_size=0;
	//char* BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
	//memcpy(Ref(related_leaf_size),BlockAddr+(addr%BlkLen), sizeof(int));
	getIdxBlockByjins(f,addr, sizeof(int),Ref(related_leaf_size));
	addr+= sizeof(int);
	//读入各个相关leaf node的id
	for(int i=0;i<related_leaf_size;i++){
		int leaf_id = -1;
		//char* BlockAddr = getFlatBlockForIdx(f, addr/BlkLen);
		//memcpy(Ref(leaf_id),BlockAddr+(addr%BlkLen), sizeof(int));
		getIdxBlockByjins(f,addr, sizeof(int),Ref(leaf_id));
		addr+= sizeof(int);
		related_leaves.push_back(leaf_id);
	}
	return related_leaves;
}



//返回object node 带有term的孩子节点的列表

//不建议
vector<int> getObjectTermRelatedEntry_slow(int term, int node){

	TreeNode n = getGIMTreeNodeData(node,OnlyO);

	vector<int> related_entry;
	int BaseAddr = getO2UTermInvPostAddr(term).obase;  //取得object term对应的列表基地址

	if(n.term_innerAddr_OMap[term].size()==0){
		return related_entry;
	}

	int innerAddr = n.term_innerAddr_OMap[term][0];   //取得该node在term列表中的地址偏移
	//int innerAddr = n.term_innerAddr_OMap[term][0];   //取得该node在term列表中的地址偏移
	FILE* f = O2UInvFile; char* data;

	int addr = BaseAddr+innerAddr;  //合成基地址

	int entry_size=0; //读取entry 的个数
	getDataBlockByjins(f,addr, sizeof(int),Ref(entry_size));
	addr+= sizeof(int);

	//读入各个相关leaf node的id
	for(int i=0;i<entry_size;i++){
		int entry_id = -1;
		getDataBlockByjins(f,addr, sizeof(int),Ref(entry_id));
		addr+= sizeof(int);
		related_entry.push_back(entry_id);
	}
	return related_entry;
}

//建议
vector<int> getObjectTermRelatedEntry(int term, int node){

	TreeNode n = getGIMTreeNodeData_fast(node,OnlyO,term);

	vector<int> related_entry;
	int BaseAddr = getO2UTermInvPostAddr(term).obase;  //取得object term对应的列表基地址

	if(n.term_innerAddr_OMap[term].size()==0){
		return related_entry;
	}

	int innerAddr = n.term_innerAddr_OMap[term][0];   //取得该node在term列表中的地址偏移
	//int innerAddr = n.term_innerAddr_OMap[term][0];   //取得该node在term列表中的地址偏移
	FILE* f = O2UInvFile; char* data;

	int addr = BaseAddr+innerAddr;  //合成基地址

	int entry_size=0; //读取entry 的个数
	getIdxBlockByjins(f,addr, sizeof(int),Ref(entry_size));
	addr+= sizeof(int);

	//读入各个相关leaf node的id
	for(int i=0;i<entry_size;i++){
		int entry_id = -1;
		getIdxBlockByjins(f,addr, sizeof(int),Ref(entry_id));
		addr+= sizeof(int);
		related_entry.push_back(entry_id);
	}
	return related_entry;
}



//返回user带有term的孩子节点的列表


vector<int> getUsrTermRelatedEntry(int term, int node){

	TreeNode n = getGIMTreeNodeData(node,OnlyU);

	vector<int> related_entry;
	int BaseAddr = getO2UTermInvPostAddr(term).ubase;;//取得该term对应的 user entry 基地址

	if(n.term_innerAddr_UMap[term].size()==0){  //gangc
		return related_entry;
	}

	int innerAddr = n.term_innerAddr_UMap[term][0]; //取得该node在term列表中的地址偏移
	FILE* f = O2UInvFile; char* data;
	int addr = BaseAddr+innerAddr;  //基地址

	//读入entry size
	int entry_size=0;
	getIdxBlockByjins(f,addr, sizeof(int),Ref(entry_size));
	addr+= sizeof(int);

	//读入各个相关leaf node的id
	for(int i=0;i<entry_size;i++){
		int entry_id = -1;
		getIdxBlockByjins(f, addr, sizeof(int), Ref(entry_id));
		addr+= sizeof(int);
		related_entry.push_back(entry_id);
	}
	return related_entry;
}



//返回user带有term的底层叶子节点的列表

vector<int> getUsrTermRelatedLeafNode(int term, int node){

	TreeNode n = getGIMTreeNodeData(node,OnlyU);

	vector<int> related_entry; vector<int> related_leaves;
	int BaseAddr = getO2UTermInvPostAddr(term).ubase;;//取得该term对应的 user entry 基地址

	if(n.term_innerAddr_UMap[term].size()==0){  //gangc
		return related_entry;
	}

	int innerAddr = n.term_innerAddr_UMap[term][0]; //取得该node在term列表中的地址偏移
	FILE* f = O2UInvFile; char* data;
	int addr = BaseAddr+innerAddr;  //基地址

	//读入entry size
	int entry_size=0;
	getIdxBlockByjins(f, addr, sizeof(int),Ref(entry_size));
	addr+= sizeof(int);
	//读入各个相关entry node的id
	for(int i=0;i<entry_size;i++){
		int entry_id = -1;
		getIdxBlockByjins(f, addr, sizeof(int),Ref(entry_id));
		addr+= sizeof(int);
		related_entry.push_back(entry_id);
	}

	//读入related leaf size
	int leaf_size=0;
	getIdxBlockByjins(f, addr, sizeof(int),Ref(leaf_size));
	addr+= sizeof(int);
	//读入各个相关leaf node的id
	for(int i=0;i<leaf_size;i++){
		int leaf_id = -1;
		getIdxBlockByjins(f, addr, sizeof(int),Ref(leaf_id));
		addr+= sizeof(int);
		related_leaves.push_back(leaf_id);
	}

	return related_leaves;
}



void OpenDiskComm(const char* fileprefix,int _cachesize) {
	char tmpFileName[255];
	BlkLen=getBlockLength();
	int header_size = sizeof(char) + sizeof(int);
	gBuffer=new char[BlkLen];
	i_capacity=(BlkLen-header_size)/(sizeof(int)+sizeof(int));
	printf("Blocklength=%d, i_cap=%d\n",BlkLen,i_capacity);

	//重要重要！
	InitFreqCache(FC_A); //邻接表文件的Frequent Cache
	InitFreqCache(FC_P); //按邻边关系组织的兴趣点数据集文件的Frequent Cache
	InitFreqCache(FC_GN); // Gtree 节点索引文件的Frequent Cache
	InitFreqCache(FC_O2UL);  //按叶节点关系组织双色体对象数据文件的Frequent Cache
	InitFreqCache(FC_O2UIdx); //(poi index for the address in poi data)//对应O2UL文件中双色体对象（用户或兴趣点）地址的索引文件的Frequent Cache
	InitFreqCache(FC_O2UIV);  // 双色体对象关键词信息倒排索引文件的Frequent Cache




	InitCache(_cachesize);

	//open POI data 1 file

	sprintf(tmpFileName,"%s.poi_d",fileprefix);
	PtFile=fopen(tmpFileName,"r");
	CheckFile(PtFile,tmpFileName);
	PtFileSize=getFileSize(PtFile);
	cout<<".p_d(poi data) file open !"<<endl;        //--------remove

	/*//open poi btree index
	sprintf(tmpFileName,"%s.poi_bt",fileprefix);
	PtTree = new BTree(tmpFileName, 128); // number of cached file entries: 128
	PtNum = PtTree->UserField;                         // UserField 代表什么？
	cout<<".p_bt(poi index) file open !"<<endl;
*/

#ifdef datanew
	//open bichromatic data file
	sprintf(tmpFileName,"%s.o2u_ld",fileprefix);
	O2ULeafFile=fopen(tmpFileName,"r");
	CheckFile(O2ULeafFile,tmpFileName);
	O2ULFileSize=getFileSize(O2ULeafFile);
	cout<<".o2u_d(bichromatic data new) file open !"<<endl;

	//open bichromatic address index file
	sprintf(tmpFileName,"%s.o2u_idx",fileprefix);
	O2UIdxFile =fopen(tmpFileName,"r");
	CheckFile(O2UIdxFile,tmpFileName);
	cout<<".o2u_idx(bichromatic addr index) file open !"<<endl;        //--------remove
	O2UIdxFileSize=getFileSize(O2UIdxFile);   //huilai
	getBichromaticIdx();
#endif


	//open o2u inverted list file
	sprintf(tmpFileName,"%s.o2u_invList",fileprefix);
	O2UInvFile=fopen(tmpFileName,"r");
	CheckFile(O2UInvFile,tmpFileName);
	cout<<"双色体(倒排索引)posting list 文件 open!"<<endl;
	O2UIVFileSize = getFileSize(O2UInvFile);



	//open gtree index file
	sprintf(tmpFileName,"%s-disk",FILE_GTREE);
	GtreeNodeFile =fopen(tmpFileName,"r");
	CheckFile(GtreeNodeFile,tmpFileName);
	GNFileSize= getFileSize(GtreeNodeFile);
	cout<<"gtree 索引文件 open!"<<endl;

	//open adjacent list file
	sprintf(tmpFileName,"%s.adj",fileprefix);
	AdjFile=fopen(tmpFileName,"r");
	CheckFile(AdjFile,tmpFileName);
	cout<<"adj(邻接表)文件 open!"<<endl;
	int num = VertexNum;
	fread(Ref(num),1,sizeof(int),AdjFile);	//	VertexNum=AdjTree->UserField;
	AdjFileSize=getFileSize(AdjFile);
	printf("PtFileSize: %d, AdjFileSize: %d, VertexNum: %d, PtNum: %d\n",PtFileSize,AdjFileSize,VertexNum,PtNum);

	//getchar();
}

void OpenDiskComm_Plus_NVD(const char* fileprefix,int _cachesize) {
	char tmpFileName[255];
	BlkLen=getBlockLength();
	int header_size = sizeof(char) + sizeof(int);
	gBuffer=new char[BlkLen];
	i_capacity=(BlkLen-header_size)/(sizeof(int)+sizeof(int));
	printf("Blocklength=%d, i_cap=%d\n",BlkLen,i_capacity);

	//重要重要！
	InitFreqCache(FC_A); //邻接表文件的Frequent Cache
	InitFreqCache(FC_P); //按邻边关系组织的兴趣点数据集文件的Frequent Cache
	InitFreqCache(FC_GN); // Gtree 节点索引文件的Frequent Cache
	InitFreqCache(FC_O2UL);  //按叶节点关系组织双色体对象数据文件的Frequent Cache
	InitFreqCache(FC_O2UIdx); //(poi index for the address in poi data)//对应O2UL文件中双色体对象（用户或兴趣点）地址的索引文件的Frequent Cache
	InitFreqCache(FC_O2UIV);  // 双色体对象关键词信息倒排索引文件的Frequent Cache

	//for nvd
	InitFreqCache(FC_v2p);
	InitFreqCache(FC_NVD_ADJGraph);
	/*InitFreqCache(FC_v2pIdx);
	InitFreqCache(FC_NVD_ADJGraphIdx);*/


	InitCache(_cachesize);

	//open POI data 1 file

	sprintf(tmpFileName,"%s.poi_d",fileprefix);
	PtFile=fopen(tmpFileName,"r");
	CheckFile(PtFile,tmpFileName);
	PtFileSize=getFileSize(PtFile);
	cout<<".p_d(poi data) file open !"<<endl;        //--------remove
	


#ifdef datanew
	//open bichromatic data file
	sprintf(tmpFileName,"%s.o2u_ld",fileprefix);
	O2ULeafFile=fopen(tmpFileName,"r");
	CheckFile(O2ULeafFile,tmpFileName);
	O2ULFileSize=getFileSize(O2ULeafFile);
	cout<<".o2u_d(bichromatic data new) file open !"<<endl;

	//open bichromatic address index file
	sprintf(tmpFileName,"%s.o2u_idx",fileprefix);
	O2UIdxFile =fopen(tmpFileName,"r");
	CheckFile(O2UIdxFile,tmpFileName);
	cout<<".o2u_idx(bichromatic addr index) file open !"<<endl;        //--------remove
	O2UIdxFileSize=getFileSize(O2UIdxFile);   //huilai
	getBichromaticIdx();
#endif


	//open o2u inverted list file
	sprintf(tmpFileName,"%s.o2u_invList",fileprefix);
	O2UInvFile=fopen(tmpFileName,"r");
	CheckFile(O2UInvFile,tmpFileName);
	cout<<"双色体(倒排索引)posting list 文件 open!"<<endl;
	O2UIVFileSize = getFileSize(O2UInvFile);



	//open gtree index file
	sprintf(tmpFileName,"%s-disk",FILE_GTREE);
	GtreeNodeFile =fopen(tmpFileName,"r");
	CheckFile(GtreeNodeFile,tmpFileName);
	GNFileSize= getFileSize(GtreeNodeFile);
	cout<<"gtree 索引文件 open!"<<endl;

	//open adjacent list file
	sprintf(tmpFileName,"%s.adj",fileprefix);
	AdjFile=fopen(tmpFileName,"r");
	CheckFile(AdjFile,tmpFileName);
	cout<<"adj(邻接表)文件 open!"<<endl;
	int num = VertexNum;
	fread(Ref(num),1,sizeof(int),AdjFile);	//	VertexNum=AdjTree->UserField;
	AdjFileSize=getFileSize(AdjFile);
	printf("PtFileSize: %d, AdjFileSize: %d, VertexNum: %d, PtNum: %d\n",PtFileSize,AdjFileSize,VertexNum,PtNum);

    //this is file address for NVD document

    char NVD_fileName[255]; char Road_fileName[255];
    char pre_fileName[255];

    sprintf(pre_fileName,"../../../data/%s/NVD/%s",dataset_name,road_data);

    //v2p address映射文件现在取消不用了
    /*sprintf(tmpFileName,"%s.v2p",pre_fileName);
	VertexHashFile =fopen(tmpFileName,"r");
	CheckFile(VertexHashFile,tmpFileName);
	V2PFileSize = getFileSize(VertexHashFile);
	cout<<"v2p(vertex poi哈希映射)文件 open!"<<endl;*/

	//open NVD graph file
	sprintf(tmpFileName,"%s.nvd",pre_fileName);
	NVD_ADJGraph = fopen(tmpFileName,"r");
	CheckFile(NVD_ADJGraph, tmpFileName);
	NVDGFileSize = getFileSize(NVD_ADJGraph);
	cout<<"nvd(Network Voronoi Diagram)文件 open!"<<endl;

	

	//getchar();
}


void OpenDiskComm_Plus_NVD_varyingO(const char* fileprefix,int _cachesize,int objectType, float rate) {
    char tmpFileName[255];
    BlkLen=getBlockLength();
    int header_size = sizeof(char) + sizeof(int);
    gBuffer=new char[BlkLen];
    i_capacity=(BlkLen-header_size)/(sizeof(int)+sizeof(int));
    printf("Blocklength=%d, i_cap=%d\n",BlkLen,i_capacity);

    //重要重要！
    InitFreqCache(FC_A); //邻接表文件的Frequent Cache
    InitFreqCache(FC_P); //按邻边关系组织的兴趣点数据集文件的Frequent Cache
    InitFreqCache(FC_GN); // Gtree 节点索引文件的Frequent Cache
    InitFreqCache(FC_O2UL);  //按叶节点关系组织双色体对象数据文件的Frequent Cache
    InitFreqCache(FC_O2UIdx); //(poi index for the address in poi data)//对应O2UL文件中双色体对象（用户或兴趣点）地址的索引文件的Frequent Cache
    InitFreqCache(FC_O2UIV);  // 双色体对象关键词信息倒排索引文件的Frequent Cache

    //for nvd
    InitFreqCache(FC_v2p);
    InitFreqCache(FC_NVD_ADJGraph);
    /*InitFreqCache(FC_v2pIdx);
    InitFreqCache(FC_NVD_ADJGraphIdx);*/


    InitCache(_cachesize);

    //open POI data 1 file

    sprintf(tmpFileName,"%s.poi_d",fileprefix);
    PtFile=fopen(tmpFileName,"r");
    CheckFile(PtFile,tmpFileName);
    PtFileSize=getFileSize(PtFile);
    cout<<".p_d(poi data) file open !"<<endl;        //--------remove



#ifdef datanew
    //open bichromatic data file
    sprintf(tmpFileName,"%s.o2u_ld",fileprefix);
    O2ULeafFile=fopen(tmpFileName,"r");
    CheckFile(O2ULeafFile,tmpFileName);
    O2ULFileSize=getFileSize(O2ULeafFile);
    cout<<".o2u_d(bichromatic data new) file open !"<<endl;

    //open bichromatic address index file
    sprintf(tmpFileName,"%s.o2u_idx",fileprefix);
    O2UIdxFile =fopen(tmpFileName,"r");
    CheckFile(O2UIdxFile,tmpFileName);
    cout<<".o2u_idx(bichromatic addr index) file open !"<<endl;        //--------remove
    O2UIdxFileSize=getFileSize(O2UIdxFile);   //huilai
    getBichromaticIdx();
#endif


    //open o2u inverted list file
    sprintf(tmpFileName,"%s.o2u_invList",fileprefix);
    O2UInvFile=fopen(tmpFileName,"r");
    CheckFile(O2UInvFile,tmpFileName);
    cout<<"双色体(倒排索引)posting list 文件 open!"<<endl;
    O2UIVFileSize = getFileSize(O2UInvFile);



    //open gtree index file
    sprintf(tmpFileName,"%s.gtree-disk",fileprefix);
    GtreeNodeFile =fopen(tmpFileName,"r");
    CheckFile(GtreeNodeFile,tmpFileName);
    GNFileSize= getFileSize(GtreeNodeFile);
    cout<<"gtree 索引文件 open!"<<endl;

    //open adjacent list file
    sprintf(tmpFileName,"%s.adj",fileprefix);
    AdjFile=fopen(tmpFileName,"r");
    CheckFile(AdjFile,tmpFileName);
    cout<<"adj(邻接表)文件 open!"<<endl;
    int num = VertexNum;
    fread(Ref(num),1,sizeof(int),AdjFile);	//	VertexNum=AdjTree->UserField;
    AdjFileSize=getFileSize(AdjFile);
    printf("PtFileSize: %d, AdjFileSize: %d, VertexNum: %d, PtNum: %d\n",PtFileSize,AdjFileSize,VertexNum,PtNum);

    //this is file address for NVD document

    char NVD_fileName[255]; char Road_fileName[255];
    char pre_fileName[255];

    if(objectType==USER){
        sprintf(pre_fileName,"../../../data/%s/NVD/%s",dataset_name,road_data);
    }
    else if(objectType==BUSINESS){
        sprintf(pre_fileName,"../../../data/%s/NVD/VaryingO/%d%/%s",dataset_name,(int)(rate*100),road_data);
    }


    //v2p address映射文件现在取消不用了
    /*sprintf(tmpFileName,"%s.v2p",pre_fileName);
	VertexHashFile =fopen(tmpFileName,"r");
	CheckFile(VertexHashFile,tmpFileName);
	V2PFileSize = getFileSize(VertexHashFile);
	cout<<"v2p(vertex poi哈希映射)文件 open!"<<endl;*/

    //open NVD graph file
    sprintf(tmpFileName,"%s.nvd",pre_fileName);
    NVD_ADJGraph = fopen(tmpFileName,"r");
    CheckFile(NVD_ADJGraph, tmpFileName);
    NVDGFileSize = getFileSize(NVD_ADJGraph);
    cout<<"nvd(Network Voronoi Diagram)文件 open!"<<endl;



    //getchar();
}



void CacheReflash(){

	BlkLen=getBlockLength();
	int header_size = sizeof(char) + sizeof(int);
	ConfigType setting;
	AddConfigFromFile(setting,"../map_exp");
	//AddConfigFromCmdLine(setting,argc,argv);
	//ListConfig(setting);
	//设定输入数据的文件目录
	string filename="../../../data/";
	filename=filename+getConfigStr("map",setting)+"/"+getConfigStr("data",setting);
	const char* fileprefix=filename.c_str();
	//获取参数信息中的cache页数
	int _cachepages_size = getConfigInt("cachepages",setting,false,DEFAULT_CACHESIZE);
	gBuffer=new char[BlkLen];
	//i_capacity=(BlkLen-header_size)/(sizeof(int)+sizeof(int));


	InitFreqCache(FC_A);
	InitFreqCache(FC_P);
	InitFreqCache(FC_O2UIdx);
	InitFreqCache(FC_O2UIV);
	InitFreqCache(FC_O2UL);
	InitCache(_cachepages_size);
	for(User u :Users){
		u.topkScore_current =-1.0;
		u.topkScore_Final =-1.0;
	}
}


void CloseDiskComm() {
	fclose(PtFile);		fclose(AdjFile);
	delete PtTree;
	delete[] gBuffer;

	//printPageAccess();
	DestroyCache();
	DestroyFreqCache(FC_A);		DestroyFreqCache(FC_P);
	DestroyFreqCache(FC_GN);
	DestroyFreqCache(FC_O2UIdx);
	DestroyFreqCache(FC_O2UIV);
	DestroyFreqCache(FC_O2UL);

	cout<<"CloseDiskComm COMPLETE!"<<endl;
}

void printSubTree(BTree* bt,NodeStruct* node,int indent) {
	char space[indent+1];
	for (int i=0;i<indent;i++) space[i]=' ';
	space[indent]='\0';

	int curLevel=(int)(node->level);
	printf("%sLevel: %d\n",space,curLevel);

	NodeStruct* cNode=createNodeStruct();
	for(int i=0;i<node->num_entries;i++) {
		printf("%s(%d,%d)\n",space,node->key[i],node->son[i]);
		if (curLevel>1) {
			ReadIndexBlock(bt,node->son[i],cNode);
			printSubTree(bt,cNode,indent+1);
		} else if (curLevel==1) {
			// fetch data page
		}
	}
	printf("\n");
}

void printTree(BTree* bt) {
	NodeStruct* root_node=createNodeStruct();
	ReadIndexBlock(bt,bt->root,root_node);
	printSubTree(bt,root_node,0);
}

// functions for visualizing clusters

FastArray<float> xcrd,ycrd;
FILE** views;	//=fopen(visualf,"w");
int Cnum;

void openVisualFiles(const char* nodename,const char* outprefix,int _Cnum) {
	char nodef[255],visualf[255];
	sprintf(nodef,"data/%s",nodename);
	FILE* cnode=fopen(nodef,"r");

	int id,Ni,Nj;
	float dist,x,y;

	// read node info with format: NodeId xcrd ycrd
	while (!feof(cnode)) {	// assume NodeID in ascending order from 0
		fscanf(cnode,"%d %f %f\n", &id, &x, &y);
		xcrd.push_back(x);	ycrd.push_back(y);
	}
	fclose(cnode);

	Cnum=_Cnum;
	views=new FILE*[Cnum+1];		// last file for outlier
	for (int i=0;i<Cnum+1;i++) {	// including file for outliers
		if (i==Cnum)
			sprintf(visualf,"visual/%s.odat",outprefix);
		else
			sprintf(visualf,"visual/%s.pdat%d",outprefix,i);
		remove(visualf);
		views[i]=fopen(visualf,"w");
	}
}

void writeVisualRecord(int Ni,int Nj,float vP,float eDist,int Cid) {
	float x1,y1,x2,y2,x,y;

	x1=xcrd[Ni];	y1=ycrd[Ni];
	x2=xcrd[Nj];	y2=ycrd[Nj];

	x=x1+(vP/eDist)*(x2-x1);
	y=y1+(vP/eDist)*(y2-y1);

	if (Cid<0)
		fprintf(views[Cnum],"%f %f\n",x,y);
	else
		fprintf(views[Cid],"%f %f\n",x,y);
}

void closeVisualFiles() {	// close files
	xcrd.clear();	ycrd.clear();
	for (int i=0;i<Cnum+1;i++) fclose(views[i]);	// including file for outliers
	delete[] views;
}



#define DISKBASED_H
#endif