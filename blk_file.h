#ifndef __BLKFILE
#define __BLKFILE

#define BFHEAD_LENGTH (sizeof(int)*2) 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef char Block[];
#define TRUE 1
#define FALSE 0

class BlockFile {
   FILE* fp;			
   char* filename;		
   int blocklength;	        
   int act_block; 	        // SEEK_CUR 
   int number;		        
   bool new_flag;		

   void put_bytes(const char* bytes,int num)  { fwrite(bytes,num,1,fp); }
   void get_bytes(char* bytes,int num)	      { fread(bytes,num,1,fp); }
   void fwrite_number(int num);	
   int fread_number();		
   void seek_block(int bnum)    
   { fseek(fp,(bnum-act_block)*blocklength,SEEK_CUR); }

public:
   BlockFile(char* name, int b_length);
   ~BlockFile();

   void read_header(char * header);	
   void set_header(char* header);	

   bool read_block(Block b,int i);
   bool write_block(const Block b,int i);
   int append_block(const Block b);	
   bool file_new()			{ return new_flag; }		
   int get_blocklength()	{ return blocklength; }
   int get_num_of_blocks()	{ return number; }
};

class CachedBlockFile : public BlockFile {
   enum uses {free,used};
   int ptr;      
   int cachesize;	
   int *cache_cont;
   uses *fuf_cont;		
   char **cache;
   int next();	
   int in_cache(int index);	

public:
   CachedBlockFile(char* name,int blength, int csize);
   ~CachedBlockFile();

   bool read_block(Block b,int i);
   bool write_block(const Block b,int i);
   void flush();	
};

#endif //__BLKFILE
