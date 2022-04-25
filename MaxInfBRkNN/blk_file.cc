#include "blk_file.h"

//===========================BlockFile============================
void BlockFile::fwrite_number(int value) {
   put_bytes((char *) &value,sizeof(int));
}

int BlockFile::fread_number() {
   char ca[sizeof(int)];
   get_bytes(ca,sizeof(int));
   return *((int *)ca);
}

BlockFile::BlockFile(char* name,int b_length) {
	char *buffer;
	int l;
	
	filename = new char[strlen(name) + 1];
	strcpy(filename,name);
	blocklength = b_length;	
	number = 0;
	
	if ((fp=fopen(name,"rb+"))!=0) {
		new_flag = FALSE;
		blocklength = fread_number();
		number = fread_number();
	} else{
     	if (blocklength < BFHEAD_LENGTH) {
     		printf("BlockFile::error in initialization\n");
			exit(1);
		}

		fp=fopen(filename,"wb+");
		if (fp==0) {
			printf("BlockFile::cannot create new file\n");
			exit(1);
		}
		new_flag = TRUE;
		fwrite_number(blocklength);
		
		// Gesamtzahl Bloecke ist 0
		fwrite_number(0);
		buffer = new char[(l=blocklength-(int)ftell(fp))];
		memset(buffer, 0, sizeof(buffer));
		put_bytes(buffer,l);
		
		delete [] buffer;
	}
	fseek(fp,0,SEEK_SET);
	act_block=0;
}

BlockFile::~BlockFile() {
   delete[] filename;
   fclose(fp);
}

void BlockFile::read_header(char* buffer) {
   fseek(fp,BFHEAD_LENGTH,SEEK_SET);
   get_bytes(buffer,blocklength-BFHEAD_LENGTH);

   if(number<1) {
       fseek(fp,0,SEEK_SET);
       act_block=0;
   } else
       act_block=1;
}

void BlockFile::set_header(char* header) {
	fseek(fp,BFHEAD_LENGTH,SEEK_SET);
	put_bytes(header,blocklength-BFHEAD_LENGTH);
	if(number<1) {
		fseek(fp,0,SEEK_SET);
		act_block=0;
	} else
		act_block=1;
}

bool BlockFile::read_block(Block b,int pos) {
   pos++;
   if (pos<=number && pos>0)
       seek_block(pos);
   else
       return FALSE;

   get_bytes(b,blocklength);
   if (pos+1>number) {
       fseek(fp,0,SEEK_SET);
       act_block=0;
   } else
       act_block=pos+1;

   return TRUE;
}

bool BlockFile::write_block(const Block block,int pos) {
   pos++;

   if (pos<=number && pos>0)
       seek_block(pos);
   else
       return FALSE;

   put_bytes(block,blocklength);
   if (pos+1>number) {
       fseek(fp,0,SEEK_SET);
       act_block=0;
   } else
       act_block=pos+1;

   return TRUE;
}

int BlockFile::append_block(const Block block) {
   fseek(fp,0,SEEK_END);
   put_bytes(block,blocklength);
   number++;
   fseek(fp,sizeof(int),SEEK_SET);
   fwrite_number(number);
   fseek(fp,-blocklength,SEEK_END);

   return (act_block=number)-1;
}

//========================CachedBlockFile=========================
int CachedBlockFile::next() {
	int ret_val, tmp;

	if (cachesize==0)
		return -1;
	else
		if (fuf_cont[ptr]==free) {
			ret_val = ptr++;
			if (ptr==cachesize) ptr=0;
			return ret_val;
		} else {
			tmp = ptr+1;
			if (tmp==cachesize) tmp=0;
			while (tmp!=ptr && fuf_cont[tmp]!=free)
				if(++tmp==cachesize) tmp=0;

			if (ptr==tmp) {
				BlockFile::write_block(cache[ptr],cache_cont[ptr]-1);
				fuf_cont[ptr]=free;
				ret_val=ptr++;			
				if (ptr==cachesize) ptr=0;
				
				return ret_val;
			} else
				return tmp;
		}
}

int CachedBlockFile::in_cache(int index) {
   for (int i = 0; i<cachesize; i++)
	   if (cache_cont[i] == index && fuf_cont[i] != free)
		   return i;
   return -1;
}

CachedBlockFile::CachedBlockFile(char* name,int blength, int csize)
   : BlockFile(name,blength) {
	ptr=0;
	if (csize>=0)
		cachesize=csize;
	else
		printf("CachedBlockFile::CachedBlockFile: neg. Cachesize\n");
		
	cache_cont = new int[cachesize];
	fuf_cont = new uses[cachesize];
	cache = new char*[cachesize];
	for (int i=0; i<cachesize; i++) {
		cache_cont[i]=0;  
		fuf_cont[i]=free;
		cache[i] = new char[get_blocklength()];
	}
}

CachedBlockFile::~CachedBlockFile() {
	flush();
	delete[] cache_cont;
	delete[] fuf_cont;

	for (int i=0;i<cachesize;i++)
		delete[] cache[i];
	delete[] cache;
}

bool CachedBlockFile::read_block(Block block,int index) {
	int c_ind;

	index++;
	if(index<=get_num_of_blocks() && index>0) {
		if((c_ind=in_cache(index))>=0)
			memcpy(block,cache[c_ind],get_blocklength());
		else {
			c_ind = next();
			if (c_ind >= 0) {
				BlockFile::read_block(cache[c_ind],index-1); 
				cache_cont[c_ind]=index;
				fuf_cont[c_ind]=used;
				memcpy(block,cache[c_ind],get_blocklength());
			}
			else
				BlockFile::read_block(block,index-1); 
		}
		return TRUE;
	} else
		return FALSE;
}

bool CachedBlockFile::write_block(const Block block,int index) {
	int c_ind;

	index++;	
	if(index <= get_num_of_blocks() && index > 0) {
		c_ind = in_cache(index);
		if(c_ind >= 0)	
			memcpy(cache[c_ind], block, get_blocklength());
		else {
			c_ind = next();
			if (c_ind >= 0) {
				memcpy(cache[c_ind],block,get_blocklength());
				cache_cont[c_ind]=index;
				fuf_cont[c_ind]=used;
			} else
				BlockFile::write_block(block,index-1);
		}
		return TRUE;
	} else
		return FALSE;
}

void CachedBlockFile::flush() {
	for (int i=0; i<cachesize; i++)
		if (fuf_cont[i]!=free)	
			BlockFile::write_block(cache[i], cache_cont[i]-1); 
}
