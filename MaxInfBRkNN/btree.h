#ifndef __BTREE


#include "blk_file.h"

class BTDirNode;

class BTree {
public:
    friend class BTDirNode;
	int UserField;
    int root;						// block # of root node
    BTDirNode  *root_ptr;           // root-node
    int num_of_inodes;	        	// # of stored directory pages
    CachedBlockFile *file;	  		// storage manager for harddisc blocks
    void load_root();            	// loads root_node into memory
    char *header;
protected:
    char *user_header;
    void read_header(char *buffer);      // reads Rtree header
    void write_header(char *buffer);     // writes Rtree header
public:
    BTree(char *fname,int _b_length, int cache_size);
    BTree(char *fname, int cache_size);
    virtual ~BTree();
};

class DirEntry {
    friend class BTDirNode ;
    friend class BTree ;

public:
    BTree  *my_tree;        // pointer to my B-tree
    BTDirNode  *son_ptr;       // pointer to son if in main mem.
    int key;              // key value
    int son;                // block # of son

    BTDirNode * get_son();  	// returns the son, loads son if necessary
    void read_from_buffer(char *buffer);// reads data from buffer
    void write_to_buffer(char *buffer); // writes data to buffer

    static const int EntrySize=sizeof(int)+sizeof(int);
    						// returns amount of needed buffer space
    
    virtual DirEntry  & operator = (DirEntry  &_d);
    DirEntry(BTree *bt = NULL);
    virtual ~DirEntry();
};

class BTDirNode {
    friend class BTree ;
public:
	BTree  *my_tree;               		// pointer to B-tree
	int capacity;                       // max. # of entries
	int num_entries;                    // # of used entries
	bool dirty;                         // TRUE, if node has to be written
	int block;                          // disc block
	char level;                         // level of the node in the tree
		
    DirEntry  *entries;            			// array of entries
    void read_from_buffer(char *buffer);	// reads data from buffer
    void write_to_buffer(char *buffer); 	// writes data to buffer

    BTDirNode(BTree  *rt);
    BTDirNode(BTree  *rt, int _block);
    virtual ~BTDirNode();
};

#define __BTREE

#endif // __BTREE
