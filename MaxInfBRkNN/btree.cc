#include "btree.h"

// DirEntry

DirEntry ::DirEntry(BTree  *bt) {
	my_tree = bt;
	son_ptr = NULL;
}

DirEntry ::~DirEntry() {
	if (son_ptr != NULL) delete son_ptr;
}

DirEntry & DirEntry ::operator = (DirEntry &_d) {
	son = _d.son;
	son_ptr = _d.son_ptr;
	key = _d.key;
	return *this;
}

void DirEntry ::read_from_buffer(char *buffer) {
	memcpy(&key, buffer, sizeof(int));
	memcpy(&son, buffer+sizeof(int), sizeof(int));
}

void DirEntry ::write_to_buffer(char *buffer) {
    memcpy(buffer, &key, sizeof(int));
    memcpy(buffer+sizeof(int), &son, sizeof(int));
}

BTDirNode *DirEntry ::get_son() {
    if (son_ptr == NULL) son_ptr = new BTDirNode (my_tree, son);

    return son_ptr;
}

// BTDirNode

// creates a brand new BT directory node
BTDirNode ::BTDirNode(BTree  *bt) {
	// from parent class
	my_tree = bt;	num_entries = 0;	block = -1;

	// header of page keeps node info.
	// level + num_entries
	int header_size = sizeof(char) + sizeof(int);
	capacity = (bt->file->get_blocklength() - header_size) / (DirEntry::EntrySize);

	// Initialize entries
	entries = new DirEntry[capacity];
	for (int iter=0;iter<capacity;iter++) {
		entries[iter].my_tree=my_tree;
		entries[iter].son_ptr=NULL;
	}

	// create new block for the node
	char* b = new char[bt->file->get_blocklength()];
	block = bt->file->append_block(b);
	delete [] b;

	bt->num_of_inodes++;
	dirty = true;	// must be written to disk before destruction
}

// reads an existing BT directory node
BTDirNode ::BTDirNode(BTree  *bt, int _block) {
	// from parent class
	my_tree = bt;	num_entries = 0;	block = -1;

	// header of page keeps node info.
	// level + num_entries
	int header_size = sizeof(char) + sizeof(int);
	capacity = (bt->file->get_blocklength() - header_size) / (DirEntry::EntrySize);


	entries = new DirEntry[capacity];		// Initialize entries
	for (int iter=0;iter<capacity;iter++) {
		entries[iter].my_tree=my_tree;
		entries[iter].son_ptr=NULL;
	}

	// now load block and read BTDirNode data from it.
	block = _block;
	char* b = new char[bt->file->get_blocklength()];
	bt->file->read_block(b, block);
	read_from_buffer(b);
	delete [] b;

	dirty = FALSE;		// not dirty yet
}

BTDirNode ::~BTDirNode() {
    if (dirty){
		// Update changes on disk
		char* b = new char[my_tree->file->get_blocklength()];
		write_to_buffer(b);
		my_tree->file->write_block(b, block);
        delete [] b;
    }
    delete [] entries;
}

void BTDirNode ::read_from_buffer(char *buffer) {
    int j=0;

    // first read header info
    memcpy(&level, &buffer[j], sizeof(char));
    j += sizeof(char);

    memcpy(&num_entries, &buffer[j], sizeof(int));
    j += sizeof(int);

    // then read entries
    int s = DirEntry::EntrySize;
    for (int i = 0; i < num_entries; i++) {
		entries[i].read_from_buffer(&buffer[j]);
		j += s;
    }
}

void BTDirNode ::write_to_buffer(char *buffer) {
    int j=0;

    // first, write header info
    memcpy(&buffer[j], &level, sizeof(char));
    j += sizeof(char);
    memcpy(&buffer[j], &num_entries, sizeof(int));
    j += sizeof(int);

    // then, write entries
    int s = DirEntry::EntrySize;
    for (int i = 0; i < num_entries; i++) {
		entries[i].write_to_buffer(&buffer[j]);
		j += s;
    }
}

// BTree

// Construction of a new BTree  //jordan
BTree ::BTree(char *fname, int _b_length, int cache_size) {
    file = new CachedBlockFile(fname, _b_length, cache_size);

    // first block is header
    header = new char [file->get_blocklength()];

    root = 0;
    root_ptr = NULL;
    num_of_inodes = 0;

    //create an (empty) root
    root_ptr = new BTDirNode (this);
    root = root_ptr->block;
}

// load an existing BTree
BTree ::BTree(char *fname, int cache_size) {
    file = new CachedBlockFile(fname, 0, cache_size);	// load file

    // read header
    header = new char [file->get_blocklength()];
    file->read_header(header);
    read_header(header);
    root_ptr = NULL;
    load_root();
}

BTree ::~BTree() {
    write_header(header);
    file->set_header(header);
    delete [] header;

    if (root_ptr != NULL) {
        delete root_ptr;
        root_ptr = NULL;
    }
    delete file;
    //printf("saved B-Tree containing %d internal nodes\n",num_of_inodes);
}

void BTree ::read_header(char *buffer) {
    int i = 0;

    memcpy(&num_of_inodes, &buffer[i], sizeof(num_of_inodes));
    i += sizeof(num_of_inodes);

    memcpy(&root, &buffer[i], sizeof(root));
    i += sizeof(root);

	memcpy(&UserField, &buffer[i], sizeof(UserField));
    i += sizeof(UserField);

    user_header = &buffer[i];
}

void BTree ::write_header(char *buffer) {
    int i = 0;

    memcpy(&buffer[i], &num_of_inodes, sizeof(num_of_inodes));
    i += sizeof(num_of_inodes);

    memcpy(&buffer[i], &root, sizeof(root));
    i += sizeof(root);

	memcpy(&buffer[i], &UserField, sizeof(UserField));
    i += sizeof(UserField);

    user_header = &buffer[i];
}

void BTree ::load_root() {
    if (root_ptr == NULL)
		root_ptr = new BTDirNode (this, root);
}
