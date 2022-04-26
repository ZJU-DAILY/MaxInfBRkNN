#ifndef HASHTREEMAP_CC
#define HASHTREEMAP_CC

#include <map>
#include <set>
#include <vector>

using namespace std;

// TODO : change function find to return the address instead of copying the
// element
//typedef unsigned int ID;
//typedef set<ID> T;

template <class ID, class T> 
class HashTreeMap {

	int _length; // number of rows in the hashtable
	int _size; // total number of elements
	vector<map<ID, T> >* ht;

public:
	HashTreeMap(int size);
	~HashTreeMap();
	void insert(ID i, T t);
	T find(ID i);
	int size() const;
	int length() const;
	map<ID, T>& getRow(ID i) const;

};

template <class ID, class T>
HashTreeMap<ID, T>::HashTreeMap(int l) {
	_length = l;
	ht = new vector<map<ID, T> >(l);
	_size = 0;
}

template <class ID, class T>
HashTreeMap<ID, T>::~HashTreeMap() {

	/*
	for (int i=0; i<length; i++) {
		map<ID, T>& curRow = ht->at(i);

		for (map<ID,T>::iterator j=curRow.begin(); j!=curRow.end(); j++) {
			delete j->second;
		}
//		curRow.clear();
//		delete curRow->second;
	}
	*/
	delete ht;
	_length = 0;
	_size = 0;
}


template <class ID, class T>
void HashTreeMap<ID, T>::insert(ID i, T t) {
	int row = i % _length;
	map<ID, T>& element = ht->at(row);

	element.insert( pair<ID, T>(i, t) );
	_size++;
}

template <class ID, class T>
T HashTreeMap<ID, T>::find(ID i) {
	int row = i % _length;
	map<ID, T>& element = ht->at(row);

//	map<ID, T>::iterator j = element.find(i);

	if (element.find(i) == element.end()) 
		return NULL;
	else 
		return element[i];
	return NULL;
}

template <class ID, class T>
int HashTreeMap<ID, T>::size() const {
	return _size;
}

template <class ID, class T>
int HashTreeMap<ID, T>::length() const {
	return _length;
}

template <class ID, class T>
map<ID, T>& HashTreeMap<ID, T>::getRow(ID i) const {
	map<ID, T>& element = ht->at(i);

	return element;
}


#endif
