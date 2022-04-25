#ifndef MAXBRGSTKNN_UTILITY_CC

#include "utility.h"
#include "drand48.h"
#include <iostream>
#include "IOcontroller.h"
using namespace std;

//int START_TIME;
clock_t START_TIME, END_TIME;
double duration;
double ZipfMaxVal;

void InitClock() {START_TIME=clock();srand(time(NULL));srand48(time(NULL));}
void PrintElapsed() { END_TIME=clock();duration=((double)(END_TIME-START_TIME)/CLOCKS_PER_SEC*1000); printf("Elapsed: %d msec\n",(int)duration);}

void CheckFile(FILE* fp,const char* filename) {
	if (fp==NULL) {
		printf("Invalid file '%s'\n",filename);
		exit(0);
	}
}



void TrimSpace(char* str) {
	if (str==NULL) return;

	char space[]={'\t','\n','\f','\r',' '};
	int pos=0;
	for (int i=0;i<strlen(str);i++) {
		bool found=false;
		for (int j=0;j<5;j++)
			if (str[i]==space[j]) found=true;

		if (!found) {
			str[pos]=str[i];
			pos++;
		}
	}
	str[pos]='\0';
}

void AddConfigFromFile(ConfigType &cr,const char* filename) {
	const int LINE_LEN=1024;
	char line[LINE_LEN],key[LINE_LEN],value[LINE_LEN];
	//printf("before\n");
	//getchar();
	ifstream br(filename);

  	if (! br.is_open())
  	{ printf("jin, Error opening file \"%s\"",filename); exit (1); }

	while (br.getline(line,LINE_LEN)){
		if (strstr(line,"//")!=NULL) continue; // remove comments
		char* chPos=strchr(line,'=');
		if (chPos!=NULL) {
			int pos=((int)(chPos-line))/sizeof(char);
			int keyLen=pos;
			int valueLen=strlen(line)-1-keyLen;
			memcpy(key,&line[0],keyLen);	key[keyLen]='\0';
			memcpy(value,&line[pos+1],valueLen);	value[valueLen]='\0';
			TrimSpace(key);	TrimSpace(value);
			cr[key]=value;
		}
	}
	br.close();
}

void AddConfigFromCmdLine(ConfigType &cr,int argc,char** argv) {
	int i=0;
	while (i<argc) {
		while ((i<argc)&&(argv[i][0]!='-')) i++;	// shortcut condition
		if (i+1<argc) {
			char* key=&(argv[i][1]);
			char* value=argv[i+1];
			TrimSpace(key);	TrimSpace(value);
			cr[key]=value;
			i+=2;
		} else
			return;
	}
}

void ListConfig(ConfigType &cr) {
	ConfigType::iterator p=cr.begin();
	while (p!=cr.end()) {
		printf("%s=%s\n",p->first.c_str(),p->second.c_str());
		p++;
	}
}

float getConfigFloat(const char* key,ConfigType &cr,bool required,float _default) {
	float value=_default;
	if (cr.count(key))
		value=atof(cr[key].c_str());
	else {
		if (required) {
			printf("Config key \"%s\" not found\n",key);
			exit(1);
		}
	}
	return value;
}

int getConfigInt(const char* key,ConfigType &cr,bool required,int _default) {
	int value=_default;
	if (cr.count(key))
		value=atoi(cr[key].c_str());
	else {
		if (required) {
			printf("Config key \"%s\" not found\n",key);
			exit(1);
		}
	}
	return value;
}

const char* getConfigStr(const char* key,ConfigType &cr) {
	const char* value;
	if (cr.count(key))
		value=cr[key].c_str();
	else {
		printf("Config key \"%s\" not found\n",key);
		exit(1);
	}
	return value;
}





#define MAXBRGSTKNN_UTILITY_CC

#endif //MAXBRGSTKNN_UTILITY_CC