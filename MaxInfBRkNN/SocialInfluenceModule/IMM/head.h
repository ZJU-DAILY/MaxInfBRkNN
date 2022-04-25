#ifndef __HEAD_H__
#define __HEAD_H__


/*
#ifdef DEBUG_INFO
this_is_deprecated
#endif
#ifdef DEBUG_TRACE
this_is_deprecated
#endif


#if defined(WIN32)
#elif defined(__CYGWIN__) // cygwin*/
#include <sys/time.h>
//#else //linux
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
//#endif


#include <iostream>
#include <set>
#include <list>
#include <sstream>
#include <cmath>
#include <queue>
#include <fstream>
#include <string>
#include <cstdio>
#include <functional>
#include <algorithm>
#include <climits>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <map>
#include <deque>
using namespace std;
typedef unsigned int uint;
typedef unsigned char uint8;
typedef long long int64;
typedef unsigned long long uint64;
typedef pair<int,int> ipair;
typedef pair<double,double> dpair;
#define MP make_pair
#define F first
#define S second

#ifndef TIMES_PER_SEC
#define TIMES_PER_SEC (2393.910e6)
#endif




typedef char int8;
typedef unsigned char uint8;
typedef long long int64;
typedef unsigned long long uint64;

#ifdef WIN32
#include < time.h >
#include <windows.h> //I've ommited this line.
#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif
struct timezone
{
    int  tz_minuteswest; /* minutes W of Greenwich */
    int  tz_dsttime;     /* type of dst correction */
};

int gettimeofday(struct timeval *tv, struct timezone *tz)
{
    FILETIME ft;
    unsigned __int64 tmpres = 0;
    static int tzflag;

    if (NULL != tv)
    {
        GetSystemTimeAsFileTime(&ft);

        tmpres |= ft.dwHighDateTime;
        tmpres <<= 32;
        tmpres |= ft.dwLowDateTime;

        /*converting file time to unix epoch*/
        tmpres -= DELTA_EPOCH_IN_MICROSECS;
        tmpres /= 10;  /*convert into microseconds*/
        tv->tv_sec = (long)(tmpres / 1000000UL);
        tv->tv_usec = (long)(tmpres % 1000000UL);
    }

    if (NULL != tz)
    {
        if (!tzflag)
        {
            _tzset();
            tzflag++;
        }
        tz->tz_minuteswest = _timezone / 60;
        tz->tz_dsttime = _daylight;
    }

    return 0;
}
#endif


#ifndef WIN32
#ifdef __CYGWIN__

//CYGWIN
uint64 rdtsc()
{
    uint64 t0;
    asm volatile("rdtsc" : "=A"(t0));
    return t0;

}
#else
//LINUX
uint64 rdtsc(void)
{
    unsigned a, d;
    //asm("cpuid");
    asm volatile("rdtsc" : "=a" (a), "=d" (d));
    return (((uint64)a) | (((uint64)d) << 32));
}
#endif
#endif


string nowStr(){
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    char str[100];
    sprintf(str,"%d_%d_%d_%d_%d_%d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    return string(str);
}


map<string, timeval> __head_h_start;
void timer_init(string arg="default"){
    timeval ts;
    gettimeofday(&ts,NULL);
    __head_h_start[arg]=ts;
}
int64 timer_elapse(string arg="default"){ // unit ms
    cout << "fuck " <<endl;
    struct timeval now;
    gettimeofday(&now,NULL);
    int64 sec=now.tv_sec-__head_h_start[arg].tv_sec;
    int64 usec=now.tv_usec-__head_h_start[arg].tv_usec;
    return sec * 1000 + usec/1000;
}
#define SIZE(t) (int)(t.size())
#define ALL(t) (t).begin(), (t).end()
#define FOR(i,n) for(int (i)=0; (i)<((int)(n)); (i)++)
#ifdef WIN32
#define FORE(i, x) for (typeid((x).begin()) i = (x).begin(); (i) != (x).end(); (i)++)
#else
#define FORE(i, x) for (__typeof((x).begin()) i = (x).begin(); (i) != (x).end(); (i)++)
#endif


static inline string &ltrim(string &s) { s.erase(s.begin(), find_if(s.begin(), s.end(), not1(ptr_fun<int, int>(isspace)))); return s; }
static inline string &rtrim(string &s) { s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end()); return s; }
static inline string &trim(string &s) { return ltrim(rtrim(s)); }
string __n_variable(string t, int n){ t=t+','; int i=0; if(n) for(; i<SIZE(t)&&n; i++) if(t[i]==',') n--; n=i; for(;t[i]!=',';i++); t=t.substr(n, i-n); trim(t);if(t[0]=='"') return ""; return t+"="; }
#define __expand_nv(x) __n_variable(t, x)<< t##x << " "
template<class T0>
void ___debug(string t,deque<T0> t0,ostream&os){os<<__n_variable(t,0);FOR(i, SIZE(t0))os<<t0[i]<<" ";}
template<class T0>
void ___debug(string t,set<T0> t0,ostream&os){os<<__n_variable(t,0);FORE(i,t0)os<<*i<<" ";}
template<class T0>
void ___debug(string t,vector<T0> t0,ostream&os){os<<__n_variable(t,0);FOR(i, SIZE(t0))os<<t0[i]<<" ";}
template<class T0,class T1>
void ___debug(string t,vector<pair<T0,T1> > t0,ostream&os){os<<__n_variable(t,0);FOR(i, SIZE(t0))os<<t0[i].F<<","<<t0[i].S<<" ";}
template<class T0>
void ___debug(string t,T0 t0,ostream&os){os<<__expand_nv(0);}
template<class T0,class T1>
void ___debug(string t,T0 t0,T1 t1,ostream&os){os<<__expand_nv(0)<<__expand_nv(1);}
template<class T0,class T1,class T2>
void ___debug(string t,T0 t0,T1 t1,T2 t2,ostream&os){os<<__expand_nv(0)<<__expand_nv(1)<<__expand_nv(2);}
template<class T0,class T1,class T2,class T3>
void ___debug(string t,T0 t0,T1 t1,T2 t2,T3 t3,ostream&os){os<<__expand_nv(0)<<__expand_nv(1)<<__expand_nv(2)<<__expand_nv(3);}
template<class T0,class T1,class T2,class T3,class T4>
void ___debug(string t,T0 t0,T1 t1,T2 t2,T3 t3, T4 t4,ostream&os){os<<__expand_nv(0)<<__expand_nv(1)<<__expand_nv(2)<<__expand_nv(3)<<__expand_nv(4);}


//#define DO_ONCE


#define RUN_TIME(...) { int64 t=rdtsc();  __VA_ARGS__; t=rdtsc()-t; cout<<  #__VA_ARGS__ << " : " << t/TIMES_PER_SEC <<"s"<<endl;  }

#ifdef HEAD_TRACE
#define TRACE(...) {{ ___debug( #__VA_ARGS__,  __VA_ARGS__,cerr); cerr<<endl;  } }
#define IF_TRACE(args) args
#define TRACE_LINE(...) { ___debug( #__VA_ARGS__,  __VA_ARGS__,cerr); cerr<<"                    \033[100D";  }
#define TRACE_SKIP(a, ...) { static int c=-1; c++; if(c%a==0)TRACE( __VA_ARGS__); }
#define TRACE_LINE_SKIP(a, ...) { static int c=-1; c++; if(c%a==0) TRACE_LINE(__VA_ARGS__);  }
#define TRACE_LINE_END(...) {cerr<<endl; }
ofstream __HEAD_H__LOG("log.txt");
#define TRACE_LOG(...) { __HEAD_H__LOG.close(); ofstream cerr("log.txt", ofstream::out|ofstream::app); ___debug( #__VA_ARGS__,  __VA_ARGS__, cerr); cerr<<endl;  }
#else
#define TRACE(...) ;
#define IF_TRACE(args) ;
#define TRACE_LINE(...) ;
#define TRACE_SKIP(a, ...) ;
#define TRACE_LINE_SKIP(a, ...) ;
#define TRACE_LINE_END(...) ;
#define TRACE_LOG(...) ;
#endif //HEAD_TRACE


#ifdef HEAD_INFO
ofstream __HEAD_H_FOUT;
void setInfoFile(string s){__HEAD_H_FOUT.open(s.c_str()); }
#define ASSERT(v) {if (!(v)) {cerr<<"ASSERT FAIL @ "<<__FILE__<<":"<<__LINE__<<endl; exit(1);}}
#define INFO(...) do { ___debug( #__VA_ARGS__,  __VA_ARGS__,cout); cout<<endl; if(__HEAD_H_FOUT.is_open()){___debug( #__VA_ARGS__,  __VA_ARGS__,__HEAD_H_FOUT); __HEAD_H_FOUT<<endl;}  } while(0)
#define ASSERTT(v, ...) {if (!(v)) {cerr<<"ASSERT FAIL @ "<<__FILE__<<":"<<__LINE__<<endl; INFO(__VA_ARGS__); exit(1);}}
#else
#define ASSERTT(v, ...) ;
#define ASSERT(v ) ;
#define INFO(...) ;
#endif

template<class T>
string toStr(T t){ stringstream ss; ss<<t; return ss.str(); }
string strcat2(const string a, const string b){return (a+b); }

string currentTimestampStr() {
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    char buf[1000];
    sprintf(buf, "%04d-%02d-%02d %02d:%02d:%02d", now->tm_year+1900, now->tm_mon+1, now->tm_mday, now->tm_hour, now->tm_min, now->tm_sec);
    return string(buf);
}



class Timer
{
    public:
        static vector<int64> timeUsed;
        static vector<string> timeUsedDesc;
        int id;
        uint64 startTime;
        bool showOnDestroy;
        Timer(int id, string desc="", bool showOnDestroy=false)
        {
            this->id=id;
            while((int)timeUsed.size() <= id)
            {
                timeUsed.push_back(0);
                timeUsedDesc.push_back("");
            }
            timeUsedDesc[id]=desc;
            startTime=rdtsc();
            this->showOnDestroy=showOnDestroy;
        }
        ~Timer()
        {
            timeUsed[id]+=(rdtsc()-startTime);
        }
        static void show(bool debug=false)
        {
            INFO("########## Timer ##########");

            for (int i=0;i<(int)timeUsed.size();i++) if (timeUsed[i]>0)
            {
                char str[100];
                sprintf(str,"%.6lf",timeUsed[i]/TIMES_PER_SEC );
                string s=str;
                if ((int)s.size()<15) s=" "+s;
                char t[100];
                memset(t, 0, sizeof t);
                sprintf(t,"Spend %s seconds on %s",s.c_str(), timeUsedDesc[i].c_str());
                cout<< t << endl;
            }
        }
        static void clearAll()
        {
            timeUsed.clear();
            timeUsedDesc.clear();
        }
};
vector<int64> Timer::timeUsed;
vector<string> Timer::timeUsedDesc;


class Counter
{
    public:
        static int cnt[1000];
        Counter(int id=0)
        {
            cnt[id]++;
        }
        ~Counter()
        {
        }
        static void show()
        {
            for (int i=0;i<1000;i++) if (cnt[i]>0)
                INFO("Counter", i,cnt[i]);
        }
};
int Counter::cnt[1000]={0};

// return the output of the command by string
string exec(const char* cmd) {
    FILE* pipe = popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
        if(fgets(buffer, 128, pipe) != NULL)
            result += buffer;
    }
    pclose(pipe);
    return result;
}
string getIpAddress(){
    string cmd=" ifconfig |grep 'inet addr' |grep -v 127.0.0.1 |awk '{print $2}' |python -c 'print raw_input().split(\":\")[-1]' ";
    return exec(cmd.c_str());
}

string __head_version = "";
class OutputInfo
{
    public:
    OutputInfo(int argn, char ** argv){
        cout<<"\e\[0;32mProgram Start at: " << currentTimestampStr()<<"\e[0m"<<endl;
        cout<<"\e\[0;32mProgram version: " << __head_version << "\e[0m"<<endl;
        cout<<"Arguments: ";
        for(int i=0; i<argn; i++){
            cout<<argv[i]<<" ";
        }
        cout<<endl;
        cout<<"--------------------------------------------------------------------------------" <<endl;
    }
    ~OutputInfo(){
        cout<<"\e\[0;31mProgram Terminate at: " << currentTimestampStr()<< "\e[0m"<<endl;
        cout<<"\e\[0;31mProgram version: " << __head_version << "\e[0m"<<endl;
    }
};
#endif //__HEAD_H__
