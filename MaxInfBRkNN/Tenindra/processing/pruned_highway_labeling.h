/*
 The MIT License (MIT)

 Copyright (c) 2015 Yuki Kawata

 Permission is hereby granted, free of charge, to any person obtaining a copy of
 this software and associated documentation files (the "Software"), to deal in
 the Software without restriction, including without limitation the rights to
 use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 the Software, and to permit persons to whom the Software is furnished to do so,
 subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
 Obtained from: https://github.com/kawatea/pruned-highway-labeling
*/

#ifndef __PRUNED_HIGHWAY_LABELING_H__
#define __PRUNED_HIGHWAY_LABELING_H__

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <queue>
#include <algorithm>
#include <malloc.h>
#include <xmmintrin.h>
#include <sys/time.h>



class PrunedHighwayLabeling {
    public:
    
    PrunedHighwayLabeling() : V(0), label(NULL), load_time(0), construct_time(0) {}
    ~PrunedHighwayLabeling() { Free(); };
    
    void ConstructLabel(const char *file);
    void LoadLabel(const char *file);
    void StoreLabel(const char *file);
    double computeIndexSize();
    double getConstructionTime();
    
    inline int Query(int v, int w) {
        //printf("phl for dist(v%d,v%d)\n",v,w);
        if (v == w) return 0;
        if (contract[v] == w || contract[w] == v || (contract[v] != -1 && contract[v] == contract[w])) return label[v].time + label[w].time;
        int time = INF;
        unsigned *vf = label[v].path;
        int *vs = label[v].cost;
        unsigned *wf = label[w].path;
        int *ws = label[w].cost;
        
        _mm_prefetch(vf, _MM_HINT_T0);
        _mm_prefetch(vs, _MM_HINT_T0);
        _mm_prefetch(wf, _MM_HINT_T0);
        _mm_prefetch(ws, _MM_HINT_T0);
        
        while (true) {
            if ((*vf & PATH_MASK) == (*wf & PATH_MASK)) {
                if (*vf == GUARD) return time + label[v].time + label[w].time;
                for (unsigned i = 0, j = 0; ; ) {
                    if (*(vs + i) == *(ws + j)) {
                        ChangeMin(time, *(vs + i + 1) + *(ws + j + 1));
                        i += 2;
                        j += 2;
                        if (i == (*vf & NUM_MASK) || j == (*wf & NUM_MASK)) break;
                    } else if (*(vs + i) < *(ws + j)) {
                        ChangeMin(time, *(vs + i + 1) + *(ws + j + 1) + *(ws + j) - *(vs + i));
                        i += 2;
                        if (i == (*vf & NUM_MASK)) break;
                    } else {
                        ChangeMin(time, *(vs + i + 1) + *(ws + j + 1) + *(vs + i) - *(ws + j));
                        j += 2;
                        if (j == (*wf & NUM_MASK)) break;
                    }
                }
                vs += *vf & NUM_MASK;
                vf++;
                ws += *wf & NUM_MASK;
                wf++;
            } else if ((*vf & PATH_MASK) < (*wf & PATH_MASK)) {
                if (*wf == GUARD) return time + label[v].time + label[w].time;
                vs += *vf & NUM_MASK;
                vf++;
            } else {
                if (*vf == GUARD) return time + label[v].time + label[w].time;
                ws += *wf & NUM_MASK;
                wf++;
            }
        }
    }
    
    void Free(void);
    
    int NumVertices(void) { return V; }
    void Statistics(void);
    
    private:
    
    static const int LEVEL = 4;
    static const int INF = 1000000000;
    static const unsigned GUARD = 0xFFFFFFFFU;
    static const unsigned PATH_MASK = 0xFFFFFC00U;
    static const unsigned NUM_MASK = 0x000003FFU;
    
    struct road {
        int from;
        int to;
        int time;
        int dist;
    };
    
    struct edge {
        int to;
        int time;
        int level;
    };
    
    struct label_t {
        int time;
        unsigned *path;
        int *cost;
    } __attribute__((aligned(64)));
    
    int V;
    label_t *label;
    std::vector <int> contract;
    double load_time, construct_time;
    
    double GetTime(void) {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return tv.tv_sec + tv.tv_usec * 1e-6;
    }
    
    int GetLevel(int time, int dist) {
        if (time == 0) return 0;
        if (dist > time * 0.9) return 0;
        if (dist > time * 0.7) return 1;
        if (dist > time * 0.5) return 2;
        return 3;
    }
    
    inline void ChangeMin(int &now, int next) { if (next < now) now = next; }
    
    inline bool Prune(std::pair<std::vector <unsigned>, std::vector <std::pair<int, int> > > &lv, std::pair<std::vector <unsigned>, std::vector <std::pair<int, int> > > &lw, int time) {
        unsigned vfnum = 0;
        unsigned vsnum = 0;
        unsigned wfnum = 0;
        unsigned wsnum = 0;
        while (true) {
            if (vfnum == lv.first.size() || wfnum == lw.first.size()) return false;
            if ((lv.first[vfnum] & PATH_MASK) == (lw.first[wfnum] & PATH_MASK)) {
                for (unsigned i = 0, j = 0; ; ) {
                    if (lv.second[vsnum + i].first == lw.second[wsnum + j].first) {
                        if (lv.second[vsnum + i].second + lw.second[wsnum + j].second <= time) return true;
                        i++;
                        j++;
                        if (i == (lv.first[vfnum] & NUM_MASK) || j == (lw.first[wfnum] & NUM_MASK)) break;
                    } else if (lv.second[vsnum + i].first < lw.second[wsnum + j].first) {
                        if (lv.second[vsnum + i].second + lw.second[wsnum + j].second + lw.second[wsnum + j].first - lv.second[vsnum + i].first <= time) return true;
                        i++;
                        if (i == (lv.first[vfnum] & NUM_MASK)) break;
                    } else {
                        if (lv.second[vsnum + i].second + lw.second[wsnum + j].second + lv.second[vsnum + i].first - lw.second[wsnum + j].first <= time) return true;
                        j++;
                        if (j == (lw.first[wfnum] & NUM_MASK)) break;
                    }
                }
                vsnum += lv.first[vfnum] & NUM_MASK;
                vfnum++;
                wsnum += lw.first[wfnum] & NUM_MASK;
                wfnum++;
            } else if ((lv.first[vfnum] & PATH_MASK) < (lw.first[wfnum] & PATH_MASK)) {
                vsnum += lv.first[vfnum] & NUM_MASK;
                vfnum++;
            } else {
                wsnum += lw.first[wfnum] & NUM_MASK;
                wfnum++;
            }
        }
    }
};

#endif
