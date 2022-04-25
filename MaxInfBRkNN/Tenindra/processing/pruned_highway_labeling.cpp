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

#include "pruned_highway_labeling.h"
#include "time.h"

const int PrunedHighwayLabeling::LEVEL;
const int PrunedHighwayLabeling::INF;
const unsigned PrunedHighwayLabeling::GUARD;
const unsigned PrunedHighwayLabeling::PATH_MASK;
const unsigned PrunedHighwayLabeling::NUM_MASK;

void PrunedHighwayLabeling::ConstructLabel(const char *file) {
    Free();
    
    // read graph
    // each line should contain two vertices, travel time and geometrical length
    // treat a graph as an undirected graph
    load_time = -GetTime();
    std::vector <std::vector <edge> > graph;
    {
        std::vector <road> edges;
        FILE *in = fopen(file, "r");
        if (in == NULL) {
            fprintf(stderr, "Can't open the graph file\n");
            return;
        }
        for (int from, to, time, dist; fscanf(in, "%d %d %d %d", &from, &to, &time, &dist) != EOF; ) {
            V = std::max(V, from + 1);
            V = std::max(V, to + 1);
            edges.push_back((road){from, to, time, dist});
        }
        fclose(in);
        
        graph.resize(V);
        for (size_t i = 0; i < edges.size(); i++) {
            int from = edges[i].from;
            int to = edges[i].to;
            int time = edges[i].time;
            int dist = edges[i].dist;
            int level = GetLevel(time, dist);
            graph[from].push_back((edge){to, time, level});
            graph[to].push_back((edge){from, time, level});
        }
    }
    load_time += GetTime();
    
    // allocate memory for labels
    label = (label_t *)memalign(64, sizeof(label_t) * V);
    if (label == NULL) {
        V = 0;
        fprintf(stderr, "Out of memory\n");
        return;
    }
    for (int v = 0; v < V; v++) {
        label[v].time = 0;
        label[v].path = NULL;
        label[v].cost = NULL;
    }
    
    // contract degree one vertices and order vertices
    printf("contract degree one vertices and order vertices...\n");
    construct_time = -GetTime();
    std::vector <std::vector <int> > order(LEVEL);
    {
        // contract degree one vertices
        contract.resize(V, -1);
        for (int v = 0; v < V; v++) {
            if (graph[v].size() == 1) contract[v] = graph[v][0].to;
        }
        for (int v = 0; v < V; v++) {
            if (contract[v] != -1) {
                int to = graph[v][0].to;
                for (std::vector <edge>::iterator it = graph[to].begin(); ; it++) {
                    if (it->to == v) {
                        graph[to].erase(it);
                        break;
                    }
                }
                label[v].time = graph[v][0].time;
                std::vector <edge>().swap(graph[v]);
            }
        }
        
        // order vertices
        printf("order vertices...\n");
        std::vector <int> vertex_level(V, LEVEL);
        std::vector <std::vector <std::pair<double, int> > > rank(LEVEL);
        for (int v = 0; v < V; v++) {
            if (contract[v] == -1) {
                for (size_t i = 0; i < graph[v].size(); i++) {
                    vertex_level[v] = std::min(vertex_level[v], graph[v][i].level);
                }
                for (int i = vertex_level[v]; i < LEVEL; i++) rank[i].push_back(std::make_pair((double)rand() / RAND_MAX, v));
            }
        }
        for (int i = 0; i < LEVEL; i++) std::sort(rank[i].begin(), rank[i].end());
        for (int i = 0; i < LEVEL; i++) {
            if (i < LEVEL - 1) {
                for (size_t j = 0; j * 1000 < rank[i].size(); j++) order[i].push_back(rank[i][j].second);
            } else {
                for (size_t j = 0; j < rank[i].size(); j++) order[i].push_back(rank[i][j].second);
            }
        }
    }
    
    // pruned highway labeling
    // select path and add labels for the path
    printf("pruned highway labeling， select path and add labels for the path...\n");
    {
        unsigned path_num = 0;
        std::vector <bool> used(V, false);
        std::vector <bool> visited(V, false);
        std::vector <int> cost(V, INF);
        std::vector <int> dist(V);
        std::vector <int> parent(V, -1);
        std::vector <int> child(V, 0);
        std::vector <int> check(V);
        std::vector <int> visit(V);
        std::vector <long long> last(V, -1);
        std::vector <std::pair<std::vector <unsigned>, std::vector <std::pair<int, int> > > > tmp_label(V);
        std::priority_queue <std::pair<int, int>, std::vector <std::pair<int, int> >, std::greater <std::pair<int, int> > > path_que;
        std::priority_queue <std::pair<int, std::pair<int, int> >, std::vector <std::pair<int, std::pair<int, int> > >, std::greater <std::pair<int, std::pair<int, int> > > > label_que;
        
        for (int i = 0; i < LEVEL; i++) {
            printf("当前level =%d\n",i);
            for (size_t j = 0; j < order[i].size(); j++) {
                // select path
                int v = order[i][j];
                if (used[v]) continue;
                int check_num = 0;
                int visit_num = 0;
                cost[v] = 0;
                check[check_num++] = v;
                path_que.push(std::make_pair(0, v));
                while (!path_que.empty()) {
                    int time = path_que.top().first;
                    int now = path_que.top().second;
                    path_que.pop();
                    if (cost[now] < time) continue;
                    visited[now] = true;
                    visit[visit_num++] = now;
                    for (size_t k = 0; k < graph[now].size(); k++) {
                        int w = graph[now][k].to;
                        if (used[w]) continue;
                        if (cost[w] > time + graph[now][k].time) {
                            if (cost[w] == INF) check[check_num++] = w;
                            cost[w] = time + graph[now][k].time;
                            parent[w] = now;
                            if (graph[now][k].level <= i) {
                                child[w] = 1;
                            } else {
                                child[w] = 0;
                            }
                            if (Prune(tmp_label[v], tmp_label[w], cost[w])) continue;
                            path_que.push(std::make_pair(cost[w], w));
                        }
                    }
                }
                for (int k = visit_num - 1; k > 0; k--) child[parent[visit[k]]] += child[visit[k]];
                int end = v;
                while (true) {
                    int next = -1, max_child = 0;
                    for (size_t k = 0; k < graph[end].size(); k++) {
                        int w = graph[end][k].to;
                        if (visited[w] && parent[w] == end && child[w] > max_child) {
                            next = w;
                            max_child = child[w];
                        }
                    }
                    if (next == -1) break;
                    end = next;
                }
                
                //add labels
                printf("add labels...\n");
                if (end != v) {
                    int child_num = 0;
                    while (end != -1) {
                        if (i == LEVEL - 1 || child_num * 20 <= child[end] * 19) {
                            used[end] = true;
                            label_que.push(std::make_pair(cost[end], std::make_pair(end, end)));
                        }
                        child_num = child[end];
                        end = parent[end];
                    }
                    printf("first while finished!\n");
                    long long num = (long long)path_num << 22;
                    while (!label_que.empty()) {
                        int time = label_que.top().first;
                        int now = label_que.top().second.first;
                        int start = label_que.top().second.second;
                        label_que.pop();
                        int _size = label_que.size();
                        if(_size%10000==0)
                            printf("after pop: label_que.size=%d\n",_size);
                        if (last[now] >= num + cost[start]) continue;
                        if ((last[now] >> 32 << 32) == num) {
                            if (dist[now] + cost[start] * 2 <= time) continue;
                            if (Prune(tmp_label[start], tmp_label[now], time - cost[start])) continue;
                        } else {
                            if (!used[now] && Prune(tmp_label[start], tmp_label[now], time - cost[start])) continue;
                        }
                        last[now] = num + cost[start];
                        dist[now] = time - cost[start] * 2;
                        if (tmp_label[now].first.size() == 0 || (tmp_label[now].first.back() & PATH_MASK) != path_num) {
                            tmp_label[now].first.push_back(path_num + 1);
                        } else {
                            tmp_label[now].first.back()++;
                        }
                        tmp_label[now].second.push_back(std::make_pair(cost[start], time - cost[start]));
                        for (size_t k = 0; k < graph[now].size(); k++) {
                            int w = graph[now][k].to;
                            if (!used[w] && last[w] < num + cost[start]){
                                label_que.push(std::make_pair(time + graph[now][k].time, std::make_pair(w, start)));
                                _size = label_que.size();
                                if(_size%5000==0)
                                    printf("after add: label_que.size=%d\n",_size);
                            }


                        }
                    }
                    printf("2th while finished!\n");
                    path_num += 1 << 10;
                }
                
                for (int k = 0; k < check_num; k++) {
                    visited[check[k]] = false;
                    cost[check[k]] = INF;
                    parent[check[k]] = -1;
                    child[check[k]] = 0;
                }
            }
        }
        
        // convert labels
        printf("convert labels...\n");
        for (int v = 0; v < V; v++) {
            if (contract[v] == -1) {
                label[v].path = (unsigned *)memalign(64, sizeof(unsigned) * (tmp_label[v].first.size() + 1));
                label[v].cost = (int *)memalign(64, sizeof(int) * tmp_label[v].second.size() * 2);
                if (label[v].path == NULL || label[v].cost == NULL) {
                    Free();
                    fprintf(stderr, "Out of memory\n");
                    return;
                }
                for (size_t i = 0; i < tmp_label[v].first.size(); i++) label[v].path[i] = (tmp_label[v].first[i] & PATH_MASK) + (tmp_label[v].first[i] & NUM_MASK) * 2;
                label[v].path[tmp_label[v].first.size()] = GUARD;
                for (size_t i = 0; i < tmp_label[v].second.size(); i++) {
                    label[v].cost[i * 2] = tmp_label[v].second[i].first;
                    label[v].cost[i * 2 + 1] = tmp_label[v].second[i].second;
                }
                std::vector <unsigned>().swap(tmp_label[v].first);
                std::vector <std::pair<int, int> >().swap(tmp_label[v].second);
            }
        }
        for (int v = 0; v < V; v++) {
            if (contract[v] != -1) {
                label[v].path = label[contract[v]].path;
                label[v].cost = label[contract[v]].cost;
            }
        }
    }
    construct_time += GetTime();
}

void PrunedHighwayLabeling::LoadLabel(const char *file) {
    printf("load PHL indexing... ");
    Free();
    clock_t  load_startTime,load_endTime;
    load_startTime = clock();


    FILE *in = fopen(file, "rb");
    if (in == NULL) {
        fprintf(stderr, "Can't open the phl label file\n");
        return;
    }
    fread(&V, sizeof(int), 1, in);
    contract.resize(V);
    fread(&contract[0], sizeof(int), V, in);
    label = (label_t *)memalign(64, sizeof(label_t) * V);
    if (label == NULL) {
        V = 0;
        fprintf(stderr, "Out of memory\n");
        return;
    }
    for (int v = 0; v < V; v++) {
        label[v].time = 0;
        label[v].path = NULL;
        label[v].cost = NULL;
    }
    for (int v = 0; v < V; v++) {
        if (contract[v] != -1) {
            fread(&label[v].time, sizeof(int), 1, in);
        } else {
            int fnum, snum;
            fread(&fnum, sizeof(int), 1, in);
            fread(&snum, sizeof(int), 1, in);
            label[v].path = (unsigned *)memalign(64, sizeof(unsigned) * fnum);
            label[v].cost = (int *)memalign(64, sizeof(int) * snum);
            if (label[v].path == NULL || label[v].cost == NULL) {
                Free();
                fprintf(stderr, "Out of memory\n");
                return;
            }
            fread(label[v].path, sizeof(unsigned), fnum, in);
            fread(label[v].cost, sizeof(int), snum, in);
        }
    }
    for (int v = 0; v < V; v++) {
        if (contract[v] != -1) {
            label[v].path = label[contract[v]].path;
            label[v].cost = label[contract[v]].cost;
        }
    }
    fclose(in);
    printf("DONE!");
    load_endTime = clock();
    double time = (double)(load_endTime-load_startTime)/CLOCKS_PER_SEC;
    printf("%f sec!\n",time);
}

void PrunedHighwayLabeling::StoreLabel(const char *file) {
    FILE *out = fopen(file, "wb");
    if (out == NULL) {
        fprintf(stderr, "Can't open the label file\n");
        return;
    }
    fwrite(&V, sizeof(int), 1, out);
    fwrite(&contract[0], sizeof(int), V, out);
    for (int v = 0; v < V; v++) {
        if (contract[v] != -1) {
            fwrite(&label[v].time, sizeof(int), 1, out);
        } else {
            int fnum = 1, snum = 0;
            for (int i = 0; ; i++) {
                if (label[v].path[i] == GUARD) break;
                fnum++;
                snum += label[v].path[i] & NUM_MASK;
            }
            fwrite(&fnum, sizeof(int), 1, out);
            fwrite(&snum, sizeof(int), 1, out);
            fwrite(label[v].path, sizeof(unsigned), fnum, out);
            fwrite(label[v].cost, sizeof(int), snum, out);
        }
    }
    fclose(out);
}

void PrunedHighwayLabeling::Free(void) {
    for (int v = 0; v < V; v++) {
        if (contract[v] == -1) {
            free(label[v].path);
            free(label[v].cost);
            label[v].path = NULL;
            label[v].cost = NULL;
        }
    }
    free(label);
    V = 0;
    label = NULL;
}

void PrunedHighwayLabeling::Statistics(void) {
    long long sum_label = 0, sum_memory = 0;
    
    for (int v = 0; v < V; v++) {
        if (contract[v] != -1) continue;
        for (int i = 0; ; i++) {
            if (label[v].path[i] == GUARD) break;
            sum_label += (label[v].path[i] & NUM_MASK) / 2;
            sum_memory += sizeof(unsigned);
        }
        sum_memory += sizeof(unsigned);
    }
    sum_memory += sizeof(int) * sum_label * 2;
    
    printf("Load Time : %lf sec\n", load_time);
    printf("Construct Time : %lf sec\n", construct_time);
    printf("Average Label Size : %lld\n", sum_label / V);
    printf("Memory Size : %lld byte\n", sum_memory);
}

double PrunedHighwayLabeling::computeIndexSize() {
    long long sum_label = 0, sum_memory = 0;
    for (int v = 0; v < V; v++) {
        if (contract[v] != -1) continue;
        for (int i = 0; ; i++) {
            if (label[v].path[i] == GUARD) break;
            sum_label += (label[v].path[i] & NUM_MASK) / 2;
            sum_memory += sizeof(unsigned);
        }
        sum_memory += sizeof(unsigned);
    }
    sum_memory += sizeof(int) * sum_label * 2;
    double memoryUsage = sum_memory;
    return memoryUsage/(1024*1024); // in MB
}

double PrunedHighwayLabeling::getConstructionTime() {
    return construct_time;
}
