#ifndef MAXINFRKGSKQ_DJ



#include <iostream>
#include <queue>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <sstream>
#include <limits> // for numeric_limits

#include <set>
#include <utility> // for pair
#include <algorithm>
#include <iterator>
#include <ctime>

using namespace std;



typedef int vertex_t;
typedef float weight_t;

const weight_t max_weight = std::numeric_limits<float>::infinity();


struct neighbor {
    vertex_t target;
    weight_t weight;

    neighbor(vertex_t arg_target, weight_t arg_weight)
            : target(arg_target), weight(arg_weight) {}
};

typedef std::vector<std::vector<neighbor> > adjacency_list_t;


void DijkstraComputePaths(vertex_t source,
                          const adjacency_list_t &adjacency_list,
                          std::vector<weight_t> &min_distance) {
    int n = adjacency_list.size();
    min_distance.clear();
    min_distance.resize(n, max_weight);
    min_distance[source] = 0;
    std::set<std::pair<weight_t, vertex_t> > vertex_queue;
    vertex_queue.insert(std::make_pair(min_distance[source], source));

    while (!vertex_queue.empty()) {
        weight_t dist = vertex_queue.begin()->first;
        vertex_t u = vertex_queue.begin()->second;
        vertex_queue.erase(vertex_queue.begin());

        // Visit each edge exiting u
        const std::vector<neighbor> &neighbors = adjacency_list[u];
        for (std::vector<neighbor>::const_iterator neighbor_iter = neighbors.begin();
             neighbor_iter != neighbors.end();
             neighbor_iter++) {
            vertex_t v = neighbor_iter->target;
            weight_t weight = neighbor_iter->weight;
            weight_t distance_through_u = dist + weight;
            if (distance_through_u < min_distance[v]) {
                vertex_queue.erase(std::make_pair(min_distance[v], v));

                min_distance[v] = distance_through_u;
                vertex_queue.insert(std::make_pair(min_distance[v], v));

            }

        }
    }
}


void DijkstraComputePaths(vertex_t source,
                          const adjacency_list_t &adjacency_list,
                          std::vector<weight_t> &min_distance,
                          std::vector<vertex_t> &previous) {
    int n = adjacency_list.size();
    min_distance.clear();
    min_distance.resize(n, max_weight);
    min_distance[source] = 0;
    previous.clear();
    previous.resize(n, -1);
    std::set<std::pair<weight_t, vertex_t> > vertex_queue;
    vertex_queue.insert(std::make_pair(min_distance[source], source));

    while (!vertex_queue.empty()) {
        weight_t dist = vertex_queue.begin()->first;
        vertex_t u = vertex_queue.begin()->second;
        vertex_queue.erase(vertex_queue.begin());

        // Visit each edge exiting u
        const std::vector<neighbor> &neighbors = adjacency_list[u];
        for (std::vector<neighbor>::const_iterator neighbor_iter = neighbors.begin();
             neighbor_iter != neighbors.end();
             neighbor_iter++) {
            vertex_t v = neighbor_iter->target;
            weight_t weight = neighbor_iter->weight;
            weight_t distance_through_u = dist + weight;
            if (distance_through_u < min_distance[v]) {
                vertex_queue.erase(std::make_pair(min_distance[v], v));

                min_distance[v] = distance_through_u;
                previous[v] = u;
                vertex_queue.insert(std::make_pair(min_distance[v], v));

            }

        }
    }
}

std::list<vertex_t> DijkstraGetShortestPathTo(
        vertex_t vertex, const std::vector<vertex_t> &previous) {
    std::list<vertex_t> path;
    for (; vertex != -1; vertex = previous[vertex])
        path.push_front(vertex);
    return path;
}

int testDijisktra()
{
    ifstream fin;
    string roadPath = "../data/cal.cedge";
    fin.open(roadPath.c_str());

    int Ni, Nj;
    float dist;

    int nodeNum = 21048;
    adjacency_list_t adjacency_list(nodeNum);
    // read
    string str;
    while (getline(fin, str)) {
        istringstream issRead(str);
        issRead >> Ni >> Nj >> dist;
        adjacency_list[Ni].push_back(neighbor(Nj, dist));
        adjacency_list[Nj].push_back(neighbor(Ni, dist));
    }

    fin.close();

    //for(int i=0;i<1000;i++)
    //{
    /*std::vector<weight_t> min_distance;
    std::vector<vertex_t> previous;
    DijkstraComputePaths(Qni_1, adjacency_list, min_distance, previous);

    //std::list<vertex_t> path = DijkstraGetShortestPathTo(j, previous);
    //std::cout << "Path : ";
    //std::copy(path.begin(), path.end(), std::ostream_iterator<vertex_t>(std::cout, " "));
    //std::cout << std::endl;

*/

    int Qni = 0, Qnj = 1000;

    std::vector<weight_t> min_distance;
    std::vector<vertex_t> previous;
    DijkstraComputePaths(Qni, adjacency_list, min_distance);
    //cout << "Distance from " << Qni << " to " << Qnj << ": " << min_distance[Qnj] << endl;

    return 0;
}

#define  MAXINFRKGSKQ_DJ
#endif