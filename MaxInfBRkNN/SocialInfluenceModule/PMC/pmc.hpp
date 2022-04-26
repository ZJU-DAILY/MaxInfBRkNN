#include <vector>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include "../../results.h"
#include "../CELF_PLUS/MC.h"


using namespace std;
using namespace _MC;

class PrunedEstimater {  //可以认为是存储DAG内容的数据结构
private:
	int n_super; //强联通图的 super Node 个数
	int n_u; //social network user node 的个数

	vector<int> weight;
	vector<int> comp;
	vector<int> sigmas;   //每个顶点的weight



	vector<int> pmoc;
	vector<int> at_p;

	vector<int> up;  //?  待更新的点（user 或 poi）

	vector<bool> memo, removed;

	vector<int> es, rs;
	vector<int> at_e, at_r;

	vector<bool> visited;

	int hub;
	vector<bool> descendant, ancestor;
	bool flag;

	// line 4 - 8
	void first() {
		hub = 0;
		//选hub
		for (int i = 0; i < n_super; i++) { 
			if ((at_e[i + 1] - at_e[i]) + (at_r[i + 1] - at_r[i])
				> (at_e[hub + 1] - at_e[hub]) + (at_r[hub + 1] - at_r[hub])) {
				hub = i;
			}
		}

		descendant.resize(n_super);
		queue<int> Q;
		Q.push(hub);
		for (; !Q.empty();) {   //标记hub的后继节点
			// forall v, !remove[v]
			const int v = Q.front();
			Q.pop();
			descendant[v] = true;
			for (int i = at_e[v]; i < at_e[v + 1]; i++) {
				const int u = es[i];
				if (!descendant[u]) {
					descendant[u] = true;
					Q.push(u);
				}
			}
		}

		ancestor.resize(n_super);
		Q.push(hub);
		for (; !Q.empty();) {  //标记hub的前继节点
			const int v = Q.front();
			Q.pop();
			ancestor[v] = true;
			for (int i = at_r[v]; i < at_r[v + 1]; i++) {
				const int u = rs[i];
				if (!ancestor[u]) {
					ancestor[u] = true;
					Q.push(u);
				}
			}
		}
		ancestor[hub] = false;

		// 计算每个super node的 weight
		for (int i = 0; i < n_super; i++) {  
			sigma(i);
		}
		ancestor.assign(n_super, false);
		descendant.assign(n_super, false);

		//把所有个体（user/poi）都加入待更新名单
		for (int i = 0; i < n_u; i++) {
			up.push_back(i);
		}
		//cout<<"n_u="<<n_u<<endl;
		int ii = n_u;
	}

	void first_POI(vector<int> stores) {
		hub = 0;
		//选hub
		for (int i = 0; i < n_super; i++) {
			if ((at_e[i + 1] - at_e[i]) + (at_r[i + 1] - at_r[i])
				> (at_e[hub + 1] - at_e[hub]) + (at_r[hub + 1] - at_r[hub])) {
				hub = i;
			}
		}

		descendant.resize(n_super);
		queue<int> Q;
		Q.push(hub);
		for (; !Q.empty();) {   //标记hub的后继节点
			// forall v, !remove[v]
			const int v = Q.front();
			Q.pop();
			descendant[v] = true;
			for (int i = at_e[v]; i < at_e[v + 1]; i++) {
				const int u = es[i];
				if (!descendant[u]) {
					descendant[u] = true;
					Q.push(u);
				}
			}
		}

		ancestor.resize(n_super);
		Q.push(hub);
		for (; !Q.empty();) {  //标记hub的前继节点
			const int v = Q.front();
			Q.pop();
			ancestor[v] = true;
			for (int i = at_r[v]; i < at_r[v + 1]; i++) {
				const int u = rs[i];
				if (!ancestor[u]) {
					ancestor[u] = true;
					Q.push(u);
				}
			}
		}
		ancestor[hub] = false;

		// 计算每个super node的 weight
		for (int i = 0; i < n_super; i++) {
			sigma(i);
		}
		ancestor.assign(n_super, false);
		descendant.assign(n_super, false);

		//首次执行，把所有 poi 都加入待更新名单
		for (int i = 0; i < stores.size(); i++) {
			up.push_back(stores[i]);
		}
		//cout<<"np="<<up.size()<<endl;
		int ii = up.size();
	}


	//计算每个super node 的 sigma
	int sigma(const int v0) {  //jordan
		if (memo[v0]) {
			return sigmas[v0];
		}
		memo[v0] = true;
		if (removed[v0]) {
			return sigmas[v0] = 0;
		} else {
			int child = unique_child(v0);
			if (child == -1) {
				return sigmas[v0] = weight[v0];
			} else if (child >= 0) {
				return sigmas[v0] = sigma(child) + weight[v0];
			} else {
				int delta = 0;
				vector<int> vec;
				visited[v0] = true;
				vec.push_back(v0);
				queue<int> Q;
				Q.push(v0);
				bool prune = ancestor[v0];

				if (prune) {
					delta += sigma(hub); // pruned the rearchable super node
				}
				//执行 pruned BFS
				for (; !Q.empty();) {
					const int v = Q.front();
					Q.pop();
					if (removed[v]) {
						continue;
					}
					if (prune && descendant[v]) {
						continue;
					}
					delta += weight[v];
					for (int i = at_e[v]; i < at_e[v + 1]; i++) {
						const int u = es[i];
						if (removed[u]) {
							continue;
						}
						if (!visited[u]) {
							visited[u] = true;
							vec.push_back(u);
							Q.push(u);
						}
					}
				}
				for (int i = 0; i < vec.size(); i++) {
					visited[vec[i]] = false;
				}
				return sigmas[v0] = delta;  //jordan
			}
		}
	}


    // ? 
	inline int unique_child(const int v) {
		int outdeg = 0, child = -1;
		for (int i = at_e[v]; i < at_e[v + 1]; i++) {
			const int u = es[i];
			if (!removed[u]) {
				outdeg++;
				child = u;
			}
		}
		if (outdeg == 0) {
			return -1;
		} else if (outdeg == 1) {
			return child;
		} else {
			return -2;
		}
	}




public:

	//根据 super Node(comp)以及 super Edge(es2) 构建 DAG
	void init(const int _n, vector<pair<int, int> > &_es, // _es: DAG node 间的边集合(border edge)
							   vector<int> &_comp) {
		flag = true;
		n_super = _n;  // n: DAG node 的个数
		n_u = _comp.size();

		visited.resize(n_super, false);

		int m = _es.size(); // m : DAG node 间边（超边）的个数
		vector<int> outdeg(n_super), indeg(n_super);

		for (int i = 0; i < m; i++) {
			int a = _es[i].first;  // DAG node 1
			int b = _es[i].second; // DAG node 2
			outdeg[a]++;
			indeg[b]++;
		}
		es.resize(m, -1);
		rs.resize(m, -1);

		at_e.resize(n_super + 1, 0);
		at_r.resize(n_super + 1, 0);

		at_e[0] = at_r[0] = 0;

		//?:
		for (int i = 1; i <= n_super; i++) {
			at_e[i] = at_e[i - 1] + outdeg[i - 1];
			at_r[i] = at_r[i - 1] + indeg[i - 1];
		}

		//?:
		for (int i = 0; i < m; i++) {
			int a = _es[i].first, b = _es[i].second;
			es[at_e[a]++] = b;
			rs[at_r[b]++] = a;
		}

		at_e[0] = at_r[0] = 0;
		//?:
		for (int i = 1; i <= n_super; i++) {
			at_e[i] = at_e[i - 1] + outdeg[i - 1];
			at_r[i] = at_r[i - 1] + indeg[i - 1];
		}

		sigmas.resize(n_super);
		comp = _comp;
		vector<pair<int, int> > ps;
		for (int i = 0; i < n_u; i++) {
			ps.push_back(make_pair(comp[i], i));
		}
		sort(ps.begin(), ps.end());
		at_p.resize(n_super + 1);
		for (int i = 0; i < n_u; i++) {
			pmoc.push_back(ps[i].second);
			at_p[ps[i].first + 1]++;
		}
		for (int i = 1; i <= n_super; i++) {
			at_p[i] += at_p[i - 1];
		}

		memo.resize(n_super);
		removed.resize(n_super);

		weight.resize(n_u, 0); //最大n_u个super Node
		for (int i = 0; i < n_u; i++) {
			weight[comp[i]]++;   //每个community weight ++
		}

		first();  //jordan
	}


	void init_POI(const int _n, vector<pair<int, int> > &_es, // _es: DAG node 间的边集合(border edge)
			  vector<int> &_comp, vector<int> stores) {
		flag = true;
		n_super = _n;  // n: DAG node 的个数
		n_u = _comp.size();

		visited.resize(n_super, false);

		int m = _es.size(); // m : DAG node 间边（超边）的个数
		vector<int> outdeg(n_super), indeg(n_super);

		for (int i = 0; i < m; i++) {
			int a = _es[i].first;  // DAG node 1
			int b = _es[i].second; // DAG node 2
			outdeg[a]++;
			indeg[b]++;
		}
		es.resize(m, -1);
		rs.resize(m, -1);

		at_e.resize(n_super + 1, 0);
		at_r.resize(n_super + 1, 0);

		at_e[0] = at_r[0] = 0;

		//?:
		for (int i = 1; i <= n_super; i++) {
			at_e[i] = at_e[i - 1] + outdeg[i - 1];
			at_r[i] = at_r[i - 1] + indeg[i - 1];
		}

		//?:
		for (int i = 0; i < m; i++) {
			int a = _es[i].first, b = _es[i].second;
			es[at_e[a]++] = b;
			rs[at_r[b]++] = a;
		}

		at_e[0] = at_r[0] = 0;
		//?:
		for (int i = 1; i <= n_super; i++) {
			at_e[i] = at_e[i - 1] + outdeg[i - 1];
			at_r[i] = at_r[i - 1] + indeg[i - 1];
		}

		sigmas.resize(n_super);
		comp = _comp;
		vector<pair<int, int> > ps;
		for (int i = 0; i < n_u; i++) {
			ps.push_back(make_pair(comp[i], i));
		}
		sort(ps.begin(), ps.end());
		at_p.resize(n_super + 1);
		for (int i = 0; i < n_u; i++) {
			pmoc.push_back(ps[i].second);
			at_p[ps[i].first + 1]++;
		}
		for (int i = 1; i <= n_super; i++) {
			at_p[i] += at_p[i - 1];
		}

		memo.resize(n_super);
		removed.resize(n_super);

		weight.resize(n_u, 0); //最大n_u个super Node
		for (int i = 0; i < n_u; i++) {
			weight[comp[i]]++;   //每个community weight ++
		}

		first_POI(stores);  //jordan
	}


		// 返回某个user的weight
	int sigma1(const int v) {
		return sigma(comp[v]);
	}




	void add(int v0) {
		v0 = comp[v0];
		queue<int> Q;
		Q.push(v0);
		removed[v0] = true;
		vector<int> rm;
		for (; !Q.empty();) {
			const int v = Q.front();
			Q.pop();
			rm.push_back(v);
			for (int i = at_e[v]; i < at_e[v + 1]; i++) {
				const int u = es[i];
				if (!removed[u]) {
					Q.push(u);
					removed[u] = true;
				}
			}
		}

		up.clear();

		vector<int> vec;
		for (int i = 0; i < (int) rm.size(); i++) {
			const int v = rm[i];
			memo[v] = false; // for update()
			for (int j = at_p[v]; j < at_p[v + 1]; j++) {
				up.push_back(pmoc[j]);
			}
			for (int j = at_r[v]; j < at_r[v + 1]; j++) {
				const int u = rs[j];
				if (!removed[u] && !visited[u]) {
					visited[u] = true;
					vec.push_back(u);
					Q.push(u);
				}
			}
		}
		// reachable to removed node
		for (; !Q.empty();) {
			const int v = Q.front();
			Q.pop();
			memo[v] = false;
			for (int j = at_p[v]; j < at_p[v + 1]; j++) {
				up.push_back(pmoc[j]);
			}
			for (int i = at_r[v]; i < at_r[v + 1]; i++) {
				const int u = rs[i];
				if (!visited[u]) {
					visited[u] = true;
					vec.push_back(u);
					Q.push(u);
				}
			}
		}
		for (int i = 0; i < vec.size(); i++) {
			visited[vec[i]] = false;
		}
	}


	void addPOI(int poi, MQueryResults& results){
		for(ResultDetail rr: results[poi]){
			int u_id = rr.usr_id;
			add(u_id);
		}
	}

	//jordan
	void update(vector<long long> &sums) {  //sum就是 gain, <long: u_id, long: weight>

		for (int i = 0; i < (int) up.size(); i++) {  //up是什么？需要update的 user node
			int v = up[i];
			if (!flag) {
				sums[v] -= sigmas[comp[v]];
			}
		}
		for (int i = 0; i < (int) up.size(); i++) {
			int v = up[i];
			sums[v] += sigma1(v);
		}
		flag = false;
	}

	//new jins: 更新兴趣点影响力的marginal gain
	void updateForPOI(vector<long long> &gain, vector<int>& stores, map<int, int>& store_idxMap, MQueryResults& results) {  //sum就是 gain, <long: u_id, long: weight>

		for (int i = 0; i < (int) up.size(); i++) {  //up是什么？需要update的 user node
			int poi = up[i];
			int poi_idx = store_idxMap[poi];
			if (!flag) {
				for(ResultDetail rr: results[poi]){
					int potential_usr = rr.usr_id;
					int v = potential_usr;
					int delta_sigma = sigmas[comp[v]];
					gain[poi_idx] -= delta_sigma;
				}

			}
		}
		for (int i = 0; i < (int) up.size(); i++) {
			int poi = up[i];
			int poi_idx = store_idxMap[poi];
			for(ResultDetail rr: results[poi]){
				int potential_usr = rr.usr_id;
				int v = potential_usr;
				int delta_sigma = sigma1(v);
				gain[poi_idx] += delta_sigma;
			}

		}
		flag = false;
	}



};

class InfluenceMaximizer_PMC : public MC {
private:
	int n, m; // |V|, |E|
	int np; //|Pc| 候选兴趣点的数量
	vector<pair<pair<int, int>, double> > es; //social network


	vector<int> es1;  	//es1: 有效边的 end vertex;
	vector<int> at_e;	// at_e: 有效边的 start vertex的（在各sampling graph instance中的）累计出度

	vector<int> rs1;   // rs1: 有效边的 start vertex;
	vector<int> at_r;  // at_r: 有效边的 end vertex d的（在各sampling graph instance中的) 累计入度

	//UserList curSeedSet;
	//POIList curPOISet;


	int scc(vector<int> &comp);


public:

	InfluenceMaximizer_PMC (AnyOption* opt);

	~InfluenceMaximizer_PMC();

	vector<int> runIM(const int k,
			const int R);

	vector<int> minPOIbyPMC(vector<int>& stores, BatchResults& results, int b, int R);



	void addSocialLinkPMC(int u1, int u2, double prob) {

		//users.insert(u1);
		//users.insert(u2);
		FriendsMap* neighbors = AM->find(u1);
		if (neighbors == NULL) {
			neighbors = new FriendsMap();
			neighbors->insert(pair<UID, float>(u2, prob));
			AM->insert(u1, neighbors);
		} else {

			neighbors->insert(pair<UID, float>(u2, prob));
		}
		//以上继承 MC的 代码
		es.push_back(make_pair(make_pair(u1, u2), prob));
	}






};

// Random Number Generator
class Xorshift {
public:
	Xorshift(int seed) {
		x = _(seed, 0);
		y = _(x, 1);
		z = _(y, 2);
		w = _(z, 3);
	}

	int _(int s, int i) {
		return 1812433253 * (s ^ (s >> 30)) + i + 1;
	}

	inline int gen_int() {
		unsigned int t = x ^ (x << 11);
		x = y;
		y = z;
		z = w;
		return w = w ^ (w >> 19) ^ t ^ (t >> 8);
	}

	inline int gen_int(int n) {
		return (int) (n * gen_double());
	}

	inline double gen_double() {
		unsigned int a = ((unsigned int) gen_int()) >> 5, b =
				((unsigned int) gen_int()) >> 6;
		return (a * 67108864.0 + b) * (1.0 / (1LL << 53));
	}

private:
	unsigned int x, y, z, w;
};
