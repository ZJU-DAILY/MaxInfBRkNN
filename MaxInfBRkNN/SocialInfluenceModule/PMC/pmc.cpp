
#include "pmc.hpp"

using namespace std;


InfluenceMaximizer_PMC::InfluenceMaximizer_PMC(AnyOption *opt) : MC(opt){
	//MC(opt);
}

InfluenceMaximizer_PMC::~InfluenceMaximizer_PMC() {

}

vector<int> InfluenceMaximizer_PMC::runIM(
		const int k, const int R) {
	n = 0;
	m = es.size();
	for (int i = 0; i < (int) es.size(); i++) {
		n = max(n, max(es[i].first.first, es[i].first.second) + 1);
	}

	sort(es.begin(), es.end());

	es1.resize(m);
	rs1.resize(m);
	at_e.resize(n + 1);
	at_r.resize(n + 1);

	vector<PrunedEstimater> infs(R);
	vector<int> seeds;

	// sampling processing and build DAG for each instance
	for (int t = 0; t < R; t++) {
		Xorshift xs = Xorshift(t);

		int mp = 0;
		at_e.assign(n + 1, 0);
		at_r.assign(n + 1, 0);
		vector<pair<int, int> > ps;
		// coin flip for each edge
		for (int i = 0; i < m; i++) {
			if (xs.gen_double() < es[i].second) {
				es1[mp++] = es[i].first.second;
				at_e[es[i].first.first + 1]++;
				ps.push_back(make_pair(es[i].first.second, es[i].first.first));
			}
		}
		at_e[0] = 0;
		sort(ps.begin(), ps.end());
		//compute SCC and obtain DAG
		for (int i = 0; i < mp; i++) {
			rs1[i] = ps[i].second;
			at_r[ps[i].first + 1]++;
		}
		for (int i = 1; i <= n; i++) {
			at_e[i] += at_e[i - 1];
			at_r[i] += at_r[i - 1];
		}

		vector<int> comp(n);

		int nscc = scc(comp);

		vector<pair<int, int> > es2;
		for (int u = 0; u < n; u++) {
			int a = comp[u];
			for (int i = at_e[u]; i < at_e[u + 1]; i++) {
				int b = comp[es1[i]];
				if (a != b) {
					es2.push_back(make_pair(a, b));   //强连通图的超边
				}
			}
		}

		sort(es2.begin(), es2.end());
		es2.erase(unique(es2.begin(), es2.end()), es2.end());

		infs[t].init(nscc, es2, comp);  //建立第t个DAG的 PrunedEstimater
	}

	vector<long long> gain(n);  // jordan (gain)
	vector<int> S;

	for (int t = 0; t < k; t++) {
		for (int j = 0; j < R; j++) {
			infs[j].update(gain);  //jordan (update)
		}
		int next = 0;
		for (int i = 0; i < n; i++) {
			if (gain[i] > gain[next]) {
				next = i;
			}
		}

		S.push_back(next);
		for (int j = 0; j < R; j++) {
			infs[j].add(next);
		}
		seeds.push_back(next);
	}
	for(int u:seeds){
		curSeedSet.insert(u);
	}
	return seeds;
}



int InfluenceMaximizer_PMC::scc(vector<int> &comp) {  //计算scc
	vector<bool> vis(n);
	stack<pair<int, int> > S;
	vector<int> lis;
	int k = 0;
	for (int i = 0; i < n; i++) {
		S.push(make_pair(i, 0));  //second 可能是 comp的标记
	}
	for (; !S.empty();) {
		int v = S.top().first, state = S.top().second;
		S.pop();
		if (state == 0) {
			if (vis[v]) {
				continue;
			}
			vis[v] = true;
			S.push(make_pair(v, 1));
			for (int i = at_e[v]; i < at_e[v + 1]; i++) {
				int u = es1[i];
				S.push(make_pair(u, 0));
			}
		} else {  // state != 0
			lis.push_back(v);
		}
	}
	for (int i = 0; i < n; i++) {
		S.push(make_pair(lis[i], -1));
	}
	vis.assign(n, false);
	for (; !S.empty();) {
		int v = S.top().first, arg = S.top().second;
		S.pop();
		if (vis[v]) {
			continue;
		}
		vis[v] = true;
		comp[v] = arg == -1 ? k++ : arg;
		for (int i = at_r[v]; i < at_r[v + 1]; i++) {
			int u = rs1[i];
			S.push(make_pair(u, comp[v]));
		}
	}
	return k;  //k是comp的数量
}


//new jins
vector<int> InfluenceMaximizer_PMC:: minPOIbyPMC(vector<int>& stores, BatchResults& results, int b, int R){
	cout<<"mining promising poi by PMC..."<<endl;
	int k = b;
	n = 0;
	m = es.size();
	for (int i = 0; i < (int) es.size(); i++) {
		n = max(n, max(es[i].first.first, es[i].first.second) + 1);
	}

	np = stores.size();  // the number of chain stores

	int _n = n;
	int _np = np;


	sort(es.begin(), es.end());

	es1.resize(m);
	rs1.resize(m);
	at_e.resize(n + 1);
	at_r.resize(n + 1);

	vector<PrunedEstimater> infs(R);
	vector<int> pois;

	map<int,int> store_idxMap;   //<id, idx>
	map<int,int> store_IDMap;   //<idx, id>
	for(int i=0;i<stores.size();i++){
		int store_id = stores[i];
		store_idxMap[store_id] = i;
	}


	// sampling stage for DAG generation
	for (int t = 0; t < R; t++) {
		Xorshift xs = Xorshift(t);

		int mp = 0;
		at_e.assign(n + 1, 0);
		at_r.assign(n + 1, 0);
		vector<pair<int, int> > ps; // 依然存在的边集合 the survived edge set
		// 对每条边进行coin flip， 保留并记录有效边的数目
		for (int i = 0; i < m; i++) {
			if (xs.gen_double() < es[i].second) {
				es1[mp++] = es[i].first.second;
				at_e[es[i].first.first + 1]++; //at_e: 有效边的 start vertex的（在各sampling graph instance中的）累计出度
				ps.push_back(make_pair(es[i].first.second, es[i].first.first));  // (end, start)
			}
		}
		at_e[0] = 0;
		sort(ps.begin(), ps.end());

		//进行计算 SCC 前的标记 ？
		for (int i = 0; i < mp; i++) {
			rs1[i] = ps[i].second;
			at_r[ps[i].first + 1]++;  //有效边的end vertex （在所有的累计入度
		}
		for (int i = 1; i <= n; i++) {
			at_e[i] += at_e[i - 1];
			at_r[i] += at_r[i - 1];
		}

		vector<int> comp(n);

		//compute SCC （得到DAG的（超） 点)，以及其个数
		int nscc = scc(comp);

		vector<pair<int, int> > es2;  // the edges between super nodes in DAG

		//构建 DAG 的（超） 边
		for (int u = 0; u < n; u++) {
			int a = comp[u];  // super node in DAG
			for (int i = at_e[u]; i < at_e[u + 1]; i++) {
				int b = comp[es1[i]];
				if (a != b) {
					es2.push_back(make_pair(a, b));   //强连通图间的边
				}
			}
		}

		sort(es2.begin(), es2.end());
		es2.erase(unique(es2.begin(), es2.end()), es2.end());

		//根据 super Node(comp)以及 super Edge(es2) 构建 DAG, 并计算其中个super node 的 weight
		infs[t].init_POI(nscc, es2, comp,stores);
	}

	//line 10, begin selection, （以上都不用改动， 以下需要改动）

	vector<long long> gain(np);  //这里要改为np jordan  <long: ,   long: >
	vector<int> S;
	int poi_size = stores.size();
	bool poi_selectFlag[poi_size];
	for(int i=0;i<stores.size();i++){
		poi_selectFlag[i]= false;
	}

	for (int t = 0; t < k; t++) {
		//cout<<"----------------------the "<<t<<" th seletion iteration----------------------"<<endl;
		for (int j = 0; j < R; j++) {
			infs[j].updateForPOI(gain,stores,store_idxMap, results);  //jordan, 第一次更新各个user(poi) 的gain, 需要改！
		}
		int next_th = 0; long current_max = -1;
		for (int i = 0; i < stores.size(); i++) {  //max{i,v}, 这里的 i<n 要改为 i<np, 选出有当前最大marginal gain 的 user(poi)
			if(poi_selectFlag[i]) {
				//cout<<"i="<<i<<"continue!"<<endl;
				continue;
			}

			else if (gain[i] > current_max) {
				next_th = i;
				current_max = gain[i];
				//cout<<"gain["<<i<<"] > current_max"<<endl;
			}
		}

		S.push_back(next_th);
		int poi_nextID = stores[next_th];
		poi_selectFlag[next_th] = true;
		//cout<<"set poi_selectFlag["<<next_th<<"]="<<poi_selectFlag[next_th]<<endl;
		//cout<<"selecting "<<"next_th="<<next_th<<" poi,poi_id="<<poi_nextID<<endl;

		for (int j = 0; j < R; j++) {
			//infs[j].add(next);   //next为选中的 user node， 需要改！
			infs[j].addPOI(poi_nextID, results);
		}

		pois.push_back(poi_nextID);
		//将poi selection 结果加入
		curPOISet.insert(poi_nextID);
		for(ResultDetail rr: results[poi_nextID]){  //将poi的potential customers 加入 seed set
			int u_id = rr.usr_id;
			curSeedSet.insert(u_id);
		}

	}
	curSeedSet.clear();
	for(int i=30;i<70;i++)
        curSeedSet.insert(i);
	return pois;
}

