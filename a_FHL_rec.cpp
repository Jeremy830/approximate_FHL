#include<cstdio>
#include<cstring>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<set>
#include<map>
#include<queue>
#include<algorithm>
//#include<ctime>
#include<sstream>
#include<chrono>
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include "Semaphore.h"


using namespace std;
using namespace chrono;

/*********************please set the following parameters*************************/
const int ism=30;//input the number of threads
const int fn = 75;//input the number of subgraphs
const int All_v = 264347;//input the number of vertices in the original graph
const double a_g=1.2;//input the approximation ratio
const int sw_method=99;/*input the method of skyline computation                                                
[sw_method: 
0 for building index, 
10 for building F1 index(optimal), 
99 for building F2-Recycle(Final), 
1 for reading index existed
100,1001 for testing others*/
/**********************************************************************************/

/*NOTE1: the input file should be in the folder "data", which include road network, 
graph partitions, boundary vertices, and skylines between any two boundaries in each 
partition. The example data file can be obtained via:
https://www.dropbox.com/scl/fi/hiuv1rn6q324qr2scdaqi/data.zip?rlkey=rwq1lzgg0glekhlpl5ql5n94j&dl=0
*/

/*Note2: crease a folder "index" in the same folder of the executable file to store the index*/




Semaphore* sm = new Semaphore(ism);
//Semaphore* sm2 = new Semaphore(1);
int tree_width = 0;
const int INF = 999999999;
vector<map<int, int>> VMap;
vector<pair<int, int>> tmp_rectangle;



struct Graph
{
	int n, m, sum;
	vector<int> V;
	vector<int> Orign_V;
	vector< map< int, int > > W;
	vector< map< int, int > > C;
	vector<pair<int, int>> es;
	vector<pair<int, int>> Orign_es;
	vector<vector<int>> r_VMap;
	pair<int, int> s;
	pair<int, int> Orign_s;
	vector < map<int, vector<pair<int, int>>>> skylines, Edge,_Edge;
	vector < map<int, vector<pair<int, int>>>> Orign_Edge;
	vector<vector<pair<int, vector<pair<int, int>>>>> Ed;
	vector<vector<int>> bound_ord;
	vector<vector<int>> lambda;
	int* DD, * DD2, * NUM;
	int* _DD, * _DD2;
	bool* changed;
	vector<int> D, Orign_D, _D;


	Graph() {
		n = m = 0;
		sum = 0;
		V.clear();
		W.clear();
		C.clear();
		Edge.clear();
		_Edge.clear();
		bound_ord.clear();
		Orign_V.clear();
		Orign_Edge.clear();
	}
	Graph(char* file, char* subFile, int sub_num, int func) {
		//	cout << "file:" << file << endl;
		vector< map< int, int > > Orign_W;
		vector< map< int, int > > Orign_C;
		Orign_W.clear();
		Orign_C.clear();
		Graph();
		ifstream ifs(file);
		if (!ifs.is_open()) {
			cerr << "open file error!" << endl;
			return;
		}
		ifs >> sum >> m;
		map<int, int> tmp;

		for (int i = 0; i <= sum; i++) {
			map< int, int > v;
			v.clear();
			Orign_W.push_back(v);
			Orign_C.push_back(v);
		}
		Orign_Edge.resize(sum + 1);


		for (int i = 0; i < m; i++) {
			int x, y, w, c = 0;
			ifs >> x >> y >> w >> c;
			if (x > sum || y > sum)
				continue;
			if (Orign_W[x].find(y) != Orign_W[x].end() || Orign_C[x].find(y) != Orign_C[x].end()) {
				if (Orign_W[x][y] > w) {
					Orign_W[x][y] = w;
					Orign_W[y][x] = w;
				}
				if (Orign_C[x][y] > c) {
					Orign_C[x][y] = c;
					Orign_C[y][x] = c;
				}
			}
			else {
				Orign_W[x].insert(make_pair(y, w));
				Orign_W[y].insert(make_pair(x, w));
				Orign_C[x].insert(make_pair(y, c));
				Orign_C[y].insert(make_pair(x, c));
				Orign_s = make_pair(w, c);
				Orign_es.push_back(Orign_s);

				Orign_Edge[x].insert(make_pair(y, Orign_es));
				Orign_Edge[y].insert(make_pair(x, Orign_es));

				Orign_es.clear();
			}
		}
		ifs.close();
		Orign_D.clear();
		Orign_D.push_back(0);
		for (int i = 1; i <= sum; i++){
			Orign_D.push_back(Orign_Edge[i].size());
		}
		bound_ord.resize(fn);

		switch (func)
		{
		case 0:
		   Sub_Graph(subFile, sub_num, Orign_W, Orign_C);
		   break;
		case 1:
		   cout << "Start to build bound tree..." << endl;
		   Bound_Graph(subFile, sub_num, Orign_W, Orign_C);
		   break;
		default:
			cerr << "Error: Wrong function number!" << endl;
			break;
		}
	}

	void Sub_Graph(char* file, int subN, vector< map< int, int > > Orign_W, vector< map< int, int > > Orign_C) {		
		ifstream ifs(file);
		int tem;
		int _tem = 0;
		if (!ifs.is_open()) {
			cerr << "open file error!" << endl;
			return;
		}
		ifs >> n;
		r_VMap.resize(fn);

		vector<int> t_l;
		t_l.clear();
		lambda.clear();
		t_l.assign(n+1,0);
		lambda.assign(n+1,t_l);

		for (int i = 0; i < n + 1; i++) {
			r_VMap[subN].push_back(_tem);
		}
		for (int i = 0; i < n; i++) {
			ifs >> tem;
			VMap[tem].insert(make_pair(subN, i + 1));
			r_VMap[subN][i + 1] = tem;
		}
		ifs.close();

		for (int i = 0; i <= n; i++) {
			map< int, int > v;
			v.clear();
			W.push_back(v);
			C.push_back(v);
		}
		bool check;
		Edge.resize(n + 1);
		int x, y, w, c = 0;
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				check = isEdgeExist_Orign(r_VMap[subN][i + 1], r_VMap[subN][j + 1]);
				if (check) {
					x = i + 1;
					y = j + 1;
					w = Orign_W[r_VMap[subN][i + 1]][r_VMap[subN][j + 1]];
					c = Orign_C[r_VMap[subN][i + 1]][r_VMap[subN][j + 1]];
					if (x > n || y > n)
						continue;
					if (W[x].find(y) != W[x].end() || C[x].find(y) != C[x].end()) {
						if (W[x][y] > w) {
							W[x][y] = w;
							W[y][x] = w;
						}
						if (C[x][y] > c) {
							C[x][y] = c;
							C[y][x] = c;
						}
					}
					else {
						W[x].insert(make_pair(y, w));
						W[y].insert(make_pair(x, w));
						C[x].insert(make_pair(y, c));
						C[y].insert(make_pair(x, c));
						s = make_pair(w, c);
						es.push_back(s);

						Edge[x].insert(make_pair(y, es));
						Edge[y].insert(make_pair(x, es));

						es.clear();
					}
				}
			}
		}
		vector<int> it_bounds;
		int tmp_w = 1;
		int tmp_c = 1;
		string bb = "";
		bb = to_string(static_cast<long long>(subN));
		ifstream in("data/sub_bound_skylines" + bb + ".txt");
		char str[255];
		int cnt = 0;
		while (!in.eof()) {
			in.getline(str, sizeof(str));
			cnt++;
		}
		in.close();
		ifstream _in("data/sub_bound_skylines" + bb + ".txt");
		//cout << "cnt:" << cnt << endl;
		for (int i = 0; i < cnt - 1; i++) {
			int x, y, w, c = 0;
			_in >> x >> y >> w >> c;
			if(x==INF||y==INF||w==INF||c==INF){
				cerr<<"Error: INF in bound file!"<<endl;
			}
			x = VMap[x][subN];
			y = VMap[y][subN];
			s = make_pair(w, c);
			Edge[x][y].push_back(s);
			Edge[y][x].push_back(s);
		}
		_in.close();
		D.clear();
		_D.clear();
		D.push_back(0);
		for (int i = 1; i <= n; i++){
			D.push_back(Edge[i].size());
		}
		_D=D;
		_Edge=Edge;
	}

	void Bound_Graph(char* file, int subN, vector< map< int, int > > Orign_W, vector< map< int, int > > Orign_C) {
		cout << "file:" << file << endl;
		ifstream ifs(file);
		int tem;
		int _tem = 0;

		if (!ifs.is_open()) {
			cout << "open file error!" << endl;
		}
		ifs >> n;
		r_VMap.resize(fn + 1);

		vector<int> t_l;
		t_l.clear();
		lambda.clear();
		t_l.assign(n+1,0);
		lambda.assign(n+1,t_l);
		for (int i = 0; i < n + 1; i++) {
			r_VMap[subN].push_back(_tem);
		}
		for (int i = 0; i < n; i++) {
			ifs >> tem;
			VMap[tem].insert(make_pair(subN, i + 1));
			r_VMap[subN][i + 1] = tem;
		}
		ifs.close();

		for (int i = 0; i <= n; i++) {
			map< int, int > v;
			v.clear();
			W.push_back(v);
			C.push_back(v);
		}
		bool check;
		Edge.resize(n + 1);
		ifs.close();
		int bn, bm = 0;
		ifstream bifs("data/boundData.txt");
		if (!bifs.is_open()) {
			cout << "open boundData error!" << endl;
		}
		bifs >> bn >> bm;
		for (int i = 0; i < bm; i++) {
			int x, y, w, c = 0;
			bifs >> x >> y >> w >> c;
			x = VMap[x][subN];
			y = VMap[y][subN];

			int rep = 0;
			if (Edge[x].find(y) != Edge[x].end() || Edge[x].find(y) != Edge[x].end()) {
				for (int size_x = 0; size_x < Edge[x][y].size(); size_x++) {
					if (Edge[x][y][size_x].first == w && Edge[x][y][size_x].second == c) {
						rep++;
					}
				}
			}
			if (rep == 0) {
				s = make_pair(w, c);
				Edge[x][y].push_back(s);
				Edge[y][x].push_back(s);
			}
		}

		D.clear();
		_D.clear();
		D.push_back(0);
		for (int i = 1; i <= n; i++){
			D.push_back(Edge[i].size());
		}
		_D=D;
		_Edge=Edge;
	}

	void EdgeInitialize() {
		Ed.clear();
		for (int i = 0; i <= n; i++) {
			vector<pair<int, vector<pair<int, int>>>>  ed;
			ed.clear();
			for (map<int, vector<pair<int, int>>>::iterator it = Edge[i].begin(); it != Edge[i].end(); it++) {
				ed.push_back(*it);
			}
			Ed.push_back(ed);
		}
	}
	bool isEdgeExist(int u, int v) {
		if (Edge[u].find(v) == Edge[u].end())
			return false;
		else return true;
	}

	bool isEdgeExist_Orign(int u, int v) {
		if (Orign_Edge[u].find(v) == Orign_Edge[u].end())
			return false;
		else return true;
	}

	void insertEdge(int u, int v, vector<pair<int, int>> e) {
		if (Edge[u].find(v) != Edge[u].end()) return;
		Edge[u].insert(make_pair(v, e));
		Edge[v].insert(make_pair(u, e));
		D[u]++;
		D[v]++;
	}

	void deleteEdge(int u, int v) {
		if (Edge[u].find(v) == Edge[u].end()) return;

		Edge[u].erase(Edge[u].find(v));
		Edge[v].erase(Edge[v].find(u));

		D[u]--;
		D[v]--;
	}

	void _insertEdge(int u, int v, vector<pair<int, int>> e) {
		if (_Edge[u].find(v) != _Edge[u].end()) return;
		_Edge[u].insert(make_pair(v, e));
		_Edge[v].insert(make_pair(u, e));
		_D[u]++;
		_D[v]++;
	}

	void _deleteEdge(int u, int v) {
		if (_Edge[u].find(v) == _Edge[u].end()) return;

		_Edge[u].erase(_Edge[u].find(v));
		_Edge[v].erase(_Edge[v].find(u));

		_D[u]--;
		_D[v]--;
	}

	bool _isEdgeExist(int u, int v) {
		if (_Edge[u].find(v) == _Edge[u].end())
			return false;
		else return true;
	}
};

struct SelEle {
	int x;
	Graph* G;
	SelEle();
	SelEle(int _x, Graph* _g) {
		x = _x;
		G = _g;
	}
	bool operator< (const SelEle se)const {
		if (G->DD[x] != G->DD[se.x])
			return G->DD[x] < G->DD[se.x];
		if (G->DD2[x] != G->DD2[se.x])
			return G->DD2[x] < G->DD2[se.x];
		return x < se.x;
	}
};

struct Node {
	vector<int> vert, pos;
	vector<int> ch;
	vector<vector<pair<int, int>>> skyToAnc, SK;
	double a_alc;
	double a_rest;
	int height;
	int pa;
	int uniqueVertex;
	Node() {
		vert.clear();
		pos.clear();
		skyToAnc.clear();
		ch.clear();
		pa = -1;
		uniqueVertex = -1;
		height = 0;
		a_alc=1;
		a_rest=a_g;
	}
};

struct Sub_Tree_Decomposition {
	Graph G, H;
	set<SelEle> deg;
	int maxSize;
	Sub_Tree_Decomposition() {}
	vector<vector<int>> neighbor;
	vector<vector<vector<pair<int, int>>>> adj_skyline;
	vector<vector<vector<pair<int, int>>>> EXL;
	vector<int> ord;
	int heightMax;
	int max_lam=0;
	vector<vector<double>> rec,alc;
	double a_g_s=pow(a_g,0.25);

	struct judge {
		bool operator()(const pair<int, int> a, const pair<int, int> b) {
			return a.first < b.first;
		};
	};

	struct judgeCost {
		bool operator()(const pair<int, int> a, const pair<int, int> b) {
			return a.second < b.second;
		};
	};

	struct hop {
		int i;
		int j;
		pair<int, int> val;
		hop(int x, int y, pair<int, int> v) : i(x), j(y), val(v) {};
	};

	struct pk {
		bool operator()(const hop& a, const hop& b) const {
			return a.val.first < b.val.first;
		}
	};
	struct pk_c {
		bool operator()(const hop& a, const hop& b) const {
			return a.val.second < b.val.second;
		}
	};

	//exact version
	void mergeSort(const vector<pair<int, int>>& s1, const vector<pair<int, int>>& s2, vector<pair<int, int>>& result) {
		int size1 = s1.size();
		int size2 = s2.size();
		result.resize(size1 + size2);
		int tmp_index, tmp_1_index, tmp_2_index;
		tmp_1_index = 0;
		tmp_2_index = 0;
		for (tmp_index = 0; (tmp_1_index < size1) && (tmp_2_index < size2); tmp_index++) {
			if (s1[tmp_1_index].first <= s2[tmp_2_index].first) {

				result[tmp_index] = s1[tmp_1_index];
				tmp_1_index++;
			}
			else {
				result[tmp_index] = s2[tmp_2_index];
				tmp_2_index++;
			}
		}
		if (tmp_1_index == size1) {
			while (tmp_2_index < size2) {
				result[tmp_index] = s2[tmp_2_index];
				tmp_index++;
				tmp_2_index++;
			}
		}
		else if (tmp_2_index == size2) {
			while (tmp_1_index < size1) {
				result[tmp_index] = s1[tmp_1_index];
				tmp_index++;
				tmp_1_index++;
			}
		}
	}

    //exact version
	vector<pair<int, int>> screenSkylines(vector<pair<int, int>>& input) {
		vector<pair<int, int>> result;
		result.push_back(input[0]);
		//cout << "input el:" << input[0].first << ' ' << input[0].second << endl;
		int size = input.size();

		if (size <= 1) {
			return result;
		}
		else {
			for (int i = 1; i < size; i++) {
				if (result[result.size() - 1].first == input[i].first && result[result.size() - 1].second > input[i].second) {
					result[result.size() - 1].second = input[i].second;
				}
				if (result[result.size() - 1].first != input[i].first && result[result.size() - 1].second > input[i].second) {
					result.push_back(input[i]);
				}
				else {
					continue;
				}
			}
			return result;
		}
	}

    //approximate version
	void a_mergeSort(const vector<pair<int, int>>& s1, const vector<pair<int, int>>& s2, vector<pair<int, int>>& result) {
		int size1 = s1.size();
		int size2 = s2.size();
		result.resize(size1 + size2);
		int tmp_index, tmp_1_index, tmp_2_index;
		tmp_1_index = 0;
		tmp_2_index = 0;
		for (tmp_index = 0; (tmp_1_index < size1) && (tmp_2_index < size2); tmp_index++) {
			if (s1[tmp_1_index].second <= s2[tmp_2_index].second) {

				result[tmp_index] = s1[tmp_1_index];
				tmp_1_index++;
			}
			else {
				result[tmp_index] = s2[tmp_2_index];
				tmp_2_index++;
			}
		}
		if (tmp_1_index == size1) {
			while (tmp_2_index < size2) {
				result[tmp_index] = s2[tmp_2_index];
				tmp_index++;
				tmp_2_index++;
			}
		}
		else if (tmp_2_index == size2) {
			while (tmp_1_index < size1) {
				result[tmp_index] = s1[tmp_1_index];
				tmp_index++;
				tmp_1_index++;
			}
		}
	}

    //approximate version
	vector<pair<int, int>> a_screenSkylines(vector<pair<int, int>>& input, double a_l) {
		vector<pair<int, int>> result;
		result.push_back(input[0]);
		int size = input.size();

		if (size <= 1) {
			return result;
		}
		else {
			for (int i = 1; i < size; i++) {
				if (result[result.size() - 1].second== input[i].second && result[result.size() - 1].first > input[i].first) {
					result[result.size() - 1].first = input[i].first;
				}
				if (result[result.size() - 1].second != input[i].second && result[result.size() - 1].first > a_l*input[i].first) {
					result.push_back(input[i]);
				}
				else {
					continue;
				}
			}
			return result;
		}
	}

	//approximate version
	vector<pair<int, int>> al_screenSkylines(vector<pair<int, int>>& input, double a_l, int p, int anc) {
		vector<pair<int, int>> result;
		result.push_back(input[0]);
		int size = input.size();

		if (size <= 1) {
			return result;
		}
		else {
			for (int i = 1; i < size; i++) {
				if (result[result.size() - 1].second== input[i].second && result[result.size() - 1].first > input[i].first) {
					result[result.size() - 1].first = input[i].first;
				}
				if (result[result.size() - 1].second != input[i].second && result[result.size() - 1].first > a_l*input[i].first) {
					result.push_back(input[i]);
				}
				else {
					if(result[result.size() - 1].first > input[i].first){
						double tmp;
						tmp=double(result[result.size()-1].first)/double(input[i].first);
						if (rec[p][anc]>tmp){
							rec[p][anc]=tmp;
						}
					}
				}
			}
			if(rec[p][anc]<a_l){
				rec[p][anc]=a_l/rec[p][anc];
			}
			else{
				rec[p][anc]=1;
			}
			return result;
		}
	}

	//approximate version
	vector<pair<int, int>> a_concatenation(const vector<pair<int, int>>& h1, const vector<pair<int, int>>& h2, const double a_l) {
		int size1 = h1.size();
		int size2 = h2.size();
		pair<int, int> wc;
		vector<pair<int, int>> result;
		if (size1 == 0 || size2 == 0) {
			result.push_back(make_pair(INF, INF));
			return result;
		}

		vector<vector<bool>> visited(size1, vector<bool>(size2, false));
		priority_queue<hop, vector<hop>, pk_c> pq;
		int w = h1[0].first + h2[0].first;
		int c = h1[0].second + h2[0].second;
		wc = make_pair(w, c);
		pq.push(hop(0, 0, wc));
		visited[0][0] = true;
		int count = 0;
		result.push_back(wc);
		while (!pq.empty()) {
			hop el = pq.top();
			pq.pop();
			int n = result.size();



			if (a_l*el.val.first < result[n - 1].first) {
				result.push_back(el.val);
			}


			int ex1 = el.i + 1;
			int ey1 = el.j;
			int ex2 = el.i;
			int ey2 = el.j + 1;
			if (ex1 < size1) {
				if (!visited[ex1][ey1]) {
					w = h1[ex1].first + h2[ey1].first;
					c = h1[ex1].second + h2[ey1].second;
					wc = make_pair(w, c);
					if (a_l*w < result[n - 1].first) {
						pq.push(hop(ex1, ey1, wc));
						visited[ex1][ey1]=true;
					}
				}
			}
			if (ey2 < size2) {
				if (!visited[ex2][ey2]) {
					w = h1[ex2].first + h2[ey2].first;
					c = h1[ex2].second + h2[ey2].second;
					wc = make_pair(w, c);
					if (a_l*w < result[n - 1].first) {
						pq.push(hop(ex2, ey2, wc));
						visited[ex2][ey2]=true;
					}
				}
			}

		}
		return result;
	}

	//exact version
	vector<pair<int, int>> concatenation(const vector<pair<int, int>>& h1, const vector<pair<int, int>>& h2) {
		int size1 = h1.size();
		int size2 = h2.size();

		pair<int, int> wc;
		vector<pair<int, int>> result;
		if (size1 == 0 || size2 == 0) {
			result.push_back(make_pair(INF, INF));
			return result;
		}

		vector<vector<bool>> visited(size1, vector<bool>(size2, false));
		priority_queue<hop, vector<hop>, pk> pq;
		int w = h1[0].first + h2[0].first;
		int c = h1[0].second + h2[0].second;
		wc = make_pair(w, c);
		pq.push(hop(0, 0, wc));
		visited[0][0] = true;
		int count = 0;
		result.push_back(wc);
		while (!pq.empty()) {
			hop el = pq.top();
			pq.pop();
			int n = result.size();



			if (el.val.second < result[n - 1].second) {
				result.push_back(el.val);
			}


			int ex1 = el.i + 1;
			int ey1 = el.j;
			int ex2 = el.i;
			int ey2 = el.j + 1;
			if (ex1 < size1) {
				if (!visited[ex1][ey1]) {
					w = h1[ex1].first + h2[ey1].first;
					c = h1[ex1].second + h2[ey1].second;
					wc = make_pair(w, c);
					if (c < result[n - 1].second) {
						pq.push(hop(ex1, ey1, wc));
						visited[ex1][ey1]=true;
					}
				}
			}
			if (ey2 < size2) {
				if (!visited[ex2][ey2]) {
					w = h1[ex2].first + h2[ey2].first;
					c = h1[ex2].second + h2[ey2].second;
					wc = make_pair(w, c);
					if (c < result[n - 1].second) {
						pq.push(hop(ex2, ey2, wc));
						visited[ex2][ey2]=true;
					}
				}
			}

		}
		return result;
	}

	void reduce(int subN) {
		deg.clear();
		neighbor.clear();
		adj_skyline.clear();
		vector<int> vectmp;
		vectmp.clear();
		for (int i = 0; i <= G.n; i++) {
			neighbor.push_back(vectmp);

		}
		adj_skyline.resize(G.n + 1);

		bool* _exist;
		_exist = (bool*)malloc(sizeof(bool) * (G.n + 1));
		for (int i = 1; i <= G.n; i++)
			_exist[i] = true;


		G.DD = (int*)malloc(sizeof(int) * (G.n + 1));
		G.DD2 = (int*)malloc(sizeof(int) * (G.n + 1));
		G._DD = (int*)malloc(sizeof(int) * (G.n + 1));
		G._DD2 = (int*)malloc(sizeof(int) * (G.n + 1));
		G.NUM = (int*)malloc(sizeof(int) * (G.n + 1));
		G.changed = (bool*)malloc(sizeof(bool) * (G.n + 1));
		for (int i = 0; i <= G.n; i++)
			G.NUM[i] = i;
		for (int i = 1; i <= G.n; i++) {
			int j = rand() % G.n + 1;
			int x = G.NUM[i];
			G.NUM[i] = G.NUM[j];
			G.NUM[j] = x;
		}
		for (int i = 1; i <= G.n; i++) {
			G.DD[i] = G.D[i];
			G.DD2[i] = G.D[i];
			G._DD[i] = G.D[i];
			G._DD2[i] = G.D[i];
			G.changed[i] = false;
			deg.insert(SelEle(i, &G));
		}
		int b_size;
		vector<int> sub_b;
		string bb = "";
		bb = to_string(static_cast<long long>(subN));
		ifstream read_bound("data/subgraphBound" + bb + ".txt");
		read_bound >> b_size;
		sub_b.resize(b_size);
		for (int i = 0; i < b_size; i++) {
			int bds;
			read_bound >> bds;
			int map_bds = VMap[bds][subN];
			sub_b[i] = map_bds;
			deg.erase(SelEle(map_bds, &G));
		}
		bool* exist;
		exist = (bool*)malloc(sizeof(bool) * (G.n + 1));
		for (int i = 1; i <= G.n; i++)
			exist[i] = true;
		ord.clear();
		int cnt = 0;
		int deg_size = deg.size();
		//cout << deg.size();
		while (!deg.empty())
		{
			vector<int> neigh;
			vector<vector<pair<int, int>>> adjsk;
			judge ju;
			cnt++;
			int bound;
			int x = (*deg.begin()).x;
			if (cnt == deg_size) {
				ord.push_back(x);
				deg.erase(SelEle(x, &G));
				exist[x] = false;
				neigh.clear();
				adjsk.clear();

				for (map<int, vector<pair<int, int>>>::iterator it = G.Edge[x].begin(); it != G.Edge[x].end(); it++) {
					int y = (*it).first;
					if (exist[y]) {
						neigh.push_back(y);
						adjsk.push_back((*it).second);
					}
				}
				int k = -1;
				for (int i = 0; i < neigh.size(); i++) {
					int y = neigh[i];
					G.deleteEdge(x, y);
					G._DD[y] = G.D[y];
					G.changed[y] = true;
				}

				for (int pu = 0; pu < neigh.size(); pu++) {
					for (int pv = 0; pv < neigh.size(); pv++) {
						if (pu != pv) {
							int u = neigh[pu];
							int v = neigh[pv];

							vector<pair<int, int>> u_sk = adjsk[pu];
							vector<pair<int, int>> v_sk = adjsk[pv];
							vector<pair<int, int>> uv_sk;
							vector<pair<int, int>> temp;
							vector<pair<int, int>> temp1;

							sort(u_sk.begin(), u_sk.end(), ju);
							sort(v_sk.begin(), v_sk.end(), ju);

							uv_sk = concatenation(u_sk, v_sk);

							if (G.isEdgeExist(u, v)) {
								mergeSort(G.Edge[u][v], uv_sk, temp);

								temp1 = screenSkylines(temp);
								G.Edge[u][v] = temp1;
								G.Edge[v][u] = temp1;
							}
							else {
								G.insertEdge(u, v, uv_sk);
								G._DD[u] = G.D[u];
								G._DD[v] = G.D[v];
								++G._DD2[u];
								++G._DD2[v];
								G.changed[u] = true;
								G.changed[v] = true;
							}
						}
					}
				}
				if (neigh.size() > tree_width)
					tree_width = neigh.size();
				neighbor[x] = neigh;
				adj_skyline[x] = adjsk;
				for (int it = 0; it < b_size; it++) {
					G.DD[sub_b[it]] = G._DD[sub_b[it]];
					G.DD2[sub_b[it]] = G._DD2[sub_b[it]];
					deg.insert(SelEle(sub_b[it], &G));
				}
				x = (*deg.begin()).x;
			}

			while (true) {
				if (G.changed[x]) {
					deg.erase(SelEle(x, &G));
					G.DD[x] = G._DD[x];
					G.DD2[x] = G._DD2[x];
					deg.insert(SelEle(x, &G));
					G.changed[x] = false;
					x = (*deg.begin()).x;
				}
				else break;
			}
			ord.push_back(x);
			deg.erase(deg.begin());
			exist[x] = false;
			if (cnt >= deg_size) {
				G.bound_ord[subN].push_back(G.r_VMap[subN][x]);
			}

			neigh.clear();
			adjsk.clear();

			for (map<int, vector<pair<int, int>>>::iterator it = G.Edge[x].begin(); it != G.Edge[x].end(); it++) {
				int y = (*it).first;
				if (exist[y]) {
					neigh.push_back(y);
					adjsk.push_back((*it).second);
				}
			}

			int k = -1;
			for (int i = 0; i < neigh.size(); i++) {
				int y = neigh[i];
				G.deleteEdge(x, y);
				G._DD[y] = G.D[y];
				G.changed[y] = true;
			}

			for (int pu = 0; pu < neigh.size(); pu++) {
				for (int pv = 0; pv < neigh.size(); pv++) {
					if (pu != pv) {
						int u = neigh[pu];
						int v = neigh[pv];

						vector<pair<int, int>> u_sk = adjsk[pu];
						vector<pair<int, int>> v_sk = adjsk[pv];
						vector<pair<int, int>> uv_sk;
						vector<pair<int, int>> temp;
						vector<pair<int, int>> temp1;

						sort(u_sk.begin(), u_sk.end(), ju);
						sort(v_sk.begin(), v_sk.end(), ju);

						uv_sk = concatenation(u_sk, v_sk);

						if (G.isEdgeExist(u, v)) {
							mergeSort(G.Edge[u][v], uv_sk, temp);

							temp1 = screenSkylines(temp);
							G.Edge[u][v] = temp1;
							G.Edge[v][u] = temp1;
						}
						else {
							G.insertEdge(u, v, uv_sk);
							G._DD[u] = G.D[u];
							G._DD[v] = G.D[v];
							++G._DD2[u];
							++G._DD2[v];
							G.changed[u] = true;
							G.changed[v] = true;
						}
					}
				}
			}
			if (neigh.size() > tree_width)
				tree_width = neigh.size();
			neighbor[x] = neigh;
			adj_skyline[x] = adjsk;
		}

		free(G.DD);
		free(G.DD2);
		free(exist);
	}

	void a_reduce_tree_build(int subN) {
		deg.clear();
		G.DD = (int*)malloc(sizeof(int) * (G.n + 1));
		G.DD2 = (int*)malloc(sizeof(int) * (G.n + 1));
		G._DD = (int*)malloc(sizeof(int) * (G.n + 1));
		G._DD2 = (int*)malloc(sizeof(int) * (G.n + 1));
		G.NUM = (int*)malloc(sizeof(int) * (G.n + 1));
		G.changed = (bool*)malloc(sizeof(bool) * (G.n + 1));
		for (int i = 0; i <= G.n; i++)
			G.NUM[i] = i;
		for (int i = 1; i <= G.n; i++) {
			int j = rand() % G.n + 1;
			int x = G.NUM[i];
			G.NUM[i] = G.NUM[j];
			G.NUM[j] = x;
		}
		for (int i = 1; i <= G.n; i++) {
			G.DD[i] = G._D[i];
			G.DD2[i] = G._D[i];
			G._DD[i] = G._D[i];
			G._DD2[i] = G._D[i];
			G.changed[i] = false;
			deg.insert(SelEle(i, &G));
		}
		int b_size;
		vector<int> sub_b;
		string bb = "";
		bb = to_string(static_cast<long long>(subN));
		ifstream read_bound("data/subgraphBound" + bb + ".txt");
		read_bound >> b_size;
		sub_b.resize(b_size);
		for (int i = 0; i < b_size; i++) {
			int bds;
			read_bound >> bds;
			int map_bds = VMap[bds][subN];
			sub_b[i] = map_bds;
			deg.erase(SelEle(map_bds, &G));
		}
		read_bound.close();
		bool* exist;
		exist = (bool*)malloc(sizeof(bool) * (G.n + 1));
		for (int i = 1; i <= G.n; i++)
			exist[i] = true;
		ord.clear();
		int cnt = 0;
		int deg_size = deg.size();
		while (!deg.empty())
		{
			vector<int> neigh;
			judge ju;
			cnt++;
			int bound;
			int x = (*deg.begin()).x;
			if (cnt == deg_size) {
				ord.push_back(x);
				deg.erase(SelEle(x, &G));
				exist[x] = false;
				neigh.clear();

				for (map<int, vector<pair<int, int>>>::iterator it = G._Edge[x].begin(); it != G._Edge[x].end(); it++) {
					int y = (*it).first;
					if (exist[y]) {
						neigh.push_back(y);
					}
				}

				int k = -1;
				for (int i = 0; i < neigh.size(); i++) {
					int y = neigh[i];
					G._deleteEdge(x, y);
					G._DD[y] = G._D[y];
					G.changed[y] = true;
				}

				for (int pu = 0; pu < neigh.size(); pu++) {
					for (int pv = 0; pv < neigh.size(); pv++) {
						if (pu != pv) {
							vector<pair<int, int>> uv_sk;	
							int u = neigh[pu];
							int v = neigh[pv];

							int tmp=max(G.lambda[x][u],G.lambda[x][v])+1;
							if(tmp>G.lambda[u][v]){
								G.lambda[u][v]=tmp;
								G.lambda[v][u]=tmp;
								if(tmp>max_lam){
									max_lam=tmp;
								}
							}			

							if (!G._isEdgeExist(u, v)) {					
								G._insertEdge(u, v, uv_sk);
								G._DD[u] = G._D[u];
								G._DD[v] = G._D[v];
								++G._DD2[u];
								++G._DD2[v];
								G.changed[u] = true;
								G.changed[v] = true;
							}
						}
					}
				}
				if (neigh.size() > tree_width)
					tree_width = neigh.size();
				for (int it = 0; it < b_size; it++) {
					G.DD[sub_b[it]] = G._DD[sub_b[it]];
					G.DD2[sub_b[it]] = G._DD2[sub_b[it]];
					deg.insert(SelEle(sub_b[it], &G));
				}
				x = (*deg.begin()).x;
			}


			while (true) {
				if (G.changed[x]) {
					deg.erase(SelEle(x, &G));
					G.DD[x] = G._DD[x];
					G.DD2[x] = G._DD2[x];
					deg.insert(SelEle(x, &G));
					G.changed[x] = false;
					x = (*deg.begin()).x;
				}
				else break;
			}
			ord.push_back(x);
			deg.erase(deg.begin());
			exist[x] = false;
			
			neigh.clear();

			for (map<int, vector<pair<int, int>>>::iterator it = G._Edge[x].begin(); it != G._Edge[x].end(); it++) {
				int y = (*it).first;
				if (exist[y]) {
					neigh.push_back(y);
				}
			}

			int k = -1;
			for (int i = 0; i < neigh.size(); i++) {
				int y = neigh[i];
				G._deleteEdge(x, y);
				G._DD[y] = G._D[y];
				G.changed[y] = true;
			}

			for (int pu = 0; pu < neigh.size(); pu++) {
				for (int pv = 0; pv < neigh.size(); pv++) {
					if (pu != pv) {
						int u = neigh[pu];
						int v = neigh[pv];
						vector<pair<int, int>> uv_sk;
						int tmp=max(G.lambda[x][u],G.lambda[x][v])+1;
						if(tmp>G.lambda[u][v]){
							G.lambda[u][v]=tmp;
							G.lambda[v][u]=tmp;
							if(tmp>max_lam){
								max_lam=tmp;
							}
						}			

						if (!G._isEdgeExist(u, v)) {					
							G._insertEdge(u, v, uv_sk);
							G._DD[u] = G._D[u];
							G._DD[v] = G._D[v];
							++G._DD2[u];
							++G._DD2[v];
							G.changed[u] = true;
							G.changed[v] = true;
						}
					}
				}
			}
			if (neigh.size() > tree_width)
				tree_width = neigh.size();
		}

		free(G.DD);
		free(G.DD2);
		free(exist);
	}

	void a_reduce(int subN) {
		deg.clear();
		neighbor.clear();
		adj_skyline.clear();
		vector<int> vectmp;
		vectmp.clear();
		for (int i = 0; i <= G.n; i++) {
			neighbor.push_back(vectmp);

		}
		adj_skyline.resize(G.n + 1);

		G.DD = (int*)malloc(sizeof(int) * (G.n + 1));
		G.DD2 = (int*)malloc(sizeof(int) * (G.n + 1));
		G._DD = (int*)malloc(sizeof(int) * (G.n + 1));
		G._DD2 = (int*)malloc(sizeof(int) * (G.n + 1));
		G.NUM = (int*)malloc(sizeof(int) * (G.n + 1));
		G.changed = (bool*)malloc(sizeof(bool) * (G.n + 1));
		for (int i = 0; i <= G.n; i++)
			G.NUM[i] = i;
		for (int i = 1; i <= G.n; i++) {
			int j = rand() % G.n + 1;
			int x = G.NUM[i];
			G.NUM[i] = G.NUM[j];
			G.NUM[j] = x;
		}
		for (int i = 1; i <= G.n; i++) {
			G.DD[i] = G.D[i];
			G.DD2[i] = G.D[i];
			G._DD[i] = G.D[i];
			G._DD2[i] = G.D[i];
			G.changed[i] = false;
			deg.insert(SelEle(i, &G));
		}
		int b_size;
		vector<int> sub_b;
		string bb = "";
		bb = to_string(static_cast<long long>(subN));
		ifstream read_bound("subgraphBound" + bb + ".txt");
		read_bound >> b_size;
		sub_b.resize(b_size);
		for (int i = 0; i < b_size; i++) {
			int bds;
			read_bound >> bds;
			int map_bds = VMap[bds][subN];
			sub_b[i] = map_bds;
			deg.erase(SelEle(map_bds, &G));
		}
		bool* exist;
		exist = (bool*)malloc(sizeof(bool) * (G.n + 1));
		for (int i = 1; i <= G.n; i++)
			exist[i] = true;
		ord.clear();
		int cnt = 0;
		int deg_size = deg.size();
		while (!deg.empty())
		{
			vector<int> neigh;
			vector<vector<pair<int, int>>> adjsk;
			judge ju;
			judgeCost juc;
			cnt++;
			int bound;
			int x = (*deg.begin()).x;
			if (cnt == deg_size) {
				ord.push_back(x);
				deg.erase(SelEle(x, &G));
				exist[x] = false;

				neigh.clear();
				adjsk.clear();

				for (map<int, vector<pair<int, int>>>::iterator it = G.Edge[x].begin(); it != G.Edge[x].end(); it++) {
					int y = (*it).first;
					if (exist[y]) {
						neigh.push_back(y);
						adjsk.push_back((*it).second);
					}
				}

				int k = -1;
				for (int i = 0; i < neigh.size(); i++) {
					int y = neigh[i];
					G.deleteEdge(x, y);
					G._DD[y] = G.D[y];
					G.changed[y] = true;
				}

				for (int pu = 0; pu < neigh.size(); pu++) {
					for (int pv = 0; pv < neigh.size(); pv++) {
						if (pu != pv) {
							int u = neigh[pu];
							int v = neigh[pv];

							vector<pair<int, int>> u_sk = adjsk[pu];
							vector<pair<int, int>> v_sk = adjsk[pv];
							vector<pair<int, int>> uv_sk;
							vector<pair<int, int>> temp;
							vector<pair<int, int>> temp1;

							double a_l;
							double kk=G.lambda[u][v];
							double rate=kk/double(max_lam);
							a_l=pow(a_g_s,rate);
							
							sort(u_sk.begin(), u_sk.end(), juc);
							sort(v_sk.begin(), v_sk.end(), juc);

							uv_sk = a_concatenation(u_sk, v_sk,a_l);

							if (G.isEdgeExist(u, v)) {
								a_mergeSort(G.Edge[u][v], uv_sk, temp);

								temp1 = a_screenSkylines(temp,a_l);
								G.Edge[u][v] = temp1;
								G.Edge[v][u] = temp1;
							}
							else {
								G.insertEdge(u, v, uv_sk);
								G._DD[u] = G.D[u];
								G._DD[v] = G.D[v];
								++G._DD2[u];
								++G._DD2[v];
								G.changed[u] = true;
								G.changed[v] = true;
							}
						}
					}
				}
				if (neigh.size() > tree_width)
					tree_width = neigh.size();
				neighbor[x] = neigh;
				adj_skyline[x] = adjsk;
				for(int i=0;i<adj_skyline[x].size();i++){
					sort(adj_skyline[x][i].begin(),adj_skyline[x][i].end(),ju);
				}
			
				for (int it = 0; it < b_size; it++) {
					G.DD[sub_b[it]] = G._DD[sub_b[it]];
					G.DD2[sub_b[it]] = G._DD2[sub_b[it]];
					deg.insert(SelEle(sub_b[it], &G));
				}
				x = (*deg.begin()).x;
			}


			while (true) {
				if (G.changed[x]) {
					deg.erase(SelEle(x, &G));
					G.DD[x] = G._DD[x];
					G.DD2[x] = G._DD2[x];
					deg.insert(SelEle(x, &G));
					G.changed[x] = false;
					x = (*deg.begin()).x;
				}
				else break;
			}
			ord.push_back(x);
			deg.erase(deg.begin());
			exist[x] = false;
			if (cnt >= deg_size) {
				G.bound_ord[subN].push_back(G.r_VMap[subN][x]);
			}

			neigh.clear();
			adjsk.clear();

			for (map<int, vector<pair<int, int>>>::iterator it = G.Edge[x].begin(); it != G.Edge[x].end(); it++) {
				int y = (*it).first;
				if (exist[y]) {
					neigh.push_back(y);
					adjsk.push_back((*it).second);
				}
			}

			int k = -1;
			for (int i = 0; i < neigh.size(); i++) {
				int y = neigh[i];
				G.deleteEdge(x, y);
				G._DD[y] = G.D[y];
				G.changed[y] = true;
			}

			for (int pu = 0; pu < neigh.size(); pu++) {
				for (int pv = 0; pv < neigh.size(); pv++) {
					if (pu != pv) {
						int u = neigh[pu];
						int v = neigh[pv];

						vector<pair<int, int>> u_sk = adjsk[pu];
						vector<pair<int, int>> v_sk = adjsk[pv];
						vector<pair<int, int>> uv_sk;
						vector<pair<int, int>> temp;
						vector<pair<int, int>> temp1;

						double a_l;
						double kk=G.lambda[u][v];
						double rate=kk/double(max_lam);
						a_l=pow(a_g_s,rate);

						sort(u_sk.begin(), u_sk.end(), juc);
						sort(v_sk.begin(), v_sk.end(), juc);

						uv_sk = a_concatenation(u_sk, v_sk,a_l);

						if (G.isEdgeExist(u, v)) {
							a_mergeSort(G.Edge[u][v], uv_sk, temp);

							temp1 = a_screenSkylines(temp,a_l);
							G.Edge[u][v] = temp1;
							G.Edge[v][u] = temp1;
						}
						else {
							G.insertEdge(u, v, uv_sk);
							G._DD[u] = G.D[u];
							G._DD[v] = G.D[v];
							++G._DD2[u];
							++G._DD2[v];
							G.changed[u] = true;
							G.changed[v] = true;
						}
					}
				}
			}
			if (neigh.size() > tree_width)
				tree_width = neigh.size();
			neighbor[x] = neigh;
			adj_skyline[x] = adjsk;
			for(int i=0;i<adj_skyline[x].size();i++){
				sort(adj_skyline[x][i].begin(),adj_skyline[x][i].end(),ju);
			}
		}

		free(G.DD);
		free(G.DD2);
		free(exist);
	}



	vector<Node> Tree;
	int root;
	int* belong, * rank;

	int match(int x, vector<int>& neigh) {
		int nearest = neigh[0];

		for (int i = 1; i < neigh.size(); i++)
			if (rank[neigh[i]] > rank[nearest])
				nearest = neigh[i];
		int p = belong[nearest];

		vector<int> a = Tree[p].vert;
		if (Tree[p].uniqueVertex >= 0) {
			a.push_back(Tree[p].uniqueVertex);
		}
		sort(a.begin(), a.end());
		int i, j;
		for (i = 0, j = 0; (i < neigh.size()) && (j < a.size()); ) {
			if (neigh[i] == a[j]) {
				i++; j++;
			}
			else if (neigh[i] < a[j])
				break;
			else j++;
		}
		if (i >= neigh.size()) {
			return p;
		}
		printf("no match!\n");
	}

	void makeTree() {
		belong = (int*)malloc(sizeof(int) * (H.n + 1));
		rank = (int*)malloc(sizeof(int) * (H.n + 1));
		int len = ord.size() - 1;
		Node rootn;
		Tree.clear();
		heightMax = 0;

		int x = ord[len];
		rootn.vert = neighbor[x];
		rootn.SK = adj_skyline[x];
		rootn.uniqueVertex = x;
		rootn.pa = -1;
		rootn.height = 1;
		rank[x] = 0;
		belong[x] = 0;
		Tree.push_back(rootn);
		len--;

		for (; len >= 0; len--) {
			int x = ord[len];
			int c = 0;
			Node nod;
			nod.vert = neighbor[x];
			nod.SK = adj_skyline[x];
			nod.uniqueVertex = x;
			int pa = match(x, neighbor[x]);
			Tree[pa].ch.push_back(Tree.size());
			nod.pa = pa;
			nod.height = Tree[pa].height + 1;
			if (nod.height > heightMax)
				heightMax = nod.height;
			rank[x] = Tree.size();
			belong[x] = Tree.size();
			Tree.push_back(nod);
		}

		root = 0;
	}

	int* toRMQ;
	vector<int> EulerSeq;
	vector<vector<int>> RMQIndex;
	void makeRMQDFS(int p, int height) {
		toRMQ[p] = EulerSeq.size();
		EulerSeq.push_back(p);
		for (int i = 0; i < Tree[p].ch.size(); i++) {
			makeRMQDFS(Tree[p].ch[i], height + 1);
			EulerSeq.push_back(p);
		}
	}
	void makeRMQ() {
		EulerSeq.clear();
		toRMQ = (int*)malloc(sizeof(int) * (G.n + 1));
		makeRMQDFS(root, 1);
		RMQIndex.clear();
		RMQIndex.push_back(EulerSeq);
		int m = EulerSeq.size();
		for (int i = 2, k = 1; i < m; i = i * 2, k++) {
			vector<int> tmp;
			tmp.clear();
			tmp.resize(EulerSeq.size());
			for (int j = 0; j < m - i; j++) {
				int x = RMQIndex[k - 1][j], y = RMQIndex[k - 1][j + i / 2];
				if (Tree[x].height < Tree[y].height)
					tmp[j] = x;
				else tmp[j] = y;
			}
			RMQIndex.push_back(tmp);
		}
	}
	int LCAQuery(int _p, int _q) {
		int p = toRMQ[_p], q = toRMQ[_q];
		if (p > q) {
			int x = p;
			p = q;
			q = x;
		}
		int len = q - p + 1;
		int i = 1, k = 0;
		while (i * 2 < len) {
			i *= 2;
			k++;
		}
		q = q - i + 1;
		if (Tree[RMQIndex[k][p]].height < Tree[RMQIndex[k][q]].height)
			return RMQIndex[k][p];
		else return RMQIndex[k][q];
	}

	vector<pair<int, int>> distanceQueryAncestorToPosterity(int p, int q) {
		vector<pair<int, int>> temp;
		temp.push_back(make_pair(0, 0));
		if (p == q) return temp;
		int x = belong[p], y = belong[q];
		//cout << "check:" << Tree[x].pos.size();
		//cout << ' ' << endl;
		return Tree[y].skyToAnc[Tree[x].pos[Tree[x].pos.size() - 1]];
	}

	void calculateIndexSizeDFS(int p, int pre, int tmp, long long& result) {
		for (int i = 0; i < Tree[p].ch.size(); i++) {
			calculateIndexSizeDFS(Tree[p].ch[i], pre + 1, (pre + 1) + tmp, result);
		}
		if (tmp + (pre + 1) > result) result = tmp + (pre + 1);
		//		result += pre;
	}
	long long calculateIndexSize() {
		long long res = Tree[root].vert.size();
		for (int i = 0; i < Tree[root].ch.size(); i++) {
			calculateIndexSizeDFS(Tree[root].ch[i], Tree[root].vert.size(), Tree[root].vert.size(), res);
		}
		return res;
	}

	void makeIndexDFS1(int p, vector<int>& list, int* toList) {
		vector<pair<int, int>> z;
		vector<pair<int, int>> _z;
		vector<pair<int, int>> zz;
		judge ju;
		zz.resize(list.size());
		Tree[p].pos.resize(Tree[p].vert.size() + 1);
		Tree[p].skyToAnc.resize(list.size());
		for (int i = 0; i < list.size(); i++) {
			Tree[p].skyToAnc[i].resize(list.size());
		}
		//	printf("step1");
		for (int i = 0; i < Tree[p].vert.size(); i++) {
			int j;
			for (j = 0; j < list.size(); j++)
				if (list[j] == Tree[p].vert[i])
					break;
			Tree[p].pos[i] = j;
		}
		Tree[p].pos[Tree[p].vert.size()] = list.size();
		//printf("Step 2");
		for (int i = 0; i < list.size(); i++) {
			for (int _j = 0; _j < Tree[p].skyToAnc[i].size(); _j++) {
				Tree[p].skyToAnc[i][_j].first = INF;
				Tree[p].skyToAnc[i][_j].second = INF;
			}
		}
		int x = Tree[p].uniqueVertex;

		vector<pair<int, int>> upperLeft;
		vector<pair<int, int>> downRight;
		vector<pair<int, int>> upperSet;
		vector<pair<int, int>> lowerSet;
		vector<vector<pair<int, int>>> hops;
		vector<int> ph;
		judgeCost juC;
		hops.resize(Tree[p].vert.size());
		int bbxMin, bbyMin, bbxMax, bbyMax;
		for (int i = 0; i < Tree[p].vert.size(); i++) {
			vector<pair<int, int>> assign;

			mergeSort(Tree[p].skyToAnc[toList[Tree[p].vert[i]]], Tree[p].SK[i], assign);
			Tree[p].skyToAnc[toList[Tree[p].vert[i]]] = screenSkylines(assign);

			int x = Tree[p].vert[i];
			int k;
			for (k = 0; k < list.size(); k++) {
				if (list[k] == x) break;
			}
			for (int j = 0; j < list.size(); j++) {
				int y = list[j];
				if (k < j) {
					z = distanceQueryAncestorToPosterity(x, y);
				}
				else
				{
					z = distanceQueryAncestorToPosterity(y, x);
				}

				int px1 = z[0].first + Tree[p].skyToAnc[toList[Tree[p].vert[i]]][0].first;
				int py1 = z[0].second + Tree[p].skyToAnc[toList[Tree[p].vert[i]]][0].second;
				upperLeft.push_back(make_pair(px1, py1));
				hops[i].push_back(make_pair(px1, py1));
				int px2 = z[z.size() - 1].first + Tree[p].skyToAnc[toList[Tree[p].vert[i]]][Tree[p].skyToAnc[toList[Tree[p].vert[i]]].size() - 1].first;
				int py2 = z[z.size() - 1].second + Tree[p].skyToAnc[toList[Tree[p].vert[i]]][Tree[p].skyToAnc[toList[Tree[p].vert[i]]].size() - 1].second;
				hops[i].push_back(make_pair(px2, py2));
				downRight.push_back(make_pair(px2, py2));
			}
		}
		upperSet = upperLeft;
		lowerSet = downRight;
		sort(upperLeft.begin(), upperLeft.end(), ju);
		sort(downRight.begin(), downRight.end(), juC);

		bbxMin = upperLeft[0].first;
		bbyMax = upperLeft[0].second;

		bbxMax = downRight[0].first;
		bbyMin = downRight[0].second;
		
		for (int i = 0; i < hops.size(); ++i) {
			if (bbxMin <= hops[i][1].first && bbxMax >= hops[i][0].first && bbyMin <= hops[i][0].second && bbyMax >= hops[i][1].second) {
				ph.push_back(i);
			}
		}
		upperSet.clear();
		lowerSet.clear();
		upperLeft.clear();
		downRight.clear();
		hops.clear();

		for (int i = 0; i < ph.size(); i++) {

			int x = Tree[p].vert[ph[i]];
			int k;
			for (k = 0; k < list.size(); k++) {
				if (list[k] == x) break;
			}
			for (int j = 0; j < list.size(); j++) {
				int y = list[j];
				if (k < j) {
					z = distanceQueryAncestorToPosterity(x, y);
				}
				else
				{
					z = distanceQueryAncestorToPosterity(y, x);
				}

				_z = concatenation(z, Tree[p].skyToAnc[toList[Tree[p].vert[ph[i]]]]);
				mergeSort(_z, Tree[p].skyToAnc[toList[y]], zz);
				Tree[p].skyToAnc[toList[y]] = screenSkylines(zz);
			}
		}
		ph.clear();
	
		toList[Tree[p].uniqueVertex] = list.size();
		list.push_back(Tree[p].uniqueVertex);
		for (int i = 0; i < Tree[p].ch.size(); i++) {
			makeIndexDFS(Tree[p].ch[i], list, toList);
		}
		list.pop_back();
		sort(Tree[p].pos.begin(), Tree[p].pos.end());
	}

	void makeIndexDFS(int p, vector<int>& list, int* toList) {
		vector<pair<int, int>> z;
		vector<pair<int, int>> _z;
		judge ju;
		Tree[p].pos.resize(Tree[p].vert.size() + 1);
		Tree[p].skyToAnc.resize(list.size());
		for (int i = 0; i < list.size(); i++) {
			Tree[p].skyToAnc[i].resize(list.size());
		}
		for (int i = 0; i < Tree[p].vert.size(); i++) {
			int j;
			for (j = 0; j < list.size(); j++)
				if (list[j] == Tree[p].vert[i])
					break;
			Tree[p].pos[i] = j;
		}
		Tree[p].pos[Tree[p].vert.size()] = list.size();

		for (int i = 0; i < list.size(); i++) {
			for (int _j = 0; _j < Tree[p].skyToAnc[i].size(); _j++) {
				Tree[p].skyToAnc[i][_j].first = INF;
				Tree[p].skyToAnc[i][_j].second = INF;
			}
		}
		int x = Tree[p].uniqueVertex;

		for (int i = 0; i < Tree[p].vert.size(); i++) {
			vector<pair<int, int>> assign;
			mergeSort(Tree[p].skyToAnc[toList[Tree[p].vert[i]]], Tree[p].SK[i], assign);
			Tree[p].skyToAnc[toList[Tree[p].vert[i]]] = screenSkylines(assign);
			
			int x = Tree[p].vert[i];
			int k;
			for (k = 0; k < list.size(); k++) {
				if (list[k] == x) break;
			}
			for (int j = 0; j < list.size(); j++) {
				int y = list[j];
				vector<pair<int, int>> zz;
				if (k < j) {
					z = distanceQueryAncestorToPosterity(x, y);
				}
				else
				{
					z = distanceQueryAncestorToPosterity(y, x);
				}
				
				_z = concatenation(z, Tree[p].skyToAnc[toList[Tree[p].vert[i]]]);
				mergeSort(_z, Tree[p].skyToAnc[toList[y]], zz);
				Tree[p].skyToAnc[toList[y]] = screenSkylines(zz);		
			}
			//}
		}
		toList[Tree[p].uniqueVertex] = list.size();
		list.push_back(Tree[p].uniqueVertex);
		for (int i = 0; i < Tree[p].ch.size(); i++) {
			makeIndexDFS(Tree[p].ch[i], list, toList);
		}
		list.pop_back();
	
		sort(Tree[p].pos.begin(), Tree[p].pos.end());
	}

	
double theta(vector<pair<int,int>>& h1, vector<pair<int,int>>& h2, double a_l){
	int size1=h1.size();
	int size2=h2.size();

	double wmax, wmin,a_i;
	a_i=1;
	wmax=h1[0].first+h2[0].first;
	wmin=h1[size1-1].first+h2[size2-1].first;
	double pre_size=log(wmax/wmin)/log(a_l);
	double cat_size=size1*size2;
	//cout<<"x: "<<double(wmax/wmin)<<"\t"<<"logx: "<<log(wmax/wmin)<<"\t"<<"base: "<<log(a_l)<<endl;
	//cout<<"cat: "<<cat_size<<'\t'<<"pre_size: "<<pre_size<<endl;
	if(cat_size>pre_size){
		/*if(cat_size>1){
		cout<<"x: "<<double(wmax/wmin)<<"\t"<<"logx: "<<log(wmax/wmin)<<"\t"<<"base: "<<log(a_l)<<endl;
	    cout<<"cat: "<<cat_size<<'\t'<<"pre_size: "<<pre_size<<endl;
		}*/
		a_i=a_l;
	}
	return a_i;
}

double theta_non(vector<pair<int,int>>& h1, vector<pair<int,int>>& h2, double a_l){
	double a_i;
	a_i=a_l;
	return a_i;
}


void a_makeIndexDFS(int p, vector<int>& list, int* toList) {
		vector<pair<int, int>> z;
		vector<pair<int, int>> _z;
		
		judge ju;
		judgeCost juc;
		double a_l;
		int depth;
		depth=Tree[p].height;
		int parent=Tree[p].pa;

		if(parent==root){
			Tree[root].a_rest=a_g_s;
		}
		
		Tree[p].pos.resize(Tree[p].vert.size() + 1);
		Tree[p].skyToAnc.resize(list.size());
		for (int i = 0; i < list.size(); i++) {
			Tree[p].skyToAnc[i].resize(list.size());
		}
		for (int i = 0; i < Tree[p].vert.size(); i++) {
			int j;
			for (j = 0; j < list.size(); j++)
				if (list[j] == Tree[p].vert[i])
					break;
			Tree[p].pos[i] = j;
		}
		Tree[p].pos[Tree[p].vert.size()] = list.size();
		for (int i = 0; i < list.size(); i++) {
			for (int _j = 0; _j < Tree[p].skyToAnc[i].size(); _j++) {
				Tree[p].skyToAnc[i][_j].first = INF;
				Tree[p].skyToAnc[i][_j].second = INF;
			}
		}
		int x = Tree[p].uniqueVertex;

		for (int i = 0; i < Tree[p].vert.size(); i++) {
			vector<pair<int, int>> assign;
			
			sort(Tree[p].SK[i].begin(),Tree[p].SK[i].end(),juc);
			a_mergeSort(Tree[p].skyToAnc[toList[Tree[p].vert[i]]], Tree[p].SK[i], assign);
			Tree[p].skyToAnc[toList[Tree[p].vert[i]]] = a_screenSkylines(assign,1);

			int x = Tree[p].vert[i];
			int k;
			for (k = 0; k < list.size(); k++) {
				if (list[k] == x) break;
			}
			for (int j = 0; j < list.size(); j++) {
				int y = list[j];
				vector<pair<int, int>> zz;

				int ay=belong[y];
				int h_ay=Tree[ay].height;
				double kk=abs(h_ay-depth);
				double nn=heightMax-h_ay;
				double rate=1/(nn-1);
				if(kk==1){
					a_l=pow(Tree[ay].a_rest,rate);
				}
				else{
					double a_l1=pow(Tree[ay].a_rest,rate);
					double rate1=(kk-1)/(heightMax-1);
					double a_l2=pow(a_g_s,rate1);
					a_l=a_l1*a_l2;
				}
				if (k < j) {
					z = distanceQueryAncestorToPosterity(x, y);
				}
				else
				{
					z = distanceQueryAncestorToPosterity(y, x);
				}
				sort(z.begin(),z.end(),juc);
				sort(Tree[p].skyToAnc[toList[Tree[p].vert[i]]].begin(),Tree[p].skyToAnc[toList[Tree[p].vert[i]]].end(),juc);

				double a_i;
				
				a_i=theta_non(z,Tree[p].skyToAnc[toList[Tree[p].vert[i]]],a_l);
				_z = a_concatenation(z,Tree[p].skyToAnc[toList[Tree[p].vert[i]]],a_i);
				a_mergeSort(_z,Tree[p].skyToAnc[toList[y]],zz);
			
				Tree[p].skyToAnc[toList[y]] = al_screenSkylines(zz,a_i,p,ay);
				if(alc[p][ay]<=a_i/rec[p][ay]){
					alc[p][ay]=a_i/rec[p][ay];	
				}					
			}
		}
		
		Tree[p].a_rest=a_g_s/alc[p][root];
		toList[Tree[p].uniqueVertex] = list.size();
		list.push_back(Tree[p].uniqueVertex);
		for (int i = 0; i < Tree[p].ch.size(); i++) {
			a_makeIndexDFS(Tree[p].ch[i], list, toList);
		}
		list.pop_back();
	
		sort(Tree[p].pos.begin(), Tree[p].pos.end());
	}

	void makeIndex() {
		makeRMQ();
		H.EdgeInitialize();
	}

	void readPositionIndex(int it, vector<vector<int>>& position) {
		string bb = "";
		bb = to_string(static_cast<long long>(it));
		//ifstream in("index\\subtree_position_index" + bb + ".txt");
		ifstream in("index/subtree_position_index" + bb + ".txt");
		string line;
		vector<int> temp;
		while (getline(in, line)) {
			stringstream ss(line);
			int x;
			while (ss >> x) {
				//cout << x << '\t';
				temp.push_back(x);
			}
			position.push_back(temp);
			temp.clear();
			//cout << endl;
		}
		in.close();
	}

	void readSkylineIndex(int it) {
		string bb = "";
		bb = to_string(static_cast<long long>(it));
		//ifstream in("index\\subtree_skyline_index" + bb + ".txt");
		ifstream in("index/subtree_skyline_index" + bb + ".txt");
		int n, ns, vert, sks;
		for (; in >> n >> ns >> vert >> sks;) {
			Tree[n].skyToAnc.resize(ns);
			for (int i = 0; i < sks; i++) {
				int x, y;
				in >> x >> y;
				Tree[n].skyToAnc[vert].push_back(make_pair(x, y));
			}
		}
		in.close();
	}

	void readSkylineIndex1(char* file) {
		ifstream in(file);
		string line;
		vector<pair<int, int>> temp1;
		vector<vector<pair<int, int>>> temp;

		while (getline(in, line)) {
			int num = 0;
			int num1;
			stringstream ss;
			ss.clear();
			ss.str("");
			ss << line;
			int i, j, k;
			int x, y;
			ss >> i;
			ss >> j;
			ss >> k;
			Tree[i].skyToAnc.resize(k);
			while (ss >> x >> y) {
				Tree[i].skyToAnc[j].push_back(make_pair(x, y));
			}
			temp.clear();
		}
		in.close();
	}


	void printIndex(int it) {
		ofstream outfile;
		string bb = "";
		bb = to_string(static_cast<long long>(it));
		outfile.open("index/subtree_position_index" + bb + ".txt");
		for (int i = 0; i < Tree.size(); i++) {
			for (int j = 0; j < Tree[i].pos.size(); j++) {
				outfile << Tree[i].pos[j] << '\t';
			}
			outfile << endl;
		}
		outfile.close();

		ofstream outfile1;
		outfile1.open("index/subtree_skyline_index" + bb + ".txt");
		for (int i = 0; i < Tree.size(); i++) {
			//outfile1 << i << '\t';
			for (int j = 0; j < Tree[i].skyToAnc.size(); j++) {
				outfile1 << i << '\t' << Tree[i].skyToAnc.size() << '\t' << j << '\t' << Tree[i].skyToAnc[j].size() << '\t';
				for (int k = 0; k < Tree[i].skyToAnc[j].size(); k++) {
					outfile1 << Tree[i].skyToAnc[j][k].first << '\t' << Tree[i].skyToAnc[j][k].second << '\t';
				}
				outfile1 << endl;
			}
		}
		outfile1.close();
	}

	void a_printIndex_reduce(int it) {
		ofstream outfile;
		string bb = "";
		bb = to_string(static_cast<long long>(it));
		outfile.open("index/reduce_subtree_position_index" + bb + ".txt");
		for (int i = 0; i < Tree.size(); i++) {
			for (int j = 0; j < Tree[i].pos.size(); j++) {
				outfile << Tree[i].pos[j] << '\t';
			}
			outfile << endl;
		}
		outfile.close();

		ofstream outfile1;
		outfile1.open("index/reduce_subtree_skyline_index" + bb + ".txt");
		for (int i = 0; i < Tree.size(); i++) {
			//outfile1 << i << '\t';
			for (int j = 0; j < Tree[i].skyToAnc.size(); j++) {
				outfile1 << i << '\t' << Tree[i].skyToAnc.size() << '\t' << j << '\t' << Tree[i].skyToAnc[j].size() << '\t';
				for (int k = 0; k < Tree[i].skyToAnc[j].size(); k++) {
					outfile1 << Tree[i].skyToAnc[j][k].first << '\t' << Tree[i].skyToAnc[j][k].second << '\t';
				}
				outfile1 << endl;
			}
		}
		outfile1.close();
	}

	void a_printIndex(int it) {
			ofstream outfile;
			string bb = "";
			bb = to_string(static_cast<long long>(it));
			outfile.open("index/a_subtree_position_index" + bb + ".txt");
			for (int i = 0; i < Tree.size(); i++) {
				for (int j = 0; j < Tree[i].pos.size(); j++) {
					outfile << Tree[i].pos[j] << '\t';
				}
				outfile << endl;
			}
			outfile.close();

			ofstream outfile1;
			outfile1.open("index/a_subtree_skyline_index" + bb + ".txt");
			for (int i = 0; i < Tree.size(); i++) {
				//outfile1 << i << '\t';
				for (int j = 0; j < Tree[i].skyToAnc.size(); j++) {
					outfile1 << i << '\t' << Tree[i].skyToAnc.size() << '\t' << j << '\t' << Tree[i].skyToAnc[j].size() << '\t';
					for (int k = 0; k < Tree[i].skyToAnc[j].size(); k++) {
						outfile1 << Tree[i].skyToAnc[j][k].first << '\t' << Tree[i].skyToAnc[j][k].second << '\t';
					}
					outfile1 << endl;
				}
			}
			outfile1.close();
	}

	void a_printIndex_theta(int it) {
			ofstream outfile;
			string bb = "";
			bb = to_string(static_cast<long long>(it));
			outfile.open("index/rec_subtree_position_index" + bb + ".txt");
			for (int i = 0; i < Tree.size(); i++) {
				for (int j = 0; j < Tree[i].pos.size(); j++) {
					outfile << Tree[i].pos[j] << '\t';
				}
				outfile << endl;
			}
			outfile.close();

			ofstream outfile1;
			outfile1.open("index/rec_subtree_skyline_index" + bb + ".txt");
			for (int i = 0; i < Tree.size(); i++) {
				//outfile1 << i << '\t';
				for (int j = 0; j < Tree[i].skyToAnc.size(); j++) {
					outfile1 << i << '\t' << Tree[i].skyToAnc.size() << '\t' << j << '\t' << Tree[i].skyToAnc[j].size() << '\t';
					for (int k = 0; k < Tree[i].skyToAnc[j].size(); k++) {
						outfile1 << Tree[i].skyToAnc[j][k].first << '\t' << Tree[i].skyToAnc[j][k].second << '\t';
					}
					outfile1 << endl;
				}
			}
			outfile1.close();
	}


	void a_resort() {
		judge ju;
		for (int i = 0; i < Tree.size(); i++) {
			for (int j = 0; j < Tree[i].skyToAnc.size(); j++) {
				sort(Tree[i].skyToAnc[j].begin(),Tree[i].skyToAnc[j].end(),ju);
			}
		}
	}

	void loadIndex(int a, int b) {
		//fstream f;
		//f.open(file);
		if (a == 0) {
			cout << "Re-building CSP-2Hop index" << endl;
			makeIndex2();
			reducePos();
			printIndex(b);
		}
		else if (a == 99) {
			cout << "Re-building Approximate-Recycle CSP-2Hop index" << endl;
			a_makeIndex2();
			reducePos();
			a_resort();
			a_printIndex_theta(b);
		}
		else if (a == 100) {
			cout << "Re-building Approximate CSP-2Hop index" << endl;
			a_makeIndex2();
			reducePos();
			a_resort();
			a_printIndex(b);
		}

		else if(a==1001){
			cout << "Re-building Approximate-Reduce CSP-2Hop index" << endl;
			makeIndex2();
			reducePos();
			a_printIndex_reduce(b);
		}
		
		else if (a == 10) {
			cout << "Re-building CSP-2Hop index(optimal)" << endl;
			makeIndex3();
			reducePos();
			printIndex(b);
		}
		else {
			cout << "Loading CSP-2Hop index" << endl;
			vector<vector<int>> ps;
			//char positionFile[255] = "position_index(BAY_simple).txt";
			//char skylineFile[255] = "skyline_index(BAY_simple).txt";
			readPositionIndex(b, ps);
			for (int i = 0; i < ps.size(); i++) {
				Tree[i].pos = ps[i];
			}
			readSkylineIndex(b);
		}
	}
	void makeIndex2() {
		vector<int> list;
		list.clear();
		int* toList;
		toList = (int*)malloc(sizeof(int) * (G.n + 1));
		Tree[root].pos.clear();

		toList[Tree[root].uniqueVertex] = 0;
		list.push_back(Tree[root].uniqueVertex);
		Tree[root].pos.push_back(0);
		for (int i = 0; i < Tree[root].ch.size(); ++i) {
			makeIndexDFS(Tree[root].ch[i], list, toList);
			break;
		}
		free(toList);

	}
	void a_makeIndex2() {
		vector<int> list;
		list.clear();
		int* toList;
		toList = (int*)malloc(sizeof(int) * (G.n + 1));
		Tree[root].pos.clear();

		toList[Tree[root].uniqueVertex] = 0;
		list.push_back(Tree[root].uniqueVertex);
		Tree[root].pos.push_back(0);
		vector<double> tmp_re;
		vector<double> tmp_alc;
	    tmp_re.assign(G.n+1,INF);
		tmp_alc.assign(G.n+1,1);
	    rec.assign(G.n+1,tmp_re);
	    alc.assign(G.n+1,tmp_alc);
		
		for (int i = 0; i < Tree[root].ch.size(); ++i) {
			a_makeIndexDFS(Tree[root].ch[i], list, toList);
			break;
		}
		free(toList);

	}


	void makeIndex3() {
		vector<int> list;
		list.clear();
		int* toList;
		toList = (int*)malloc(sizeof(int) * (G.n + 1));
		Tree[root].pos.clear();

		toList[Tree[root].uniqueVertex] = 0;
		list.push_back(Tree[root].uniqueVertex);
		Tree[root].pos.push_back(0);

		for (int i = 0; i < Tree[root].ch.size(); i++) {
			makeIndexDFS1(Tree[root].ch[i], list, toList);
			break;
		}
		free(toList);
	}

	void reducePosDFS(int p) {
		if (Tree[p].ch.size() == 2) {
			int t = Tree[p].ch[0];
			if (Tree[Tree[p].ch[0]].pos.size() > Tree[Tree[p].ch[1]].pos.size())
				t = Tree[p].ch[1];
			Tree[p].pos = Tree[t].pos;
			if (Tree[p].pos.size() != 0) {
				Tree[p].pos.erase(Tree[p].pos.begin() + (Tree[p].pos.size() - 1));
			}

		}
		for (int i = 0; i < Tree[p].ch.size(); i++)
			reducePosDFS(Tree[p].ch[i]);
	}
	void reducePos() {
		reducePosDFS(root);
	}
	void cntSize() {
		long long s_nonroot = 0;
		long long s_size = 0;

		long long s_dis = 0;
		for (int i = 0; i < Tree.size(); ++i) {
			s_nonroot += Tree[i].height - 1;
			s_size += Tree[i].vert.size();
			s_dis += Tree[i].height;
		}
		long long s_root = (long long)Tree[0].vert.size() * (long long)G.n;
		printf("tree width: %d\n", tree_width);
		printf("nonroot idx size = %0.3lfGB, avg node size=%0.3lf, avg label size=%0.3lf\n",
			s_nonroot * 4.0 / 1000000000.0, s_size * 1.0 / G.n, s_dis * 1.0 / G.n);
	}

	inline pair<int, int> concatenation_C(const vector<pair<int, int>>& h1, const vector<pair<int, int>>& h2, const int& C) {
		int size1 = h1.size();
		int size2 = h2.size();

		pair<int, int> wc;
		vector<pair<int, int>> result;
		if (size1 == 0 || size2 == 0) {
			result.push_back(make_pair(INF, INF));
			return result[0];
		}

		vector<vector<bool>> visited(size1, vector<bool>(size2, false));
		priority_queue<hop, vector<hop>, pk> pq;
		int w = h1[0].first + h2[0].first;
		int c = h1[0].second + h2[0].second;

		wc = make_pair(w, c);
		pq.push(hop(0, 0, wc));
		visited[0][0] = true;
		int count = 0;
		result.push_back(wc);
		while (!pq.empty()) {
			hop el = pq.top();
			pq.pop();

			int n = result.size();

			if (n > 0) {
				if (el.val.second <= result[n - 1].second) {
					result.push_back(el.val);
				}
			}
			if (el.val.second <= C) break;

			int ex1 = el.i + 1;
			int ey1 = el.j;
			int ex2 = el.i;
			int ey2 = el.j + 1;
			if (ex1 < size1) {
				if (!visited[ex1][ey1]) {
					w = h1[ex1].first + h2[ey1].first;
					c = h1[ex1].second + h2[ey1].second;
					wc = make_pair(w, c);
					if (c <= result[n - 1].second) {
						pq.push(hop(ex1, ey1, wc));
						visited[ex1][ey1]=true;
					}
				}
			}
			if (ey2 < size2) {
				if (!visited[ex2][ey2]) {
					w = h1[ex2].first + h2[ey2].first;
					c = h1[ex2].second + h2[ey2].second;
					wc = make_pair(w, c);
					if (c <= result[n - 1].second) {
						pq.push(hop(ex2, ey2, wc));
						visited[ex2][ey2]=true;
					}
				}
			}

		}
		if (result[result.size() - 1].second > C) {
			result.push_back(make_pair(INF, INF));
			return result[result.size() - 1];
		}

		return result[result.size() - 1];
	}


	vector<pair<int, int>> upper;
	vector<pair<int, int>> lower;
	inline vector<int> multiHop(vector<vector<pair<int, int>>>& x, vector<vector<pair<int, int>>>& y, vector<int>& p) {
		vector<int> res;
		vector<vector<pair<int, int>>> hops;
		vector<pair<int, int>> upperLeft;
		vector<pair<int, int>> downRight;
		int bbxMin, bbyMin, bbxMax, bbyMax;
		judge ju;
		judgeCost juC;
		upper.clear();
		lower.clear();

		hops.resize(p.size());
		if (p.size() == 1) {
			res.push_back(0);
			return res;
		}
		for (int i = 0; i < p.size(); i++) {
			int px1 = x[p[i]][0].first + y[p[i]][0].first;
			int py1 = x[p[i]][0].second + y[p[i]][0].second;
			upperLeft.push_back(make_pair(px1, py1));
			hops[i].push_back(make_pair(px1, py1));
			int px2 = x[p[i]][x[p[i]].size() - 1].first + y[p[i]][y[p[i]].size() - 1].first;
			int py2 = x[p[i]][x[p[i]].size() - 1].second + y[p[i]][y[p[i]].size() - 1].second;
			hops[i].push_back(make_pair(px2, py2));
			downRight.push_back(make_pair(px2, py2));
		}
		upper = upperLeft;
		lower = downRight;
		sort(upperLeft.begin(), upperLeft.end(), ju);
		sort(downRight.begin(), downRight.end(), juC);

		bbxMin = upperLeft[0].first;
		bbyMax = upperLeft[0].second;

		bbxMax = downRight[0].first;
		bbyMin = downRight[0].second;
		for (int i = 0; i < p.size(); i++) {
			if (bbxMin <= hops[i][1].first && bbxMax >= hops[i][0].first && bbyMin <= hops[i][0].second && bbyMax >= hops[i][1].second) {
				res.push_back(i);
			}
		}
		return res;
	}

	inline pair<int, int> multiHops_skylineQuery(int p, int q, int C) {
		pair<int, int> temp;
		pair<int, int> candi;
		vector<pair<int, int>> s;
		vector<vector<pair<int, int>>> sky_x, sky_y;
		temp = make_pair(0, 0);
		if (p == q) return temp;
		int x = belong[p], y = belong[q];
		int lca = LCAQuery(x, y);
		//cout << "lca:" << lca << endl;
		if (lca == x || lca == y) {
			//queryCnt++;
			if (lca == y) {
				int v = y;
				y = x;
				x = v;
				v = p;
				p = q;
				q = v;
			}
			if (lca == 0) {
				s = Tree[y].skyToAnc[0];
			}
			else {
				s = Tree[y].skyToAnc[Tree[x].pos[Tree[x].pos.size() - 1]];
			}
			
			if (s[s.size() - 1].second > C) {
				return temp;
			}
			else {
				for (int i = 0; i < s.size(); i++) {
					if (s[i].second <= C) {
						temp = s[i];
						break;
					}
				}
				return temp;
			}
		}
		else {
			int w = INF;
			int c = INF;

			pair<int, int> wc;
			sky_x = Tree[x].skyToAnc;
			sky_y = Tree[y].skyToAnc;
			vector<int> hp;
			vector<int> _p = Tree[lca].pos;
			int ps = Tree[lca].pos.size();
			hp = multiHop(sky_x, sky_y, _p);
			vector<pair<int, int>> h;
			for (int i = 0; i < hp.size(); i++) {
				int pos2 = hp[i];
				if (upper[pos2].second <= C) {
					if (upper[pos2].first < w) {
						w = upper[pos2].first;
						c = upper[pos2].second;
					}
				}
				else if (lower[pos2].second > C) {
					continue;
				}
				else {
					candi = concatenation_C(sky_x[_p[pos2]], sky_y[_p[pos2]], C);
					if (w > candi.first) {
						w = candi.first;
						c = candi.second;
					}
				}
			}
			wc = make_pair(w, c);
			return wc;
		}
	}

	inline vector<pair<int, int>> multiHops_skylines(int p, int q) {
		vector<pair<int, int>> temp;
		vector<pair<int, int>> candi;
		vector<pair<int, int>> temp1;
		vector<pair<int, int>> s;
		vector<vector<pair<int, int>>> sky_x, sky_y;
		s.clear();
		temp.push_back(make_pair(0, 0));
		temp1.push_back(make_pair(INF, INF));
		if (p == q) return temp;
		int x = belong[p], y = belong[q];
		int lca = LCAQuery(x, y);
		if (lca == x || lca == y) {
			if (lca == y) {
				int v = y;
				y = x;
				x = v;
				v = p;
				p = q;
				q = v;
			}
			if (lca == 0) {
				if (Tree[y].skyToAnc.size() == 0) {
					return temp1;
				}
				x = x + 1;
				s = Tree[y].skyToAnc[Tree[x].pos[Tree[x].pos.size() - 1]];
			}
			else {
				s = Tree[y].skyToAnc[Tree[x].pos[Tree[x].pos.size() - 1]];
			}
			return s;
		}
		else {
			int w = INF;
			int c = INF;
			pair<int, int> wc;
			sky_x = Tree[x].skyToAnc;
			sky_y = Tree[y].skyToAnc;
			if (lca == 0) {
				if (sky_x.size() == 0 && sky_y.size() != 0) {
					s = Tree[y].skyToAnc[0];
					return s;
				}
				if (sky_x.size() != 0 && sky_y.size() == 0) {
					s = Tree[x].skyToAnc[0];
					return s;
				}
			}

			vector<int> hp;
			vector<int> _p = Tree[lca].pos;
			int ps = Tree[lca].pos.size();
			hp = multiHop(sky_x, sky_y, _p);
			for (int i = 0; i < hp.size(); i++) {
				int pos2 = hp[i];
				candi = concatenation(sky_x[_p[pos2]], sky_y[_p[pos2]]);
				mergeSort(temp1, candi, s);
				temp1 = screenSkylines(s);
			}
			return temp1;
		}
	}

	inline vector<pair<int, int>> skylines_bound(int p, int q) {
		vector<pair<int, int>> temp, temp1;
		vector<pair<int, int>> candi;
		vector<pair<int, int>> s;
		vector<vector<pair<int, int>>> sky_x, sky_y;
		temp.push_back(make_pair(0, 0));
		temp1.push_back(make_pair(INF, INF));
		if (p == q) return temp;
		int x = belong[p], y = belong[q];
		int lca = LCAQuery(x, y);
		if (lca == x || lca == y) {
			if (lca == y) {
				int v = y;
				y = x;
				x = v;
				v = p;
				p = q;
				q = v;
			}
			int s_size = Tree[y].skyToAnc[Tree[x].pos[Tree[x].pos.size() - 1]].size();
			if (s_size <= 1) {
				s = Tree[y].skyToAnc[Tree[x].pos[Tree[x].pos.size() - 1]];
			}
			else
			{
				vector<pair<int, int>> tmp_s = Tree[y].skyToAnc[Tree[x].pos[Tree[x].pos.size() - 1]];
				s.push_back(tmp_s[0]);
				s.push_back(tmp_s[s_size - 1]);
			}
			return s;
		}
		else {
			int w = INF;
			int c = INF;

			pair<int, int> wc;
			sky_x = Tree[x].skyToAnc;
			sky_y = Tree[y].skyToAnc;
			vector<int> p = Tree[lca].pos;
			int ps = Tree[lca].pos.size();
			vector<pair<int, int>> h;
			s = rectangle(sky_x, sky_y, p);
		}
		return s;
	}

	inline vector<pair<int, int>> rectangle(vector<vector<pair<int, int>>>& x, vector<vector<pair<int, int>>>& y, vector<int>& p) {
		vector<pair<int, int>> res;
		vector<vector<pair<int, int>>> hops;
		vector<pair<int, int>> upperLeft;
		vector<pair<int, int>> downRight;
		int bbxMin, bbyMin, bbxMax, bbyMax;
		judge ju;
		judgeCost juC;
		hops.resize(p.size());

		for (int i = 0; i < p.size(); ++i) {
			int px1 = x[p[i]][0].first + y[p[i]][0].first;
			int py1 = x[p[i]][0].second + y[p[i]][0].second;
			upperLeft.push_back(make_pair(px1, py1));
			hops[i].push_back(make_pair(px1, py1));
			int px2 = x[p[i]][x[p[i]].size() - 1].first + y[p[i]][y[p[i]].size() - 1].first;
			int py2 = x[p[i]][x[p[i]].size() - 1].second + y[p[i]][y[p[i]].size() - 1].second;
			hops[i].push_back(make_pair(px2, py2));
			downRight.push_back(make_pair(px2, py2));
		}

		sort(upperLeft.begin(), upperLeft.end(), ju);
		sort(downRight.begin(), downRight.end(), juC);

		bbxMin = upperLeft[0].first;
		bbyMax = upperLeft[0].second;
		res.push_back(make_pair(bbxMin, bbyMax));
		bbxMax = downRight[0].first;
		bbyMin = downRight[0].second;
		res.push_back(make_pair(bbxMax, bbyMin));
		return res;
	}

	inline pair<int, int> simple_skylineQuery(int p, int q, int C) {
		pair<int, int> temp;
		pair<int, int> candi;
		vector<pair<int, int>> s;
		vector<vector<pair<int, int>>> sky_x, sky_y;
		temp = make_pair(0, 0);
		if (p == q) return temp;
		int x = belong[p], y = belong[q];
		int lca = LCAQuery(x, y);
		if (lca == x || lca == y) {
			if (lca == y) {
				int v = y;
				y = x;
				x = v;
				v = p;
				p = q;
				q = v;
			}
			s = Tree[y].skyToAnc[Tree[x].pos[Tree[x].pos.size() - 1]];
			if (s[s.size() - 1].second > C) {
				return temp;
			}
			else {
				for (int i = 0; i < s.size(); i++) {
					if (s[i].second <= C) {
						temp = s[i];
						break;
					}
				}
				return temp;
			}
		}
		else {
			int w = INF;
			int c = INF;

			pair<int, int> wc;
			sky_x = Tree[x].skyToAnc;
			sky_y = Tree[y].skyToAnc;
			vector<int> p = Tree[lca].pos;
			int ps = Tree[lca].pos.size();
			vector<pair<int, int>> h;
			for (int i = 0; i < ps; i++) {
				candi = concatenation_C(sky_x[p[i]], sky_y[p[i]], C);
				if (w > candi.first) {
					w = candi.first;
					c = candi.second;
				}

			}
			wc = make_pair(w, c);
			return wc;
		}
	}

	struct judge2 {
		bool operator()(const pair<int, int> a, const pair<int, int> b) {
			return a.first <= b.first;
		};
	};

	inline pair<int, int> base_skylineQuery(int p, int q, int C) {
		pair<int, int> temp;
		pair<int, int> candi;
		vector<pair<int, int>> s;
		vector<vector<pair<int, int>>> sky_x, sky_y;
		judge ju2;
		temp = make_pair(0, 0);
		if (p == q) return temp;
		int x = belong[p], y = belong[q];
		int lca = LCAQuery(x, y);
		if (lca == x || lca == y) {
			if (lca == y) {
				int v = y;
				y = x;
				x = v;
				v = p;
				p = q;
				q = v;
			}
			s = Tree[y].skyToAnc[Tree[x].pos[Tree[x].pos.size() - 1]];
			if (s[s.size() - 1].second > C) {
				return temp;
			}
			else {
				for (int i = 0; i < s.size(); i++) {
					if (s[i].second <= C) {
						temp = s[i];
						break;
					}
				}
				return temp;
			}
		}
		else {
			int w = INF;
			int c = INF;

			pair<int, int> wc;
			sky_x = Tree[x].skyToAnc;
			sky_y = Tree[y].skyToAnc;
			vector<int> p = Tree[lca].pos;
			int ps = Tree[lca].pos.size();
			vector<pair<int, int>> h;
			vector<pair<int, int>> candiSet;
			candiSet.clear();
			for (int i = 0; i < ps; i++) {
				for (int x = 0; x < sky_x[p[i]].size(); x++) {
					for (int y = 0; y < sky_y[p[i]].size(); y++) {
						int w_d = sky_x[p[i]][x].first + sky_y[p[i]][y].first;
						int w_c = sky_x[p[i]][x].second + sky_y[p[i]][y].second;
						if (w_c <= C) {
							candiSet.push_back(make_pair(w_d, w_c));
						}
					}
				}
				sort(candiSet.begin(), candiSet.end(), ju2);

				if (candiSet.size() > 0) {
					candi = candiSet[0];
					if (w > candi.first) {
						w = candi.first;
						c = candi.second;
					}
				}
				candiSet.clear();
			}
			wc = make_pair(w, c);
			return wc;
		}
	}

	inline vector<pair<int, int>> skylines(int p, int q) {
		vector<pair<int, int>> temp, temp1;
		vector<pair<int, int>> candi;
		vector<pair<int, int>> s;
		vector<vector<pair<int, int>>> sky_x, sky_y;
		temp.push_back(make_pair(0, 0));
		temp1.push_back(make_pair(INF, INF));
		if (p == q) return temp;
		int x = belong[p], y = belong[q];
		int lca = LCAQuery(x, y);
		if (lca == x || lca == y) {
			if (lca == y) {
				int v = y;
				y = x;
				x = v;
				v = p;
				p = q;
				q = v;
			}
			s = Tree[y].skyToAnc[Tree[x].pos[Tree[x].pos.size() - 1]];
			return s;
		}
		else {
			int w = INF;
			int c = INF;

			pair<int, int> wc;
			sky_x = Tree[x].skyToAnc;
			sky_y = Tree[y].skyToAnc;
			vector<int> p = Tree[lca].pos;
			int ps = Tree[lca].pos.size();
			vector<pair<int, int>> h;
			for (int i = 0; i < ps; i++) {
				candi = concatenation(sky_x[p[i]], sky_y[p[i]]);
				mergeSort(temp1, candi, s);
				temp1 = screenSkylines(s);
			}
			return temp1;
		}
	}
};

struct Bound_Tree_Decomposition {
	Graph G, H;
	vector<Sub_Tree_Decomposition> TDs;
	set<SelEle> deg;
	int maxSize;
	Bound_Tree_Decomposition() {}
	vector<vector<int>> neighbor;
	vector<vector<vector<pair<int, int>>>> adj_skyline;
	vector<int> ord;
	int heightMax;
	int max_lam=0;
	vector<vector<double>> rec,alc;
	double a_g_d=pow(a_g,0.5);

	struct judge {
		bool operator()(const pair<int, int> a, const pair<int, int> b) {
			return a.first < b.first;
		};
	};

	struct judgeCost {
		bool operator()(const pair<int, int> a, const pair<int, int> b) {
			return a.second < b.second;
		};
	};

	struct hop {
		int i;
		int j;
		pair<int, int> val;
		hop(int x, int y, pair<int, int> v) : i(x), j(y), val(v) {};
	};

	struct pk {
		bool operator()(const hop& a, const hop& b) const {
			return a.val.first < b.val.first;
		}
	};

	struct pk_c {
		bool operator()(const hop& a, const hop& b) const {
			return a.val.second < b.val.second;
		}
	};

	void mergeSort(const vector<pair<int, int>>& s1, const vector<pair<int, int>>& s2, vector<pair<int, int>>& result) {
		int size1 = s1.size();
		int size2 = s2.size();
		result.resize(size1 + size2);
		int tmp_index, tmp_1_index, tmp_2_index;
		tmp_1_index = 0;
		tmp_2_index = 0;
		for (tmp_index = 0; (tmp_1_index < size1) && (tmp_2_index < size2); tmp_index++) {
			if (s1[tmp_1_index].first <= s2[tmp_2_index].first) {

				result[tmp_index] = s1[tmp_1_index];
				tmp_1_index++;
			}
			else {
				result[tmp_index] = s2[tmp_2_index];
				tmp_2_index++;
			}
		}
		if (tmp_1_index == size1) {
			while (tmp_2_index < size2) {
				result[tmp_index] = s2[tmp_2_index];
				tmp_index++;
				tmp_2_index++;
			}
		}
		else if (tmp_2_index == size2) {
			while (tmp_1_index < size1) {
				result[tmp_index] = s1[tmp_1_index];
				tmp_index++;
				tmp_1_index++;
			}
		}
	}

	vector<pair<int, int>> screenSkylines(vector<pair<int, int>>& input) {
		vector<pair<int, int>> result;
		result.push_back(input[0]);
		int size = input.size();

		if (size <= 1) {
			return result;
		}
		else {
			for (int i = 1; i < size; i++) {
				if (result[result.size() - 1].first == input[i].first && result[result.size() - 1].second > input[i].second) {
					result[result.size() - 1].second = input[i].second;
				}
				if (result[result.size() - 1].first != input[i].first && result[result.size() - 1].second > input[i].second) {
					result.push_back(input[i]);
				}
				else {
					continue;
				}
			}
			return result;
		}
	}

	vector<pair<int, int>> concatenation(const vector<pair<int, int>>& h1, const vector<pair<int, int>>& h2) {
		int size1 = h1.size();
		int size2 = h2.size();
		pair<int, int> wc;
		vector<pair<int, int>> result;
		if (size1 == 0 || size2 == 0) {
			result.push_back(make_pair(INF, INF));
			return result;
		}
		vector<vector<bool>> visited(size1, vector<bool>(size2, false));
		priority_queue<hop, vector<hop>, pk> pq;
		int w = h1[0].first + h2[0].first;
		int c = h1[0].second + h2[0].second;
		wc = make_pair(w, c);
		pq.push(hop(0, 0, wc));
		visited[0][0] = true;
		result.push_back(wc);
		while (!pq.empty()) {
			hop el = pq.top();
			pq.pop();

			int ex1 = el.i + 1;
			int ey1 = el.j;
			int ex2 = el.i;
			int ey2 = el.j + 1;
			if (ex1 < size1) {
				if (!visited[ex1][ey1]) {
					w = h1[ex1].first + h2[ey1].first;
					c = h1[ex1].second + h2[ey1].second;
					wc = make_pair(w, c);
					if (c < result[result.size() - 1].second) {
						pq.push(hop(ex1, ey1, wc));
						visited[ex1][ey1]=true;
					}
				}
			}
			if (ey2 < size2) {
				if (!visited[ex2][ey2]) {
					w = h1[ex2].first + h2[ey2].first;
					c = h1[ex2].second + h2[ey2].second;
					wc = make_pair(w, c);
					if (c < result[result.size() - 1].second) {
						pq.push(hop(ex2, ey2, wc));
						visited[ex2][ey2]=true;
					}
				}
			}
		}
		return result;
	}

	void reduce() {
		neighbor.clear();
		adj_skyline.clear();
		vector<int> vectmp;
		vectmp.clear();
		for (int i = 0; i <= G.n; i++) {
			neighbor.push_back(vectmp);

		}
		adj_skyline.resize(G.n + 1);
		bool* _exist;
		_exist = (bool*)malloc(sizeof(bool) * (G.n + 1));
		for (int i = 1; i <= G.n; i++)
			_exist[i] = true;

		vector<vector<int>> bo;
		vector<pair<int, int>> bo_size;
		vector<int> sub_b;
		bo.resize(fn);
		for (int i = 0; i < fn; i++) {
			bo[i] = TDs[i].G.bound_ord[i];
			bo_size.push_back(make_pair(bo[i].size(), i));
		}
		judge2 ju;
		sort(bo_size.begin(), bo_size.end(), ju);
		for (int i = 0; i < bo.size(); i++) {
			int tmp = bo_size[i].second;
			for (int j = 0; j < bo[tmp].size(); j++) {
				int map_bds = VMap[bo[tmp][j]][fn];
				if (map_bds != 0) {
					sub_b.push_back(map_bds);
				}
			}
		}
		bool* exist;
		exist = (bool*)malloc(sizeof(bool) * (G.n + 1));
		for (int i = 1; i <= G.n; i++)
			exist[i] = true;
		ord.clear();
		int cnt = 0;
		while (!sub_b.empty())
		{
			vector<int> neigh;
			vector<vector<pair<int, int>>> adjsk;
			judge ju;
			cnt++;
			int bound;
			int x = sub_b[0];
			ord.push_back(x);
			sub_b.erase(sub_b.begin());
			exist[x] = false;
			neigh.clear();
			adjsk.clear();

			for (map<int, vector<pair<int, int>>>::iterator it = G.Edge[x].begin(); it != G.Edge[x].end(); it++) {
				int y = (*it).first;
				if (exist[y]) {
					neigh.push_back(y);
					adjsk.push_back((*it).second);
				}
			}
			int k = -1;
			for (int i = 0; i < neigh.size(); i++) {
				int y = neigh[i];
				G.deleteEdge(x, y);
			}

			for (int pu = 0; pu < neigh.size(); pu++) {
				for (int pv = 0; pv < neigh.size(); pv++) {
					if (pu != pv) {
						int u = neigh[pu];
						int v = neigh[pv];

						vector<pair<int, int>> u_sk = adjsk[pu];
						vector<pair<int, int>> v_sk = adjsk[pv];
						vector<pair<int, int>> uv_sk;
						vector<pair<int, int>> temp;
						vector<pair<int, int>> temp1;
						if (G.isEdgeExist(u, v)) {
							continue;
						}
						else {
							sort(u_sk.begin(), u_sk.end(), ju);
							sort(v_sk.begin(), v_sk.end(), ju);
							uv_sk = concatenation(u_sk, v_sk);
							G.insertEdge(u, v, uv_sk);
						}
					}
				}
			}
			if (neigh.size() > tree_width)
				tree_width = neigh.size();
			neighbor[x] = neigh;
			adj_skyline[x] = adjsk;
		}
	}

	void a_reduce() {
		neighbor.clear();
		adj_skyline.clear();
		vector<int> vectmp;
		vectmp.clear();
		for (int i = 0; i <= G.n; i++) {
			neighbor.push_back(vectmp);

		}
		adj_skyline.resize(G.n + 1);
		vector<vector<int>> bo;
		vector<pair<int, int>> bo_size;
		vector<int> sub_b;

		bo.resize(fn);
		for (int i = 0; i < fn; i++) {
			bo[i] = TDs[i].G.bound_ord[i];
			bo_size.push_back(make_pair(bo[i].size(), i));
		}
		judge2 ju;
		judgeCost juc;
		sort(bo_size.begin(), bo_size.end(), ju);
		for (int i = 0; i < bo.size(); i++) {
			int tmp = bo_size[i].second;
			for (int j = 0; j < bo[tmp].size(); j++) {
				int map_bds = VMap[bo[tmp][j]][fn];
				if (map_bds != 0) {
					sub_b.push_back(map_bds);
				}
			}
		}
		bool* exist;
		exist = (bool*)malloc(sizeof(bool) * (G.n + 1));
		for (int i = 1; i <= G.n; i++)
			exist[i] = true;
		ord.clear();
		int cnt = 0;
		while (!sub_b.empty())
		{
			vector<int> neigh;
			vector<vector<pair<int, int>>> adjsk;
			judge ju;
			cnt++;
			int bound;
			int x = sub_b[0];
			ord.push_back(x);
			sub_b.erase(sub_b.begin());
			exist[x] = false;
			neigh.clear();
			adjsk.clear();

			for (map<int, vector<pair<int, int>>>::iterator it = G.Edge[x].begin(); it != G.Edge[x].end(); it++) {
				int y = (*it).first;
				if (exist[y]) {
					neigh.push_back(y);
					adjsk.push_back((*it).second);
				}
			}

			int k = -1;
			for (int i = 0; i < neigh.size(); i++) {
				int y = neigh[i];
				G.deleteEdge(x, y);
			}

			for (int pu = 0; pu < neigh.size(); pu++) {
				for (int pv = 0; pv < neigh.size(); pv++) {
					if (pu != pv) {
						int u = neigh[pu];
						int v = neigh[pv];

						vector<pair<int, int>> u_sk = adjsk[pu];
						vector<pair<int, int>> v_sk = adjsk[pv];
						vector<pair<int, int>> uv_sk;
						vector<pair<int, int>> temp;
						vector<pair<int, int>> temp1;

						if (G.isEdgeExist(u, v)) {
							continue;
						}
						else {
							double a_l;
							double kk=G.lambda[u][v];
							double rate=kk/double(max_lam);
							a_l=pow(a_g_d,rate);
							
							sort(u_sk.begin(), u_sk.end(), juc);
							sort(v_sk.begin(), v_sk.end(), juc);

							uv_sk = a_concatenation(u_sk, v_sk,a_l);
							G.insertEdge(u, v, uv_sk);
						}
					}
				}
			}
			if (neigh.size() > tree_width)
				tree_width = neigh.size();
			neighbor[x] = neigh;
			adj_skyline[x] = adjsk;
			for(int i=0;i<adj_skyline[x].size();i++){
				sort(adj_skyline[x][i].begin(),adj_skyline[x][i].end(),ju);
			}

		}
	}

	void a_reduce_tree_bulit() {
		//deg.clear();
		neighbor.clear();

		vector<int> vectmp;
		vectmp.clear();
		for (int i = 0; i <= G.n; i++) {
			neighbor.push_back(vectmp);

		}

		vector<vector<int>> bo;
		vector<pair<int, int>> bo_size;
		vector<int> sub_b;
		//sub_b.resize(2542);
		bo.resize(fn);
		for (int i = 0; i < fn; i++) {
			bo[i] = TDs[i].G.bound_ord[i];
			bo_size.push_back(make_pair(bo[i].size(), i));
		}
		judge2 ju;
		sort(bo_size.begin(), bo_size.end(), ju);
		for (int i = 0; i < bo.size(); i++) {
			int tmp = bo_size[i].second;
			for (int j = 0; j < bo[tmp].size(); j++) {
				int map_bds = VMap[bo[tmp][j]][fn];
				if (map_bds != 0) {
					sub_b.push_back(map_bds);
				}
			}
		}
		bool* exist;
		exist = (bool*)malloc(sizeof(bool) * (G.n + 1));
		for (int i = 1; i <= G.n; i++)
			exist[i] = true;
		ord.clear();
		int cnt = 0;
		while (!sub_b.empty())
		{
			vector<int> neigh;
			judge ju;
			cnt++;
			int bound;
			int x = sub_b[0];
			ord.push_back(x);
			sub_b.erase(sub_b.begin());
			exist[x] = false;
			neigh.clear();
			for (map<int, vector<pair<int, int>>>::iterator it = G._Edge[x].begin(); it != G._Edge[x].end(); it++) {
				int y = (*it).first;
				if (exist[y]) {
					neigh.push_back(y);
				}
			}

			int k = -1;
			for (int i = 0; i < neigh.size(); i++) {
				int y = neigh[i];
				G._deleteEdge(x, y);
			}

			for (int pu = 0; pu < neigh.size(); pu++) {
				for (int pv = 0; pv < neigh.size(); pv++) {
					if (pu != pv) {
						int u = neigh[pu];
						int v = neigh[pv];
						vector<pair<int, int>> uv_sk;
						if (G._isEdgeExist(u, v)) {
							continue;
						}
						else {
							int tmp=max(G.lambda[x][u],G.lambda[x][v])+1;
							if(tmp>G.lambda[u][v]){
								G.lambda[u][v]=tmp;
								G.lambda[v][u]=tmp;
								if(tmp>max_lam){
									max_lam=tmp;
								}
							}			
							G._insertEdge(u, v, uv_sk);
						}
					}
				}
			}
			if (neigh.size() > tree_width)
				tree_width = neigh.size();

		}
		free(exist);
	}

	vector<Node> Tree;
	int root;
	int* belong, * rank;
	int match(int x, vector<int>& neigh) {
		int nearest = neigh[0];
		for (int i = 1; i < neigh.size(); i++)
			if (rank[neigh[i]] > rank[nearest])
				nearest = neigh[i];
		int p = belong[nearest];
		vector<int> a = Tree[p].vert;
		if (Tree[p].uniqueVertex >= 0) {
			a.push_back(Tree[p].uniqueVertex);
		}
		sort(a.begin(), a.end());
		int i, j;
		for (i = 0, j = 0; (i < neigh.size()) && (j < a.size()); ) {
			if (neigh[i] == a[j]) {
				i++; j++;
			}
			else if (neigh[i] < a[j])
				break;
			else j++;
		}
		if (i >= neigh.size()) {
			return p;
		}
		printf("no match!\n");
	}

	void makeTree() {
		belong = (int*)malloc(sizeof(int) * (H.n + 1));
		rank = (int*)malloc(sizeof(int) * (H.n + 1));
		int len = ord.size() - 1;
		Node rootn;
		Tree.clear();
		heightMax = 0;

		int x = ord[len];
		rootn.vert = neighbor[x];
		rootn.SK = adj_skyline[x];
		rootn.uniqueVertex = x;
		rootn.pa = -1;
		rootn.height = 1;
		rank[x] = 0;
		belong[x] = 0;
		Tree.push_back(rootn);
		len--;

		for (; len >= 0; len--) {
			int x = ord[len];
			int c = 0;
			Node nod;
			nod.vert = neighbor[x];
			nod.SK = adj_skyline[x];
			nod.uniqueVertex = x;
			int pa = match(x, neighbor[x]);
			Tree[pa].ch.push_back(Tree.size());
			nod.pa = pa;
			nod.height = Tree[pa].height + 1;
			if (nod.height > heightMax)
				heightMax = nod.height;
			rank[x] = Tree.size();
			belong[x] = Tree.size();
			Tree.push_back(nod);
		}

		root = 0;
	}

	int* toRMQ;
	vector<int> EulerSeq;
	vector<vector<int>> RMQIndex;
	void makeRMQDFS(int p, int height) {
		toRMQ[p] = EulerSeq.size();
		EulerSeq.push_back(p);
		for (int i = 0; i < Tree[p].ch.size(); i++) {
			makeRMQDFS(Tree[p].ch[i], height + 1);
			EulerSeq.push_back(p);
		}
	}
	void makeRMQ() {
		EulerSeq.clear();
		toRMQ = (int*)malloc(sizeof(int) * (G.n + 1));
		makeRMQDFS(root, 1);
		RMQIndex.clear();
		RMQIndex.push_back(EulerSeq);
		int m = EulerSeq.size();
		for (int i = 2, k = 1; i < m; i = i * 2, k++) {
			vector<int> tmp;
			tmp.clear();
			tmp.resize(EulerSeq.size());
			for (int j = 0; j < m - i; j++) {
				int x = RMQIndex[k - 1][j], y = RMQIndex[k - 1][j + i / 2];
				if (Tree[x].height < Tree[y].height)
					tmp[j] = x;
				else tmp[j] = y;
			}
			RMQIndex.push_back(tmp);
		}
	}
	int LCAQuery(int _p, int _q) {
		int p = toRMQ[_p], q = toRMQ[_q];
		if (p > q) {
			int x = p;
			p = q;
			q = x;
		}
		int len = q - p + 1;
		int i = 1, k = 0;
		while (i * 2 < len) {
			i *= 2;
			k++;
		}
		q = q - i + 1;
		if (Tree[RMQIndex[k][p]].height < Tree[RMQIndex[k][q]].height)
			return RMQIndex[k][p];
		else return RMQIndex[k][q];
	}

	vector<pair<int, int>> distanceQueryAncestorToPosterity(int p, int q) {
		vector<pair<int, int>> temp;
		temp.push_back(make_pair(0, 0));
		if (p == q) return temp;
		int x = belong[p], y = belong[q];
		return Tree[y].skyToAnc[Tree[x].pos[Tree[x].pos.size() - 1]];
	}

	void calculateIndexSizeDFS(int p, int pre, int tmp, long long& result) {
		for (int i = 0; i < Tree[p].ch.size(); i++) {
			calculateIndexSizeDFS(Tree[p].ch[i], pre + 1, (pre + 1) + tmp, result);
		}
		if (tmp + (pre + 1) > result) result = tmp + (pre + 1);
	}
	long long calculateIndexSize() {
		long long res = Tree[root].vert.size();
		for (int i = 0; i < Tree[root].ch.size(); i++) {
			calculateIndexSizeDFS(Tree[root].ch[i], Tree[root].vert.size(), Tree[root].vert.size(), res);
		}
		return res;
	}

	void makeIndexDFS1(int p, vector<int>& list, int* toList) {
		vector<pair<int, int>> z;
		vector<pair<int, int>> _z;
		vector<pair<int, int>> zz;
		judge ju;
		zz.resize(list.size());
		Tree[p].pos.resize(Tree[p].vert.size() + 1);
		Tree[p].skyToAnc.resize(list.size());
		for (int i = 0; i < list.size(); i++) {
			Tree[p].skyToAnc[i].resize(list.size());
		}
		for (int i = 0; i < Tree[p].vert.size(); i++) {
			int j;
			for (j = 0; j < list.size(); j++)
				if (list[j] == Tree[p].vert[i])
					break;
			Tree[p].pos[i] = j;
		}
		Tree[p].pos[Tree[p].vert.size()] = list.size();
		for (int i = 0; i < list.size(); i++) {
			for (int _j = 0; _j < Tree[p].skyToAnc[i].size(); _j++) {
				Tree[p].skyToAnc[i][_j].first = INF;
				Tree[p].skyToAnc[i][_j].second = INF;
			}
		}
		int x = Tree[p].uniqueVertex;

		vector<pair<int, int>> upperLeft;
		vector<pair<int, int>> downRight;
		vector<pair<int, int>> upperSet;
		vector<pair<int, int>> lowerSet;
		vector<vector<pair<int, int>>> hops;
		vector<int> ph;
		judgeCost juC;
		hops.resize(Tree[p].vert.size());
		int bbxMin, bbyMin, bbxMax, bbyMax;
		for (int i = 0; i < Tree[p].vert.size(); i++) {
			vector<pair<int, int>> assign;
			mergeSort(Tree[p].skyToAnc[toList[Tree[p].vert[i]]], Tree[p].SK[i], assign);
			Tree[p].skyToAnc[toList[Tree[p].vert[i]]] = screenSkylines(assign);
			int x = Tree[p].vert[i];
			int k;
			for (k = 0; k < list.size(); k++) {
				if (list[k] == x) break;
			}
			for (int j = 0; j < list.size(); j++) {
				int y = list[j];
				if (k < j) {
					z = distanceQueryAncestorToPosterity(x, y);
				}
				else
				{
					z = distanceQueryAncestorToPosterity(y, x);
				}
				int px1 = z[0].first + Tree[p].skyToAnc[toList[Tree[p].vert[i]]][0].first;
				int py1 = z[0].second + Tree[p].skyToAnc[toList[Tree[p].vert[i]]][0].second;
				upperLeft.push_back(make_pair(px1, py1));
				hops[i].push_back(make_pair(px1, py1));
				int px2 = z[z.size() - 1].first + Tree[p].skyToAnc[toList[Tree[p].vert[i]]][Tree[p].skyToAnc[toList[Tree[p].vert[i]]].size() - 1].first;
				int py2 = z[z.size() - 1].second + Tree[p].skyToAnc[toList[Tree[p].vert[i]]][Tree[p].skyToAnc[toList[Tree[p].vert[i]]].size() - 1].second;
				hops[i].push_back(make_pair(px2, py2));
				downRight.push_back(make_pair(px2, py2));
			}
		}
		upperSet = upperLeft;
		lowerSet = downRight;
		sort(upperLeft.begin(), upperLeft.end(), ju);
		sort(downRight.begin(), downRight.end(), juC);

		bbxMin = upperLeft[0].first;
		bbyMax = upperLeft[0].second;

		bbxMax = downRight[0].first;
		bbyMin = downRight[0].second;
		
		for (int i = 0; i < hops.size(); ++i) {
			if (bbxMin <= hops[i][1].first && bbxMax >= hops[i][0].first && bbyMin <= hops[i][0].second && bbyMax >= hops[i][1].second) {
				ph.push_back(i);
			}
		}
		upperSet.clear();
		lowerSet.clear();
		upperLeft.clear();
		downRight.clear();
		hops.clear();

		for (int i = 0; i < ph.size(); i++) {

			int x = Tree[p].vert[ph[i]];
			int k;
			for (k = 0; k < list.size(); k++) {
				if (list[k] == x) break;
			}
			for (int j = 0; j < list.size(); j++) {
				int y = list[j];
				if (k < j) {
					z = distanceQueryAncestorToPosterity(x, y);
				}
				else
				{
					z = distanceQueryAncestorToPosterity(y, x);
				}

				_z = concatenation(z, Tree[p].skyToAnc[toList[Tree[p].vert[ph[i]]]]);
				mergeSort(_z, Tree[p].skyToAnc[toList[y]], zz);
				Tree[p].skyToAnc[toList[y]] = screenSkylines(zz);

			}
		}
		ph.clear();
		toList[Tree[p].uniqueVertex] = list.size();
		list.push_back(Tree[p].uniqueVertex);
		for (int i = 0; i < Tree[p].ch.size(); i++) {
			makeIndexDFS(Tree[p].ch[i], list, toList);
		}
		list.pop_back();
		sort(Tree[p].pos.begin(), Tree[p].pos.end());
	}

	void makeIndexDFS(int p, vector<int>& list, int* toList) {
		vector<pair<int, int>> z;
		vector<pair<int, int>> _z;
		vector<pair<int, int>> zz;
		judge ju;
		zz.resize(list.size());
		Tree[p].pos.resize(Tree[p].vert.size() + 1);
		Tree[p].skyToAnc.resize(list.size());
		for (int i = 0; i < list.size(); i++) {
			Tree[p].skyToAnc[i].resize(list.size());
		}
		for (int i = 0; i < Tree[p].vert.size(); i++) {
			int j;
			for (j = 0; j < list.size(); j++)
				if (list[j] == Tree[p].vert[i])
					break;
			Tree[p].pos[i] = j;
		}
		Tree[p].pos[Tree[p].vert.size()] = list.size();
		for (int i = 0; i < list.size(); i++) {
			for (int _j = 0; _j < Tree[p].skyToAnc[i].size(); _j++) {
				Tree[p].skyToAnc[i][_j].first = INF;
				Tree[p].skyToAnc[i][_j].second = INF;
			}
		}
		int x = Tree[p].uniqueVertex;

		for (int i = 0; i < Tree[p].vert.size(); i++) {
			vector<pair<int, int>> assign;
			mergeSort(Tree[p].skyToAnc[toList[Tree[p].vert[i]]], Tree[p].SK[i], assign);
			Tree[p].skyToAnc[toList[Tree[p].vert[i]]] = screenSkylines(assign);
			int x = Tree[p].vert[i];
			int k;
			for (k = 0; k < list.size(); k++) {
				if (list[k] == x) break;
			}
			for (int j = 0; j < list.size(); j++) {
				int y = list[j];
				if (k < j) {
					z = distanceQueryAncestorToPosterity(x, y);
				}
				else
				{
					z = distanceQueryAncestorToPosterity(y, x);
				}
				_z = concatenation(z, Tree[p].skyToAnc[toList[Tree[p].vert[i]]]);

				mergeSort(_z, Tree[p].skyToAnc[toList[y]], zz);
				Tree[p].skyToAnc[toList[y]] = screenSkylines(zz);
			}
		}
		toList[Tree[p].uniqueVertex] = list.size();
		list.push_back(Tree[p].uniqueVertex);
		for (int i = 0; i < Tree[p].ch.size(); i++) {
			makeIndexDFS(Tree[p].ch[i], list, toList);
		}
		list.pop_back();
		sort(Tree[p].pos.begin(), Tree[p].pos.end());
	}

	void makeIndexDFS_new(int p, vector<int>& list, int* toList) {
		vector<pair<int, int>> z;
		vector<pair<int, int>> _z;
		judge ju;
		Tree[p].pos.resize(Tree[p].vert.size() + 1);
		Tree[p].skyToAnc.resize(list.size());
		for (int i = 0; i < list.size(); i++) {
			Tree[p].skyToAnc[i].resize(list.size());
		}
		for (int i = 0; i < Tree[p].vert.size(); i++) {
			int j;
			for (j = 0; j < list.size(); j++)
				if (list[j] == Tree[p].vert[i])
					break;
			Tree[p].pos[i] = j;
		}
		Tree[p].pos[Tree[p].vert.size()] = list.size();
		pair<int, int> pp = make_pair(INF, INF);
		vector<pair<int, int> > vp(1, pp);
		Tree[p].skyToAnc.assign(list.size(), vp);

		int x = Tree[p].uniqueVertex;

		vector<vector<int> > vvInput;
		vector<int> vT;
		vvInput.assign(ism, vT);
		for (int i = 0; i < list.size(); i++) {
			vvInput[i % ism].push_back(list[i]);
		}

		cout << "Tree " << p << "\t width:" << Tree[p].pos.size() << endl;
		vector<vector<pair<int, int> > > vvpSkyline;
		for (int i = 0; i < Tree[p].vert.size(); i++)
		{
			int x = Tree[p].vert[i];
			vvpSkyline.push_back(Tree[p].SK[i]);
		}
		boost::thread_group threads;
		for (int i = 0; i < ism; i++)
		{
			if (!vvInput[i].empty())
			{
				threads.add_thread(new boost::thread(&Bound_Tree_Decomposition::propagation, this, p, boost::ref(vvInput[i]), boost::ref(toList), list.size(), boost::ref(list), boost::ref(vvpSkyline)));
			}
		}
		threads.join_all();

		toList[Tree[p].uniqueVertex] = list.size();
		list.push_back(Tree[p].uniqueVertex);
		for (int i = 0; i < Tree[p].ch.size(); i++) {
			makeIndexDFS_new(Tree[p].ch[i], list, toList);
		}
		list.pop_back();
		
		sort(Tree[p].pos.begin(), Tree[p].pos.end());
	}

	void propagation(int p, vector<int>& list, int* toList, int ancNum, vector<int>& ancList, vector<vector<pair<int, int>>>& vvpSkyline) {
		vector<pair<int, int>> z;
		vector<pair<int, int>> _z;
		for (int i = 0; i < list.size(); i++)
		{
			int y = list[i];
			for (int j = 0; j < Tree[p].vert.size(); j++)
			{
				int x = Tree[p].vert[j];
				vector<pair<int, int>> zz;
				if (x == y)
				{
					vector<pair<int, int>> assign;
					mergeSort(vvpSkyline[j], Tree[p].SK[j], assign);
					Tree[p].skyToAnc[toList[x]] = screenSkylines(assign);
					continue;
				}
				int ix, iy;
				for (int k = 0; k < ancList.size(); k++)
				{
					if (ancList[k] == x)
						ix = k;
					if (ancList[k] == y)
						iy = k;
				}

				if (ix < iy)
				{
					z = distanceQueryAncestorToPosterity(x, y);
				}
				else
				{
					z = distanceQueryAncestorToPosterity(y, x);
				}
				_z = concatenation(z, vvpSkyline[j]);
				mergeSort(_z, Tree[p].skyToAnc[toList[y]], zz);

				Tree[p].skyToAnc[toList[y]] = screenSkylines(zz);
			}
		}
	}


void a_makeIndexDFS_new(int p, vector<int>& list, int* toList) {
		vector<pair<int, int>> z;
		vector<pair<int, int>> _z;
		judge ju;
		judgeCost juc;
		Tree[p].pos.resize(Tree[p].vert.size() + 1);
		Tree[p].skyToAnc.resize(list.size());

		for (int i = 0; i < list.size(); i++) {
			Tree[p].skyToAnc[i].resize(list.size());
		}
		for (int i = 0; i < Tree[p].vert.size(); i++) {
			int j;
			for (j = 0; j < list.size(); j++)
				if (list[j] == Tree[p].vert[i])
					break;
			Tree[p].pos[i] = j;
		}
		Tree[p].pos[Tree[p].vert.size()] = list.size();
		pair<int, int> pp = make_pair(INF, INF);
		vector<pair<int, int> > vp(1, pp);
		Tree[p].skyToAnc.assign(list.size(), vp);
		int x = Tree[p].uniqueVertex;
		vector<vector<int> > vvInput;
		vector<int> vT;
		vvInput.assign(ism, vT);
		for (int i = 0; i < list.size(); i++) {
			vvInput[i % ism].push_back(list[i]);
		}

		cout << "Tree " << p << "\t height:" << Tree[p].height << endl;
		vector<vector<pair<int, int> > > vvpSkyline;
		for (int i = 0; i < Tree[p].vert.size(); i++)
		{
			int x = Tree[p].vert[i];
			sort(Tree[p].SK[i].begin(),Tree[p].SK[i].end(),juc);
			vvpSkyline.push_back(Tree[p].SK[i]);
		}
		boost::thread_group threads;
		for (int i = 0; i < ism; i++)
		{
			if (!vvInput[i].empty())
			{
				threads.add_thread(new boost::thread(&Bound_Tree_Decomposition::a_propagation, this, p, boost::ref(vvInput[i]), boost::ref(toList), list.size(), boost::ref(list), boost::ref(vvpSkyline)));
			}
		}
		threads.join_all();
		Tree[p].a_rest=a_g_d/alc[p][root];

		toList[Tree[p].uniqueVertex] = list.size();
		list.push_back(Tree[p].uniqueVertex);
		for (int i = 0; i < Tree[p].ch.size(); i++) {
			a_makeIndexDFS_new(Tree[p].ch[i], list, toList);
		}
		list.pop_back();
		sort(Tree[p].pos.begin(), Tree[p].pos.end());
	}

	double theta(vector<pair<int,int>>& h1, vector<pair<int,int>>& h2, double a_l){
		int size1=h1.size();
		int size2=h2.size();

		double wmax, wmin,a_i;
		a_i=1;
		wmax=h1[0].first+h2[0].first;
		wmin=h1[size1-1].first+h2[size2-1].first;
		double pre_size=log(wmax/wmin)/log(a_l);
		double cat_size=size1*size2;
		
		if(cat_size>pre_size){
			a_i=a_l;
		}
		return a_i;
	}



	double theta_non(vector<pair<int,int>>& h1, vector<pair<int,int>>& h2, double a_l){
		double a_i;
		a_i=a_l;
		return a_i;
	}
	void a_propagation(int p, vector<int>& list, int* toList, int ancNum, vector<int>& ancList, vector<vector<pair<int, int>>>& vvpSkyline) {
		judgeCost juc;
		vector<pair<int, int>> z;
		vector<pair<int, int>> _z;
		double a_l;
		int depth;
		depth=Tree[p].height;
		int parent=Tree[p].pa;

		if(parent==root){
			Tree[root].a_rest=a_g_d;
		}

		for (int i = 0; i < list.size(); i++)
		{
			int y = list[i];
			int ay=belong[y];
			int h_ay=Tree[ay].height;
			double kk=abs(h_ay-depth);
			double nn=heightMax-h_ay;
			double rate=1/(nn-1);
			if(kk==1){
				a_l=pow(Tree[ay].a_rest,rate);
			}
			else{
				double a_l1=pow(Tree[ay].a_rest,rate);
				double rate1=(kk-1)/(heightMax-1);
				double a_l2=pow(a_g_d,rate1);
				a_l=a_l1*a_l2;
			}
			double a_i;
			for (int j = 0; j < Tree[p].vert.size(); j++)
			{
				int x = Tree[p].vert[j];
				vector<pair<int, int>> zz;
				
				if (x == y)
				{
					vector<pair<int, int>> assign;
					a_mergeSort(vvpSkyline[j], Tree[p].SK[j], assign);
					Tree[p].skyToAnc[toList[x]] = a_screenSkylines(assign,1);
					continue;
				}
				int ix, iy;
				for (int k = 0; k < ancList.size(); k++)
				{
					if (ancList[k] == x)
						ix = k;
					if (ancList[k] == y)
						iy = k;
				}

				if (ix < iy)
				{
					z = distanceQueryAncestorToPosterity(x, y);
				}
				else
				{
					z = distanceQueryAncestorToPosterity(y, x);
				}
				sort(z.begin(),z.end(),juc);
				sort(Tree[p].skyToAnc[toList[y]].begin(),Tree[p].skyToAnc[toList[y]].end(),juc);
			
				a_i=theta_non(z,vvpSkyline[j],a_l);
				_z = a_concatenation(z, vvpSkyline[j],a_i);
				a_mergeSort(_z, Tree[p].skyToAnc[toList[y]], zz);
				Tree[p].skyToAnc[toList[y]] = al_screenSkylines(zz,a_i,p,ay);
			}
			if(alc[p][ay]<=a_i/rec[p][ay]){
				alc[p][ay]=a_i/rec[p][ay];		
			}
		}
	}


void a_mergeSort(const vector<pair<int, int>>& s1, const vector<pair<int, int>>& s2, vector<pair<int, int>>& result) {
		int size1 = s1.size();
		int size2 = s2.size();
		result.resize(size1 + size2);
		int tmp_index, tmp_1_index, tmp_2_index;
		tmp_1_index = 0;
		tmp_2_index = 0;
		for (tmp_index = 0; (tmp_1_index < size1) && (tmp_2_index < size2); tmp_index++) {
			if (s1[tmp_1_index].second <= s2[tmp_2_index].second) {

				result[tmp_index] = s1[tmp_1_index];
				tmp_1_index++;
			}
			else {
				result[tmp_index] = s2[tmp_2_index];
				tmp_2_index++;
			}
		}
		if (tmp_1_index == size1) {
			while (tmp_2_index < size2) {
				result[tmp_index] = s2[tmp_2_index];
				tmp_index++;
				tmp_2_index++;
			}
		}
		else if (tmp_2_index == size2) {
			while (tmp_1_index < size1) {
				result[tmp_index] = s1[tmp_1_index];
				tmp_index++;
				tmp_1_index++;
			}
		}
	}

	vector<pair<int, int>> a_screenSkylines(vector<pair<int, int>>& input, double a_l) {
		vector<pair<int, int>> result;
		result.push_back(input[0]);
		int size = input.size();

		if (size <= 1) {
			return result;
		}
		else {
			for (int i = 1; i < size; i++) {
				if (result[result.size() - 1].second== input[i].second && result[result.size() - 1].first > input[i].first) {
					result[result.size() - 1].first = input[i].first;
				}
				if (result[result.size() - 1].second != input[i].second && result[result.size() - 1].first > a_l*input[i].first) {
					result.push_back(input[i]);
				}
				else {
					continue;
				}
			}
			return result;
		}
	}

	vector<pair<int, int>> al_screenSkylines(vector<pair<int, int>>& input, double a_l, int p, int anc) {
		vector<pair<int, int>> result;
		result.push_back(input[0]);
		int size = input.size();

		if (size <= 1) {
			return result;
		}
		else {
			for (int i = 1; i < size; i++) {
				if (result[result.size() - 1].second== input[i].second && result[result.size() - 1].first > input[i].first) {
					result[result.size() - 1].first = input[i].first;
				}
				if (result[result.size() - 1].second != input[i].second && result[result.size() - 1].first > a_l*input[i].first) {
					result.push_back(input[i]);
				}
				else {
					if(result[result.size() - 1].first > input[i].first){
						double tmp;
						tmp=double(result[result.size()-1].first)/double(input[i].first);
						if (rec[p][anc]>tmp){
							rec[p][anc]=tmp;
						}
					}
				}
			}
			if(rec[p][anc]<a_l){
				rec[p][anc]=a_l/rec[p][anc];
			}
			else{
				rec[p][anc]=1;
			}
			return result;
		}
	}

	vector<pair<int, int>> a_concatenation(const vector<pair<int, int>>& h1, const vector<pair<int, int>>& h2, const double a_l) {
		int size1 = h1.size();
		int size2 = h2.size();
		pair<int, int> wc;
		vector<pair<int, int>> result;
		if (size1 == 0 || size2 == 0) {
			result.push_back(make_pair(INF, INF));
			return result;
		}

		vector<vector<bool>> visited(size1, vector<bool>(size2, false));
		priority_queue<hop, vector<hop>, pk_c> pq;
		int w = h1[0].first + h2[0].first;
		int c = h1[0].second + h2[0].second;
		wc = make_pair(w, c);
		pq.push(hop(0, 0, wc));
		visited[0][0] = true;
		int count = 0;
		result.push_back(wc);
		while (!pq.empty()) {
			hop el = pq.top();
			pq.pop();
			int n = result.size();
			if (a_l*el.val.first < result[n - 1].first) {
				result.push_back(el.val);
			}
			int ex1 = el.i + 1;
			int ey1 = el.j;
			int ex2 = el.i;
			int ey2 = el.j + 1;
			if (ex1 < size1) {
				if (!visited[ex1][ey1]) {
					w = h1[ex1].first + h2[ey1].first;
					c = h1[ex1].second + h2[ey1].second;
					wc = make_pair(w, c);
					if (a_l*w < result[n - 1].first) {
						pq.push(hop(ex1, ey1, wc));
						visited[ex1][ey1]=true;
					}
				}
			}
			if (ey2 < size2) {
				if (!visited[ex2][ey2]) {
					w = h1[ex2].first + h2[ey2].first;
					c = h1[ex2].second + h2[ey2].second;
					wc = make_pair(w, c);
					if (a_l*w < result[n - 1].first) {
						pq.push(hop(ex2, ey2, wc));
						visited[ex2][ey2]=true;
					}
				}
			}

		}
		return result;
	}

	void makeIndex() {
		makeRMQ();
		H.EdgeInitialize();
	}

	void readPositionIndex(char* file, vector<vector<int>>& position) {
		ifstream in(file);
		string line;
		vector<int> temp;
		while (getline(in, line)) {
			stringstream ss(line);
			int x;
			while (ss >> x) {
				temp.push_back(x);
			}
			position.push_back(temp);
			temp.clear();
		}
		in.close();
	}

	void readSkylineIndex(char* file) {
		ifstream in(file);
		int n, ns, vert, sks;
		for (; in >> n >> ns >> vert >> sks;) {
			Tree[n].skyToAnc.resize(ns);
			for (int i = 0; i < sks; i++) {
				int x, y;
				in >> x >> y;
				Tree[n].skyToAnc[vert].push_back(make_pair(x, y));
			}
		}
		in.close();
	}

	void readSkylineIndex1(char* file) {
		ifstream in(file);
		string line;
		vector<pair<int, int>> temp1;
		vector<vector<pair<int, int>>> temp;

		while (getline(in, line)) {
			int num = 0;
			int num1;
			stringstream ss;
			ss.clear();
			ss.str("");
			ss << line;
			int i, j, k;
			int x, y;
			ss >> i;
			ss >> j;
			ss >> k;
			Tree[i].skyToAnc.resize(k);
			while (ss >> x >> y) {
				Tree[i].skyToAnc[j].push_back(make_pair(x, y));
			}
			temp.clear();
		}
		in.close();
	}

	void printIndex() {
		ofstream outfile;
		outfile.open("index/position_index(bound).txt");
		for (int i = 0; i < Tree.size(); i++) {
			for (int j = 0; j < Tree[i].pos.size(); j++) {
				outfile << Tree[i].pos[j] << '\t';
			}
			outfile << endl;
		}
		outfile.close();

		ofstream outfile1;
		outfile1.open("index/skyline_index(bound).txt");
		for (int i = 0; i < Tree.size(); i++) {
			for (int j = 0; j < Tree[i].skyToAnc.size(); j++) {
				outfile1 << i << '\t' << Tree[i].skyToAnc.size() << '\t' << j << '\t' << Tree[i].skyToAnc[j].size() << '\t';
				for (int k = 0; k < Tree[i].skyToAnc[j].size(); k++) {
					outfile1 << Tree[i].skyToAnc[j][k].first << '\t' << Tree[i].skyToAnc[j][k].second << '\t';
				}
				outfile1 << endl;
			}
		}
		outfile1.close();
	}

	void a_printIndex_reduce(){
		ofstream outfile;
		outfile.open("index/reduce_position_index(bound).txt");
		for (int i = 0; i < Tree.size(); i++) {
			for (int j = 0; j < Tree[i].pos.size(); j++) {
				outfile << Tree[i].pos[j] << '\t';
			}
			outfile << endl;
		}
		outfile.close();

		ofstream outfile1;
		outfile1.open("index/reduce_skyline_index(bound).txt");
		for (int i = 0; i < Tree.size(); i++) {
			//outfile1 << i << '\t';
			for (int j = 0; j < Tree[i].skyToAnc.size(); j++) {
				outfile1 << i << '\t' << Tree[i].skyToAnc.size() << '\t' << j << '\t' << Tree[i].skyToAnc[j].size() << '\t';
				for (int k = 0; k < Tree[i].skyToAnc[j].size(); k++) {
					outfile1 << Tree[i].skyToAnc[j][k].first << '\t' << Tree[i].skyToAnc[j][k].second << '\t';
				}
				outfile1 << endl;
			}
		}
		outfile1.close();
	}

	void a_printIndex() {
		ofstream outfile;
		outfile.open("index/a_position_index(bound).txt");
		for (int i = 0; i < Tree.size(); i++) {
			for (int j = 0; j < Tree[i].pos.size(); j++) {
				outfile << Tree[i].pos[j] << '\t';
			}
			outfile << endl;
		}
		outfile.close();

		ofstream outfile1;
		outfile1.open("index/a_skyline_index(bound).txt");
		for (int i = 0; i < Tree.size(); i++) {
			//outfile1 << i << '\t';
			for (int j = 0; j < Tree[i].skyToAnc.size(); j++) {
				outfile1 << i << '\t' << Tree[i].skyToAnc.size() << '\t' << j << '\t' << Tree[i].skyToAnc[j].size() << '\t';
				for (int k = 0; k < Tree[i].skyToAnc[j].size(); k++) {
					outfile1 << Tree[i].skyToAnc[j][k].first << '\t' << Tree[i].skyToAnc[j][k].second << '\t';
				}
				outfile1 << endl;
			}
		}
		outfile1.close();
	}

	void a_printIndex_theta() {
		ofstream outfile;
		outfile.open("index/rec_position_index(bound).txt");
		for (int i = 0; i < Tree.size(); i++) {
			for (int j = 0; j < Tree[i].pos.size(); j++) {
				outfile << Tree[i].pos[j] << '\t';
			}
			outfile << endl;
		}
		outfile.close();

		ofstream outfile1;
		outfile1.open("index/rec_skyline_index(bound).txt");
		for (int i = 0; i < Tree.size(); i++) {
			//outfile1 << i << '\t';
			for (int j = 0; j < Tree[i].skyToAnc.size(); j++) {
				outfile1 << i << '\t' << Tree[i].skyToAnc.size() << '\t' << j << '\t' << Tree[i].skyToAnc[j].size() << '\t';
				for (int k = 0; k < Tree[i].skyToAnc[j].size(); k++) {
					outfile1 << Tree[i].skyToAnc[j][k].first << '\t' << Tree[i].skyToAnc[j][k].second << '\t';
				}
				outfile1 << endl;
			}
		}
		outfile1.close();
	}

	void a_resort() {
		judge ju;
		for (int i = 0; i < Tree.size(); i++) {
			for (int j = 0; j < Tree[i].skyToAnc.size(); j++) {
				sort(Tree[i].skyToAnc[j].begin(),Tree[i].skyToAnc[j].end(),ju);
			}
		}
	}

	void loadIndex(int a) {
		if (a == 0) {
			cout << "Re-building CSP-2Hop index" << endl;

			std::chrono::high_resolution_clock::time_point t1;
			std::chrono::high_resolution_clock::time_point t2;
			std::chrono::duration<double> time_span;
			t1 = std::chrono::high_resolution_clock::now();
			makeIndex2();
			reducePos();
			t2 = std::chrono::high_resolution_clock::now();
			time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
			cout << "Boundary Tree Time:" << time_span.count() << endl;
			printIndex();
		}
		else if (a == 99) {
			cout << "Re-building Approximate-Recycle CSP-2Hop index" << endl;
			std::chrono::high_resolution_clock::time_point t1;
			std::chrono::high_resolution_clock::time_point t2;
			std::chrono::duration<double> time_span;
			t1 = std::chrono::high_resolution_clock::now();
			a_makeIndex2();
			reducePos();
			t2 = std::chrono::high_resolution_clock::now();
			time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
			cout << "Boundary Tree Time:" << time_span.count() << endl;
			a_resort();
			a_printIndex_theta();
		}
		else if(a==1001){
			cout << "Re-building Approximate-Reduce CSP-2Hop index" << endl;
			std::chrono::high_resolution_clock::time_point t1;
			std::chrono::high_resolution_clock::time_point t2;
			std::chrono::duration<double> time_span;
			t1 = std::chrono::high_resolution_clock::now();
			makeIndex2();
			reducePos();
			t2 = std::chrono::high_resolution_clock::now();
			time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
			cout << "Boundary Tree Time:" << time_span.count() << endl;
			a_printIndex_reduce();
		}
		else if (a == 100) {
			cout << "Re-building Approximate CSP-2Hop index" << endl;
			std::chrono::high_resolution_clock::time_point t1;
			std::chrono::high_resolution_clock::time_point t2;
			std::chrono::duration<double> time_span;
			t1 = std::chrono::high_resolution_clock::now();
			a_makeIndex2();
			reducePos();
			t2 = std::chrono::high_resolution_clock::now();
			time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
			cout << "Boundary Tree Time:" << time_span.count() << endl;
			a_resort();
			a_printIndex();
		}
		else if (a == 10) {
			cout << "Re-building CSP-2Hop index(optimal)" << endl;
			makeIndex3();
			reducePos();
			printIndex();
		}
		else {
			cout << "Loading CSP-2Hop index" << endl;
			vector<vector<int>> ps;
			char positionFile[255] = "index/position_index(bound).txt";
			char skylineFile[255] = "index/skyline_index(bound).txt";
			readPositionIndex(positionFile, ps);
			for (int i = 0; i < ps.size(); i++) {
				Tree[i].pos = ps[i];
			}
			readSkylineIndex(skylineFile);
		}
	}
	void makeIndex2() {
		vector<int> list;
		list.clear();
		int* toList;
		toList = (int*)malloc(sizeof(int) * (G.n + 1));
		Tree[root].pos.clear();

		toList[Tree[root].uniqueVertex] = 0;
		list.push_back(Tree[root].uniqueVertex);
		Tree[root].pos.push_back(0);

		for (int i = 0; i < Tree[root].ch.size(); ++i) {
			//makeIndexDFS(Tree[root].ch[i], list, toList);
			makeIndexDFS_new(Tree[root].ch[i], list, toList);
			break;
		}
		free(toList);

	}
	void a_makeIndex2() {
		vector<int> list;
		list.clear();
		int* toList;
		toList = (int*)malloc(sizeof(int) * (G.n + 1));
		Tree[root].pos.clear();

		toList[Tree[root].uniqueVertex] = 0;
		list.push_back(Tree[root].uniqueVertex);
		Tree[root].pos.push_back(0);
		vector<double> tmp_re;
		vector<double> tmp_alc;
	    tmp_re.assign(G.n+1,INF);
		tmp_alc.assign(G.n+1,1);
	    rec.assign(G.n+1,tmp_re);
	    alc.assign(G.n+1,tmp_alc);
		for (int i = 0; i < Tree[root].ch.size(); ++i) {
			a_makeIndexDFS_new(Tree[root].ch[i], list, toList);
			break;
		}
		free(toList);

	}

	void makeIndex3() {
		vector<int> list;
		list.clear();
		int* toList;
		toList = (int*)malloc(sizeof(int) * (G.n + 1));
		Tree[root].pos.clear();

		toList[Tree[root].uniqueVertex] = 0;
		list.push_back(Tree[root].uniqueVertex);
		Tree[root].pos.push_back(0);

		for (int i = 0; i < Tree[root].ch.size(); i++) {
			makeIndexDFS1(Tree[root].ch[i], list, toList);
			break;
		}
		free(toList);
	}

	void reducePosDFS(int p) {
		if (Tree[p].ch.size() == 2) {
			int t = Tree[p].ch[0];
			if (Tree[Tree[p].ch[0]].pos.size() > Tree[Tree[p].ch[1]].pos.size())
				t = Tree[p].ch[1];
			Tree[p].pos = Tree[t].pos;
			if (Tree[p].pos.size() != 0) {
				Tree[p].pos.erase(Tree[p].pos.begin() + (Tree[p].pos.size() - 1));
			}

		}
		for (int i = 0; i < Tree[p].ch.size(); i++)
			reducePosDFS(Tree[p].ch[i]);
	}
	void reducePos() {
		reducePosDFS(root);
	}
	void cntSize() {
		long long s_nonroot = 0;
		long long s_size = 0;

		long long s_dis = 0;
		for (int i = 0; i < Tree.size(); ++i) {
			s_nonroot += Tree[i].height - 1;
			s_size += Tree[i].vert.size();
			s_dis += Tree[i].height;
		}
		long long s_root = (long long)Tree[0].vert.size() * (long long)G.n;
		printf("tree width: %d\n", tree_width);
		printf("nonroot idx size = %0.3lfGB, avg node size=%0.3lf, avg label size=%0.3lf\n",
			s_nonroot * 4.0 / 1000000000.0, s_size * 1.0 / G.n, s_dis * 1.0 / G.n);
	}

	inline pair<int, int> concatenation_C(const vector<pair<int, int>>& h1, const vector<pair<int, int>>& h2, const int& C) {
		int size1 = h1.size();
		int size2 = h2.size();
		pair<int, int> wc;
		vector<pair<int, int>> result;
		/*if (size1 == 0 || size2 == 0) {
			result.push_back(make_pair(INF, INF));
			return result[0];
		}*/
		vector<vector<bool>> visited(size1, vector<bool>(size2, false));
		priority_queue<hop, vector<hop>, pk> pq;
		int w = h1[0].first + h2[0].first;
		int c = h1[0].second + h2[0].second;
		wc = make_pair(w, c);
		pq.push(hop(0, 0, wc));
		visited[0][0] = true;
		int count = 0;
		result.push_back(wc);
		while (!pq.empty()) {
			hop el = pq.top();
			pq.pop();
			int n = result.size();
			if (n > 0) {
				if (el.val.second <= result[n - 1].second) {
					result.push_back(el.val);
				}
			}
			if (el.val.second <= C){ 
				break;
			}
			int ex1 = el.i + 1;
			int ey1 = el.j;
			int ex2 = el.i;
			int ey2 = el.j + 1;
			if (ex1 < size1) {
				if (!visited[ex1][ey1]) {
					w = h1[ex1].first + h2[ey1].first;
					c = h1[ex1].second + h2[ey1].second;
					wc = make_pair(w, c);
					if (c <= result[n - 1].second) {
						pq.push(hop(ex1, ey1, wc));
						visited[ex1][ey1]=true;
					}
				}
			}
			if (ey2 < size2) {
				if (!visited[ex2][ey2]) {
					w = h1[ex2].first + h2[ey2].first;
					c = h1[ex2].second + h2[ey2].second;
					wc = make_pair(w, c);
					if (c <= result[n - 1].second) {
						pq.push(hop(ex2, ey2, wc));
						visited[ex2][ey2]=true;
					}
				}
			}

		}
		/*if (result[result.size() - 1].second > C) {
			result.push_back(make_pair(INF, INF));
			return result[result.size() - 1];
		}*/
		//cout<<"result size:"<<result.size()<<endl;
		return result[result.size() - 1];
	}


	vector<pair<int, int>> upper;
	vector<pair<int, int>> lower;
	inline vector<int> multiHop(vector<vector<pair<int, int>>>& x, vector<vector<pair<int, int>>>& y, vector<int>& p) {
		int bbxMin, bbyMin, bbxMax, bbyMax;
		vector<int> res;
		vector<vector<pair<int, int>>> hops;
		vector<pair<int, int>> upperLeft;
		vector<pair<int, int>> downRight;
		judge ju;
		judgeCost juC;
		upper.clear();
		lower.clear();
		tmp_rectangle.clear();
		hops.resize(p.size());
		if (p.size() == 1) {
			res.push_back(0);
			return res;
		}
		for (int i = 0; i < p.size(); ++i) {
				int px1 = x[p[i]][0].first + y[p[i]][0].first;
				int py1 = x[p[i]][0].second + y[p[i]][0].second;
				upperLeft.push_back(make_pair(px1, py1));
				hops[i].push_back(make_pair(px1, py1));
				int px2 = x[p[i]][x[p[i]].size() - 1].first + y[p[i]][y[p[i]].size() - 1].first;
				int py2 = x[p[i]][x[p[i]].size() - 1].second + y[p[i]][y[p[i]].size() - 1].second;
				hops[i].push_back(make_pair(px2, py2));
				downRight.push_back(make_pair(px2, py2));
		}
		upper = upperLeft;
		lower = downRight;
		sort(upperLeft.begin(), upperLeft.end(), ju);
		sort(downRight.begin(), downRight.end(), juC);

		bbxMin = upperLeft[0].first;
		bbyMax = upperLeft[0].second;
		tmp_rectangle.push_back(make_pair(bbxMin, bbyMax));

		bbxMax = downRight[0].first;
		bbyMin = downRight[0].second;
		tmp_rectangle.push_back(make_pair(bbxMax, bbyMin));

		for (int i = 0; i < p.size(); ++i) {
			if (bbxMin <= hops[i][1].first && bbxMax >= hops[i][0].first && bbyMin <= hops[i][0].second && bbyMax >= hops[i][1].second) {
				res.push_back(i);
			}
		}
		return res;
	}

	inline vector<pair<int, int>> skylines_bound(int p, int q) {
		vector<pair<int, int>> temp, temp1;
		vector<pair<int, int>> candi;
		vector<pair<int, int>> s;
		vector<vector<pair<int, int>>> sky_x, sky_y;
		temp.push_back(make_pair(0, 0));
		temp1.push_back(make_pair(INF, INF));
		if (p == q) return temp;
		int x = belong[p], y = belong[q];
		int lca = LCAQuery(x, y);
		if (lca == x || lca == y) {
			if (lca == y) {
				int v = y;
				y = x;
				x = v;
				v = p;
				p = q;
				q = v;
			}
			int s_size = Tree[y].skyToAnc[Tree[x].pos[Tree[x].pos.size() - 1]].size();
			if (s_size <= 1) {
				s = Tree[y].skyToAnc[Tree[x].pos[Tree[x].pos.size() - 1]];
			}
			else
			{
				vector<pair<int, int>> tmp_s = Tree[y].skyToAnc[Tree[x].pos[Tree[x].pos.size() - 1]];
				s.push_back(tmp_s[0]);
				s.push_back(tmp_s[s_size - 1]);
			}
			return s;
		}
		else {
			int w = INF;
			int c = INF;

			pair<int, int> wc;
			sky_x = Tree[x].skyToAnc;
			sky_y = Tree[y].skyToAnc;
			vector<int> p;
			p.clear();
			p=Tree[lca].pos;
			int ps = Tree[lca].pos.size();
			vector<pair<int, int>> h;
			s = rectangle(sky_x, sky_y, p);
		}
		return s;
	}

	inline vector<pair<int, int>> rectangle(vector<vector<pair<int, int>>>& x, vector<vector<pair<int, int>>>& y, vector<int>& p) {
		vector<pair<int, int>> res;
		vector<vector<pair<int, int>>> hops;
		vector<pair<int, int>> upperLeft;
		vector<pair<int, int>> downRight;
		int bbxMin, bbyMin, bbxMax, bbyMax;
		judge ju;
		judgeCost juC;
		hops.resize(p.size());

		for (int i = 0; i < p.size(); ++i) {
			int px1 = x[p[i]][0].first + y[p[i]][0].first;
			int py1 = x[p[i]][0].second + y[p[i]][0].second;
			upperLeft.push_back(make_pair(px1, py1));
			hops[i].push_back(make_pair(px1, py1));
			int px2 = x[p[i]][x[p[i]].size() - 1].first + y[p[i]][y[p[i]].size() - 1].first;
			int py2 = x[p[i]][x[p[i]].size() - 1].second + y[p[i]][y[p[i]].size() - 1].second;
			hops[i].push_back(make_pair(px2, py2));
			downRight.push_back(make_pair(px2, py2));
		}

		sort(upperLeft.begin(), upperLeft.end(), ju);
		sort(downRight.begin(), downRight.end(), juC);

		bbxMin = upperLeft[0].first;
		bbyMax = upperLeft[0].second;
		res.push_back(make_pair(bbxMin, bbyMax));
		bbxMax = downRight[0].first;
		bbyMin = downRight[0].second;
		res.push_back(make_pair(bbxMax, bbyMin));
		return res;
	}

	inline pair<int, int> multiHops_skylineQuery(int p, int q, int C) {
		pair<int, int> temp;
		pair<int, int> candi;
		vector<pair<int, int>> s;
		vector<vector<pair<int, int>>> sky_x, sky_y;
		temp = make_pair(0, 0);
		if (p == q) return temp;
		int x = belong[p], y = belong[q];
		int lca = LCAQuery(x, y);
		if (lca == x || lca == y) {
			if (lca == y) {
				int v = y;
				y = x;
				x = v;
				v = p;
				p = q;
				q = v;
			}
			s = Tree[y].skyToAnc[Tree[x].pos[Tree[x].pos.size() - 1]];
			if (s[s.size() - 1].second > C) {
				return temp;
			}
			else {
				for (int i = 0; i < s.size(); i++) {
					if (s[i].second <= C) {
						temp = s[i];
						break;
					}
				}
				return temp;
			}
		}
		else {
			int w = INF;
			int c = INF;

			pair<int, int> wc;
			sky_x = Tree[x].skyToAnc;
			sky_y = Tree[y].skyToAnc;
			vector<int> hp;
			vector<int> p;
			p.clear();
			p = Tree[lca].pos;
			int ps = Tree[lca].pos.size();
			hp = multiHop(sky_x, sky_y, p);
			for (int i = 0; i < hp.size(); i++) {
				int pos2 = hp[i];
				if (upper[pos2].second <= C) {
					if (upper[pos2].first < w) {
						w = upper[pos2].first;
						c = upper[pos2].second;
					}
				}
				else if (lower[pos2].second > C) {
					continue;
				}
				else {
					candi = concatenation_C(sky_x[p[pos2]], sky_y[p[pos2]], C);
					if (w > candi.first) {
						w = candi.first;
						c = candi.second;
					}
				}
			}
			wc = make_pair(w, c);
			return wc;
		}
	}

	inline pair<int, int> Forest_CSP(const vector<vector<pair<int, int>>>& h1, const vector<vector<pair<int, int>>>& h2, const vector<int>& B, int C) {
		pair<int, int> candi;
		vector<pair<int, int>> s;
		vector<vector<pair<int, int>>> sky_x, sky_y;

		int w = INF;
		int c = INF;

		pair<int, int> wc;
		sky_x = h1;
		sky_y = h2;
		vector<int> hp;
		vector<int> p;
		int ps = B.size();
		for (int i = 0; i < ps; i++) {
			p.push_back(i);
		}
		hp = multiHop(sky_x, sky_y, p);
		for (int i = 0; i < hp.size(); i++) {
			int pos2 = hp[i];
			if (upper[pos2].second <= C) {
				if (upper[pos2].first < w) {
					w = upper[pos2].first;
					c = upper[pos2].second;
				}
			}
			else if (lower[pos2].second > C) {
				continue;
			}
			else {
				candi = concatenation_C(sky_x[p[pos2]], sky_y[p[pos2]], C);
				if (w > candi.first) {
					w = candi.first;
					c = candi.second;
				}
			}
		}
		wc = make_pair(w, c);
		return wc;
	}

	inline pair<int, int> Updated_Forest_CSP(const vector<vector<pair<int, int>>>& h1, const vector<vector<pair<int, int>>>& h2, const vector<int>& B, const vector<int>& B_b, int C) {
		pair<int, int> temp;
		pair<int, int> candi;
		vector<pair<int, int>> s;
		vector<vector<pair<int, int>>> sky_x, sky_y;
		temp = make_pair(0, 0);

		int w = INF;
		int c = INF;

		pair<int, int> wc;
		sky_x = h1;
		sky_y = h2;
		vector<int> hp;
		vector<int> p;
		int ps = B.size();
		for (int i = 0; i < ps; i++) {
			int pos2 = B[B_b[i]];
			int size_x = sky_x[pos2].size();
			int size_y = sky_y[pos2].size();
			int tmp_w_u = sky_x[pos2][0].first + sky_y[pos2][0].first;
			int tmp_c_u = sky_x[pos2][0].second + sky_y[pos2][0].second;
			int tmp_w_l = sky_x[pos2][size_x - 1].first + sky_y[pos2][size_y - 1].first;
			int tmp_c_l = sky_x[pos2][size_x - 1].second + sky_y[pos2][size_y - 1].second;
			if (tmp_c_u <= C) {
				if (tmp_w_u < w) {
					w = tmp_w_u;
					c = tmp_c_u;
				}
			}
			else if (tmp_c_l > C) {
				continue;
			}
			else {
				candi = concatenation_C(sky_x[pos2], sky_y[pos2], C);
				if (w > candi.first) {
					w = candi.first;
					c = candi.second;
				}
			}
		}
		wc = make_pair(w, c);
		return wc;
	}

	inline vector<pair<int, int>> multiHops_skylines(int p, int q) {
		vector<pair<int, int>> temp;
		vector<pair<int, int>> candi;
		vector<pair<int, int>> temp1;
		vector<pair<int, int>> s;
		vector<vector<pair<int, int>>> sky_x, sky_y;
		//temp.push_back(make_pair(0, 0));
		//temp1.push_back(make_pair(INF, INF));
		//if (p == q) return temp;
		int x = belong[p], y = belong[q];
		int lca = LCAQuery(x, y);
		if (lca == x || lca == y) {
			if (lca == y) {
				int v = y;
				y = x;
				x = v;
				v = p;
				p = q;
				q = v;
			}
			/*if (lca == 0) {
				s = Tree[y].skyToAnc[0];
			}*/
			//else {
			s = Tree[y].skyToAnc[Tree[x].pos[Tree[x].pos.size() - 1]];
			//}
			return s;
		}
		else {
			//int w = INF;
			//int c = INF;

			//pair<int, int> wc;
			sky_x = Tree[x].skyToAnc;
			sky_y = Tree[y].skyToAnc;
			vector<int> hp;
			vector<int> _p;
			_p.clear();
			_p = Tree[lca].pos;
			int ps = Tree[lca].pos.size();
			hp = multiHop(sky_x, sky_y, _p);
			for (int i = 0; i < hp.size(); i++) {
				int pos2 = hp[i];
				candi = concatenation(sky_x[_p[pos2]], sky_y[_p[pos2]]);
				//mergeSort(temp1, candi, s);
				//temp1 = screenSkylines(s);
			}
			//return temp1;
			return candi;
		}
	}

	inline vector<pair<int, int>>  Forest_Hops_computing(vector<vector<pair<int, int>>>& h1, vector<vector<pair<int, int>>>& h2, vector<int>& B1) {
		vector<pair<int, int>> temp;
		vector<pair<int, int>> candi;
		vector<pair<int, int>> temp1;
		vector<pair<int, int>> s;
		vector<vector<pair<int, int>>> sky_x, sky_y;
		//temp.push_back(make_pair(0, 0));
		//temp1.push_back(make_pair(INF, INF));

		//int w = INF;
		//int c = INF;

		pair<int, int> wc;
		sky_x = h1;
		sky_y = h2;
		vector<int> hp;
		vector<int> p;
		int ps = B1.size();
		for (int i = 0; i < ps; i++) {
			p.push_back(i);
		}
		hp = multiHop(sky_x, sky_y, p);
		if(hp.size()==0){
			cerr<<"error: no hp"<<endl;
		}
		for (int i = 0; i < hp.size(); i++) {
			int pos2 = hp[i];
			candi = concatenation(sky_x[p[pos2]], sky_y[p[pos2]]);
			//mergeSort(temp1, candi, s);
			//temp1 = screenSkylines(s);
		}

		//return temp1;
		return candi;
	}

	inline vector<int>  Step1_Hops_bound(vector<vector<pair<int, int>>>& h1, vector<vector<pair<int, int>>>& h2, vector<int>& B1) {
		vector<pair<int, int>> temp;
		vector<pair<int, int>> candi;
		vector<pair<int, int>> s;
		vector<vector<pair<int, int>>> sky_x, sky_y;
		//temp.push_back(make_pair(0, 0));


		//int w = INF;
		//int c = INF;

		//pair<int, int> wc;
		sky_x = h1;
		sky_y = h2;
		vector<int> hp;
		vector<int> p;
		int ps = B1.size();
		for (int i = 0; i < ps; i++) {
			p.push_back(i);
		}
		hp = multiHop(sky_x, sky_y, p);
		return hp;
	}


	inline pair<int, int> simple_skylineQuery(int p, int q, int C) {
		pair<int, int> temp;
		pair<int, int> candi;
		vector<pair<int, int>> s;
		vector<vector<pair<int, int>>> sky_x, sky_y;
		temp = make_pair(0, 0);
		if (p == q) return temp;
		int x = belong[p], y = belong[q];
		int lca = LCAQuery(x, y);
		if (lca == x || lca == y) {
			if (lca == y) {
				int v = y;
				y = x;
				x = v;
				v = p;
				p = q;
				q = v;
			}
			s = Tree[y].skyToAnc[Tree[x].pos[Tree[x].pos.size() - 1]];
			if (s[s.size() - 1].second > C) {
				return temp;
			}
			else {
				for (int i = 0; i < s.size(); i++) {
					if (s[i].second <= C) {
						temp = s[i];
						break;
					}
				}
				return temp;
			}
		}
		else {
			int w = INF;
			int c = INF;
			pair<int, int> wc;
			sky_x = Tree[x].skyToAnc;
			sky_y = Tree[y].skyToAnc;
			vector<int> p = Tree[lca].pos;
			int ps = Tree[lca].pos.size();
			vector<pair<int, int>> h;
			for (int i = 0; i < ps; i++) {
				candi = concatenation_C(sky_x[p[i]], sky_y[p[i]], C);
				if (w > candi.first) {
					w = candi.first;
					c = candi.second;
				}

			}
			wc = make_pair(w, c);
			return wc;
		}
	}

	struct judge2 {
		bool operator()(const pair<int, int> a, const pair<int, int> b) {
			return a.first <= b.first;
		};
	};

	inline pair<int, int> base_skylineQuery(int p, int q, int C) {
		pair<int, int> temp;
		pair<int, int> candi;
		vector<pair<int, int>> s;
		vector<vector<pair<int, int>>> sky_x, sky_y;
		judge ju2;
		temp = make_pair(0, 0);
		if (p == q) return temp;
		int x = belong[p], y = belong[q];
		int lca = LCAQuery(x, y);
		if (lca == x || lca == y) {
			if (lca == y) {
				int v = y;
				y = x;
				x = v;
				v = p;
				p = q;
				q = v;
			}
			s = Tree[y].skyToAnc[Tree[x].pos[Tree[x].pos.size() - 1]];
			if (s[s.size() - 1].second > C) {
				return temp;
			}
			else {
				for (int i = 0; i < s.size(); i++) {
					if (s[i].second <= C) {
						temp = s[i];
						break;
					}
				}
				return temp;
			}
		}
		else {
			int w = INF;
			int c = INF;
			pair<int, int> wc;
			sky_x = Tree[x].skyToAnc;
			sky_y = Tree[y].skyToAnc;
			vector<int> p = Tree[lca].pos;
			int ps = Tree[lca].pos.size();
			vector<pair<int, int>> h;
			vector<pair<int, int>> candiSet;
			candiSet.clear();
			for (int i = 0; i < ps; i++) {
				for (int x = 0; x < sky_x[p[i]].size(); x++) {
					for (int y = 0; y < sky_y[p[i]].size(); y++) {
						int w_d = sky_x[p[i]][x].first + sky_y[p[i]][y].first;
						int w_c = sky_x[p[i]][x].second + sky_y[p[i]][y].second;
						if (w_c <= C) {
							candiSet.push_back(make_pair(w_d, w_c));
						}
					}
				}
				sort(candiSet.begin(), candiSet.end(), ju2);

				if (candiSet.size() > 0) {
					candi = candiSet[0];
					if (w > candi.first) {
						w = candi.first;
						c = candi.second;
					}
				}
				candiSet.clear();
			}
			wc = make_pair(w, c);
			return wc;
		}
	}

	inline vector<pair<int, int>> skylines(int p, int q) {
		vector<pair<int, int>> temp, temp1;
		vector<pair<int, int>> candi;
		vector<pair<int, int>> s;
		vector<vector<pair<int, int>>> sky_x, sky_y;
		temp.push_back(make_pair(0, 0));
		temp1.push_back(make_pair(INF, INF));
		if (p == q) return temp;
		int x = belong[p], y = belong[q];
		int lca = LCAQuery(x, y);
		if (lca == x || lca == y) {
			if (lca == y) {
				int v = y;
				y = x;
				x = v;
				v = p;
				p = q;
				q = v;
			}
			s = Tree[y].skyToAnc[Tree[x].pos[Tree[x].pos.size() - 1]];
			return s;
		}
		else {
			int w = INF;
			int c = INF;

			pair<int, int> wc;
			sky_x = Tree[x].skyToAnc;
			sky_y = Tree[y].skyToAnc;
			vector<int> p = Tree[lca].pos;
			int ps = Tree[lca].pos.size();
			vector<pair<int, int>> h;
			for (int i = 0; i < ps; i++) {
				candi = concatenation(sky_x[p[i]], sky_y[p[i]]);
				mergeSort(temp1, candi, s);
				temp1 = screenSkylines(s);
			}
			return temp1;
		}
	}

};

void constructSubTree(int i, vector<Sub_Tree_Decomposition>& TDs)
{
	sm->wait();
	cout << "Sub_Tree: " << i << endl;
	Sub_Tree_Decomposition td;
	string bb = "";
	bb = to_string(static_cast<long long>(i));
	string str = "data/sp.txt" + bb;
	str = str + ".txt";
	cout << str << endl;
	char graphFileName[255] = "data/NY.txt";
	char* data_name = (char*)str.c_str();
	Graph G = Graph(graphFileName, data_name, i, 0);
	td.G = G;
	td.H = td.G;

	td.reduce(i);
	td.makeTree();
	td.makeIndex();
	td.loadIndex(sw_method, i);//input 0 for re-building index, input 10 for re-building index (optimal), input 1 for reading index existed

	td.cntSize();
	TDs[i] = td;
	//cout << VMap.size() << endl;
	sm->notify();
}


int main(int argc, char* argv[])
{
	VMap.resize(All_v);
	vector<Sub_Tree_Decomposition> TDs;
	vector<Bound_Tree_Decomposition> BTDs;
	ifstream queries("USA-NY-CSPQuery_AVG.txt");
	char graphFileName[255] = "data/NY.txt";
#if 1
	Sub_Tree_Decomposition td;
	TDs.assign(fn, td);
	boost::thread_group threads;
	for (int i = 0; i < fn; i++)
	{

		threads.add_thread(new boost::thread(&constructSubTree, i, boost::ref(TDs)));

	}
	threads.join_all();

	//cout << TDs.size() << endl;
#endif
	cout << "Bound_Tree" << endl;
	Bound_Tree_Decomposition btd;
	string str_bound = "data/bound.txt";
	char* data_name_bound = (char*)str_bound.c_str();
	Graph bG = Graph(graphFileName, data_name_bound, fn, 1);
	cout << "Graph Readed" << endl;
	btd.G = bG;
	btd.H = btd.G;
	btd.TDs = TDs;
	auto start = system_clock::now();
	cout << "reduce..." << endl;
	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;
	std::chrono::duration<double> time_span;
	t1 = std::chrono::high_resolution_clock::now();
	btd.reduce();
	cout << "make tree..." << endl;
	btd.makeTree();
	cout << "make index..." << endl;
	btd.makeIndex();
	t2 = std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	cout << "Boundary Tree building: " << time_span.count() << endl;
	btd.loadIndex(sw_method);//input 0 for re-building index, input 10 for re-building index (optimal), input 1 for reading index existed

	btd.cntSize();
	auto end = system_clock::now();
	auto duration = duration_cast<microseconds>(end - start);
	double time = double(duration.count()) * microseconds::period::num / microseconds::period::den;
	cout << "indexing time:" << time << endl;
	cout << "indexing complete!" << endl;
#if 0
	ifstream in(data_name_bound);
	int num;
	vector<int> bound;
	in >> num;
	for (int i = 0; i < num; i++) {
		int bs;
		in >> bs;
		bound.push_back(bs);
	}
	in.close();

	int s, t, C;

	ofstream out;
	out.open("TestQueries2.txt");
	pair<int, int> res;
	ofstream ofsr("result.txt");
	while (queries >> s >> t >> C) {
		//queries >> s >> t >> C;
		cout << "query: " << s << " " << t << " " << C << endl;
		C=INF;
		int iner_s, outer_s = 0;
		int iner_t, outer_t = 0;
		vector<vector<pair<int, int>>> x;
		vector<vector<pair<int, int>>> y;
		vector<int> p;
		bool same_area = false;
		for (int it = 0; it < bound.size(); it++) {
			if (s == bound[it]) {
				outer_s = 1;
			}
			if (t == bound[it]) {
				outer_t = 1;
			}
		}
		if (VMap[s].begin()->first == VMap[t].begin()->first) {
			//cout << "VMap[s].begin()->first:" << VMap[s].begin()->first << " " << "VMap[t].begin()->first:" << VMap[t].begin()->first << endl;
			same_area = true;
		}
		auto start3 = system_clock::now();
		if (outer_s == 1 && outer_t == 1 && !same_area) {
			//cout << "---bound---" << endl;
			int s_id = VMap[s][fn];
			int t_id = VMap[t][fn];
			//cout << "bound-tree search: " << s_id << " " << t_id << endl;
			res = btd.multiHops_skylineQuery(s_id, t_id, C);
			//cout << res.first << '\t' << res.second << endl;
		}

		if (same_area) {
			//cout << "---sub---" << endl;
			int sub_id_s = VMap[s].begin()->first;
			int s_id = VMap[s][sub_id_s];
			int t_id = VMap[t][sub_id_s];
			//cout << "single sub-tree search: " << s_id << " " << t_id << endl;
			res = TDs[sub_id_s].multiHops_skylineQuery(s_id, t_id, C);
			//cout << res.first << '\t' << res.second << endl;
		}
		if (!same_area && outer_s == 1 && outer_t == 0) {

			//cout << "bound-sub" << endl;
			int sub_id_t = VMap[t].begin()->first;
			int t_id = VMap[t][sub_id_t];
			vector<int> B2;
			B2=TDs[sub_id_t].G.bound_ord[sub_id_t];
			vector<vector<pair<int, int>>> t_skylines;
			//t_skylines.resize(B2.size());
			for (int it = 0; it < B2.size(); it++) {
				vector<pair<int, int>> tmp_skylines_s;
				int bound2_id = VMap[B2[it]][sub_id_t];
				tmp_skylines_s = TDs[sub_id_t].multiHops_skylines(t_id, bound2_id);
				t_skylines.push_back(tmp_skylines_s);
			}

			int s_id = VMap[s][fn];
			vector<vector<pair<int, int>>> bound_skylines;
			//bound_skylines.resize(B2.size());
			for (int it2 = 0; it2 < B2.size(); it2++) {
				vector<pair<int, int>> tmp_skylines_b;
				int b2_id = VMap[B2[it2]][fn];
				tmp_skylines_b = btd.multiHops_skylines(b2_id, s_id);
				bound_skylines.push_back(tmp_skylines_b);
			}

			res = btd.Forest_CSP(bound_skylines, t_skylines, B2, C);
			//cout << res.first << '\t' << res.second << endl;
		}

		if (!same_area && outer_s == 0 && outer_t == 1) {
			//cout << "sub-bound" << endl;
			int sub_id_s = VMap[s].begin()->first;
			int s_id = VMap[s][sub_id_s];
			vector<int> B1;
			B1=TDs[sub_id_s].G.bound_ord[sub_id_s];
			vector<vector<pair<int, int>>> s_skylines;
			for (int it = 0; it < B1.size(); it++) {
				vector<pair<int, int>> tmp_skylines_s;
				int bound1_id = VMap[B1[it]][sub_id_s];
				tmp_skylines_s = TDs[sub_id_s].multiHops_skylines(s_id, bound1_id);
				s_skylines.push_back(tmp_skylines_s);
			}


			int t_id = VMap[t][fn];


			vector<vector<pair<int, int>>> bound_skylines;;

			for (int it1 = 0; it1 < B1.size(); it1++) {
				vector<pair<int, int>> tmp_skylines_b;
				int b1_id = VMap[B1[it1]][fn];
				tmp_skylines_b = btd.multiHops_skylines(t_id, b1_id);
				bound_skylines.push_back(tmp_skylines_b);
			}


			res = btd.Forest_CSP(s_skylines, bound_skylines, B1, C);
			//cout << res.first << '\t' << res.second << endl;
		}

 //original 
		if (outer_s == 0 && outer_t == 0 && !same_area) {
			//cout << "sub-sub" << endl;
			int sub_id_s = VMap[s].begin()->first;
			//cout << "sub_id_s:" << sub_id_s << endl;

			int s_id = VMap[s][sub_id_s];
			//cout << "s_id:" << s_id << endl;
			vector<int> B1 = TDs[sub_id_s].G.bound_ord[sub_id_s];
			vector<vector<pair<int, int>>> s_skylines;
			//cout << "B1 size:" << B1.size() << " " << B1[0] << endl;
			for (int it = 0; it < B1.size(); it++) {
				vector<pair<int, int>> tmp_skylines_s;
				int bound1_id = VMap[B1[it]][sub_id_s];
				tmp_skylines_s = TDs[sub_id_s].multiHops_skylines(s_id, bound1_id);
				s_skylines.push_back(tmp_skylines_s);
				//cout << "s_skylines size:" << s_skylines[it].size() << " " << s_skylines[it][s_skylines[it].size() - 1].second << endl;
			}
			//cout << "s_skylines size:" << s_skylines.size() <<" "<<s_skylines[0][0].first<<" "<<s_skylines[0][0].second<<endl;
			int sub_id_t = VMap[t].begin()->first;
			//cout << "sub_id_t:" << sub_id_t << endl;
			int t_id = VMap[t][sub_id_t];
			//cout << "t_id:" << t_id << endl;
			vector<int> B2;
			B2 = TDs[sub_id_t].G.bound_ord[sub_id_t];
			//cout << "B2:" << B2.size() << endl;
			vector<vector<pair<int, int>>> t_skylines;
			for (int it = 0; it < B2.size(); it++) {
				vector<pair<int, int>> tmp_skylines_t;
				int bound2_id = VMap[B2[it]][sub_id_t];
				//cout << "t_id: " << t_id << " " << "bound2_id: " << bound2_id << endl;
				tmp_skylines_t = TDs[sub_id_t].multiHops_skylines(t_id, bound2_id);
				t_skylines.push_back(tmp_skylines_t);
				//cout << "t_skylines size:" << t_skylines.size() << " " << t_skylines[it][t_skylines[it].size() - 1].second << endl;
			}
			//cout << "t_skylines size:" << t_skylines.size() << " " << t_skylines[0][0].first << " " << t_skylines[0][0].second << " " << t_skylines[0][t_skylines[0].size() - 1].first << " " << t_skylines[0][t_skylines[0].size() - 1].second << endl;

			vector<vector<vector<pair<int, int>>>> bound_skylines;
			bound_skylines.resize(B2.size());
			for (int it2 = 0; it2 < B2.size(); it2++) {
				vector<pair<int, int>> tmp_skylines_b;
				int b2_id = VMap[B2[it2]][fn];
				for (int it1 = 0; it1 < B1.size(); it1++) {
					int b1_id = VMap[B1[it1]][fn];
					//cout << b2_id << " " << b1_id << endl;
					//cout << B2[it2] << " " << B1[it1] << endl;
					tmp_skylines_b = btd.multiHops_skylines(b2_id, b1_id);
					bound_skylines[it2].push_back(tmp_skylines_b);
					//cout << "bound_skylines size:" << it2 << " " << bound_skylines[it2].size() << " " << bound_skylines[it2][it1][bound_skylines[it2][it1].size() - 1].second << endl;
				}

			}

			vector<vector<pair<int, int>>> step1_skylines;
			//step1_skylines.resize(B2.size());
			for (int p1 = 0; p1 < B2.size(); p1++) {
				vector<pair<int, int>> tmp_skylines;
				tmp_skylines = btd.Forest_Hops_computing(s_skylines, bound_skylines[p1], B1);
				//cout << s_skylines.size() << " " << bound_skylines[p1].size() << " " << B1.size() << endl;
				step1_skylines.push_back(tmp_skylines);
				//cout << "step1_skylines size:" << step1_skylines.size() << " " << step1_skylines[p1][step1_skylines[p1].size() - 1].second << endl;
			}
			//cout << "step1_skylines size:" << step1_skylines.size() << " " << step1_skylines[0][0].first << " " << step1_skylines[0][0].second << " " << step1_skylines[0][step1_skylines[0].size() - 1].first << " " << step1_skylines[0][step1_skylines[0].size() - 1].second << endl;

			res = btd.Forest_CSP(step1_skylines, t_skylines, B2, C);
			//cout << res.first << '\t' << res.second << endl;
		}
#endif
#if 0//updated
		if (outer_s == 0 && outer_t == 0 && !same_area) {
			//cout << "sub-sub" << endl;
			int sub_id_s = VMap[s].begin()->first;
			//cout << "sub_id_s:" << sub_id_s << endl;

			int s_id = VMap[s][sub_id_s];
			//cout << "s_id:" << s_id << endl;
			vector<int> B1 = TDs[sub_id_s].G.bound_ord[sub_id_s];
			vector<vector<pair<int, int>>> s_skylines;
			//cout << "B1 size:" << B1.size() << " " << B1[0] << endl;
			for (int it = 0; it < B1.size(); it++) {
				vector<pair<int, int>> tmp_skylines_s;
				int bound1_id = VMap[B1[it]][sub_id_s];
				tmp_skylines_s = TDs[sub_id_s].multiHops_skylines(s_id, bound1_id);
				s_skylines.push_back(tmp_skylines_s);
				//cout << "s_skylines size:" << s_skylines[it].size() << " " << s_skylines[it][s_skylines[it].size() - 1].second << endl;
			}
			//cout << "s_skylines size:" << s_skylines.size() <<" "<<s_skylines[0][0].first<<" "<<s_skylines[0][0].second<<endl;
			int sub_id_t = VMap[t].begin()->first;
			//cout << "sub_id_t:" << sub_id_t << endl;
			int t_id = VMap[t][sub_id_t];
			//cout << "t_id:" << t_id << endl;
			vector<int> B2;
			B2 = TDs[sub_id_t].G.bound_ord[sub_id_t];
			//cout << "B2:" << B2.size() << endl;
			vector<vector<pair<int, int>>> t_skylines;
			for (int it = 0; it < B2.size(); it++) {
				vector<pair<int, int>> tmp_skylines_t;
				int bound2_id = VMap[B2[it]][sub_id_t];
				//cout << "t_id: " << t_id << " " << "bound2_id: " << bound2_id << endl;
				tmp_skylines_t = TDs[sub_id_t].multiHops_skylines(t_id, bound2_id);
				t_skylines.push_back(tmp_skylines_t);
				//cout << "t_skylines size:" << t_skylines.size() << " " << t_skylines[it][t_skylines[it].size() - 1].second << endl;
			}
			//cout << "t_skylines size:" << t_skylines.size() << " " << t_skylines[0][0].first << " " << t_skylines[0][0].second << " " << t_skylines[0][t_skylines[0].size() - 1].first << " " << t_skylines[0][t_skylines[0].size() - 1].second << endl;

			vector<vector<vector<pair<int, int>>>> bound_skylines;
			bound_skylines.resize(B2.size());
			for (int it2 = 0; it2 < B2.size(); it2++) {
				vector<pair<int, int>> tmp_skylines_b;
				int b2_id = VMap[B2[it2]][fn];
				for (int it1 = 0; it1 < B1.size(); it1++) {
					int b1_id = VMap[B1[it1]][fn];
					//cout << b2_id << " " << b1_id << endl;
					//cout << B2[it2] << " " << B1[it1] << endl;
					tmp_skylines_b = btd.skylines_bound(b2_id, b1_id);
					bound_skylines[it2].push_back(tmp_skylines_b);
					//cout << "bound_skylines size:" << it2 << " " << bound_skylines[it2].size() << " " << bound_skylines[it2][it1][bound_skylines[it2][it1].size() - 1].second << endl;
				}

			}

			vector<vector<pair<int, int>>> h1_skylines;
			vector<vector<int>> h1_hops;
			//step1_skylines.resize(B2.size());
			for (int p1 = 0; p1 < B2.size(); p1++) {
				vector<pair<int, int>> tmp_skylines;
				vector<int> tmp_hp;
				tmp_hp = btd.Step1_Hops_bound(s_skylines, bound_skylines[p1], B1);
				h1_hops.push_back(tmp_hp);
				tmp_skylines = tmp_rectangle;
				h1_skylines.push_back(tmp_skylines);
				//cout << s_skylines.size() << " " << bound_skylines[p1].size() << " " << B1.size() << endl;
			}
			//cout << "step1_skylines size:" << step1_skylines.size() << " " << step1_skylines[0][0].first << " " << step1_skylines[0][0].second << " " << step1_skylines[0][step1_skylines[0].size() - 1].first << " " << step1_skylines[0][step1_skylines[0].size() - 1].second << endl;
			vector<vector<pair<int, int>>> h2_skylines;
			vector<int> h2_hops;
			h2_hops = btd.Step1_Hops_bound(h1_skylines, t_skylines, B2);
			//cout << res.first << '\t' << res.second << endl;

			//computing paths
			vector<vector<pair<int, int>>> t_skylines_final;
			for (int it = 0; it < h2_hops.size(); it++) {
				vector<pair<int, int>> tmp_skylines_t;
				tmp_skylines_t.clear();
				int bound2_id = VMap[B2[h2_hops[it]]][sub_id_t];
				//cout << "t_id: " << t_id << " " << "bound2_id: " << bound2_id << endl;
				tmp_skylines_t = TDs[sub_id_t].multiHops_skylines(t_id, bound2_id);
				t_skylines_final.push_back(tmp_skylines_t);
				//cout << "t_skylines size:" << t_skylines.size() << " " << t_skylines[it][t_skylines[it].size() - 1].second << endl;
			}
			vector<vector<int>> H1;
			H1.resize(B1.size());
			for (int it2 = 0; it2 < h2_hops.size(); it2++) {
				for (int it1 = 0; it1 < h1_hops[h2_hops[it2]].size(); it1++) {
					vector<int> tmp_hp1;
					tmp_hp1 = h1_hops[h2_hops[it2]];
					H1[tmp_hp1[it1]].push_back(h2_hops[it2]);
				}
			}

			vector<vector<vector<pair<int, int>>>> bound_skylines;
			bound_skylines.resize(B1.size());
			for (int it2 = 0; it2 < B1.size(); it2++) {
				if (H1[it2].size() != 0) {
					vector<pair<int, int>> tmp_skylines_b;
					int b1_id = VMap[B1[it2]][fn];

					for (int it1 = 0; it1 < H1[it2].size(); it1++) {
						int b2_id = VMap[B2[H1[it2][it1]]][fn];
						//cout << b2_id << " " << b1_id << endl;
						//cout << B2[it2] << " " << B1[it1] << endl;
						tmp_skylines_b = btd.multiHops_skylines(b1_id, b2_id);
						bound_skylines[it2].push_back(tmp_skylines_b);
						//cout << "bound_skylines size:" << it2 << " " << bound_skylines[it2].size() << " " << bound_skylines[it2][it1][bound_skylines[it2][it1].size() - 1].second << endl;
					}
				}
			}
			


			vector<vector<pair<int, int>>> step1_skylines;
			step1_skylines.resize(B2.size());
			for (int p1 = 0; p1 < B2.size(); p1++) {
				if (bound_skylines[p1].size() >= 1) {
					vector<pair<int, int>> tmp_skylines;
					vector<pair<int, int>> _s;
					vector<pair<int, int>> _candi;

					tmp_skylines.push_back(make_pair(INF, INF));
					
					
					_candi = btd.concatenation(s_skylines[p1], bound_skylines[p1]);
					//cout << "candi" << candi.first << " " << candi.second << endl;
					btd.mergeSort(tmp_skylines, _candi, _s);
					tmp_skylines = btd.screenSkylines(_s);

					step1_skylines[p1]=tmp_skylines;
				}
				
				//cout << s_skylines.size() << " " << bound_skylines[p1].size() << " " << B1.size() << endl;
				
				//cout << "step1_skylines size:" << step1_skylines.size() << " " << step1_skylines[p1][step1_skylines[p1].size() - 1].second << endl;
			}

			res = btd.Updated_Forest_CSP(step1_skylines, t_skylines, B2, h2_hops, C);
		}
#endif
		//res = td.simple_skylineQuery(s, t, C);
		//res = td.multiHops_skylineQuery(s, t, C);
	
		auto end3 = system_clock::now();
		auto duration3 = duration_cast<microseconds>(end3 - start3);
		double time3 = double(duration3.count()) * microseconds::period::num / microseconds::period::den;
		//ofsr << res.first << '\t' << res.second <<'\t'<<time3<<endl;
		ofsr << time3 << endl;
		
	}
	ofsr.close();
	cout << "indexing time:" << time << endl;
	cout << "Done!" << endl;

#endif
#if 0
	for (int i = 0; i < 500; ++i) {
		queries >> s >> t >> C;
		int iner_s, outer_s = 0;
		int iner_t, outer_t = 0;
		//vector<vector<pair<int, int>>> x;
		//vector<vector<pair<int, int>>> y;
		vector<int> p;
		bool same_area = false;
		for (int it = 0; it < bound.size(); it++) {
			if (s == bound[it]) {
				outer_s = 1;
			}
			if (t == bound[it]) {
				outer_t = 1;
			}
		}
		if (VMap[s].begin()->first == VMap[t].begin()->first) {
			same_area = true;
		}
		auto start3 = system_clock::now();
		if (outer_s == 1 && outer_t == 1) {
			int s_id = VMap[s][fn];
			int t_id = VMap[t][fn];
			res = btd.multiHops_skylineQuery(s_id, t_id, C);
		}

		if (same_area) {
			int sub_id_s = VMap[s].begin()->first;
			int s_id = VMap[s][sub_id_s];
			int t_id = VMap[t][sub_id_s];
			res = TDs[sub_id_s].multiHops_skylineQuery(s_id, t_id, C);
		}
		if (!same_area && outer_s == 1 && outer_t == 0) {

			int sub_id_t = VMap[t].begin()->first;
			int t_id = VMap[t][sub_id_t];
			vector<int> B2 = TDs[sub_id_t].G.bound_ord[sub_id_t];
			vector<vector<pair<int, int>>> t_skylines;
			for (int it = 0; it < B2.size(); it++) {
				vector<pair<int, int>> tmp_skylines;
				int bound2_id = VMap[B2[it]][sub_id_t];
				tmp_skylines = TDs[sub_id_t].multiHops_skylines(t_id, bound2_id);
				t_skylines.push_back(tmp_skylines);
			}

			int s_id = VMap[s][fn];
			vector<vector<pair<int, int>>> bound_skylines;
			bound_skylines.resize(B2.size());
			for (int it2 = 0; it2 < B2.size(); it2++) {
				vector<pair<int, int>> tmp_skylines;
				int b2_id = VMap[B2[it2]][fn];
				tmp_skylines = btd.multiHops_skylines(b2_id, s_id);
				bound_skylines.push_back(tmp_skylines);
			}

			res = btd.Forest_CSP(bound_skylines, t_skylines, B2, C);
		}

		if (!same_area && outer_s == 0 && outer_t == 1) {
			int sub_id_s = VMap[s].begin()->first;
			int s_id = VMap[s][sub_id_s];
			vector<int> B1 = TDs[sub_id_s].G.bound_ord[sub_id_s];
			vector<vector<pair<int, int>>> s_skylines;
			for (int it = 0; it < B1.size(); it++) {
				vector<pair<int, int>> tmp_skylines;
				int bound1_id = VMap[B1[it]][sub_id_s];
				tmp_skylines = TDs[sub_id_s].multiHops_skylines(s_id, bound1_id);
				s_skylines.push_back(tmp_skylines);
			}


			int t_id = VMap[t][fn];


			vector<vector<pair<int, int>>> bound_skylines;;

			for (int it1 = 0; it1 < B1.size(); it1++) {
				vector<pair<int, int>> tmp_skylines;
				int b1_id = VMap[B1[it1]][fn];
				tmp_skylines = btd.multiHops_skylines(t_id, b1_id);
				bound_skylines.push_back(tmp_skylines);
			}


			res = btd.Forest_CSP(s_skylines, bound_skylines, B1, C);
		}


		if (outer_s == 0 && outer_t == 0 && !same_area) {
			int sub_id_s = VMap[s].begin()->first;
			int s_id = VMap[s][sub_id_s];
			vector<int> B1 = TDs[sub_id_s].G.bound_ord[sub_id_s];
			vector<vector<pair<int, int>>> s_skylines;
			for (int it = 0; it < B1.size(); it++) {
				vector<pair<int, int>> tmp_skylines;
				int bound1_id = VMap[B1[it]][sub_id_s];
				tmp_skylines = TDs[sub_id_s].skylines_bound(s_id, bound1_id);
				s_skylines.push_back(tmp_skylines);
			}

			int sub_id_t = VMap[t].begin()->first;
			int t_id = VMap[t][sub_id_t];
			vector<int> B2 = TDs[sub_id_t].G.bound_ord[sub_id_t];
			vector<vector<pair<int, int>>> t_skylines;
			for (int it = 0; it < B2.size(); it++) {
				vector<pair<int, int>> tmp_skylines;
				int bound2_id = VMap[B2[it]][sub_id_t];
				tmp_skylines = TDs[sub_id_t].skylines_bound(t_id, bound2_id);
				t_skylines.push_back(tmp_skylines);
			}

			vector<vector<vector<pair<int, int>>>> bound_skylines;
			bound_skylines.resize(B2.size());
			for (int it2 = 0; it2 < B2.size(); it2++) {
				vector<pair<int, int>> tmp_skylines;
				int b2_id = VMap[B2[it2]][fn];
				for (int it1 = 0; it1 < B1.size(); it1++) {
					int b1_id = VMap[B1[it1]][fn];
					tmp_skylines = btd.skylines_bound(b2_id, b1_id);
					bound_skylines[it2].push_back(tmp_skylines);
				}
			}

			vector<vector<pair<int, int>>> step1_skylines;
			step1_skylines.resize(B2.size());
			for (int p1 = 0; p1 < B2.size(); p1++) {
				step1_skylines[p1] = btd.Forest_Hops_computing(s_skylines, bound_skylines[p1], B1);
			}

			res = btd.Forest_CSP(step1_skylines, t_skylines, B2, C);
		}

		//res = td.simple_skylineQuery(s, t, C);
		//res = td.multiHops_skylineQuery(s, t, C);
		auto end3 = system_clock::now();
		auto duration3 = duration_cast<microseconds>(end3 - start3);
		double time3 = double(duration3.count()) * microseconds::period::num / microseconds::period::den;
		cout << res.first << '\t' << res.second << endl;
		out << "time:" << time3 << endl;
	}
#endif
#if 0
	pair<int, int> res;
	vector<pair<int, int>> res_skylines;
	auto start1 = system_clock::now();
	res = td.simple_skylineQuery(1, 13, 30000);
	//res = td.multiHops_skylineQuery(1, 13, 100);
	//res = td.base_skylineQuery(1, 13, 5);
	auto end1 = system_clock::now();
	auto duration1 = duration_cast<microseconds>(end1 - start1);
	double time1 = double(duration1.count()) * microseconds::period::num / microseconds::period::den;
	cout << res.first << '\t' << res.second << endl;
	cout << "time:" << time1 << endl;

	res_skylines = td.skylines(1, 13);
	for (int i = 0; i < res_skylines.size(); i++) {
		cout << res_skylines[i].first << '\t' << res_skylines[i].second << endl;
	}
#endif
	return 0;
}


