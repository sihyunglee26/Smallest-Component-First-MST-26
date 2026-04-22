/*
Authored by Sihyung Lee
Compile this file with -O3 option (e.g., g++ smallest_component_first.cpp -O3)
*/

#include <iostream>
#include <vector>
#include <numeric>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <queue>
#include <limits>
#include <chrono>
#include <iomanip>
#include <random>
#include <unordered_set>
#include <deque>
#include <stdexcept>
#include <format>
#include <optional>
#include <ostream>
#include <unordered_map>
#include <optional>

using namespace std;
 
struct Edge {
    int u, v;
    double w;
};

// pack two 32-bit integers into a single 64-bit integer, which will serve as a unique key of an undirected edge
static inline uint64_t edge_key_undirected(int a, int b) {
    if (a > b) swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | static_cast<uint32_t>(b);
}

struct Graph {
    int V = 0;
    vector<Edge> edges;
    vector<vector<pair<int, double>>> adj;
};

bool is_graph_connected(Graph g) {    
    vector<bool> visited(g.V, false);
    int num_visited = 0;
    deque<int> q;

    visited[0] = true;
    num_visited++;
    q.push_back(0);

    while(q.size() > 0) {
        int u = q.front();
        q.pop_front();
        for (const auto& [v, w] : g.adj[u]) {
            if (!visited[v]) {
                visited[v] = true;
                num_visited++;
                q.push_back(v);
            }
        }
    }    

    return num_visited >= g.V;
}

Graph build_graph(int V, const vector<Edge>& edges, bool undirected = true) {
    Graph g;
    g.V = V;
    g.edges = edges;
    g.adj.assign(V, {});
    for (const auto& e : edges) {
        g.adj[e.u].push_back({e.v, e.w});
        if (undirected) g.adj[e.v].push_back({e.u, e.w});
    }
    return g;
}

Graph read_graph_from_file(const string& filename, bool undirected = true) {
    Graph g;
    ifstream fin(filename);
    if (!fin) {
        throw runtime_error("Cannot open file: " + filename);
    }

    // Expected format per line: u v w
    // Vertices assumed 0-indexed. If not, preprocess or adapt below.
    vector<Edge> edges;
    int max_vertex = -1;

    string line;
    while (getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        istringstream iss(line);
        int u, v;
        double w;
        if (!(iss >> u >> v >> w)) continue;
        edges.push_back({u, v, w});
        max_vertex = max(max_vertex, max(u, v));
    }

    return build_graph(max_vertex + 1, edges, undirected);
}

// Standard BA generator: start with complete graph of size m+1,
// then each new node attaches to m existing nodes by preferential attachment.
Graph generate_ba_graph(int V, int m, uint64_t seed = 1) {
    if (m < 1 || V <= m + 1) {
        throw runtime_error("Invalid BA parameters");
    }

    mt19937_64 rng(seed);   // creates a high-quality, 64-bit random number generator (which uses Mersenne Twister algorithm)
    uniform_real_distribution<double> wdist(0.0, 1.0);  // create a distribution for real numbers between 0 and 1, used for assigning edge weights

    vector<Edge> edges;
    //unordered_set<uint64_t> used;   // store 64-bit key of already added edges, so that duplicate edges are not added later
    vector<int> degree_pool; // repeated node IDs proportional to degree

    int initial = m + 1;

    // create an initial clique, a complete graph of size m+1
    for (int i = 0; i < initial; ++i) {
        for (int j = i + 1; j < initial; ++j) {
            edges.push_back({i, j, wdist(rng)});
            //used.insert(edge_key_undirected(i, j));
            degree_pool.push_back(i);
            degree_pool.push_back(j);
        }
    }

    for (int new_v = initial; new_v < V; ++new_v) {
        unordered_set<int> chosen;
        uniform_int_distribution<size_t> pick(0, degree_pool.size() - 1);

        // select m previous, distinct vertices with probability proportional to their degrees
        while ((int)chosen.size() < m) {
            int candidate = degree_pool[pick(rng)];
            if (candidate == new_v) continue;
            chosen.insert(candidate);
        }

        for (int old_v : chosen) {
            edges.push_back({new_v, old_v, wdist(rng)});
            //used.insert(edge_key_undirected(new_v, old_v));
            degree_pool.push_back(new_v);
            degree_pool.push_back(old_v);
        }
    }

    return build_graph(V, edges, true);
}

Graph generate_er_graph(int V, double p, uint64_t seed = 1) {
    mt19937_64 rng(seed);   // creates a high-quality, 64-bit random number generator (which uses Mersenne Twister algorithm)
    uniform_real_distribution<double> prob(0.0, 1.0);   // create a distribution for real numbers between 0 and 1, used for adding edges
    uniform_real_distribution<double> wdist(0.0, 1.0);  // create a distribution for real numbers between 0 and 1, used for assigning edge weights

    vector<Edge> edges;
    edges.reserve((size_t)(p * V * (V - 1) / 2));

    for (int i = 0; i < V; ++i) {
        for (int j = i + 1; j < V; ++j) {
            if (prob(rng) < p) {
                edges.push_back({i, j, wdist(rng)});
            }
        }
    }
    return build_graph(V, edges, true);
}

class Stats {
public:
    long long edge_scans = 0;
    long long weight_comparisons = 0;
    long long uf_root_ops = 0;
    long long pq_ops = 0; // count only O(log N)-style ops
    long long cache_check_for_worst = 0;
    long long cache_check_for_best = 0;
    long long component_udates = 0;
    double time = 0.0;
    double rank = 0;    

    void add(const Stats& stats) {
        edge_scans += stats.edge_scans;
        weight_comparisons += stats.weight_comparisons;
        uf_root_ops += stats.uf_root_ops;
        pq_ops += stats.pq_ops;
        cache_check_for_best += stats.cache_check_for_best;
        cache_check_for_worst += stats.cache_check_for_worst;        
        component_udates += stats.component_udates;
        time += stats.time;
        rank += stats.rank;        
    }

    void divide(int length) {
        edge_scans /= length;
        weight_comparisons /= length;
        uf_root_ops /= length;
        pq_ops /= length;
        cache_check_for_best /= length;
        cache_check_for_worst /= length;
        component_udates /= length;
        time /= length;
        rank /= length;        
    }
};

struct MSTResult {
    string algo_name;
    double total_weight = 0.0;
    vector<Edge> mst_edges;
    Stats stats;
    int V, E;
    string graph_type;
    double graph_param;
};

class UnionFind {
public:
    explicit UnionFind(int n, Stats* stats = nullptr)   // "explicit" to allow only explicit use of the constructor (e.g., UnionFind uf(10)) but to forbid type conversion (e.g., UnionFind uf = 10;)
        : parent(n), sz(n, 1), stats(stats) {           // initializes class members, which is faster than doing it in the body 
        iota(parent.begin(), parent.end(), 0);          // fill the entire parent vector, starting from value 0 with a +1 step
    }

    int root(int x) {
        if (stats) stats->uf_root_ops++;
        while (x != parent[x]) {
            parent[x] = parent[parent[x]];              // path compression
            x = parent[x];
        }
        return x;
    }

    bool unite(int a, int b) {
        int ra = root(a);
        int rb = root(b);
        if (ra == rb) return false;
        if (sz[ra] < sz[rb]) swap(ra, rb);
        parent[rb] = ra;
        sz[ra] += sz[rb];
        return true;
    }

    bool connected(int a, int b) {
        return root(a) == root(b);
    }

    int size_of_root(int r) const {                     // why do we need this?
        return sz[r];
    }

private:
    vector<int> parent;
    vector<int> sz;
    Stats* stats;
};

template <typename T>
class IndexedMinHeap {
public:
    explicit IndexedMinHeap(int N, Stats* s = nullptr, bool smallest_component_first = false)   // "explicit" to allow only explicit use of the constructor (e.g., IndexedMinHeap pq(10)) but to forbid type conversion (e.g., IndexedMinHeap pq = 10;)
        : max_n(N), heap(N+1, -1), pos(N, -1), keys(N+1), stats(s) {
        if (smallest_component_first) {
            n = N;
            iota(heap.begin(), heap.end(), -1); // fill the entire parent vector, starting from value -1 with a +1 step
            iota(pos.begin(), pos.end(), 1);     // fill the entire parent vector, starting from value 1 with a +1 step
            //fill(keys.begin(), keys.end(), 1);  // all keys (component size) are initially 1
        } else {
            n = 0;
            fill(keys.begin(), keys.end(), T());    
        }
    }

    bool empty() const { return n == 0; }   // "const" after method name means that this method doesn't change any member variables
    bool contains(int id) const { return pos[id] != -1; }
    int size() const { return n; }
    bool greater(int i, int j) const { return keys[heap[i]] > keys[heap[j]]; }
    void exch(int i, int j) {
        swap(heap[i], heap[j]);
        pos[heap[i]] = i; pos[heap[j]] = j;        
    }
    void sift_up(int i) {
        int i_parent;
        while (i > 1) {            
            i_parent = i / 2;
            if (!greater(i_parent, i)) break;
            exch(i, i_parent);
            i = i_parent;
        }
    }
    void sift_down(int i) {
        int i_child;        
        while (true) {
            i_child = 2 * i;    // left child
            if (i_child > n) break; // if a child does not exist, break
            if (i_child < n && greater(i_child, i_child + 1)) { i_child++; }    // find a smaller child
            if (!greater(i, i_child)) break;
            exch(i, i_child);
            i = i_child;
        }
    }
    void insert(int id, T key) {
        if (contains(id)) { throw runtime_error("index " + to_string(id) + " is already in PQ"); }
        n++;
        pos[id] = n;
        heap[n] = id;
        keys[id] = key;
        if (stats) stats->pq_ops++;
        sift_up(n);        
    }
    int min_index() const {
        if (n == 0) { throw runtime_error("PQ has no element, so no min index exists"); }
        return heap[1];
    }
    T min_key() const {
        if (n == 0) { throw runtime_error("PQ has no element, so no min key exists"); }
        return keys[heap[1]];
    }
    pair<int, T> min() const {
        if (n == 0) { throw runtime_error("PQ has no element, so no min element exists"); }
        return {heap[1], keys[heap[1]]};
    }
    pair<int, T> del_min() {
        if (n == 0) { throw runtime_error("PQ has no element, so no element to delete"); }
        int min_id = heap[1];
        T min_key = keys[min_id];
        exch(1, n);
        n--;
        if (stats) stats->pq_ops++;
        sift_down(1);
        pos[min_id] = -1;        
        return {min_id, min_key};
    }
    T key_of(int id) const {
        if (!contains(id)) { throw runtime_error("index " + to_string(id) + " is not in PQ"); }
        return keys[id];
    }
    void update_key(int id, T new_key) {
        if (!contains(id)) { throw runtime_error("index " + to_string(id) + " is not in PQ"); }
        keys[id] = new_key;
        if (stats) stats->pq_ops++;
        sift_up(pos[id]);
        if (stats) stats->pq_ops++;
        sift_down(pos[id]);
    }
    void decrease_key(int id, T new_key) {
        if (!contains(id)) { throw runtime_error("index " + to_string(id) + " is not in PQ"); }
        if (keys[id] <= new_key) { throw runtime_error("calling decrease_key with a new key >= the previous key"); }
        keys[id] = new_key;
        if (stats) stats->pq_ops++;
        sift_up(pos[id]);
    }
    void increase_key(int id, T new_key) {
        if (!contains(id)) { throw runtime_error("index " + to_string(id) + " is not in PQ"); }
        if (keys[id] >= new_key) { throw runtime_error("calling increase_key with a new key <= the previous key"); }
        keys[id] = new_key;
        if (stats) stats->pq_ops++;
        sift_down(pos[id]);
    }
    void del(int id) {
        if (!contains(id)) { throw runtime_error("index " + to_string(id) + " is not in PQ"); }
        int i = pos[id];
        exch(i, n);
        n--;
        if (stats) stats->pq_ops++;
        sift_up(i);
        if (stats) stats->pq_ops++;
        sift_down(i);
        pos[id] = -1;        
    }

    vector<T> keys;      // key[id]: key of index id

private:
    int n, max_n;
    vector<int> heap;   // heap[i]: the id of element at heap position i
    vector<int> pos;    // pos[id]: the heap position i of the element with index id      
    Stats* stats;
};

struct Cmp {
    long long* counter;

    bool operator()(const Edge& a, const Edge& b) const {
        (*counter)++;
        return a.w < b.w;
    }
};

MSTResult mst_kruskal_with_sort(const Graph& g) {
    MSTResult res;
    auto edges = g.edges;   // creates a complete copy (instead, 'auto&' would make an alias)
    sort(edges.begin(), edges.end(), Cmp{&res.stats.weight_comparisons});
    /*sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {                
        return a.w < b.w;
    });*/

    UnionFind uf(g.V, &res.stats);

    for (const auto& e : edges) {
        res.stats.edge_scans++;
        if (uf.unite(e.u, e.v)) {
            res.total_weight += e.w;
            res.mst_edges.push_back(e);
            if ((int)res.mst_edges.size() == g.V - 1) break;
        }
    }
    return res;
}

MSTResult mst_prim_eager_with_indexed_heap(const Graph& g) {
    MSTResult res;
    const double INF = numeric_limits<double>::infinity();

    vector<bool> in_mst(g.V, false);

    using Node = pair<double, int>;     // {weight, from_vertex}
    IndexedMinHeap<Node> pq(g.V, &res.stats, false);
    
    pq.insert(0, {0.0, 0});    
    //res.stats.pq_ops++;
    
    while (!pq.empty() && (int)res.mst_edges.size() < g.V - 1) {
        auto [min_id, min_pair] = pq.del_min();
        auto [min_weight, from_id] = min_pair;
        //res.stats.pq_ops++;

        in_mst[min_id] = true;
        if (min_id != from_id) {
            res.mst_edges.push_back({from_id, min_id, min_weight});
            res.total_weight += min_weight;
        }

        for (const auto& [v, w] : g.adj[min_id]) {
            res.stats.edge_scans++;
            if (!in_mst[v]) {          
                if (!pq.contains(v)) {
                    pq.insert(v, {w, min_id});
                } else {
                    res.stats.weight_comparisons++;
                    auto [w_prev, from_id_prev] = pq.key_of(v);
                    if (w < w_prev) {
                        pq.decrease_key(v, {w, min_id});
                    }
                }                   
            }
        }
    }

    return res;
}

MSTResult mst_boruvka(const Graph& g) {
    MSTResult res;
    UnionFind uf(g.V, &res.stats);
    int num_components = g.V;

    while (num_components > 1) {
        vector<int> best_edge_idx(g.V, -1);

        for (int i = 0; i < (int)g.edges.size(); ++i) {
            const auto& e = g.edges[i];
            res.stats.edge_scans++;

            int ru = uf.root(e.u);
            int rv = uf.root(e.v);
            if (ru == rv) continue;

            if (best_edge_idx[ru] == -1) {
                best_edge_idx[ru] = i;
            } else {
                res.stats.weight_comparisons++;
                if (e.w < g.edges[best_edge_idx[ru]].w) {
                    best_edge_idx[ru] = i;
                }
            }

            if (best_edge_idx[rv] == -1) {
                best_edge_idx[rv] = i;
            } else {
                res.stats.weight_comparisons++;
                if (e.w < g.edges[best_edge_idx[rv]].w) {
                    best_edge_idx[rv] = i;
                }
            }
        }

        bool merged = false;
        for (int r = 0; r < g.V; ++r) {
            int idx = best_edge_idx[r];
            if (idx == -1) continue;
            const auto& e = g.edges[idx];
            if (uf.unite(e.u, e.v)) {
                res.mst_edges.push_back(e);
                res.total_weight += e.w;
                num_components--;
                merged = true;
            }
        }

        if (!merged) break; // disconnected graph guard
    }

    return res;
}

const int K = 2;  // cache top-k best outgoing edges per component
template<int K>
struct VertexCache {
    int size = 0;
    int to[K];
    double w[K];
};

inline void consider_for_topk_cached_edge(VertexCache<K>& cache, int to, double w, Stats& stats) {
    if (cache.size < K) {
        cache.to[cache.size] = to; cache.w[cache.size] = w;
        cache.size++;
    } else {
        // find the greatest existing entry
        int worst_i = 0;
        double worst_w = cache.w[0];
        for(int i = 1; i < cache.size; i++) {
            //stats.edge_scans++;
            stats.cache_check_for_worst++;                      
            stats.weight_comparisons++;
            if (cache.w[i] > worst_w) {
                worst_w = cache.w[i]; worst_i = i;
            }
        }
        stats.weight_comparisons++;
        if (w < worst_w) { // replace the greatest existing entry with the new one if the new one is better
            cache.to[worst_i] = to; cache.w[worst_i] = w;
        }
    }
}

inline pair<int, double> get_valid_cached_edge(const VertexCache<K>& cache, int component_id, const vector<int>& vertex_to_component, Stats& stats) {
    int best_to = -1;
    double best_w = numeric_limits<double>::infinity();
    
    for(int i = 0; i < cache.size; i++) {
        //stats.edge_scans++;
        stats.cache_check_for_best++;
        if (vertex_to_component[cache.to[i]] == component_id) continue;
        stats.weight_comparisons++;
        if (cache.w[i] < best_w) {
            best_w = cache.w[i]; best_to = cache.to[i];
        }
    }
    return {best_to, best_w};    
}

MSTResult mst_smallest_component_first_with_best_k_edge_caching_per_vertex(const Graph& g) {
    MSTResult res;
    
    vector<vector<int>> component_vertices(g.V);    // component_vertices[c] = list of vertices currently in component c
    vector<int> vertex_to_component(g.V);    
    vector<VertexCache<K>> vertex_edge_cache(g.V);    

    for (int i = 0; i < g.V; ++i) {
        component_vertices[i].push_back(i);
        vertex_to_component[i] = i;        
    }
    
    IndexedMinHeap<int> pq(g.V, &res.stats, true);
    fill(pq.keys.begin(), pq.keys.end(), 1);  // all keys (component size) are initially 1
    
    int mst_count = 0;

    while (mst_count < g.V - 1) {
        auto [min_id, min_size] = pq.del_min();        
        
        double best_w = numeric_limits<double>::infinity();
        int best_u = -1, best_v = -1, other_id = -1;

        for (int u : component_vertices[min_id]) {            
            auto [best_to_u, best_w_u] = get_valid_cached_edge(vertex_edge_cache[u], min_id, vertex_to_component, res.stats);
            
            if (best_to_u != -1) {                
                // valid cached edge found and thus use it
                if (best_w_u < best_w) {                    
                    best_w = best_w_u; best_u = u; best_v = best_to_u; 
                    other_id = vertex_to_component[best_v];
                }
            } else {                
                // do full rescan of outgoing edges and rebuild top-k cache
                vertex_edge_cache[u].size = 0;
                for (const auto& [v, w] : g.adj[u]) {
                    res.stats.edge_scans++;                    
                    if (vertex_to_component[v] == min_id) continue;      // skip an internal edge
                    consider_for_topk_cached_edge(vertex_edge_cache[u], v, w, res.stats);                    
                }                
                if (vertex_edge_cache[u].size > 0) {                    
                    auto [best_to_u, best_w_u] = get_valid_cached_edge(vertex_edge_cache[u], min_id, vertex_to_component, res.stats);
                    if (best_w_u < best_w) {                        
                        best_w = best_w_u; best_u = u; best_v = best_to_u;                         
                        other_id = vertex_to_component[best_v];                        
                    }
                }                
            }
        }          

        if (other_id == -1) break; // disconnected graph guard

        res.mst_edges.push_back({best_u, best_v, best_w});
        res.total_weight += best_w;
        mst_count++;
        
        // Merge min_id into other_id
        pq.increase_key(other_id, pq.key_of(other_id) + min_size);          
        for (int x : component_vertices[min_id]) {
            res.stats.component_udates++;
            vertex_to_component[x] = other_id;            
            component_vertices[other_id].push_back(x);
        }
        component_vertices[min_id].clear();
    }

    return res;
}

template <typename Func>        // "Func" will be some type, and the compiler will figure it out later
void benchmark(Func f, const Graph& g, int runs, MSTResult& last_result) {
    double total = 0.0;
    for (int i = 0; i < runs; ++i) {
        auto start = chrono::high_resolution_clock::now();
        last_result = f(g);        
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> diff = end - start;        // in seconds
        total += diff.count();
    }    
    last_result.stats.time = total / runs;    
}

void print_result(vector<ostream*>& outs, vector<MSTResult>& v_r_sum) {
    sort(v_r_sum.begin(), v_r_sum.end(), [](const MSTResult& a, const MSTResult& b) {
              return a.stats.time < b.stats.time;
          });   // sort in inreasing order of stats.time
    int rank = 1;
    for (MSTResult& r_sum : v_r_sum) { r_sum.stats.rank = rank++; }

    for (auto out: outs) {
        for (const MSTResult& r : v_r_sum) {
            (*out) << fixed << setprecision(9);       // print up to 9 digits below 1
            (*out) << r.algo_name << ": "
                << r.graph_type << " graph with V=" << r.V << ", E=" << r.E << ", "
                << r.stats.time << " secs, "
                //<< (int) r.stats.rank << " rank, "
                << r.stats.edge_scans << " edge scans, "
                << r.stats.weight_comparisons << " w comps, "
                << r.stats.uf_root_ops << " uf ops, "
                << r.stats.pq_ops << " pq ops, "
                << r.stats.component_udates << " component updates, "
                << r.stats.cache_check_for_best << "/" << r.stats.cache_check_for_worst << " cache checks for best/worst, "
                << "with MST weight " << r.total_weight << "\n"; 
        }
        (*out) << endl;
    }   
}

void print_result_summary(vector<ostream*>& outs, const unordered_map<string, vector<MSTResult>>& results) {    
    vector<pair<string, Stats>> final_result;
    for (auto &[name, vec] : results) {        
        Stats grand_sum;
        for (const MSTResult& mstr : vec) { grand_sum.add(mstr.stats); }
        grand_sum.divide(vec.size());
        final_result.push_back({name, grand_sum});        
    }

    sort(final_result.begin(), final_result.end(), [](const pair<string, Stats>& a, const pair<string, Stats>& b) {
              return a.second.time < b.second.time;
          });   // sort in inreasing order of time

    for (auto out: outs) {
        (*out) << "average:" << endl;
        for (auto &[name, stats] : final_result) {
            (*out) << name << ": "
                << fixed << setprecision(9)
                << stats.time << " secs, "
                << fixed << setprecision(3)
                << stats.rank << " rank, "
                << stats.edge_scans << " edge scans, "
                << stats.weight_comparisons << " w comps, "
                << stats.uf_root_ops << " uf ops, "
                << stats.pq_ops << " pq ops, "        
                << stats.component_udates << " component updates, "        
                << stats.cache_check_for_best << "/" << stats.cache_check_for_worst << " cache checks for best/worst, "
                << endl;
        }
    }
}

struct GraphType {
    string type;
    int v;
    double param;
};

int main() {
    int v = 1000;
    // Note that ER graphs can become disconnected if p is too small, thus leading to a substantially long loops until a connected graph is generated. 
    vector<GraphType> g_type_vector_er = {
        {"er", v, (double)5*2/(v-1)}, {"er", v, (double)10*2/(v-1)}, {"er", v, (double)20*2/(v-1)}, 
        {"er", v, (double)30*2/(v-1)}, {"er", v, (double)40*2/(v-1)}, {"er", v, (double)50*2/(v-1)}, {"er", v, (double)60*2/(v-1)}, 
        {"er", v, (double)80*2/(v-1)}, {"er", v, (double)100*2/(v-1)}, {"er", v, (double)120*2/(v-1)}, {"er", v, (double)140*2/(v-1)}, 
        {"er", v, (double)160*2/(v-1)}, {"er", v, (double)180*2/(v-1)}, {"er", v, (double)200*2/(v-1)}, {"er", v, (double)220*2/(v-1)}
    };
    
    vector<GraphType> g_type_vector_ba = {
        {"ba", v, 2},  {"ba", v, 5}, {"ba", v, 10}, {"ba", v, 20}, {"ba", v, 30}, {"ba", v, 40}, {"ba", v, 50}, 
        {"ba", v, 60}, {"ba", v, 80}, {"ba", v, 100}, {"ba", v, 120}, {"ba", v, 140}, {"ba", v, 160}, 
        {"ba", v, 180}, {"ba", v, 200}, {"ba", v, 220}
    };    

    vector<GraphType> g_type_vector;
    g_type_vector.insert(g_type_vector.end(), g_type_vector_er.begin(), g_type_vector_er.end());
    g_type_vector.insert(g_type_vector.end(), g_type_vector_ba.begin(), g_type_vector_ba.end());

    int num_graphs = 10; // number of different graphs generated for the same parameter set
    int runs = 5;   // number of runs for the same graph (to measure the avg. performance)
    using AlgoFunc = MSTResult (*)(const Graph&);
    vector<pair<string, AlgoFunc>> algos = {
        {"Kruskal", mst_kruskal_with_sort}, {" Prim  ", mst_prim_eager_with_indexed_heap}, {"Boruvka", mst_boruvka},         
        {"SCFcach(K=" + to_string(K) + ")", mst_smallest_component_first_with_best_k_edge_caching_per_vertex}        
    };    

    unordered_map<string, vector<MSTResult>> results;
    for (auto &[name, algo] : algos) { results[name] = {}; } 

    ofstream fout("smallest_component_first.log"); // open a log file
    vector<ostream*> outs = {&cout, &fout};
   
    for (const auto& [t, num_v, param] : g_type_vector) {
        vector<Graph> g_vector;
        if (t == "er" || t == "ba") {
            int seed = 1;
            for(int i = 0; i < num_graphs; i++) {
                Graph g;
                while (true) {
                    if (t == "er") { g = generate_er_graph(num_v, param, seed); }
                    else if (t == "ba") { g = generate_ba_graph(num_v, param, seed); }

                    if (is_graph_connected(g)) { g_vector.push_back(g); break; } 
                    seed++;
                }
            }
        } 

        vector<MSTResult> v_r_sum;
        for (auto &[name, algo] : algos) {
            double t_sum = 0;
            MSTResult r_sum;
            for (Graph g : g_vector) {
                benchmark(algo, g, runs, r_sum);
                t_sum += r_sum.stats.time;
            }
            r_sum.algo_name = name; 
            r_sum.V = g_vector[0].V; r_sum.E = g_vector[0].edges.size(); 
            r_sum.graph_type = t; r_sum.graph_param = param;
            r_sum.stats.time = t_sum / g_vector.size();  
            v_r_sum.push_back(r_sum);                       
        }
        
        print_result(outs, v_r_sum);
        for (MSTResult& r_sum : v_r_sum) { results[r_sum.algo_name].push_back(r_sum); }
    }
    print_result_summary(outs, results);    

    return 0;
}