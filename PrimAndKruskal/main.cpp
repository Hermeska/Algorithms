// 84428810
#include <algorithm>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
using std::vector;
const int kInf = 1000005;

template <typename Vertextype = int,
          typename Edge = std::pair<Vertextype, Vertextype>,
          typename EdgeWeight = double>
class Graph {
 public:
  using EdgeType = std::pair<Edge, EdgeWeight>;
  using VertexEdge = Vertextype;

  virtual void AddVertex(Vertextype vert) = 0;

  virtual void AddEdge(const Edge& edge, const EdgeWeight& weight) = 0;

  virtual size_t VertexCount() const = 0;

  virtual size_t EdgeCount() const = 0;

  virtual vector<std::pair<Vertextype, EdgeWeight>> Neighbours(
      Vertextype ver) const = 0;

  virtual typename std::vector<
      typename std::pair<Vertextype, EdgeWeight>>::const_iterator
  BeginNeighbours(Vertextype ver) const = 0;

  virtual typename std::unordered_map<
      Vertextype,
      std::vector<std::pair<Vertextype, EdgeWeight>>>::const_iterator
  Begin() const = 0;

  virtual typename std::unordered_map<
      Vertextype,
      std::vector<std::pair<Vertextype, EdgeWeight>>>::const_iterator
  End() const = 0;

  virtual typename std::vector<
      typename std::pair<Vertextype, EdgeWeight>>::const_iterator
  EndNeighbours(Vertextype ver) const = 0;

  virtual size_t VerCount() const = 0;
};
template <typename Vertextype = int,
          typename Edge = std::pair<Vertextype, Vertextype>,
          typename EdgeWeight = double>
class Tree;

template <typename Vertextype = int,
          typename Edge = std::pair<Vertextype, Vertextype>,
          typename EdgeWeight = double>
class ListGraph : public Graph<Vertextype, Edge, EdgeWeight> {
  using EdgeType = typename Graph<Vertextype, Edge, EdgeWeight>::EdgeType;

 public:
  void AddVertex(Vertextype ver) override {
    if (adj_.find(ver) == adj_.end()) {
      adj_[ver] = vector<std::pair<Vertextype, EdgeWeight>>();
      vertexes_.push_back(ver);
    }
  }

  void AddEdge(const Edge& edge, const EdgeWeight& weight) override {
    AddVertex(edge.first);
    AddVertex(edge.second);
    Vertextype v_vert = edge.first;
    Vertextype u_vert = edge.second;
    adj_.at(v_vert).push_back({u_vert, weight});
    adj_.at(u_vert).push_back({v_vert, weight});
  }

  size_t VertexCount() const override { return adj_.size(); }

  size_t EdgeCount() const override {
    size_t cnt = 0;
    for (const auto& vert : adj_) {
      cnt += vert.second.size();
    }
    return cnt;
  }
  vector<Vertextype> GetVertexes() { return vertexes_; }
  vector<std::pair<Vertextype, EdgeWeight>> Neighbours(
      Vertextype ver) const override {
    return adj_.at(ver);
  }
  std::unordered_map<Vertextype, vector<std::pair<Vertextype, EdgeWeight>>>
  GetGraph() {
    return adj_;
  }
  typename std::vector<
      typename std::pair<Vertextype, EdgeWeight>>::const_iterator
  BeginNeighbours(Vertextype ver) const override {
    return adj_.at(ver).begin();
  }

  typename std::vector<
      typename std::pair<Vertextype, EdgeWeight>>::const_iterator
  EndNeighbours(Vertextype ver) const override {
    return adj_.at(ver).end();
  }

  size_t VerCount() const override { return adj_.size(); }

  typename std::unordered_map<
      Vertextype,
      std::vector<std::pair<Vertextype, EdgeWeight>>>::const_iterator
  Begin() const override {
    return adj_.begin();
  }

  typename std::unordered_map<
      Vertextype,
      std::vector<std::pair<Vertextype, EdgeWeight>>>::const_iterator
  End() const override {
    return adj_.end();
  };

 private:
  std::unordered_map<Vertextype, vector<std::pair<Vertextype, EdgeWeight>>>
      adj_;
  std::vector<Vertextype> vertexes_;
};
template <typename Vertextype = int,
          typename Edge = std::pair<Vertextype, Vertextype>,
          typename EdgeWeight = double>
class MatrixGraph : public Graph<Vertextype, Edge, EdgeWeight> {
 public:
  using EdgeType = std::pair<Edge, EdgeWeight>;
  using VertexEdge = Vertextype;

  void AddVertex(Vertextype vert) override { adj_list_[vert]; }

  void AddEdge(const Edge& edge, const EdgeWeight& weight) override {
    adj_list_[edge.first][edge.second] = weight;
  }

  size_t VertexCount() const override { return adj_list_.size(); }

  size_t EdgeCount() const override {
    size_t count = 0;
    for (const auto& pair : adj_list_) {
      count += pair.second.size();
    }
    return count;
  }

  vector<std::pair<Vertextype, EdgeWeight>> Neighbours(
      Vertextype vert) const override {
    vector<std::pair<Vertextype, EdgeWeight>> neighbors;
    auto iter = adj_list_.find(vert);
    if (iter != adj_list_.end()) {
      for (const auto& pair : iter->second) {
        neighbors.push_back(pair);
      }
    }
    return neighbors;
  }

  typename std::unordered_map<
      Vertextype, std::unordered_map<Vertextype, EdgeWeight>>::const_iterator
  Begin() const override {
    return adj_list_.cbegin();
  }

  typename std::unordered_map<
      Vertextype, std::unordered_map<Vertextype, EdgeWeight>>::const_iterator
  End() const override {
    return adj_list_.cend();
  }

  typename std::vector<
      typename std::pair<Vertextype, EdgeWeight>>::const_iterator
  BeginNeighbours(Vertextype vert) const override {
    auto iter = adj_list_.find(vert);
    if (iter != adj_list_.end()) {
      return iter->second.cbegin();
    }
    return std::vector<std::pair<Vertextype, EdgeWeight>>::const_iterator();
  }

  typename std::vector<
      typename std::pair<Vertextype, EdgeWeight>>::const_iterator
  EndNeighbours(Vertextype vert) const override {
    auto iter = adj_list_.find(vert);
    if (iter != adj_list_.end()) {
      return iter->second.cend();
    }
  }

 private:
  std::unordered_map<Vertextype, std::unordered_map<Vertextype, EdgeWeight>>
      adj_list_;
};

template <typename Type>
class DisjointSet {
  std::unordered_map<Type, Type> parent_;
  std::unordered_map<Type, int> rank_;

 public:
  DisjointSet() {}
  void MakeSet(Type element) {
    parent_[element] = element;
    rank_[element] = 0;
  }
  Type Find(Type elemnt) {
    if (parent_[elemnt] != elemnt) {
      parent_[elemnt] = Find(parent_[elemnt]);
    }
    return parent_[elemnt];
  }

  void Merge(Type first_elem, Type second_elem) {
    Type first_elem_parent = Find(first_elem);
    Type seconde_elem_parent = Find(second_elem);
    if (first_elem_parent == seconde_elem_parent) {
      return;
    }
    if (rank_[first_elem_parent] < rank_[seconde_elem_parent]) {
      parent_[first_elem_parent] = seconde_elem_parent;
    } else if (rank_[first_elem_parent] > rank_[seconde_elem_parent]) {
      parent_[seconde_elem_parent] = first_elem_parent;
    } else {
      parent_[seconde_elem_parent] = first_elem_parent;
      rank_[first_elem_parent]++;
    }
  }
};

template <typename Vertextype, typename Edge, typename EdgeWeight>
class Tree : public ListGraph<Vertextype, Edge, EdgeWeight> {
 public:
  explicit Tree()
      : ListGraph<Vertextype, Edge, EdgeWeight>(), global_weight(0) {}
  int global_weight;

 private:
};

template <typename Vertextype = int,
          typename Edge = std::pair<Vertextype, Vertextype>,
          typename EdgeWeight = double>
Tree<Vertextype, Edge, EdgeWeight> ElPrimo(
    ListGraph<Vertextype, Edge, EdgeWeight>& graph, Vertextype start) {
  std::priority_queue<std::pair<EdgeWeight, Vertextype>,
                      vector<std::pair<EdgeWeight, Vertextype>>, std::greater<>>

      prior_queue;
  std::unordered_map<Vertextype, EdgeWeight> key;
  std::unordered_map<Vertextype, Vertextype> parents;
  std::unordered_set<Vertextype> visited;
  prior_queue.emplace(0, start);
  for (auto vertex : graph.GetVertexes()) {
    key[vertex] = kInf;
  }
  key[start] = EdgeWeight();
  while (!prior_queue.empty()) {
    auto [key_from_queue, vertex] = prior_queue.top();
    prior_queue.pop();
    visited.insert(vertex);
    for (auto& [neighbor, weight] : graph.Neighbours(vertex)) {
      if (!visited.contains(neighbor) && weight < key[neighbor]) {
        key[neighbor] = weight;
        parents[neighbor] = vertex;
        prior_queue.push({key[neighbor], neighbor});
      }
    }
  }

  Tree<Vertextype, Edge, EdgeWeight> mst;
  for (auto& vertex : graph.GetVertexes()) {
    Vertextype elem_of_parent = parents[vertex];
    EdgeWeight weight = key[vertex];
    mst.global_weight += weight;
    mst.AddEdge({elem_of_parent, vertex}, weight);
  }
  return mst;
}

template <typename Vertextype = int,
          typename Edge = std::pair<Vertextype, Vertextype>,
          typename EdgeWeight = double>
Tree<Vertextype, Edge, EdgeWeight> Kruskal(
    ListGraph<Vertextype, Edge, EdgeWeight>& graph) {
  using EdgeType = std::pair<Edge, EdgeWeight>;
  Tree<Vertextype, Edge, EdgeWeight> mst;
  DisjointSet<Vertextype> dsu;
  for (const auto& vertex : graph.GetVertexes()) {
    dsu.MakeSet(vertex);
  }
  auto copygraph = graph.GetGraph();
  std::vector<EdgeType> edges;
  for (const auto& vertex : graph) {
    for (const auto& neighbor : vertex.second) {
      edges.emplace_back({vertex.first, neighbor.first}, neighbor.second);
    }
  }
  std::sort(edges.begin(), edges.end(),
            [](const EdgeType& first_edge, const EdgeType& second_edge) {
              return first_edge.second < second_edge.second;
            });
  for (const auto& edge : edges) {
    Vertextype source = edge.first.first;
    Vertextype destination = edge.first.second;
    EdgeWeight weight = edge.second;
    if (dsu.Find(source) != dsu.Find(destination)) {
      dsu.Merge(source, destination);
      mst.AddEdge({source, destination}, weight);
      mst.global_weight += weight;
    }
  }

  return mst;
}
int main() {
  int amount_of_vertexes;
  int amount_of_edges;
  std::cin >> amount_of_vertexes >> amount_of_edges;
  ListGraph<int, std::pair<int, int>, int> graph;
  for (int i = 0; i < amount_of_vertexes; ++i) {
    graph.AddVertex(i);
  }
  int vert1;
  int vert2;
  int weight;
  for (int i = 0; i < amount_of_edges; i++) {
    std::cin >> vert1 >> vert2 >> weight;
    graph.AddEdge({vert1 - 1, vert2 - 1}, weight);
  }

  Tree<int, std::pair<int, int>, int> mst = ElPrimo(graph, 0);
  std::cout << mst.global_weight;

  return 0;
}
