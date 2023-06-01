// 83110962
#include <algorithm>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>

using std::vector;
const int kInf = 2009000999;

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
class ListGraph : public Graph<Vertextype, Edge, EdgeWeight> {
  using EdgeType = typename Graph<Vertextype, Edge, EdgeWeight>::EdgeType;

 public:
  void AddVertex(Vertextype ver) override {
    if (adj_.find(ver) == adj_.end()) {
      adj_[ver] = vector<std::pair<Vertextype, EdgeWeight>>();
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

  vector<std::pair<Vertextype, EdgeWeight>> Neighbours(
      Vertextype ver) const override {
    return adj_.at(ver);
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
template <typename Vertextype = int,
          typename Edge = std::pair<Vertextype, Vertextype>,
          typename EdgeWeight = double>

class DijkstraVisitor {
 public:
  DijkstraVisitor(ListGraph<Vertextype, Edge, EdgeWeight>& graph,
                  int start_point)
      : graph_(graph), kStartVertex(start_point) {
    distances_.reserve(graph.VertexCount());
    parents_.reserve(graph.VertexCount());
    for (auto iter = graph.Begin(); iter != graph.End(); ++iter) {
      distances_[iter->first] = kInf;
      parents_[iter->first] = -1;
    }
    distances_[start_point] = 0;
  };

  void Visit(Vertextype vertex) {
    for (const auto& [neighbor, cost] : graph_.Neighbours(vertex)) {
      EdgeWeight new_dist = distances_[vertex] + cost;
      if (new_dist < distances_[neighbor]) {
        distances_[neighbor] = new_dist;
        parents_[neighbor] = vertex;
      }
    }
  }

  void Dijkstra() {
    std::priority_queue<std::pair<EdgeWeight, Vertextype>,
                        vector<std::pair<EdgeWeight, Vertextype>>,
                        std::greater<>>
        queue;
    queue.emplace(0, kStartVertex);
    while (!queue.empty()) {
      auto [distance, vertex] = queue.top();
      queue.pop();
      if (visited_.contains(vertex)) {
        continue;
      }
      visited_.insert(vertex);
      Visit(vertex);
      for (const auto& [neighbor, cost] : graph_.Neighbours(vertex)) {
        if (!visited_.contains(neighbor)) {
          queue.push({distances_[neighbor], neighbor});
        }
      }
    }
  }

  const std::unordered_map<Vertextype, Vertextype>& Distances() const {
    return distances_;
  }

  const std::unordered_map<Vertextype, Vertextype>& Parents() const {
    return parents_;
  }

 private:
  ListGraph<Vertextype, Edge, EdgeWeight> graph_;
  const Vertextype kStartVertex;
  std::unordered_map<Vertextype, EdgeWeight> distances_;
  std::unordered_map<Vertextype, Vertextype> parents_;
  std::unordered_set<Vertextype> visited_;
};

int main() {
  int amount_of_requests;
  std::cin >> amount_of_requests;
  while (amount_of_requests != 0) {
    ListGraph<int, std::pair<int, int>, int> graph;
    int amount_of_vertexes;
    int amount_of_edges;
    std::cin >> amount_of_vertexes >> amount_of_edges;
    for (int i = 0; i < amount_of_vertexes; ++i) {
      graph.AddVertex(i);
    }
    for (int i = 0; i < amount_of_edges; ++i) {
      int temp1;
      int temp2;
      int temp_weight;
      std::cin >> temp1 >> temp2 >> temp_weight;
      graph.AddEdge({temp1, temp2}, temp_weight);
    }
    int start_point;
    std::cin >> start_point;
    DijkstraVisitor<int, std::pair<int, int>, int> dijkstra_visitor(
        graph, start_point);
    dijkstra_visitor.Dijkstra();
    std::unordered_map<int, int> output = dijkstra_visitor.Distances();
    for (int i = 0; i < static_cast<int>(output.size()); ++i) {
      std::cout << output[i] << " ";
    }
    std::cout << "\n";
    --amount_of_requests;
  }
}
