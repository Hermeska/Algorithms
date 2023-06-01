// 87071969
#include "iostream"
#include "map"
#include "queue"
#include "string"
#include "vector"

const int kAlphabetSize = 26;
const int kNull = -1;
using std::cin;
using std::cout;
using std::queue;
using std::string;
using std::vector;
class CorAss {
 public:
  CorAss();

  void Add(const string& str);

  void Update();

  vector<vector<int>> Find(const string& text);

 private:
  struct Node {
    explicit Node(int depth);
    vector<int> transitions;
    bool terminal;
    vector<int> patterns_nums;
    int depth;
  };

  void SetLinks();

  void FindCompressedLink(int vertex);

  vector<int> link_;
  vector<int> compressed_;
  vector<vector<int>> move_;
  int size_of_dictionary_;
  vector<Node> bor_;
};

CorAss::Node::Node(int depth)
    : transitions(kAlphabetSize, kNull), terminal(false), depth(depth) {}

CorAss::CorAss() : size_of_dictionary_(0) {}

void CorAss::Add(const string& str) {
  if (bor_.empty()) {
    bor_.emplace_back(0);
  }
  int vertex = 0;
  for (char symbol : str) {
    if (bor_[vertex].transitions[symbol - 'a'] == kNull) {
      bor_[vertex].transitions[symbol - 'a'] = static_cast<int>(bor_.size());
      bor_.emplace_back(bor_[vertex].depth + 1);
    }
    vertex = bor_[vertex].transitions[symbol - 'a'];
  }
  bor_[vertex].terminal = true;
  bor_[vertex].patterns_nums.push_back(size_of_dictionary_);
  ++size_of_dictionary_;
}
void CorAss::SetLinks() {
  queue<int> queue;
  queue.push(0);
  while (!queue.empty()) {
    auto vertex = queue.front();
    queue.pop();
    for (int sym = 0; sym < kAlphabetSize; ++sym) {
      auto to_vertex = bor_[vertex].transitions[sym];
      if (to_vertex == kNull) {
        continue;
      }
      link_[to_vertex] = ((vertex == 0) ? 0 : move_[link_[vertex]][sym]);
      FindCompressedLink(to_vertex);
      for (int sym_of_to = 0; sym_of_to < kAlphabetSize; ++sym_of_to) {
        if (bor_[to_vertex].transitions[sym_of_to] == kNull) {
          move_[to_vertex][sym_of_to] = move_[link_[to_vertex]][sym_of_to];
        } else {
          move_[to_vertex][sym_of_to] = bor_[to_vertex].transitions[sym_of_to];
        }
      }
      queue.push(to_vertex);
    }
  }
}

void CorAss::Update() {
  move_ = vector<vector<int>>(bor_.size(), vector<int>(kAlphabetSize));
  link_ = vector<int>(bor_.size());
  compressed_ = vector<int>(bor_.size(), 0);
  link_[0] = 0;
  for (int iter = 0; iter < kAlphabetSize; ++iter) {
    if (bor_[0].transitions[iter] != kNull) {
      move_[0][iter] = bor_[0].transitions[iter];
    } else {
      move_[0][iter] = 0;
    }
  }
  SetLinks();
}

vector<vector<int>> CorAss::Find(const string& text) {
  int cur_vertex = 0;
  vector<vector<int>> answer(size_of_dictionary_);
  for (size_t i = 0; i < text.size(); ++i) {
    cur_vertex = move_[cur_vertex][text[i] - 'a'];
    for (int link_vertex = cur_vertex; link_vertex != 0;
         link_vertex = compressed_[link_vertex]) {
      if (!bor_[link_vertex].terminal) {
        continue;
      }
      for (auto pattern : bor_[link_vertex].patterns_nums) {
        answer[pattern].push_back(i + 1 + 1 - bor_[link_vertex].depth);
      }
    }
  }
  return answer;
}

void CorAss::FindCompressedLink(int vertex) {
  auto link_vertex = link_[vertex];
  compressed_[vertex] =
      (bor_[link_vertex].terminal) ? link_vertex : compressed_[link_vertex];
}

class Solution {
 public:
  Solution();

  void Solve();

 private:
  string some_text_;
  int size_of_dict_;
  CorAss borrov_;
};

Solution::Solution() : size_of_dict_(0) {}

void Solution::Solve() {
  cin >> some_text_;
  cin >> size_of_dict_;
  size_t temp = size_of_dict_;
  while (temp > 0) {
    string some_string;
    cin >> some_string;
    borrov_.Add(some_string);
    --temp;
  }
  borrov_.Update();
  auto answer = borrov_.Find(some_text_);
  for (const auto& positions : answer) {
    cout << positions.size() << ' ';
    for (auto val : positions) {
      cout << val << ' ';
    }
    cout << std::endl;
  }
}
int main() {
  Solution solver;
  solver.Solve();
}
