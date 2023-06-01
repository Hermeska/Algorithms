// 87284286
#include "iostream"
#include "map"
#include "string"
#include "vector"
const int kStart = 96;
const int kFinish = 123;
const int kAscii = 32;

class Suff {
 public:
  Suff();
  void AddString(const std::string& string);
  bool IsInAutomaton(const std::string& string);

 private:
  static const long long kNull = -1;
  long long size_;
  std::vector<long long> len_;
  std::vector<long long> link_;
  std::vector<std::map<char, long long>> transitive_;
  long long terminal_;
  void AddChar(char symbol);
  void AddVertex();
};

Suff::Suff() : size_(1), len_(1), link_(1), transitive_(1), terminal_(0) {
  len_[terminal_] = 0;
  link_[terminal_] = kNull;
}

void Suff::AddString(const std::string& string) {
  for (auto symbol : string) {
    AddChar(symbol);
  }
}

bool Suff::IsInAutomaton(const std::string& string) {
  long long curr = 0;
  for (char symbol : string) {
    if (curr == static_cast<long long>(transitive_.size()) ||
        static_cast<int>(transitive_[curr].contains(symbol)) == 0U) {
      return false;
    }
    curr = transitive_[curr][symbol];
  }
  return true;
}

void Suff::AddChar(char symbol) {
  AddVertex();
  long long added_vertex = size_ - 1;
  len_[added_vertex] = len_[terminal_] + 1;
  long long updated_vertex = terminal_;
  while (updated_vertex != kNull &&
         !(transitive_[updated_vertex].contains(symbol))) {
    transitive_[updated_vertex][symbol] = added_vertex;
    updated_vertex = link_[updated_vertex];
  }
  if (updated_vertex == kNull) {
    link_[added_vertex] = 0;
    terminal_ = added_vertex;
    return;
  }
  long long suspicious_vertex = transitive_[updated_vertex][symbol];
  if (len_[suspicious_vertex] == len_[updated_vertex] + 1) {
    link_[added_vertex] = suspicious_vertex;
    terminal_ = added_vertex;
    return;
  }
  AddVertex();
  long long clone_of_suspicious_vertex = size_ - 1;
  len_[clone_of_suspicious_vertex] = len_[updated_vertex] + 1;
  link_[clone_of_suspicious_vertex] = link_[suspicious_vertex];
  link_[suspicious_vertex] = clone_of_suspicious_vertex;
  link_[added_vertex] = clone_of_suspicious_vertex;
  transitive_[clone_of_suspicious_vertex] = transitive_[suspicious_vertex];
  while (updated_vertex != kNull &&
         transitive_[updated_vertex].find(symbol)->second ==
             suspicious_vertex) {
    transitive_[updated_vertex][symbol] = clone_of_suspicious_vertex;
    updated_vertex = link_[updated_vertex];
  }
  terminal_ = added_vertex;
}

void Suff::AddVertex() {
  ++size_;
  link_.emplace_back(0);
  transitive_.emplace_back();
  len_.emplace_back(0);
}

void LowerCase(std::string& string) {
  for (char& symbol : string) {
    if (symbol > kStart and symbol < kFinish) {
      symbol -= kAscii;
    }
  }
}

std::string GetQuery(const std::string& query) {
  return query.substr(2, query.size() - 2);
}
void Run() {
  Suff automaton;
  std::string query;

  while (getline(std::cin, query)) {
    LowerCase(query);
    if (query[0] == '?') {
      if (automaton.IsInAutomaton(GetQuery(query))) {
        std::cout << "YES\n";
      } else {
        std::cout << "NO\n";
      }
    } else {
      automaton.AddString(GetQuery(query));
    }
  }
}
int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);
  std::cout.tie(nullptr);
  Run();
}
