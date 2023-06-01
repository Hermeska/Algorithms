// 87294347
#include "algorithm"
#include "complex"
#include "iostream"
#include "vector"

template <typename T>
class Field {
 public:
  Field();

  Field(const T& value);

  Field(const Field& other);

  Field& operator=(const Field& other);

  Field& operator+=(const Field& other);

  Field& operator-=(const Field& other);

  Field& operator*=(const Field& other);

  Field& operator/=(const Field& other);

  Field Inverse();

  bool operator==(const Field& other);

  bool operator!=(const Field& other);

  static Field GetElementaryUnit(int power);

  T& GetValue();

 private:
  constexpr static const long double kPi = 3.14159265358979323846;
  T value_;
};

template <typename Type>
Field<Type>::Field() = default;

template <typename Type>
Field<Type>::Field(const Type& value) : value_(value) {}

template <typename Type>
Field<Type>::Field(const Field& other) = default;

template <typename Type>
Field<Type>& Field<Type>::operator=(const Field<Type>& other) = default;

template <typename Type>
Type& Field<Type>::GetValue() {
  return value_;
}

template <typename Type>
Field<Type>& Field<Type>::operator+=(const Field& other) {
  value_ += other.value_;
  return *this;
}

template <typename Type>
Field<Type>& Field<Type>::operator-=(const Field& other) {
  value_ -= other.value_;
  return *this;
}

template <typename Type>
Field<Type>& Field<Type>::operator*=(const Field& other) {
  value_ *= other.value_;
  return *this;
}

template <typename Type>
Field<Type>& Field<Type>::operator/=(const Field& other) {
  value_ /= other.value_;
  return *this;
}

template <>
Field<std::complex<long double>>
Field<std::complex<long double>>::GetElementaryUnit(int power) {
  long double angle = 2 * kPi / power;
  return std::polar<long double>(1, angle);
}

template <typename Type>
bool Field<Type>::operator==(const Field& other) {
  return value_ == other.value_;
}

template <typename Type>
bool Field<Type>::operator!=(const Field& other) {
  return !(value_ == other.value_);
}

template <>
Field<std::complex<long double>> Field<std::complex<long double>>::Inverse() {
  return std::conj(value_);
}

using ComplexField = Field<std::complex<long double>>;

template <typename Type>
Field<Type> operator+(const Field<Type>& r_value, const Field<Type>& l_value) {
  Field<Type> copy = r_value;
  copy += l_value;
  return copy;
}

template <typename Type>
Field<Type> operator-(const Field<Type>& r_value, const Field<Type>& l_value) {
  Field<Type> copy = r_value;
  copy -= l_value;
  return copy;
}

template <typename Type>
Field<Type> operator*(const Field<Type>& r_value, const Field<Type>& l_value) {
  Field<Type> copy = r_value;
  copy *= l_value;
  return copy;
}

template <typename Type>
Field<Type> operator/(const Field<Type>& r_value, const Field<Type>& l_value) {
  Field<Type> copy = r_value;
  copy /= l_value;
  return copy;
}

template <typename Field>
class FFTConverter {
 public:
  void Convert(std::vector<Field>& values);
  void Invert(std::vector<Field>& values);
};

template <typename Field>
void FFTConverter<Field>::Convert(std::vector<Field>& values) {
  uint64_t size = values.size();
  if (size <= 1) {
    return;
  }

  std::vector<Field> even_part(size / 2);
  std::vector<Field> odd_part(size / 2);
  for (uint64_t i = 0; i < size / 2; i++) {
    even_part[i] = values[2 * i];
    odd_part[i] = values[2 * i + 1];
  }

  Convert(even_part);
  Convert(odd_part);

  Field root(1);
  Field elementary_unit = Field().GetElementaryUnit(size);
  for (uint64_t i = 0; i < size / 2; ++i) {
    auto adder = root * odd_part[i];
    values[i] = even_part[i] + adder;
    values[i + size / 2] = even_part[i] - adder;
    root *= elementary_unit;
  }
}

template <typename Field>
void FFTConverter<Field>::Invert(std::vector<Field>& values) {
  for (auto& element : values) {
    element = element.Inverse();
  }
  Convert(values);
  for (auto& element : values) {
    element = element.Inverse();
    element /= Field(values.size());
  }
}

template <typename T>
class Polynomial {
 public:
  Polynomial();

  Polynomial(uint64_t size);

  Polynomial(const std::vector<T>& coefficients);

  Polynomial(const Polynomial& other_polynomial);

  Polynomial& operator=(const Polynomial& other_polynomial);

  T& operator[](unsigned long long index);

  const T& operator[](unsigned long long index) const;

  Polynomial& operator*=(const Polynomial& other_polynomial);

  void WriteCoefficients();

 private:
  void AddZeroCoefficients(size_t new_degree);

  std::vector<T> coefficients_;
  uint64_t degree_;
  FFTConverter<T> converter_;
};

template <typename T>
Polynomial<T>::Polynomial() = default;

template <typename T>
Polynomial<T>::Polynomial(uint64_t size) : degree_(size) {}

template <typename T>
Polynomial<T>::Polynomial(const std::vector<T>& coefficients)
    : coefficients_(coefficients), degree_(coefficients_.size()) {}

template <typename T>
Polynomial<T>::Polynomial(const Polynomial& other_polynomial)
    : coefficients_(other_polynomial.coefficients_),
      degree_(other_polynomial.degree_) {}

template <typename T>
Polynomial<T>& Polynomial<T>::operator=(const Polynomial& other_polynomial) {
  coefficients_ = other_polynomial.coefficients_;
  degree_ = other_polynomial.degree_;
  return *this;
}

template <typename T>
T& Polynomial<T>::operator[](unsigned long long index) {
  return coefficients_[index];
}

template <typename T>
const T& Polynomial<T>::operator[](unsigned long long index) const {
  return coefficients_[index];
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator*=(const Polynomial& other_polynomial) {
  Polynomial copy = other_polynomial;
  uint64_t end_degree = degree_ + other_polynomial.degree_ - 1;
  uint64_t new_degree = 1;
  while (new_degree < degree_ + other_polynomial.degree_) {
    new_degree <<= 1;
  }
  copy.AddZeroCoefficients(new_degree);
  AddZeroCoefficients(new_degree);
  converter_.Convert(coefficients_);
  converter_.Convert(copy.coefficients_);
  for (uint64_t i = 0; i < new_degree; ++i) {
    coefficients_[i] *= copy.coefficients_[i];
  }
  converter_.Invert(coefficients_);
  AddZeroCoefficients(end_degree);
  return *this;
}

template <typename T>
void Polynomial<T>::AddZeroCoefficients(size_t new_degree) {
  coefficients_.resize(new_degree);
  degree_ = new_degree;
}

template <typename T>
void Polynomial<T>::WriteCoefficients() {
  std::cout << degree_ - 1 << ' ';
  for (uint64_t i = degree_ - 1; i > 0; --i) {
    std::cout << llround(coefficients_[i].GetValue().real()) << ' ';
  }
  std::cout << llround(coefficients_[0].GetValue().real()) << ' ';
}

class Solver {
 public:
  Solver();

  void Solve();

 private:
  uint64_t size_1_;
  uint64_t size_2_;
  std::vector<ComplexField> coefficients_1_;
  std::vector<ComplexField> coefficients_2_;

  void Read();
};

Solver::Solver() : size_1_(0), size_2_(0) { Read(); }

void Solver::Solve() {
  Polynomial<ComplexField> polynomial_1(coefficients_1_);
  Polynomial<ComplexField> polynomial_2(coefficients_2_);
  polynomial_1 *= polynomial_2;
  polynomial_1.WriteCoefficients();
}

void Solver::Read() {
  std::cin >> size_1_;
  coefficients_1_.resize(size_1_ + 1);
  for (auto& element : coefficients_1_) {
    long double val;
    std::cin >> val;
    element = ComplexField(val);
  }
  std::reverse(coefficients_1_.begin(), coefficients_1_.end());
  std::cin >> size_2_;
  coefficients_2_.resize(size_2_ + 1);
  for (auto& element : coefficients_2_) {
    long double val;
    std::cin >> val;
    element = ComplexField(val);
  }
  std::reverse(coefficients_2_.begin(), coefficients_2_.end());
}

int main() {
  Solver solver;
  solver.Solve();
}
