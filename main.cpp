#include <iostream>
#include <concepts>

// Only supports non-negative exponents.
template<std::integral T>
[[nodiscard]] auto FastPow(T base, T exp) noexcept -> T {
  if (exp <= 0) {
    return 1;
  }

  T mul{1};

  while (exp > 1) {
    if (exp % 2 == 1) {
      mul *= base;
      exp -= 1;
    }

    base *= base;
    exp /= 2;
  }

  return base * mul;
}

int main() {
  for (auto i{0}; i < 50; i++) {
    std::cout << FastPow<long long>(2, i) << '\n';
  }
}