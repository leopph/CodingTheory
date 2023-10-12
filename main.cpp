#include <iostream>
#include <concepts>

// Only supports non-negative exponents.
template <std::integral T>
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

template <std::integral T>
[[nodiscard]] auto MillerRabin(T const p) noexcept -> bool {
  if (p <= 1) {
    return false;
  }

  if (p == 2) {
    return true;
  }

  if (p % 2 == 0) {
    return false;
  }

  auto const r{
    [tmp = p - 1]() mutable {
      T ret{0};

      while (tmp % 2 == 0) {
        ++ret;
        tmp /= 2;
      }

      return ret;
    }()
  };

  auto const m{(p - 1) / FastPow<T>(2, r)};

  for (T a{2}; a < p; ++a) {
    bool pass{false};

    auto prev{p - 1};

    for (T i{0}; i <= r; ++i) {
      auto cur{FastPow(a, m * FastPow(2, i)) % p};

      if (cur == 1 && prev == p - 1) {
        pass = true;
        break;
      }

      prev = cur;
    }

    if (!pass) {
      return false;
    }
  }

  return true;
}

auto main() -> int {
  for (auto i{3}; i <= 7; i++) {
    std::cout << i << ": " << MillerRabin(i) << '\n';
  }
}
