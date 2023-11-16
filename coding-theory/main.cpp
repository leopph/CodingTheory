#include <iostream>
#include <random>
#include <cstdint>

// Only supports non-negative exponents.
[[nodiscard]] static auto FastPow(std::int64_t base, std::int64_t exp) noexcept -> std::int64_t {
  if (exp <= 0) {
    return 1;
  }

  std::int64_t mul{1};

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

// Only supports non-negative exponents.
[[nodiscard]] static auto FastModPow(std::int64_t base, std::int64_t exp, std::int64_t mod) noexcept -> std::int64_t {
  if (mod == 1) {
    return 0;
  }

  std::int64_t ret{1};
  base %= mod;

  while (exp > 0) {
    if (exp % 2 == 1) {
      ret = ret * base % mod;
    }
    exp >>= 1;
    base = base * base % mod;
  }

  return ret;
}

enum class MillerRabinMode {
  TestRandom = 0,
  // Test a random set of potential witnesses.
  TestAll = 1 // Test using all potential witnesses.
};

/**
 * \brief Tests for (probable) primality using the Miller-Rabin algorithm.
 * The number of tested potential witnesses can be changed using mode and randomWitnessCount.
 * If mode is TestAll, all numbers in [2, p - 1] are tested.
 * If mode is TestRandom, min(p - 1, randomWitnessCount) number of random
 * numbers are tested in [2, p - 1].
 * \param p The number to test the primality of
 * \param mode Specifies how witnesses are tested.
 * \param randomWitnessCount Specifies how many random witnesses are tested.
 * \return Whether the number is a (potential) prime.
 */
[[nodiscard]] static auto MillerRabin(std::int64_t const p, MillerRabinMode const mode = MillerRabinMode::TestAll, std::int64_t const randomWitnessCount = 10) noexcept -> bool {
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
      std::int64_t ret{0};

      while (tmp % 2 == 0) {
        ++ret;
        tmp /= 2;
      }

      return ret;
    }()
  };

  auto const m{(p - 1) / FastPow(2, r)};

  auto const testWitness{
    [p, r, m](std::int64_t const a) {
      auto prev{p - 1};

      for (std::int64_t i{0}; i <= r; ++i) {
        auto cur{FastModPow(a, m * FastPow(2, i), p)};

        if (cur == 1 && prev == p - 1) {
          return true;
        }

        prev = cur;
      }

      return false;
    }
  };

  if (mode == MillerRabinMode::TestAll) {
    for (std::int64_t a{2}; a < p; ++a) {
      if (!testWitness(a)) {
        return false;
      }
    }
  } else {
    std::uniform_int_distribution<std::int64_t> dist{2, p - 1};
    std::default_random_engine gen{std::random_device{}()};

    auto const actualWitnessCount{std::min(randomWitnessCount, p - 1)};

    for (std::int64_t i{0}; i < actualWitnessCount; ++i) {
      if (!testWitness(dist(gen))) {
        return false;
      }
    }
  }


  return true;
}

auto main() -> int {
  for (auto i{3ll}; i <= 50; i++) {
    std::cout << i << ": " << MillerRabin(i, MillerRabinMode::TestAll) << '\n';
  }
}
