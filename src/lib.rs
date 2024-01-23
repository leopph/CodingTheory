use num::bigint::RandBigInt;
use num::range;
use num::traits::One;
use num::traits::Pow;
use num::traits::Zero;
use num::BigInt;
use num::BigUint;
use num::Integer;
use rand::thread_rng;

pub fn fast_pow(mut base: BigUint, mut exp: BigUint) -> BigUint {
    let one: BigUint = BigUint::one();

    if exp.is_zero() {
        return one;
    }

    let mut mul = one.clone();

    while exp > one {
        if exp.is_odd() {
            mul *= &base;
            exp -= 1u8;
        }

        base = base.pow(2u8);
        exp >>= 2;
    }

    base * mul
}

pub fn fast_mod_pow(mut base: BigUint, mut exp: BigUint, modulus: &BigUint) -> BigUint {
    if modulus.is_zero() {
        return BigUint::zero();
    }

    base %= modulus;
    let mut ret = BigUint::one();

    while !exp.is_zero() {
        if exp.is_odd() {
            ret = ret * &base % modulus;
        }

        exp >>= 1;
        base = base.pow(2u8) % modulus;
    }

    ret
}

pub fn miller_rabin(p: &BigUint, test_count: BigUint) -> bool {
    let p_minus_one = p - 1u8;

    let (m, r) = {
        let mut m = p_minus_one.clone();
        let mut r: u64 = 0;

        while m.is_even() {
            m /= 2u8;
            r += 1;
        }

        (m, r)
    };

    let one = BigUint::one();

    'outer: for _ in range(BigUint::zero(), test_count) {
        let a = thread_rng().gen_biguint_range(&one, &p_minus_one);
        // Initializing to true because a^m mod p == 1 is a pass
        let mut prev_was_minus_one = true;

        for i in 0..=r {
            let x = fast_mod_pow(a.clone(), &m * fast_pow(BigUint::from(2u8), i.into()), p);

            if prev_was_minus_one && x.is_one() {
                continue 'outer;
            }

            prev_was_minus_one = x == p_minus_one;
        }

        return false;
    }

    true
}

fn get_miller_rabin_test_count(p: BigUint) -> BigUint {
    p / 4u8 + 1u8
}

// Extended euclidean algorithm in Z_n
pub fn get_modular_inverse(a: BigUint, modulus: BigUint) -> BigUint {
    let modulus = BigInt::from(modulus);

    let mut x_prev = BigInt::zero();
    let mut x_curr = BigInt::one();

    let mut r_prev = modulus.clone();
    let mut r_curr = BigInt::from(a);

    while !r_curr.is_zero() {
        let q = &r_prev / &r_curr;

        let tmp = r_curr;
        r_curr = r_prev - &q * &tmp;
        r_prev = tmp;

        let tmp = x_curr;
        x_curr = x_prev - &q * &tmp;
        x_prev = tmp;
    }

    if r_prev > BigInt::one() {
        panic!("No modular inverse.");
    }

    if x_prev < BigInt::zero() {
        x_prev += modulus;
    }

    BigUint::try_from(x_prev).unwrap()
}

pub fn gen_rand_prime(bit_size: u64) -> BigUint {
    loop {
        let tmp = thread_rng().gen_biguint(bit_size);
        if miller_rabin(&tmp, get_miller_rabin_test_count(tmp.clone())) {
            break tmp;
        }
    }
}

pub struct RSAKeys {
    pub n: BigUint,
    pub e: BigUint,
    pub d: BigUint,
}

impl RSAKeys {
    pub fn gen() -> RSAKeys {
        const PRIME_BIT_SIZE: u64 = 16;
        let p = gen_rand_prime(PRIME_BIT_SIZE);
        let q = gen_rand_prime(PRIME_BIT_SIZE);
        let n = &p * &q;
        let fi_n = (p - 1u8) * (q - 1u8);
        let e = BigUint::from(65537u64);
        let d = get_modular_inverse(e.clone(), fi_n);
        RSAKeys { n, e, d }
    }
}

pub fn rsa_encode(msg: BigUint, e: BigUint, n: &BigUint) -> BigUint {
    fast_mod_pow(msg, e, n)
}

pub fn rsa_decode(cyp: BigUint, d: BigUint, n: &BigUint) -> BigUint {
    fast_mod_pow(cyp, d, n)
}

fn calc_jacobi_symbol(num: BigInt, denom: BigUint) -> i8 {
    let mut p = BigInt::from(denom);
    let mut a = num % &p;

    let mut ret = 1;

    let three = BigInt::from(3u8);
    let four = BigInt::from(4u8);
    let five = BigInt::from(5u8);

    while !a.is_zero() {
        let p_mod_8 = &p % 8u8;

        while a.is_even() {
            a /= 2u8;

            if p_mod_8 == three || p_mod_8 == five {
                ret = -ret;
            }
        }

        std::mem::swap(&mut a, &mut p);

        if &a % &four == three && &p % &four == three {
            ret = -ret;
        }

        a %= &p;
    }

    if p == BigInt::one() {
        ret
    } else {
        0
    }
}

pub fn solovay_strassen(n: BigUint, test_count: BigUint) -> bool {
    let n = BigInt::from(n);
    let zero = BigInt::zero();
    let one = BigInt::one();
    let n_minus_one_over_two = (n.clone() - 1) / 2;

    for _ in range(BigUint::zero(), test_count) {
        let a = loop {
            let a = thread_rng().gen_bigint_range(&one, &n);
            if a.gcd(&n).is_one() {
                break a;
            }
        };

        let mut expected = a.modpow(&n_minus_one_over_two, &n);
        let mut ls = BigInt::from(calc_jacobi_symbol(a, BigUint::try_from(n.clone()).unwrap()));

        if ls < zero {
            ls += &n;
        }

        if expected < zero {
            expected += &n;
        }

        if ls != expected {
            return false;
        }
    }

    true
}

pub fn get_solovay_strassen_test_count(p: BigUint) -> BigUint {
    p / 2u8 + 1u8
}

pub fn fermat_factorize(n: &BigInt) -> (BigInt, BigInt) {
    let mut a = {
        let mut root = n.sqrt();
        if root.clone().pow(2u8) != *n {
            root += 1;
        }
        root
    };

    let mut b = a.clone().pow(2u8) - n;

    while b.sqrt().pow(2u8) != b {
        a += 1;
        b = a.clone().pow(2u8) - n;
    }

    (a, b.sqrt())
}

pub fn gcd_euclid(a: BigUint, b: BigUint) -> BigUint {
    let mut prev = a;
    let mut curr = b;

    while !curr.is_zero() {
        let next = prev % &curr;
        prev = curr;
        curr = next;
    }

    prev
}

pub fn pollard_rho_factorize(n: &BigUint) -> BigUint {
    let fun = |x: BigUint| (x.clone().pow(2u8) + 1u8) % n;

    let mut v = vec![BigUint::one(), fun(BigUint::one())];

    loop {
        let last = v.last().unwrap();

        for elem in v.iter().take(v.len() - 1) {
            let diff = last - elem;

            if !diff.is_zero() {
                let gcd = gcd_euclid(n.clone(), diff);

                if !gcd.is_one() {
                    return gcd;
                }
            }
        }

        v.push(fun(last.clone()));
    }
}

#[cfg(test)]
mod tests {
    use crate::*;

    #[test]
    fn fast_pow_test() {
        assert_eq!(
            fast_pow(BigUint::from(5u8), BigUint::from(2u8)),
            BigUint::from(25u8)
        );
    }

    #[test]
    fn fast_mod_pow_test() {
        assert_eq!(
            fast_mod_pow(BigUint::from(5u8), BigUint::from(2u8), &BigUint::from(3u8)),
            BigUint::from(1u8)
        );
    }

    fn get_real_primes() -> &'static [u16] {
        static PRIMES: &[u16] = &[3, 5, 7, 23, 383, 1031, 2087, 3359, 4447, 5519, 6329, 7919];
        PRIMES
    }

    fn get_carmichaels() -> &'static [u32] {
        static CARMICHAELS: &[u32] = &[
            561, 1105, 2465, 6601, 8911, 10585, 15841, 46657, 62745, 75361,
        ];
        CARMICHAELS
    }

    #[test]
    fn miller_rabin_primes() {
        for p in get_real_primes() {
            let p = BigUint::from(*p);
            assert!(miller_rabin(&p, get_miller_rabin_test_count(p.clone())));
        }
    }

    #[test]
    fn miller_rabin_carmichaels() {
        for p in get_carmichaels() {
            let p: BigUint = BigUint::from(*p);
            assert!(!miller_rabin(&p, get_miller_rabin_test_count(p.clone())));
        }
    }

    #[test]
    fn rsa_key_test() {
        let keys = RSAKeys::gen();
        let m = BigUint::from(69u8);
        assert_eq!(m, fast_mod_pow(m.clone(), keys.e * keys.d, &keys.n));
    }

    #[test]
    fn rsa_encode_decode_test() {
        let keys = RSAKeys::gen();
        let msg = BigUint::from(123u8);
        let cyp = rsa_encode(msg.clone(), keys.e, &keys.n);
        let decoded = rsa_decode(cyp, keys.d, &keys.n);
        assert_eq!(msg, decoded);
    }

    #[test]
    fn jacobi_symbol_test() {
        assert_eq!(calc_jacobi_symbol(BigInt::from(30), BigUint::from(37u8)), 1);
    }

    #[test]
    fn jacobi_symbol_test1() {
        assert_eq!(calc_jacobi_symbol(BigInt::from(1), BigUint::from(3u8)), 1);
    }

    #[test]
    fn solovay_strassen_primes() {
        for p in get_real_primes() {
            let p = BigUint::from(*p);
            assert!(solovay_strassen(
                p.clone(),
                get_solovay_strassen_test_count(p)
            ));
        }
    }

    #[test]
    fn solovay_strassen_carmichaels() {
        for p in get_carmichaels() {
            let p = BigUint::from(*p);
            assert!(!solovay_strassen(
                p.clone(),
                get_solovay_strassen_test_count(p)
            ));
        }
    }

    #[test]
    fn fermat_factorization() {
        let (a, b) = fermat_factorize(&BigInt::from(517));
        assert_eq!(a, BigInt::from(29));
        assert_eq!(b, BigInt::from(18));
    }

    #[test]
    fn gcd_euclid_test() {
        for (a, b, gcd) in [
            (2u32, 3u32, 1u32),
            (4, 5, 1),
            (6, 9, 3),
            (15, 105, 15),
            (42, 56, 14),
            (24826148, 45296490, 526),
        ] {
            assert_eq!(gcd_euclid(a.into(), b.into()), gcd.into());
        }
    }

    #[test]
    fn pollard_rho_factorization() {
        assert_eq!(
            BigUint::from(7u8),
            pollard_rho_factorize(&BigUint::from(91u8))
        );
    }
}
