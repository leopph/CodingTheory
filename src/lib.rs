use num::bigint::RandBigInt;
use num::bigint::ToBigInt;
use num::traits::One;
use num::traits::Pow;
use num::traits::Zero;
use num::BigInt;
use num::BigUint;
use rand::thread_rng;

pub fn fast_pow(base: &BigUint, exp: &BigUint) -> BigUint {
    let zero: BigUint = BigUint::zero();
    let one: BigUint = BigUint::one();

    if *exp == zero {
        return BigUint::one();
    }

    let mut mul = BigUint::from(1u8);
    let mut base = base.clone();
    let mut exp = exp.clone();

    while exp > one {
        if &exp % 2u8 == one {
            mul *= &base;
            exp -= &one;
        }

        base = Pow::pow(base, 2u8);
        exp >>= 2;
    }

    base * mul
}

pub fn fast_mod_pow(base: &BigUint, exp: &BigUint, modulus: &BigUint) -> BigUint {
    let zero: BigUint = BigUint::zero();
    let one: BigUint = BigUint::one();

    if *modulus == zero {
        return zero;
    }

    let mut ret = one.clone();
    let mut base = base % modulus;
    let mut exp = exp.clone();

    while exp > zero {
        if &exp % 2u8 == one {
            ret = ret * &base % modulus;
        }

        exp >>= 1;
        base = Pow::pow(base, 2u8) % modulus;
    }

    ret
}

pub fn miller_rabin(p: &BigUint, test_count: &BigUint) -> bool {
    let two = BigUint::from(2u8);

    if p == &two {
        return true;
    }

    let zero = BigUint::zero();

    let p_minus_one = p - 1u8;
    let mut m = p_minus_one.clone();
    let mut r: u64 = 0;

    while &m % 2u8 == zero {
        m /= 2u8;
        r += 1;
    }

    let m = m;
    let r = r;

    let one = BigUint::one();

    let mut test_counter = zero.clone();

    'outer: while test_counter < *test_count {
        let a = thread_rng().gen_biguint_range(&one, &p_minus_one);
        // Initializing to true because a^m mod p == 1 is a pass
        let mut prev_was_minus_one = true;

        for i in 0..=r {
            let x = fast_mod_pow(&a, &(&m * &fast_pow(&two, &i.into())), p);

            if prev_was_minus_one && x == one {
                test_counter += 1u8;
                continue 'outer;
            }

            prev_was_minus_one = x == p_minus_one;
        }

        return false;
    }

    true
}

fn get_required_miller_rabin_test_count(p: &BigUint) -> BigUint {
    p / 4u8 + 1u8
}

// Extended euclidean algorithm in Z_n
pub fn get_modular_inverse(a: &BigUint, modulus: &BigUint) -> BigUint {
    let mut x_prev = BigInt::zero();
    let mut x_curr = BigInt::one();

    let mut r_prev = modulus.to_bigint().unwrap();
    let mut r_curr = a.to_bigint().unwrap();

    let zero = BigInt::zero();

    while r_curr != zero {
        let q = &r_prev / &r_curr;

        let tmp = r_curr;
        r_curr = r_prev - &q * &tmp;
        r_prev = tmp;

        let tmp = x_curr;
        x_curr = x_prev - &q * &tmp;
        x_prev = tmp;
    }

    if r_prev > BigInt::one() {
        panic!("{} has no modular inverse mod {}.", a, modulus);
    }

    if x_prev < zero {
        x_prev += modulus.to_bigint().unwrap();
    }

    x_prev.to_biguint().unwrap()
}

pub fn gen_rand_prime(bit_size: u64) -> BigUint {
    loop {
        let tmp = thread_rng().gen_biguint(bit_size);
        if miller_rabin(&tmp, &get_required_miller_rabin_test_count(&tmp)) {
            break tmp;
        }
    }
}

pub struct RSAKeys {
    pub n: BigUint,
    pub e: BigUint,
    pub d: BigUint,
}

pub fn gen_rsa_keys() -> RSAKeys {
    const PRIME_BIT_SIZE: u64 = 16;
    let p = gen_rand_prime(PRIME_BIT_SIZE);
    let q = gen_rand_prime(PRIME_BIT_SIZE);
    let n = &p * &q;
    let fi_n = (p - 1u8) * (q - 1u8);
    let e = BigUint::from(65537u64);
    let d = get_modular_inverse(&e, &fi_n);
    RSAKeys { n, e, d }
}

pub fn rsa_encode(msg: &BigUint, e: &BigUint, n: &BigUint) -> BigUint {
    fast_mod_pow(msg, e, n)
}

pub fn rsa_decode(cyp: &BigUint, d: &BigUint, n: &BigUint) -> BigUint {
    fast_mod_pow(cyp, d, n)
}

#[cfg(test)]
mod tests {
    use crate::*;

    #[test]
    fn fast_pow_test() {
        assert_eq!(
            fast_pow(&BigUint::from(5u8), &BigUint::from(2u8)),
            BigUint::from(25u8)
        );
    }

    #[test]
    fn fast_mod_pow_test() {
        assert_eq!(
            fast_mod_pow(
                &BigUint::from(5u8),
                &BigUint::from(2u8),
                &BigUint::from(3u8)
            ),
            BigUint::from(1u8)
        );
    }

    #[test]
    fn miller_rabin_primes() {
        for p in [
            2u16, 3, 5, 7, 23, 383, 1031, 2087, 3359, 4447, 5519, 6329, 7919,
        ] {
            let p: BigUint = p.into();
            assert!(miller_rabin(&p, &get_required_miller_rabin_test_count(&p)));
        }
    }

    #[test]
    fn miller_rabin_carmichaels() {
        for p in [
            561u32, 1105, 2465, 6601, 8911, 10585, 15841, 46657, 62745, 75361,
        ] {
            let p: BigUint = p.into();
            assert!(!miller_rabin(&p, &get_required_miller_rabin_test_count(&p)));
        }
    }

    #[test]
    fn rsa_key_test() {
        let keys = gen_rsa_keys();
        let m = BigUint::from(69u8);
        assert_eq!(m, fast_mod_pow(&m, &(keys.e * keys.d), &keys.n));
    }

    #[test]
    fn rsa_encode_decode_test() {
        let keys = gen_rsa_keys();
        let msg = BigUint::from(123u8);
        let cyp = rsa_encode(&msg, &keys.e, &keys.n);
        let decoded = rsa_decode(&cyp, &keys.d, &keys.n);
        assert_eq!(msg, decoded);
    }
}
