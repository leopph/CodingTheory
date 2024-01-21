use num::bigint::RandBigInt;
use num::traits::One;
use num::traits::Pow;
use num::traits::Zero;
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
            let x = fast_mod_pow(&a, &(&m * &fast_pow(&two, &i.into())), &p);

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
    return p / 4u8 + 1u8;
}

pub fn gcd_euclid(a: &BigUint, b: &BigUint) -> BigUint {
    let zero = BigUint::zero();

    let mut prev = a.clone();
    let mut curr = b.clone();

    loop {
        let prev_prev = prev;
        prev = curr;
        curr = prev_prev % &prev;

        if curr == zero {
            return prev;
        }
    }
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
    fn gcd_euclid_test() {
        for (a, b, gcd) in [
            (2u32, 3u32, 1u32),
            (4, 5, 1),
            (6, 9, 3),
            (15, 105, 15),
            (42, 56, 14),
            (24826148, 45296490, 526),
        ] {
            assert_eq!(gcd_euclid(&a.into(), &b.into()), gcd.into());
        }
    }
}
