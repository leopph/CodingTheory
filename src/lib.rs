use num::traits::One;
use num::traits::Pow;
use num::traits::Zero;
use num::BigUint;

fn fast_pow(base: &BigUint, exp: &BigUint) -> BigUint {
    let zero: BigUint = BigUint::zero();
    let one: BigUint = BigUint::one();

    if *exp <= zero {
        return BigUint::zero();
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

fn fast_mod_pow(base: &BigUint, exp: &BigUint, modulus: &BigUint) -> BigUint {
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

fn miller_rabin(p: &BigUint) -> bool {
    let two = BigUint::from(2u8);

    if p == &two {
        return true;
    }

    let zero = BigUint::zero();
    let one = BigUint::one();

    let p_minus_one = p - &one;
    let mut m = p_minus_one.clone();
    let mut r: u64 = 0;

    while &m % &two == zero {
        m /= &two;
        r += 1;
    }

    let m = m;
    let r = r;
    let mut a = one.clone();

    // Initializing to true because a^m mod p == 1 results in a pass
    let mut prev_was_minus_one = true;

    while &a < &p {
        for i in 0..=r {
            let x = fast_mod_pow(&a, &(&m * &fast_pow(&two, &i.into())), &p);

            if prev_was_minus_one && x == one {
                return true;
            }

            prev_was_minus_one = x == p_minus_one;
        }

        a += 1u8;
    }

    false
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
    fn miller_rabin_test() {
        for p in [
            2u16, 3, 5, 7, 23, 383, 1031, 2087, 3359, 4447, 5519, 6329, 7919,
        ] {
            assert!(miller_rabin(&p.into()));
        }
    }
}
