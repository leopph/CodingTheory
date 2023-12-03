use num::traits::One;
use num::traits::Pow;
use num::traits::Zero;
use num::BigInt;

fn fast_pow(base: &BigInt, exp: &BigInt) -> BigInt {
    let zero: BigInt = BigInt::zero();
    let one: BigInt = BigInt::one();

    if *exp <= zero {
        return BigInt::zero();
    }

    let mut mul = BigInt::from(1);
    let mut base = base.clone();
    let mut exp = exp.clone();

    while exp > one {
        if &exp % 2 == one {
            mul *= &base;
            exp -= &one;
        }

        base = Pow::pow(base, 2u8);
        exp >>= 2;
    }

    base * mul
}

fn fast_mod_pow(base: &BigInt, exp: &BigInt, modulus: &BigInt) -> BigInt {
    let zero: BigInt = BigInt::zero();
    let one: BigInt = BigInt::one();

    if *modulus == zero {
        return zero;
    }

    let mut ret = one.clone();
    let mut base = base % modulus;
    let mut exp = exp.clone();

    while exp > zero {
        if &exp % 2 == one {
            ret = ret * &base % modulus;
        }

        exp >>= 1;
        base = Pow::pow(base, 2u8) % modulus;
    }

    ret
}

#[cfg(test)]
mod tests {
    use crate::*;

    #[test]
    fn fast_pow_test() {
        assert_eq!(
            fast_pow(&BigInt::from(5), &BigInt::from(2)),
            BigInt::from(25)
        );
    }

    #[test]
    fn fast_mod_pow_test() {
        assert_eq!(
            fast_mod_pow(&BigInt::from(5), &BigInt::from(2), &BigInt::from(3)),
            BigInt::from(1)
        );
    }
}
