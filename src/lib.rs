use num_bigint::BigInt;

fn fast_pow(base: &BigInt, exp: &BigInt) -> BigInt {
    let zero: BigInt = BigInt::from(0);
    let one: BigInt = BigInt::from(1);

    if exp <= &zero {
        return BigInt::from(0);
    }

    let mut mul = BigInt::from(1);
    let mut base = base.clone();
    let mut exp = exp.clone();

    while exp > one {
        if &exp % 2 == one {
            mul *= &base;
            exp -= &one;
        }

        base = base.pow(2);
        exp >>= 2;
    }

    base * mul
}

fn fast_mod_pow(base: &BigInt, exp: &BigInt, modulus: &BigInt) -> BigInt {
    let zero: BigInt = BigInt::from(0);
    let one: BigInt = BigInt::from(1);

    if modulus == &zero {
        return zero;
    }

    let mut ret = BigInt::from(1);
    let mut base = base % modulus;
    let mut exp = exp.clone();

    while exp > zero {
        if &exp % 2 == one {
            ret = ret * &base % modulus;
        }

        exp >>= 1;
        base = base.pow(2) % modulus;
    }

    ret
}

#[cfg(test)]
mod tests {
    use crate::*;
    use num_bigint::BigInt;

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
