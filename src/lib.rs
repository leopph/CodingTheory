rozsda::rozsda! {
használj num::bigint::RandBigInt;
használj num::range;
használj num::traits::Euclid;
használj num::traits::One;
használj num::traits::Pow;
használj num::traits::Zero;
használj num::BigInt;
használj num::BigUint;
használj num::Integer;
használj num::Signed;
használj rand::thread_rng;

nyilvános függvény fast_pow(megváltoztatható base: BigUint, megváltoztatható exp: BigUint) -> BigUint {
    legyen megváltoztatható ret = BigUint::one();

    amíg !exp.is_zero() {
        ha exp.is_odd() {
            ret *= &base;
        }

        base = base.pow(2u8);
        exp >>= 1;
    }

    ret
}

nyilvános függvény fast_mod_pow(megváltoztatható base: BigUint, megváltoztatható exp: BigUint, modulus: &BigUint) -> BigUint {
    ha modulus.is_zero() {
        pánik!("Modulus wmint zero.");
    }

    ha modulus.is_one() {
        visszatérít BigUint::zero();
    }

    legyen megváltoztatható ret = BigUint::one();
    base %= modulus;

    amíg !exp.is_zero() {
        ha exp.is_odd() {
            ret = ret * &base % modulus;
        }

        base = base.pow(2u8) % modulus;
        exp >>= 1;
    }

    ret
}

nyilvános függvény miller_rabin(p: &BigUint, test_count: BigUint) -> bool {
    legyen p_minus_one = p - 1u8;

    legyen (m, r) = {
        legyen megváltoztatható m = p_minus_one.clone();
        legyen megváltoztatható r: u64 = 0;

        amíg m.is_even() {
            m /= 2u8;
            r += 1;
        }

        (m, r)
    };

    legyen one = BigUint::one();

    'outer: minden _ ebben range(BigUint::zero(), test_count) {
        legyen a = thread_rng().gen_biguint_range(&one, &p_minus_one);
        // Initializing to igaz because a^m mod p == 1 is a pass
        legyen megváltoztatható prev_was_minus_one = igaz;

        minden i ebben 0..=r {
            legyen x = fast_mod_pow(a.clone(), &m * fast_pow(BigUint::ebből(2u8), i.ebbe()), p);

            ha prev_was_minus_one && x.is_one() {
                continue 'outer;
            }

            prev_was_minus_one = x == p_minus_one;
        }

        visszatérít hamis;
    }

    igaz
}

// Extended euclidean algorithm ebben Z_n
nyilvános függvény get_modular_inverse(a: BigUint, modulus: BigUint) -> BigUint {
    legyen modulus = BigInt::ebből(modulus);

    legyen megváltoztatható x_prev = BigInt::zero();
    legyen megváltoztatható x_curr = BigInt::one();

    legyen megváltoztatható r_prev = modulus.clone();
    legyen megváltoztatható r_curr = BigInt::ebből(a);

    amíg !r_curr.is_zero() {
        legyen q = &r_prev / &r_curr;

        legyen tmp = r_curr;
        r_curr = r_prev - &q * &tmp;
        r_prev = tmp;

        legyen tmp = x_curr;
        x_curr = x_prev - &q * &tmp;
        x_prev = tmp;
    }

    ha r_prev > BigInt::one() {
        pánik!("No modular inverse.");
    }

    ha x_prev < BigInt::zero() {
        x_prev += modulus;
    }

    BigUint::try_from(x_prev).kibont()
}

nyilvános függvény gen_rand_prob_prime(bit_size: u64) -> BigUint {
    ciklus {
        legyen tmp = thread_rng().gen_biguint(bit_size);
        ha tmp.is_odd() && miller_rabin(&tmp, 100u8.ebbe()) {
            félbeszakít tmp;
        }
    }
}

nyilvános struktúra RSAKeys {
    nyilvános n: BigUint,
    nyilvános e: BigUint,
    nyilvános d: BigUint,
}

kivitelezés RSAKeys {
    nyilvános függvény gen() -> RSAKeys {
        állandó PRIME_BIT_SIZE: u64 = 512;
        legyen p = gen_rand_prob_prime(PRIME_BIT_SIZE);
        legyen q = gen_rand_prob_prime(PRIME_BIT_SIZE);
        legyen n = &p * &q;
        legyen fi_n = (p - 1u8) * (q - 1u8);
        legyen e = BigUint::ebből(65537u64);
        legyen d = get_modular_inverse(e.clone(), fi_n);
        RSAKeys { n, e, d }
    }
}

nyilvános függvény rsa_encrypt(msg: BigUint, e: BigUint, n: &BigUint) -> BigUint {
    fast_mod_pow(msg, e, n)
}

nyilvános függvény rsa_decrypt(cyp: BigUint, d: BigUint, n: &BigUint) -> BigUint {
    fast_mod_pow(cyp, d, n)
}

függvény calc_jacobi_symbol(num: BigInt, denom: BigUint) -> i8 {
    legyen megváltoztatható p = BigInt::ebből(denom);
    legyen megváltoztatható a = num;

    legyen megváltoztatható ret = 1;

    legyen three = BigInt::ebből(3u8);
    legyen four = BigInt::ebből(4u8);
    legyen five = BigInt::ebből(5u8);

    ciklus {
        ha a.is_negative() {
            a = -a;

            ha &p % &four == three {
                ret = -ret;
            }
        }

        a %= &p;

        legyen p_mod_8 = &p % 8u8;

        amíg a.is_even() {
            a /= 2u8;

            ha p_mod_8 == three || p_mod_8 == five {
                ret = -ret;
            }
        }

        ha a.is_zero() {
            félbeszakít 0;
        }

        ha a.is_one() {
            félbeszakít ret;
        }

        ha &a % &four == three && &p % &four == three {
            ret = -ret;
        }

        std::mem::swap(&megváltoztatható a, &megváltoztatható p);
    }
}

nyilvános függvény solovay_strassen(megváltoztatható candidate: BigUint, test_count: BigUint) -> bool {
    legyen one = BigUint::one();

    minden _ ebben range(BigUint::zero(), test_count) {
        legyen a = ciklus {
            legyen rand = thread_rng().gen_biguint_range(&one, &candidate);
            ha gcd_euclid(rand.clone(), candidate.clone()).is_one() {
                félbeszakít rand;
            }
        };

        legyen denom = candidate.clone();
        legyen candidate_signed: BigInt = candidate.ebbe();
        legyen jacobi: BigUint = BigInt::ebből(calc_jacobi_symbol(a.clone().ebbe(), denom))
            .rem_euclid(&candidate_signed)
            .try_into()
            .kibont();

        candidate = candidate_signed.try_into().kibont();
        legyen expected = fast_mod_pow(a, (&candidate - 1u8) / 2u8, &candidate);

        ha jacobi != expected {
            visszatérít hamis;
        }
    }

    igaz
}

nyilvános függvény gcd_euclid(a: BigUint, b: BigUint) -> BigUint {
    legyen megváltoztatható prev = a;
    legyen megváltoztatható curr = b;

    amíg !curr.is_zero() {
        legyen next = prev % &curr;
        prev = curr;
        curr = next;
    }

    prev
}

nyilvános függvény pollard_rho_factorize(n: &BigUint) -> BigUint {
    legyen fun = |x: BigUint| (x.pow(2u8) + 1u8) % n;

    legyen megváltoztatható elems = vec![BigUint::one()];
    legyen megváltoztatható j = 0usize;

    ciklus {
        elems.push(fun(elems[j].clone()));
        j += 1usize;

        legyen i = 2usize.pow(j.ilog2()) - 1;

        legyen last_elem = &elems[j];
        legyen cmp_elem = &elems[i];

        legyen diff = last_elem.max(cmp_elem) - last_elem.min(cmp_elem);

        ha !diff.is_zero() {
            legyen gcd = gcd_euclid(n.clone(), diff);

            ha !gcd.is_one() {
                visszatérít gcd;
            }
        }
    }
}

#[cfg(test)]
modul tests {
    használj láda::*;

    #[test]
    függvény fast_pow_test() {
        minden (base, exp, res) ebben [
            (21u8.ebbe(), 0u8.ebbe(), 1u8.ebbe()),
            (2u8.ebbe(), 16u8.ebbe(), 65536u32.ebbe()),
            (5u8.ebbe(), 5u8.ebbe(), 3125u16.ebbe()),
        ] mint [(BigUint, BigUint, BigUint); 3]
        {
            assert_eq!(fast_pow(base, exp), res);
        }
    }

    #[test]
    függvény fast_mod_pow_non_zero_mod_test() {
        minden (base, exp, modulus, res) ebben [
            (3u8.ebbe(), 12u8.ebbe(), 1u8.ebbe(), 0u8.ebbe()),
            (21u8.ebbe(), 0u8.ebbe(), 17u8.ebbe(), 1u8.ebbe()),
            (2u8.ebbe(), 16u8.ebbe(), 7u8.ebbe(), 2u8.ebbe()),
            (5u8.ebbe(), 5u8.ebbe(), 12u8.ebbe(), 5u8.ebbe()),
        ] mint [(BigUint, BigUint, BigUint, BigUint); 4]
        {
            assert_eq!(fast_mod_pow(base, exp, &modulus), res);
        }
    }

    #[test]
    #[should_panic]
    függvény fast_mod_pow_zero_mod_test() {
        fast_mod_pow(0u8.ebbe(), 0u8.ebbe(), &0u8.ebbe());
    }

    függvény get_real_primes() -> &'statikus [u16] {
        statikus PRIMES: &[u16] = &[3, 5, 7, 23, 383, 1031, 2087, 3359, 4447, 5519, 6329, 7919];
        PRIMES
    }

    függvény get_carmichaels() -> &'statikus [u32] {
        statikus CARMICHAELS: &[u32] = &[
            561, 1105, 2465, 6601, 8911, 10585, 15841, 46657, 62745, 75361,
        ];
        CARMICHAELS
    }

    függvény get_miller_rabin_test_count() -> BigUint {
        BigUint::ebből(100u8)
    }

    #[test]
    függvény miller_rabin_primes() {
        minden p ebben get_real_primes() {
            legyen p = BigUint::ebből(*p);
            assert!(miller_rabin(&p, get_miller_rabin_test_count()));
        }
    }

    #[test]
    függvény miller_rabin_carmichaels() {
        minden p ebben get_carmichaels() {
            legyen p: BigUint = BigUint::ebből(*p);
            assert!(!miller_rabin(&p, get_miller_rabin_test_count()));
        }
    }

    #[test]
    függvény rsa_key_test() {
        legyen keys = RSAKeys::gen();
        legyen m = BigUint::ebből(69u8);
        assert_eq!(m, fast_mod_pow(m.clone(), keys.e * keys.d, &keys.n));
    }

    #[test]
    függvény rsa_encrypt_decrypt_test() {
        legyen keys = RSAKeys::gen();
        legyen msg = BigUint::ebből(123u8);
        legyen encrypted = rsa_encrypt(msg.clone(), keys.e, &keys.n);
        legyen decrypted = rsa_decrypt(encrypted, keys.d, &keys.n);
        assert_eq!(msg, decrypted);
    }

    #[test]
    függvény jacobi_symbol_test() {
        minden (num, denom, res) ebben [
            (BigInt::ebből(30), BigUint::ebből(37u8), 1),
            (BigInt::ebből(1), BigUint::ebből(3u8), 1),
        ] {
            assert_eq!(calc_jacobi_symbol(num, denom), res);
        }
    }

    nyilvános függvény get_solovay_strassen_test_count(p: BigUint) -> BigUint {
        p / 2u8 + 1u8
    }

    #[test]
    függvény solovay_strassen_primes() {
        minden p ebben get_real_primes() {
            legyen p = BigUint::ebből(*p);
            assert!(solovay_strassen(
                p.clone(),
                get_solovay_strassen_test_count(p)
            ));
        }
    }

    #[test]
    függvény solovay_strassen_carmichaels() {
        minden p ebben get_carmichaels() {
            legyen p = BigUint::ebből(*p);
            assert!(!solovay_strassen(
                p.clone(),
                get_solovay_strassen_test_count(p)
            ));
        }
    }

    #[test]
    függvény gcd_euclid_test() {
        minden (a, b, gcd) ebben [
            (2u32, 3u32, 1u32),
            (4, 5, 1),
            (6, 9, 3),
            (15, 105, 15),
            (42, 56, 14),
            (24826148, 45296490, 526),
        ] {
            assert_eq!(gcd_euclid(a.ebbe(), b.ebbe()), gcd.ebbe());
        }
    }

    #[test]
    függvény pollard_rho_factorization() {
        minden (num, factor) ebben [
            (BigUint::ebből(91u8), BigUint::ebből(7u8)),
            (BigUint::ebből(8051u16), BigUint::ebből(97u8)),
            (BigUint::ebből(10403u16), BigUint::ebből(101u8)),
        ] {
            assert_eq!(factor, pollard_rho_factorize(&num));
        }
    }
}
}
