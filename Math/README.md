# Math

Every file is self-contained: paste the struct/function you need into your
template and it compiles. Each one opens with a `WHAT / WHEN`, `API`,
`COMPLEXITY` and `PITFALLS` block — read the pitfalls, they're where the
wrong-answer verdicts live.

## Index

| File | Contents |
|---|---|
| `basics.cpp` | `fdiv`/`cdiv`, `pmod`, exact `isqrt`/`icbrt`, `mulmod`/`powmod`, `fib`, series, `floorSum`, digits, `josephus` |
| `modular.cpp` | `extgcd`, `modinv`, `invTable`, `modLinear`, `diophantine`, `crt`, `eulerPhi`, `discreteLog`, `sqrtMod`, `primitiveRoot`, `discreteRoot` |
| `modint.cpp` | `Mint<MOD>` — an int that reduces itself; `DMint` for a runtime modulus |
| `combinatorics.cpp` | `Comb` (nCr, Catalan, stars & bars, derangements, surjections), `Lucas`, `Cbig`, `pascal`, Stirling, `partitions` |
| `sieve.cpp` | `Sieve` (spf, primes, phi, mobius, divisors), `segmentedSieve`, `countPrimes`, `DivTable` (all divisors of every i), divisor tables, divisor/multiple transforms |
| `factor.cpp` | Miller-Rabin + Pollard rho: `isPrimeMR`, `factor`, `phi`; plus `divisorsSqrt` (all divisors of one n, O(sqrt n)) |
| `matrix.cpp` | `Matrix` mod p: multiply, `pow`, `det`, `rank`, `inverse`, `solve`, `linearRecMatrix` |
| `gauss.cpp` | `solveLinear` (reals), `detReal`, `solveLinearBinary` (GF(2) bitsets) |
| `xorbasis.cpp` | `XorBasis`: max/min/k-th XOR, membership, count |
| `fft.cpp` | `conv` (double), `convNTT` (mod 998244353), `convMod<M>`, `bitConv` (xor/or/and) |
| `linrec.cpp` | `berlekampMassey` + `linearRec` (Kitamasa) |
| `frac.cpp` | `Frac` exact rationals, `contFrac`, `bestApprox` |
| `bigint.cpp` | `BigInt` arbitrary precision, base 1e9 |
| `game.cpp` | `mex`, `Grundy`, `subtractionGame`, nim / misere nim |
| `theorems.md` | formula sheet: identities, theorems, and the standard counting arguments |

## Which file do I want?

**"answer mod 1e9+7"** → `modint.cpp` for the arithmetic, `combinatorics.cpp`
if binomials appear. Don't hand-roll `% MOD`.

**"all divisors of n (not just prime)"**
- one n up to ~1e14 → `divisorsSqrt(n)` in `factor.cpp` — pairs `(i, n/i)` around sqrt
- one n beyond that → `divisors(n)` in `factor.cpp` (factorise, then expand)
- many x, all ≤ 1e7 → `Sieve::divisors(x)` in `sieve.cpp`
- every i ≤ N at once → `DivTable` in `sieve.cpp` — sieve them, never loop per-number

**"is it prime / factor it"**
- many numbers, all ≤ 1e7 → `Sieve` in `sieve.cpp` (spf makes each query O(log n))
- one number up to 1e18 → `isPrimeMR` / `factor` in `factor.cpp`
- primes in a range `[L, R]` with R up to 1e12 → `segmentedSieve`

**"n is up to 1e18 and there's a recurrence"**
- you know the recurrence, order ≤ 10 → `linearRecMatrix` in `matrix.cpp`
- you know it, order is large → `linearRec` in `linrec.cpp`
- you *don't* know it → brute-force ~2d terms, then `berlekampMassey`

**"solve these equations"**
- mod a prime → `Matrix::solve`
- real numbers (expected value with cycles, physics) → `solveLinear` in `gauss.cpp`
- XOR equations → `solveLinearBinary`, or `XorBasis` if you only need the span

**"for every k, count pairs summing to k"** → `fft.cpp`. Read the precision note
before using the `double` version.

**"count things"** → `combinatorics.cpp`, and `theorems.md` for the argument
(inclusion-exclusion, Burnside, stars & bars, Catalan).

**"it overflows"** → `__int128` first, then mod arithmetic, then `frac.cpp` if
it's a ratio, and only then `bigint.cpp`.

## Conventions

- `#define int long long` is on in most files, so `int` means 64-bit.
  **Exception:** the big arrays in `sieve.cpp` are explicitly `int32_t`/`int8_t`
  so a 1e7 sieve doesn't MLE. Don't "clean that up".
- `fft.cpp`, `linrec.cpp`, `gauss.cpp` and `bigint.cpp` use explicit `ll`/`double`
  instead, so they behave the same whether or not the macro is present.
- Use `std::gcd` / `std::lcm` (C++17), **not** `__gcd` — Apple's libc++
  static-asserts on signed types, so `__gcd` won't compile locally.
- Ranges in these files are half-open or inclusive as documented per function;
  `Sieve` and `Comb` are sized by the largest value you will query, not by n.

Everything here was checked against brute force (including `Comb` vs Pascal's
triangle, `Sieve` vs trial division, `XorBasis` vs all 2^n subsets, `det` under
both prime and composite moduli, and `BigInt` against `long long`).
