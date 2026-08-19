# Formula sheet

The results worth recognising mid-contest. Code for most of these is in the
neighbouring files; this is the "what am I even looking at" reference.

---

## Modular arithmetic

**Fermat.** `p` prime, `gcd(a, p) = 1` → `a^(p-1) = 1 (mod p)`, so
`a^-1 = a^(p-2)`.

**Euler.** `gcd(a, m) = 1` → `a^phi(m) = 1 (mod m)`, so `a^-1 = a^(phi(m)-1)`.

**Reducing exponents.** `a^e = a^(e mod phi(m)) (mod m)` — **only** when
`gcd(a, m) = 1`. Otherwise, for `e >= log2(m)`:
`a^e = a^(phi(m) + (e mod phi(m))) (mod m)`.

**Wilson.** `n` is prime ⟺ `(n-1)! = -1 (mod n)`. Beautiful, and useless
computationally (it's O(n) to test one number) — but it is why `(p-1)!` shows up
in Lucas-style derivations.

**CRT.** Congruences `x = r_i (mod m_i)` with pairwise coprime `m_i` have exactly
one solution mod `prod m_i`. With non-coprime moduli a solution exists iff
`r_i = r_j (mod gcd(m_i, m_j))` for every pair — `crt()` handles this.

**Lifting the exponent (LTE).** `p` odd prime, `p | a-b`, `p ∤ a, b`:
`v_p(a^n - b^n) = v_p(a-b) + v_p(n)`, where `v_p` is the exponent of `p`.
For `p = 2` and even `n`: `v_2(a^n - b^n) = v_2(a-b) + v_2(a+b) + v_2(n) - 1`.
Use it for "how many trailing zeros / what power of p divides this".

**Quadratic residues.** `a` is a QR mod odd prime `p` ⟺ `a^((p-1)/2) = 1`.
Exactly `(p-1)/2` of the non-zero residues are QRs. Square roots come from
`sqrtMod` (Tonelli-Shanks); when `p % 4 == 3` it is just `a^((p+1)/4)`.

---

## Number theory

**Divisor function.** With `n = prod p_i^e_i`:
- `d(n) = prod (e_i + 1)` — number of divisors
- `sigma(n) = prod (p_i^(e_i+1) - 1) / (p_i - 1)` — sum of divisors
- `phi(n) = n * prod (1 - 1/p_i)` — count of `x <= n` coprime to `n`

**Bounds.** `d(n) <= 1344` for `n <= 1e9`; `<= 103680` for `n <= 1e18`.
`pi(1e6) = 78498`, `pi(1e9) = 50847534`. Sum of `d(i)` for `i <= n` is
`~ n ln n`, which is why the harmonic-loop divisor table is O(n log n).

**Totient identities.** `sum_{d | n} phi(d) = n`. `phi` is multiplicative.
`phi(p^k) = p^k - p^(k-1)`.

**Mobius.** `mu(n) = 0` if `n` has a squared prime factor, else `(-1)^(#distinct primes)`.
`sum_{d | n} mu(d) = [n == 1]` — the indicator that makes inclusion-exclusion
over divisors work.

**Mobius inversion.** `g(n) = sum_{d|n} f(d)` ⟺ `f(n) = sum_{d|n} mu(d) g(n/d)`.
The "count pairs with gcd exactly k" recipe is this in disguise: count pairs
with gcd *divisible by* k first (easy), then invert.

**Counting by gcd.** `# pairs (i,j) with gcd = 1` in `[1,n]^2` is
`sum_d mu(d) floor(n/d)^2`; pairs with `gcd = k` is the same with `n -> n/k`.
`floor(n/d)` takes only O(sqrt n) distinct values, but block-walking a *weighted*
sum needs prefix sums of the weight — for `mu` that is Mertens. See
`CoprimeCounter` in `sieve.cpp` (O(n) build, O(sqrt n) per query);
`divisorBlockSum` is the unweighted skeleton.

**Bezout.** `ax + by = c` has integer solutions ⟺ `gcd(a,b) | c`. The
solutions are `x + k(b/g)`, `y - k(a/g)`. The largest integer NOT representable
as a non-negative combination of coprime `a, b` is `ab - a - b` (Chicken McNugget).

**Fibonacci.** `gcd(F(a), F(b)) = F(gcd(a,b))`. `F(n) | F(kn)`.
`sum_{i<=n} F(i) = F(n+2) - 1`. `F(n)` is even ⟺ `3 | n`.
Zeckendorf: every positive integer is a unique sum of non-consecutive Fibonaccis.

---

## Combinatorics

**Binomial identities.**
```
C(n,k) = C(n-1,k-1) + C(n-1,k)                Pascal
sum_k C(n,k) = 2^n            sum_k (-1)^k C(n,k) = [n == 0]
sum_k C(n,k)^2 = C(2n,n)
sum_{i<=n} C(i,k) = C(n+1,k+1)                hockey stick
sum_i C(a,i) C(b,k-i) = C(a+b,k)              Vandermonde
C(n,k) C(k,j) = C(n,j) C(n-j,k-j)             subset of a subset
k C(n,k) = n C(n-1,k-1)                       absorption
```

**Stars and bars.** `k` identical items into `b` ordered boxes:
`C(k+b-1, b-1)` allowing empties, `C(k-1, b-1)` if every box must be non-empty.
With per-box caps, inclusion–exclude over which boxes overflow.

**Catalan.** `C_n = C(2n,n)/(n+1) = C(2n,n) - C(2n,n+1)`; `1, 1, 2, 5, 14, 42, ...`
Counts: balanced bracket sequences of length `2n`; binary trees with `n` nodes;
triangulations of an `(n+2)`-gon; monotone lattice paths staying weakly below the
diagonal; ways to stack-sort a permutation.

**Ballot / reflection.** Paths from `(0,0)` to `(a,b)` that touch a forbidden
line are counted by reflecting the start across that line — subtract
`C(a+b, a')` for the reflected endpoint. This is where Catalan comes from, and
it generalises to "paths avoiding a diagonal band".

**Derangements.** `D(n) = (n-1)(D(n-1) + D(n-2))`, `D(0)=1, D(1)=0`;
equivalently `n! sum (-1)^k / k!` = `round(n!/e)`.

**Stirling.** `S2(n,k)` = partitions of an `n`-set into `k` non-empty unlabelled
blocks, `S2(n,k) = k S2(n-1,k) + S2(n-1,k-1)`. `S1(n,k)` = permutations of `n`
with `k` cycles. `Bell(n) = sum_k S2(n,k)`.
Surjections from `n` onto `k` labelled boxes = `k! S2(n,k)`.

**Trees.** Labelled trees on `n` vertices: `n^(n-2)` (Cayley). With prescribed
degrees `d_i`: `(n-2)! / prod (d_i - 1)!` (Prüfer sequences).
Spanning trees of any graph: any cofactor of the Laplacian (Matrix-Tree,
`matrix.cpp`).

**Inclusion-exclusion.**
`|A_1 ∪ ... ∪ A_n| = sum |A_i| - sum |A_i ∩ A_j| + ...`
With `n <= 20` properties, iterate over bitmasks weighting by `(-1)^popcount`.

**Burnside / Polya.** Distinct objects under a symmetry group `G`:
`(1/|G|) sum_{g in G} fix(g)`.
- necklaces, `n` beads, `k` colours, rotations only: `(1/n) sum_{d|n} phi(n/d) k^d`
- bracelets (add reflections), `n` odd: `+ n k^((n+1)/2)`, all over `2n`
- `n` even: `+ (n/2)(k^(n/2) + k^(n/2+1))`, all over `2n`

**Pigeonhole, extremal.** `n+1` items in `n` boxes force a repeat. Among any
`n+1` numbers from `[1, 2n]` two are coprime, and two divide each other.

---

## Geometry-adjacent

**Pick's theorem.** A simple lattice polygon with area `A`, `I` interior lattice
points and `B` boundary points: `A = I - B/2 + 1`.
Lattice points on the segment `(x1,y1)-(x2,y2)`: `gcd(|dx|, |dy|) + 1` endpoints
included.

**Shoelace.** `2A = |sum (x_i y_{i+1} - x_{i+1} y_i)|`. Keep it doubled and
integral — never divide by 2 until the end.

---

## Games

**Sprague-Grundy.** Any impartial normal-play game = a nim pile of size
`g(state) = mex{ g(next) }`. Independent sub-games XOR together. Losing (for the
player to move) ⟺ `g = 0`.

**Nim.** Win ⟺ XOR of pile sizes ≠ 0. Misère nim: if some pile ≥ 2, the normal
rule holds; if all piles are 1, win ⟺ the number of piles is even.

**Staircase nim.** Tokens on steps, move any number down one step: XOR only the
odd-indexed steps.

**Wythoff.** Two piles, take from one or equally from both. Losing positions are
`(floor(k φ), floor(k φ²))` for `k >= 0`, φ the golden ratio.

---

## Estimation

```
n!    ~ sqrt(2 pi n) (n/e)^n            log10(1e6!) ~ 5.5e6
2^10  ~ 1e3      2^20 ~ 1e6      2^30 ~ 1e9      2^63 ~ 9.2e18
20! ~ 2.4e18 (last one in int64)        13! ~ 6.2e9 (first past int32)
harmonic sum_{i<=n} 1/i ~ ln n + 0.577
sum_{i<=n} n/i ~ n ln n                 (divisor-table cost)
```
