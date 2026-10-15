Remark (on arXiv:2408.02591). S. De claims Ramanujan's inequality for all x >= exp(43.5102147). The argument combines
(27) pi(x) < e pi(x/e) (x >= exp(43.5102147)) with (28) pi(x) > x/log x, and concludes in (29) that
  G(x) = pi(x)^2 - (ex/log x) pi(x/e) < pi(x)^2 - pi(x) pi(x) = 0.
This step requires (x/log x) * e pi(x/e) >= pi(x)^2. From (27) one only gets (x/log x) e pi(x/e) > (x/log x) pi(x),
and (x/log x) pi(x) >= pi(x)^2 would need pi(x) <= x/log x, the reverse of (28) (indeed pi(x) > x/log x for x >= 17).
So (27)-(28) do not imply G(x) < 0. This is not a fixable slip: G(x) ~ -c x^2/log^6 x, while (27) and the bounds used
(Panaitopol-type, relative accuracy ~ (log x)^{-3/2}) lose terms of order x^2/log^{3.5} x, far larger than |G|.
The tabulated values of G on [exp(43), exp(3159)] (Sect. 3) cannot come from exact values of pi there and do not
constitute a proof. The best unconditional threshold therefore remains Axler's x >= exp(3158.442).

Remark (on arXiv:2407.12052). The same author claims the inequality for all x >= exp(547). The argument does not stand:
(i) Section 6 (monotonicity of G, Prop. 6.1) assumes Sterneck's conjecture |M(x)| < (1/2) sqrt(x) for x >= 201, which is
false (Neubauer, 1963; the weaker Mertens conjecture |M(x)| < sqrt(x) was disproved by Odlyzko and te Riele, 1985).
(ii) Prop. 2.11(iii) asserts M(x) = O(sqrt x); this is an open problem (widely believed false), and the proof given yields
at most O(x^{1/2+eps}). (iii) Even granting these, Section 4 obtains G(x) = -(2x/log^2 x) sum_{n<=x} M(x/n) + O(x^2/log^2 x)
with the retained term itself only O(x^2/log^2 x), so the sign of G is not determined; the true size of G is ~ x^2/log^6 x.
(iv) The range up to exp(3158.442) is "covered" by evaluating G at 26 sample points e^547, e^647, ..., which is not a proof.
