# Closed-form Laplace transform of Heath-Brown's w (BTY Definition 1), no quadrature:
#   w(u) = s2[ s2(th ct - u/2) cos(k u) + 2 th ct - u + sin(2th - k u)/sin(2th) - 2 - 2 sin(th - k u)/sin(th) ],  k = tan th.
import cmath, math
th = 1.1338; k = math.tan(th); ct = 1/k; s2 = 1/math.cos(th)**2; D = 2*th*ct
def J0(s): return D if s == 0 else (1 - cmath.exp(-s*D))/s
def J1(s): return D*D/2 if s == 0 else (1 - cmath.exp(-s*D)*(1 + s*D))/(s*s)
def W(z):
    z = complex(z); ik = 1j*k
    al, be = th*ct, -0.5                                  # (th ct - u/2) = al + be u
    cosp = 0.5*(al*J0(z-ik) + be*J1(z-ik) + al*J0(z+ik) + be*J1(z+ik))
    lin = (2*th*ct - 2)*J0(z) - J1(z)
    def sinpart(phi): return (cmath.exp(1j*phi)*J0(z+ik) - cmath.exp(-1j*phi)*J0(z-ik))/(2j)
    return s2*(s2*cosp + lin + sinpart(2*th)/math.sin(2*th) - 2*sinpart(th)/math.sin(th))
if __name__ == "__main__":
    import mpmath as mp
    import bty2
    for z in [1.0, -1.0, 0.3+2j, 140+5j, 830-40j, -0.9+0.05j, 3+80j]:
        ref = mp.quad(lambda u: mp.e**(-mp.mpc(z)*u)*bty2.w(u), mp.linspace(0, bty2.D, 17))
        print(z, W(z), complex(ref), abs(W(z)-complex(ref)))
