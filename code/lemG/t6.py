from resid import *
mp.mp.dps=20
def Ix(w,zi,zj,c=mp.mpf('0.3')):
    f=lambda t: (lambda s2: zi**(-w-s2)*zj**s2*mp.zeta(1+w)*A(w,-w-s2,s2)/(mp.zeta(1-s2)*mp.zeta(1+w+s2)*(w+s2)**2*s2**2))(c+1j*t)
    return mp.re(mp.quad(f,[-mp.inf,-20,-3,-0.5,0,0.5,3,20,mp.inf]))/(2*mp.pi)
for zi,zj in ((60,60),(60,200),(300,300)):
    w=mp.mpf('0.05')
    li,lj=mp.log(zi),mp.log(zj)
    A12=mp.diff(lambda x,y: A(w,x,y),(0,0),(1,1)); dzl=mp.diff(lambda s: zl(s),1+w)
    main=li*lj+dzl+A12-zj**(-w)/w**2+Ix(w,zi,zj)
    ex=exact(w,zi,zj)
    print(zi,zj,"exact",float(ex),"main(with Ix)",float(main),"diff",float(ex-main),flush=True)
