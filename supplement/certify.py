#!/usr/bin/env python3
"""Outward-rounded dyadic interval certificate for the quadratic refinement.

All certification arithmetic below is integer arithmetic. Floating point is used
only to propose a bisection bracket and to display progress, never for acceptance.
"""
from __future__ import annotations
from math import isqrt, log, ceil
from dataclasses import dataclass
import json, time, sys, argparse

P = 100
S = 1 << P

def cd(a: int, b: int) -> int:
    return -((-a)//b)

class I:
    __slots__=('l','h')
    def __init__(self, a=0, b=1):
        self.l=(int(a)*S)//int(b); self.h=cd(int(a)*S,int(b))
    @classmethod
    def raw(cls,l,h=None):
        x=object.__new__(cls);x.l=int(l);x.h=int(l if h is None else h);return x
    def __add__(self,other):
        if not isinstance(other,I):other=I(other)
        return I.raw(self.l+other.l,self.h+other.h)
    __radd__=__add__
    def __neg__(self):return I.raw(-self.h,-self.l)
    def __sub__(self,other):return self+-cv(other)
    def __rsub__(self,other):return cv(other)+-self
    def __mul__(self,other):
        other=cv(other)
        if self.l>=0 and other.l>=0:
            return I.raw(self.l*other.l//S,cd(self.h*other.h,S))
        z=(self.l*other.l,self.l*other.h,self.h*other.l,self.h*other.h)
        return I.raw(min(z)//S,cd(max(z),S))
    __rmul__=__mul__
    def __truediv__(self,other):
        other=cv(other)
        if other.l<=0:raise ZeroDivisionError((other.l,other.h))
        if self.l>=0:return I.raw(self.l*S//other.h,cd(self.h*S,other.l))
        z=(self.l*S//other.l,self.l*S//other.h,self.h*S//other.l,self.h*S//other.h)
        # Generic version includes a unit of outward slack at the upper end.
        return I.raw(min(z),max(z)+1)
    def __rtruediv__(self,other):return cv(other)/self
    def pow(self,k):
        if k<0: return I(1)/self.pow(-k)
        ans=I(1);a=self
        while k:
            if k&1:ans=ans*a
            k>>=1
            if k:a=a*a
        return ans
    def sqrt(self):
        assert self.l>=0
        lo=isqrt(self.l*S); hh=isqrt(self.h*S)
        return I.raw(lo,hh+(hh*hh<self.h*S))
    def halfpow(self,k):
        """Outward enclosure of self**(k/2), k a nonnegative integer."""
        ans=self.pow(k//2)
        return ans*self.sqrt() if k&1 else ans
    def low(self):return I.raw(self.l)
    def high(self):return I.raw(self.h)
    def __repr__(self):return f'[{self.l/S:.12g},{self.h/S:.12g}]'

def cv(x):return x if isinstance(x,I) else I(x)
ONE=I(1);GAM=I(5,3)

_grid_cache={}
def grid(mhi,sub=2):
    key=(mhi,sub)
    if key in _grid_cache:return _grid_cache[key]
    K=mhi.bit_length()+3
    breaks=[0]+[S>>k for k in range(K,0,-1)]+[S]
    points={0,S}
    for l,h in zip(breaks,breaks[1:]):
        for j in range(sub+1):
            t=l+(h-l)*j//sub
            points.add(t);points.add(S-t)
    ts=sorted(points)
    ti=[I.raw(t) for t in ts]
    weights=[]
    for l,h in zip(ti,ti[1:]):
        z=h-l
        w10=z*(l/2+z/6);w11=z*(l/2+z/3)
        w20=z*(l*l/2+l*z/3+z*z/12)
        w21=z*(l*l/2+2*l*z/3+z*z/4)
        weights.append((w10,w11,w20,w21))
    _grid_cache[key]=(ti,weights)
    return ti,weights


def polar_h(ml,mh,al,ah,improve=True):
    """A certified common lower bound on 2*a*x-a*a over the box."""
    dl=(ONE-ah*ah).low();du=(ONE-al*al).high()
    alpha=(I(mh)*du/2).high()
    h=(1/(1+alpha)).low()
    if not improve or dl.l<=0:return h
    nl=I(ml+2,2);nu=I(mh+2,2)
    apow=al.pow(mh+2).low()
    def f(z):
        return (ONE+du*z).halfpow(mh+2)-apow-nl*dl*(ONE+z)
    if f(I(0)).h>=0:return h
    # Bracket proposal only; all updates are certified by exact interval tests.
    av=(nu*du).h/S
    hh=min(S,max(h.h+1,int(S*min(1.0,2*log(max(2.0,2*av))/max(av,1e-99)))))
    lo=0;hi=hh
    while hi<S and f(I.raw(hi)).h<0:
        lo=hi;hi=min(S,2*hi)
    if f(I.raw(hi)).h<0:return I.raw(max(h.l,hi))
    for _ in range(34):
        mid=(lo+hi)//2
        if mid==lo:break
        if f(I.raw(mid)).h<0:lo=mid
        else:hi=mid
    return I.raw(max(h.l,lo))


def integrals(ml,mh,al,ah,h,sub=2):
    cl=(al*al).low();cu=(ah*ah).high();d=(cu+h)/2
    if d.h>=S:return None
    ts,ww=grid(mh,sub)
    bs=[];phis=[]
    for t in ts:
        b=ONE-h*t-cl*t*(ONE-t)
        b=I.raw(max(0,b.l),max(0,b.h))
        q=b.halfpow(ml-1)
        den=ONE-d*t
        if den.l<=0:return None
        bs.append(q.high());phis.append((q/den).high())
    k1=I(0);k2=I(0)
    for j,(w10,w11,w20,w21) in enumerate(ww):
        k1=k1+w10*bs[j]+w11*bs[j+1]
        k2=k2+w20*phis[j]+w21*phis[j+1]
    return k1.high(),k2.high()


def check(ml,mh,al,ah,h=None,sub=2):
    assert 5<=ml<=mh and 0<=al.l<=ah.h<=S
    if h is None:h=polar_h(ml,mh,al,ah)
    # Since a^2+h<=2ax<2a, this condition excludes the whole box.
    if h.l>(2*ah-ah*ah).h:return 'polar',h,{}
    vals=integrals(ml,mh,al,ah,h,sub)
    if vals is None:return None,h,{}
    k1,k2=vals
    b=(ONE-h).high()
    F=b.halfpow(ml)
    R=GAM*ah*ah*mh*k1
    direct=F+R+I(mh,mh+1)*ah
    if direct.h<S:return 'direct',h,{'slack':(S-direct.h)/S}
    E=b.halfpow(ml+1).high()
    C=(GAM*ah.pow(3)*mh*(mh+1)/2*k2).high()
    # Uniform treatment of the otherwise singular boundary a=1.
    if 2*C.h<=al.l:
        alpha=(I(mh)*(ONE-al*al)/2).high()
        if 0<alpha.l and alpha.h<=I(ml-1,2).l:
            ratio=(I(mh)/(2*alpha)*(alpha/(1+alpha)).halfpow(ml+1)).high()
            if ratio.h<I(1,4).l:
                return 'boundary',h,{'C':C.h/S,'ratio':ratio.h/S}
    if 2*C.h<=ah.l:
        if (ah+E).h<S:return 'center-monotone',h,{'slack':(S-(ah+E).h)/S}
    elif E.h<S:
        lhs=4*(ah+C).pow(3)
        rhs=27*C*(ONE-E).pow(2)
        if lhs.h<rhs.l:return 'center-cubic',h,{'slack':(rhs.l-lhs.h)/S}
    return None,h,{'direct':direct.h/S,'C':C.h/S,'E':E.h/S}


def mblocks(stop):
    lo=5
    while lo<stop:
        hi=lo if lo<45 else min(stop-1,lo+max(1,lo//40))
        yield lo,hi
        lo=hi+1


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stop',type=int,default=1000000)
    ap.add_argument('--out',default='certificate.json')
    ap.add_argument('--max-depth',type=int,default=30)
    ap.add_argument('--start',type=int,default=5)
    args=ap.parse_args()
    records=[];counts={};start=time.time();tests=0
    for ml,mh in mblocks(args.stop):
        if mh<args.start:continue
        stack=[(0,S,0)]
        block=0
        while stack:
            al,ah,dep=stack.pop();tests+=1
            crit,h,info=check(ml,mh,I.raw(al),I.raw(ah))
            if crit is None and dep>=args.max_depth:
                for sub in [4,8,16]:
                    crit,h,info=check(ml,mh,I.raw(al),I.raw(ah),h,sub)
                    if crit is not None:break
                if crit is None:
                    print('FAILED',ml,mh,al/S,ah/S,dep,info,flush=True)
                    return 1
            else:sub=2
            if crit is not None:
                records.append([ml,mh,str(al),str(ah),str(h.l),crit,sub])
                counts[crit]=counts.get(crit,0)+1;block+=1
            else:
                mid=(al+ah)//2
                stack.append((mid,ah,dep+1));stack.append((al,mid,dep+1))
        print('PASS',ml,mh,'boxes',block,'total',len(records),'tests',tests,'sec',round(time.time()-start,1),flush=True)
        payload={'precision_bits':P,'m_start':args.start,'m_stop_exclusive':args.stop,'records':records,'counts':counts,'tests':tests}
        with open(args.out,'w') as f:json.dump(payload,f,separators=(',',':'))
    print('ALL PASS',counts,'boxes',len(records),'seconds',time.time()-start,flush=True)
    return 0

if __name__=='__main__':sys.exit(main())
