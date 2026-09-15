#!/usr/bin/env python3
"""Exact-rational checks for the numerical ledger in the analytic tail proof.
The monotonicity arguments extending these checks to m >= 10^6 are in the proof.
"""
from fractions import Fraction as Q
from math import factorial
import json

M=10**6
checks=[]
def test(name,condition):
    if not condition:raise AssertionError(name)
    checks.append(name)

elo=Q(27,10);ehi=Q(11,4)
s5=sum((Q(1,factorial(k)) for k in range(6)),Q(0))
test('27/10 < e < 11/4 via the exponential series',elo<s5 and s5+Q(7,4320)<ehi)
test('sqrt(e) < 5/3',ehi<Q(25,9))
test('13 < log(10^6) < 14',ehi**13<M<elo**14)
test('log(3*10^6/8) > 12',ehi**12<Q(3*M,8))
test('log(10^6/20) > 4',ehi**4<Q(M,20))
test('decreasing-interval coefficient < 7/m',Q(20,3)*Q(M,M-1)**2<7)
near_constant=Q(2048*5,3*81)*Q(10,9)**3*Q(M*M*(M+1),(M-1)**3)
test('near-boundary integral coefficient < 60/m',near_constant<60)
test('near-boundary exponent dominates 4 log m',Q(1000,6)>56 and Q(9*M,160)>56)
test('monotonicity thresholds in the two exponential estimates',1000>48 and Q(9,160)>Q(4,M))
test('polar logarithmic exponent exceeds nine',10*(1-Q(3,M))>9)
test('C < 61/m < a/2 in the near region',Q(60,M)+Q(7,M*M)<Q(61,M)<Q(9,20))
test('small-alpha endpoint estimate at m=16; sequence decreases afterwards',16*16**2<2**17 and Q(17,16)**2<2)
test('large-alpha endpoint estimate',Q(1,2*M**3)<Q(1,4))
test('finite-degree correction c >= 99999/100000',Q(M-1,M+2)>Q(99999,100000))
test('m^(4/(m+2)) < 1001/1000',1/(1-Q(56,M+2))<Q(1001,1000))
test('(8/3)^(4/3) < 4',Q(8,3)**4<4**3)
test('small-a increasing-interval coefficient < 5/(2 log m)',4*Q(5,3)*Q(1001,1000)/elo<Q(5,2))
test('small-a decreasing-interval bound < .292',7*Q(M+2,2*M)/12<Q(292,1000))
test('small-a total error < 1/2',Q(292,1000)+Q(5,26)+Q(1,10000)<Q(1,2))
test('small-a endpoint exponent and base',Q(M,M+2)>Q(4,5) and Q(3*M,8)>10**5)
test('monotonicity in delta of (1-delta)*(m*delta/2)^(-c/delta)',Q(9,10)*Q(16,9)*3>4)
test('cube-root bound for A0=375000',72**3<Q(3*M,8))
test('A0^(1/75000) < 1001/1000',1/(1-Q(14,75000))<Q(1001,1000))
test('large-a increasing-interval error < .016',Q(5,12)*Q(1,27)*Q(1001,1000)<Q(16,1000))
test('large-a total error < .017 < .05',Q(112,M)+Q(8,3*M)+Q(16,1000)<Q(17,1000)<Q(1,20))
test('1-sqrt(.9) > .05',Q(19,20)**2>Q(9,10))
print(json.dumps({'status':'PASS','checks':len(checks),'m0':M,'verified':checks},indent=2))
