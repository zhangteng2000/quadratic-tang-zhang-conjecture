#!/usr/bin/env python3
"""Replay the certificate, checking all bounds and its complete domain cover.

Usage: python verify_certificate.py certificate.json
The companion certify.py contains the small, integer-only dyadic arithmetic
kernel. No floating-point value participates in an acceptance decision.
"""
import sys,json,hashlib,time
from collections import defaultdict,Counter
from certify import I,P,S,ONE,check

def fail(message):
    raise RuntimeError(message)

def valid_h(ml,mh,al,ah,h):
    du=(ONE-al*al).high();dl=(ONE-ah*ah).low()
    alpha=(I(mh)*du/2).high()
    elementary=(ONE/(ONE+alpha)).low()
    if h.h<=elementary.l:return True
    if dl.l<=0:return False
    upper=(ONE+du*h).halfpow(mh+2)-al.pow(mh+2).low()-I(ml+2,2)*dl*(ONE+h)
    return upper.h<0

def main(path):
    raw=open(path,'rb').read();doc=json.loads(raw)
    if doc['precision_bits']!=P:fail('Precision mismatch')
    if doc['m_start']!=5 or doc['m_stop_exclusive']!=1000000:fail('Wrong domain')
    cover=defaultdict(list);counts=Counter();start=time.time()
    for idx,rec in enumerate(doc['records']):
        ml,mh,als,ahs,hs,stated,sub=rec
        al,ah,h=I.raw(int(als)),I.raw(int(ahs)),I.raw(int(hs))
        if not (5<=ml<=mh<1000000 and 0<=al.l<ah.l<=S and 0<=h.l<=S):
            fail(f'Invalid box {idx}')
        if not valid_h(ml,mh,al,ah,h):fail(f'Uncertified polar bound {idx}')
        actual,_,_=check(ml,mh,al,ah,h,sub)
        if actual is None:fail(f'Failed inequality {idx}')
        cover[(ml,mh)].append((al.l,ah.l))
        counts[actual]+=1
    previous=4
    for (ml,mh),intervals in sorted(cover.items()):
        if ml!=previous+1:fail(f'Degree gap or overlap before {ml}')
        previous=mh
        endpoint=0
        for l,h in sorted(intervals):
            if l!=endpoint:fail(f'Parameter gap or overlap for degree box {ml,mh}')
            endpoint=h
        if endpoint!=S:fail(f'Missing a=1 endpoint for degree box {ml,mh}')
    if previous!=999999:fail('Missing final degree range')
    report={'status':'PASS','arithmetic':'100-bit outward-rounded dyadic integer intervals',
            'm_range':[5,999999],'a_range':[0,1],
            'degree_blocks':len(cover),'parameter_boxes':len(doc['records']),
            'criteria':dict(counts),'certificate_sha256':hashlib.sha256(raw).hexdigest(),
            'verification_seconds':round(time.time()-start,3)}
    print(json.dumps(report,indent=2))
    return report
if __name__=='__main__':main(sys.argv[1] if len(sys.argv)>1 else 'certificate.json')
