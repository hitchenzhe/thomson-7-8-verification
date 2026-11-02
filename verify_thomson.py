#!/usr/bin/env python3
"""Independent checks for Thomson-26-10.30.tex. See README-verification.md.

Exact algebra and finite enumeration are distinguished from numerical evidence.
No SDP dual certificate is supplied or asserted by this program.
"""
import argparse
import hashlib
import itertools as it
import json
import math
import platform
from pathlib import Path
import time
from fractions import Fraction

import mpmath as mp
import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.special import eval_legendre
import sympy as sp
from sympy.polys.matrices import DomainMatrix

CHECKS = []
DETAILS = {}
EXPECTED_SOURCE_SHA256 = 'd6ea2dd5c0c0e6d62b12449a1b6cee95dfbbbbaeb6e41cd116909c4b6804b11c'


def check(name, condition, evidence=None, kind="exact"):
    ok = bool(condition)
    CHECKS.append(dict(name=name, passed=ok, kind=kind, evidence=evidence))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}", flush=True)
    if not ok:
        raise AssertionError(f"{name}: {evidence}")


def zero(expr):
    return sp.simplify(sp.expand(expr)) == 0


def num(expr):
    return mp.mpf(str(sp.N(expr, mp.mp.dps)))


def rational_enclosure(expr, digits=65):
    """Rigorous rational interval for radical expressions; no float arithmetic."""
    expr=sp.sympify(expr)
    if expr.is_Rational:
        x=Fraction(int(expr.p),int(expr.q))
        return x,x
    def multiply(A,B):
        products=[x*y for x in A for y in B]
        return min(products),max(products)
    if expr.is_Add:
        intervals=[rational_enclosure(x,digits) for x in expr.args]
        return sum(x[0] for x in intervals),sum(x[1] for x in intervals)
    if expr.is_Mul:
        result=(Fraction(1),Fraction(1))
        for x in expr.args:
            result=multiply(result,rational_enclosure(x,digits))
        return result
    if expr.is_Pow:
        base,power=expr.args
        A=rational_enclosure(base,digits)
        if power in (sp.Rational(1,2),sp.Rational(-1,2)):
            if A[0]<0:
                raise ValueError('negative radicand enclosure')
            scale=10**digits
            def floor_sqrt(x):
                return math.isqrt(x.numerator*scale*scale//x.denominator)
            result=(Fraction(floor_sqrt(A[0]),scale),Fraction(floor_sqrt(A[1])+1,scale))
            if power<0:
                if result[0]<=0:
                    raise ValueError('reciprocal interval contains zero')
                result=(1/result[1],1/result[0])
            return result
        if power.is_Integer:
            exponent=int(power)
            if exponent<0:
                if A[0]<=0<=A[1]:
                    raise ValueError('reciprocal interval contains zero')
                A=(1/A[1],1/A[0]); exponent=-exponent
            result=(Fraction(1),Fraction(1))
            for _ in range(exponent):
                result=multiply(result,A)
            return result
    raise ValueError(f'Unsupported exact interval expression: {expr}')


def n7_exact():
    s = sp.sqrt(5)
    a, b = (s-1)/4, -(s+1)/4
    u, z = sp.symbols("u z")
    S = [-sp.S.One, sp.S.Zero, a, b]
    check("N7 basic contact identities", all(zero(x) for x in
          [a+b+sp.Rational(1,2), a*b+sp.Rational(1,4), 2*b*b-1-a]))
    check("N7 S intersect (-S) is {0}",
          [x for x in S if any(zero(x+y) for y in S)] == [0])
    rows = []
    for counts in it.product(range(7), repeat=4):
        m,n0,na,nb = counts
        if sum(counts)==6 and 4*(1-m)-na-nb==0 and na==nb:
            rows.append(counts)
    check("N7 exhaustive integer row multiplicities", rows==[(0,2,2,2),(1,5,0,0)], rows,
          "exact finite enumeration")
    det = sp.Matrix([[1,b,b],[b,1,u],[b,u,1]]).det()
    check("N7 3x3 Gram determinant factorization", zero(det-(1-u)*(1+u-2*b*b)))
    feasible = [x for x in S if sp.simplify(det.subs(u,x)).is_nonnegative]
    check("N7 PSD forces mutual neighbor contact a", feasible==[a])
    check("N7 vector recurrence has zero squared residual", zero(2+2*a-4*b*b))
    C = sp.Matrix([[0,1],[-1,2*b]])
    periods = [n for n in range(3,8) if all(zero(x) for x in C**n-sp.eye(2))]
    check("N7 recurrence closes only at length 5 among 3..7", periods==[5], periods)
    # Enumerate every labelled graph with all degrees 0 or 2, without graph libraries.
    def cycle_graphs(vertices):
        if not vertices:
            yield (), ()
            return
        first, *rest = vertices
        for cycles,isolated in cycle_graphs(rest):
            yield cycles, (first,)+isolated
        for length in range(3,len(vertices)+1):
            for subset in it.combinations(rest,length-1):
                remaining = [x for x in rest if x not in subset]
                for order in it.permutations(subset):
                    if order[0] > order[-1]:
                        continue
                    for cycles,isolated in cycle_graphs(remaining):
                        yield ((first,)+order,)+cycles, isolated
    total, survivors = 0, []
    for cycles,isolated in cycle_graphs(list(range(7))):
        total += 1
        # Isolates have unique antipodes, and every nontrivial cycle must close.
        if len(isolated)%2 or any(len(c) not in periods for c in cycles):
            continue
        survivors.append((cycles,isolated))
    check("N7 full degree-0/2 graph enumeration leaves C5 plus two poles",
          len(survivors)==252 and all(len(c)==1 and len(i)==2 for c,i in survivors),
          dict(graphs_enumerated=total, surviving_labelled_graphs=len(survivors)),
          "exact finite enumeration")
    G = sp.eye(7)
    for i,j in it.combinations(range(5),2):
        G[i,j]=G[j,i]= b if (i-j)%5 in (1,4) else a
    G[5,6]=G[6,5]=-1
    # A fixed permutation identifies all 252 survivors with this matrix.
    poly = sp.Poly(G.charpoly(z).as_expr(),z)
    check("N7 recovered Gram matrix is PSD, rank 3", zero(poly.as_expr()-z**4*(z-2)*(z-sp.Rational(5,2))**2)
          and G*sp.ones(7,1)==sp.zeros(7,1),
          "characteristic polynomial = z^4 (z-2) (z-5/2)^2")
    energy = sum(1/sp.sqrt(2-2*G[i,j]) for i,j in it.combinations(range(7),2))
    expected = sp.Rational(1,2)+5*sp.sqrt(2)+sp.sqrt(50-10*s)/2+sp.sqrt(50+10*s)/2
    K = sp.QQ.algebraic_field(sp.sqrt(2),sp.sqrt(5),sp.sqrt(10-2*s),sp.sqrt(10+2*s))
    check("N7 candidate energy equals paper's radical expression", K.from_sympy(energy-expected)==K.zero)
    check("N7 pair contact multiplicities", [sum(zero(G[i,j]-x) for i,j in it.combinations(range(7),2)) for x in S]==[1,10,5,5])
    f = (2-2*u)**(-sp.Rational(1,2))
    for order in range(1,10):
        check(f"Coulomb derivative formula order {order}",
              zero(sp.diff(f,u,order)-sp.factorial2(2*order-1)*(2-2*u)**(-sp.Rational(2*order+1,2))))
    # Solve the confluent Vandermonde system in an exact algebraic number field.
    A, rhs = [], []
    for x in S:
        for d in (0,1):
            A.append([sp.diff(u**j,u,d).subs(u,x) for j in range(8)])
            rhs.append(sp.diff(f,u,d).subs(u,x))
    Am = DomainMatrix.from_Matrix(sp.Matrix(A)).convert_to(K)
    bm = DomainMatrix.from_Matrix(sp.Matrix(rhs)).convert_to(K)
    coeff = Am.inv().matmul(bm)
    check("N7 exact Hermite interpolation: all eight value/derivative equations", Am.matmul(coeff)==bm)
    hcoeff = list(coeff.to_Matrix())
    hsecond = sum(j*(j-1)*hcoeff[j]*(-1)**(j-2) for j in range(2,8))
    gamma = K.from_sympy(8*(sp.Rational(3,32)-hsecond))
    omega = sp.Poly(sp.expand((u+1)**2*u**2*(u-a)**2*(u-b)**2),u)
    check("N7 omega has rational coefficients and endpoint second derivative 1/8",
          all(x.is_Rational for x in omega.all_coeffs()) and sp.diff(omega.as_expr(),u,2).subs(u,-1)==sp.Rational(1,8))
    Hcoeff = [K.from_sympy(hcoeff[j]) if j<8 else K.zero for j in range(9)]
    Hcoeff = [Hcoeff[j]+gamma*K.from_sympy(omega.nth(j)) for j in range(9)]
    Hend2 = sum((j*(j-1)*(-1)**(j-2)*Hcoeff[j] for j in range(2,9)), K.zero)
    check("N7 enhanced polynomial matches f''(-1) exactly", Hend2==K.from_sympy(sp.Rational(3,32)))
    gn = num(K.to_sympy(gamma))
    gamma_interval=rational_enclosure(K.to_sympy(gamma))
    check("N7 gamma has a strictly positive rational enclosure", gamma_interval[0]>0,
          dict(decimal=mp.nstr(gn,45),lower=str(gamma_interval[0]),upper=str(gamma_interval[1])),
          "rigorous rational interval arithmetic")
    hn = [num(x) for x in hcoeff]
    Hn = [num(K.to_sympy(x)) for x in Hcoeff]
    contacts = [num(x) for x in S]
    points = [mp.mpf(-1)+mp.mpf(2)*i/2400 for i in range(2400)]
    for x in contacts:
        points += [x+sign*mp.mpf(10)**(-j) for j in range(2,19,2) for sign in (-1,1)
                   if -1<x+sign*mp.mpf(10)**(-j)<1]
    evaluate = lambda coeffs,x: mp.polyval(coeffs[::-1],x)
    errors, enhancements = [], []
    for x in points:
        if any(abs(x-c)<mp.mpf('1e-70') for c in contacts):
            continue
        errors.append(1/mp.sqrt(2-2*x)-evaluate(Hn,x))
        enhancements.append(evaluate(Hn,x)-evaluate(hn,x))
    check("N7 90-digit minorant grid and near-contact checks", min(errors)>0 and min(enhancements)>0,
          dict(points=len(errors), min_f_minus_H=mp.nstr(min(errors),12),
               min_H_minus_h=mp.nstr(min(enhancements),12)), "numerical sampling, not interval proof")
    eps = mp.mpf('1e-8')
    x = -1+eps
    w = mp.fprod((x-c)**2 for c in contacts)
    excess = 1/mp.sqrt(2-2*x)-evaluate(hn,x)-(gn+mp.mpf('0.0001'))*w
    check("N7 increasing gamma violates the endpoint bound", excess<0, mp.nstr(excess,15), "90-digit numerical check")
    # Exact remainder products establish the sign for every point of [-1,1).
    check("N7 exact ninth-order Hermite remainder product",
          sp.expand((u+1)*omega.as_expr())==sp.expand((u+1)**3*u**2*(u-a)**2*(u-b)**2))
    DETAILS['n7'] = dict(energy=mp.nstr(num(expected),65), gamma=mp.nstr(gn,65),
                        h7_base_coefficients=[mp.nstr(x,65) for x in hn],
                        h7_enhanced_coefficients=[mp.nstr(x,65) for x in Hn],
                        graph_count=total, surviving_labelled_graphs=252)
    return float(num(expected))


def n8_checks():
    t = sp.Symbol('t')
    s=sp.sqrt(2)
    c=[t,2*t-1,-t+(1-t)/s,-t-(1-t)/s]
    w=[8,4,8,8]
    F=t*t+6*t+1
    H=t**3+3*t*t+15*t+29
    V=(9+4*s)*F**3-32*(1-t)**3*H
    Q=sp.expand(V**2-512*(1-t)**6*F**3)
    X=(1+2*s)**2/(1-t)**3
    Y=16*(2+s)**2/(2-s+(2+s)*t)**3
    Z=16*(2-s)**2/(2+s+(2-s)*t)**3
    check("N8 elimination rational-function identity", sp.cancel((X-Y-Z)**2-4*Y*Z-Q/((1-t)**6*F**6),extension=s)==0)
    check("N8 elimination degree and leading coefficient", sp.degree(Q,t)==12 and zero(sp.Poly(Q,t).LC()-1201-328*s))
    qmod=sp.Poly(sp.expand(Q).subs(s,18),t,modulus=23).mul_ground(pow(21,-1,23))
    reduced=[int(x)%23 for x in qmod.all_coeffs()]
    expected=[1,8,22,9,17,7,20,17,14,6,22,7,18]
    values=[int(qmod.eval(i))%23 for i in range(23)]
    check("N8 polynomial reduction modulo 23", reduced==expected and (18*18)%23==2 and (1201+328*18)%23==21, reduced)
    check("N8 all 23 residue values match the manuscript and are nonzero",
          values==[18,7,18,9,20,8,2,19,1,21,15,2,3,6,13,8,9,17,6,13,6,18,14],values)
    rows=[]
    for n in it.product(range(8),repeat=4):
        if sum(n)==7 and n[1]==1 and n[2]==n[3] and n[0]+2*n[1]-n[2]-n[3]==0:
            rows.append(n)
    check("N8 unique integer row counts",rows==[(2,1,2,2)],rows,"exact finite enumeration")
    q=sp.expand(sum(ww*cc**6 for ww,cc in zip(w,c))/28)
    q_expected=(231*t**6-630*t**5+735*t**4-460*t**3+165*t*t-30*t+3)/14
    check("N8 sixth moment polynomial",zero(q-q_expected))
    qL=sp.diff(q,t).subs(t,sp.Rational(39,125))
    check("N8 sixth moment derivative rational bound",qL==sp.Rational(33046048032,213623046875) and qL>sp.Rational(3,20),str(qL))
    for power in (1,3):
        check(f"N8 weighted contact moment {power}",zero(sum(ww*cc**power for ww,cc in zip(w,c))+4))
    a,b,cc,d=c
    check("N8 algebraic inequalities factorizations",all(zero(v) for v in
        [2*d*d-1-a-(3-2*s)*t*(t-1), 1+a+2*d-(s-1)*(t-1),
         (1+b)*(1+d)-(a+d)**2-(t-1)*(1-(5-2*s)*t)/2]))
    # Verify the manuscript's six rational endpoint comparisons, exactly.
    lo,hi=sp.Rational(7071,5000),sp.Rational(14143,10000)
    tt,ss=sp.symbols('tt ss')
    XX=(1+2*ss)**2/(1-tt)**3
    YY=16*(2+ss)**2/(2-ss+(2+ss)*tt)**3
    ZZ=16*(2-ss)**2/(2+ss+(2-ss)*tt)**3
    L,U=sp.Rational(39,125),sp.Rational(8,25)
    comparisons=[XX.subs({tt:L,ss:hi})<sp.Rational(671,100)**2,
                 YY.subs({tt:L,ss:lo})>sp.Rational(643,100)**2,
                 ZZ.subs({tt:L,ss:hi})>sp.Rational(17,50)**2,
                 XX.subs({tt:U,ss:lo})>sp.Rational(341,50)**2,
                 YY.subs({tt:U,ss:hi})<sp.Rational(629,100)**2,
                 ZZ.subs({tt:U,ss:lo})<sp.Rational(7,20)**2]
    check("N8 six rational endpoint comparisons",all(comparisons))
    sn=mp.sqrt(2)
    phi=lambda v:(2+4*sn)/mp.sqrt(1-v)+8/mp.sqrt(2-sn+(2+sn)*v)+8/mp.sqrt(2+sn+(2-sn)*v)
    deriv=lambda v:(2+4*sn)/(2*(1-v)**mp.mpf('1.5'))-4*(2+sn)/(2-sn+(2+sn)*v)**mp.mpf('1.5')-4*(2-sn)/(2+sn+(2-sn)*v)**mp.mpf('1.5')
    root=mp.findroot(deriv,(mp.mpf('.312'),mp.mpf('.32')))
    check("N8 stationary parameter interval",mp.mpf(39)/125<root<mp.mpf(8)/25 and abs(deriv(root))<mp.mpf('1e-75'),mp.nstr(root,65),"90-digit numerical check")
    check("N8 f^(8)(-1)/8! equals 6435/16777216",sp.factorial2(15)/sp.Integer(4)**sp.Rational(17,2)/sp.factorial(8)==sp.Rational(6435,16777216))
    cn=[num(x.subs(t,sp.Float(str(root),90))) for x in c]
    mat=mp.matrix([[x**j if k==0 else (j*x**(j-1) if j else 0) for j in range(8)] for x in cn for k in (0,1)])
    rhs=mp.matrix([1/mp.sqrt(2-2*x) if k==0 else (2-2*x)**mp.mpf('-1.5') for x in cn for k in (0,1)])
    coeff=mp.lu_solve(mat,rhs)
    gam=mp.mpf(6435)/16777216
    gaps=[]
    for i in range(2000):
        x=mp.mpf(-1)+mp.mpf(2)*i/2000
        gaps.append(1/mp.sqrt(2-2*x)-mp.polyval(list(coeff)[::-1],x)-gam*mp.fprod((x-y)**2 for y in cn))
    check("N8 enhanced Hermite minorant sampled on 2000 points",min(gaps)>0,mp.nstr(min(gaps),15),"90-digit numerical sampling, not interval proof")
    # The geometry is reconstructed directly from the circulant first row.
    row=sp.Matrix([1,d,a,cc,b,cc,a,d])
    G=sp.Matrix(8,8,lambda i,j:row[(j-i)%8])
    check("N8 circulant Gram eigenvalue identity",all(zero(x) for x in G*G-4*(1-t)*G-4*t*sp.Matrix(8,8,lambda i,j:(-1)**(i-j)*(3*t-1))))
    yn=np.array([[math.sqrt(1-float(root))*math.cos(3*math.pi*k/4),math.sqrt(1-float(root))*math.sin(3*math.pi*k/4),(-1)**k*math.sqrt(float(root))] for k in range(8)])
    Gn=np.array(G.subs(t,sp.Float(str(root),40))).astype(float)
    eig=np.linalg.eigvalsh(Gn)
    check("N8 explicit coordinates, circulant Gram and rank",np.max(abs(yn@yn.T-Gn))<2e-14 and np.sum(eig>1e-9)==3 and min(eig)>-1e-12,
          eig.tolist(),"floating-point numerical check")
    DETAILS['n8']=dict(t_star=mp.nstr(root,65),energy=mp.nstr(phi(root),65),
                       contacts=[mp.nstr(x,45) for x in cn],qprime=mp.nstr(num(sp.diff(q,t).subs(t,sp.Float(str(root),90))),45))
    return yn,float(phi(root))


def moment_checks(points,label):
    N=len(points)
    G=points@points.T
    injections=np.array(list(it.permutations(range(N),5)))
    samples=G[injections[:,:,None],injections[:,None,:]]
    minors=np.stack([np.linalg.det(samples[:,J][:,:,J]) for J in it.combinations(range(5),4)])
    check(f"{label} all sampled 4x4 Gram determinants vanish",np.max(abs(minors))<2e-12,
          dict(injections=len(injections),max_abs=float(np.max(abs(minors)))),"floating-point numerical check")
    for r in range(9):
        direct=np.mean(G[np.triu_indices(N,1)]**r)
        check(f"{label} edge moment normalization r={r}",abs(np.mean(samples[:,0,1]**r)-direct)<1e-13,kind="floating-point numerical check")
    def multi3(degree):
        return [x for x in it.product(range(degree+1),repeat=3) if sum(x)<=degree]
    for ell in range(9):
        alpha=np.array(multi3((8-ell)//2))
        P3=np.prod(samples[:,:3,3,None]**alpha.T[None,:,:],axis=1)
        P4=np.prod(samples[:,:3,4,None]**alpha.T[None,:,:],axis=1)
        expansion=(N-3)*(P3.T@P3)/len(samples)+(N-3)*(N-4)*(P3.T@(eval_legendre(ell,samples[:,3,4])[:,None]*P4))/len(samples)
        direct=np.zeros_like(expansion)
        for triple in it.permutations(range(N),3):
            rem=[j for j in range(N) if j not in triple]
            features=np.prod(G[np.ix_(triple,rem)][:,:,None]**alpha.T[:,None,:],axis=0)
            direct+=features.T@eval_legendre(ell,G[np.ix_(rem,rem)])@features
        direct/=N*(N-1)*(N-2)
        discrepancy=np.max(abs(expansion-direct))
        mineig=np.linalg.eigvalsh((direct+direct.T)/2).min()
        check(f"{label} F_{ell}: full matrix expansion and PSD",discrepancy<2e-11 and mineig>-2e-11,
              dict(dimension=len(alpha),max_abs_difference=float(discrepancy),min_eigenvalue=float(mineig)),"floating-point numerical check")
    # Degree-4, degree-3 and degree-2 compressed monomial bases for M,D,T.
    edges=list(it.combinations(range(5),2))
    values=np.array([samples[:,i,j] for i,j in edges]).T
    def features(degree):
        out=[np.ones(len(samples))]
        for j in range(10):
            out.append(values[:,j])
            out.append(values[:,j]**degree)
        return np.array(out).T
    v4,v3,v2=(features(k) for k in (4,3,2))
    M=v4.T@v4/len(samples)
    D=v3.T@((1-samples[:,0,1]**2)[:,None]*v3)/len(samples)
    T=np.einsum('nab,ni,nj->aibj',samples,v2,v2,optimize=True).reshape(5*v2.shape[1],5*v2.shape[1])/len(samples)
    eigs=[float(np.linalg.eigvalsh((B+B.T)/2).min()) for B in (M,D,T)]
    check(f"{label} selected monomial compressions of M,D,T are PSD",min(eigs)>-2e-11,eigs,
          "floating-point numerical check; compressed bases only")


def energy_grad(raw):
    raw=raw.reshape(-1,3)
    norms=np.linalg.norm(raw,axis=1)
    if np.min(norms)<1e-14:
        return 1e50,np.zeros(raw.size)
    X=raw/norms[:,None]
    delta=X[:,None,:]-X[None,:,:]
    d2=np.sum(delta*delta,axis=2)
    np.fill_diagonal(d2,np.inf)
    if d2.min()<1e-24:
        return 1e50,np.zeros(raw.size)
    energy=np.sum(d2**-.5)/2
    gradient=-np.sum(delta*d2[:,:,None]**-1.5,axis=1)
    gradient=(gradient-X*np.sum(gradient*X,axis=1)[:,None])/norms[:,None]
    return float(energy),gradient.ravel()


def numerical_search(n,reference,starts,seed):
    rng=np.random.default_rng(seed)
    raw=rng.normal(size=(n,3))
    _,grad=energy_grad(raw.ravel())
    direction=rng.normal(size=raw.size); direction/=np.linalg.norm(direction)
    step=1e-6
    fd=(energy_grad(raw.ravel()+step*direction)[0]-energy_grad(raw.ravel()-step*direction)[0])/(2*step)
    check(f"N{n} optimization gradient independent finite difference",abs(fd-grad@direction)<2e-7,
          dict(analytic=float(grad@direction),finite_difference=fd),"floating-point numerical check")
    records=[]
    for trial in range(starts):
        initial=rng.normal(size=(n,3))
        result=minimize(energy_grad,initial.ravel(),jac=True,method='BFGS',
                        options=dict(gtol=2e-9,maxiter=1500))
        X=result.x.reshape(n,3); X/=np.linalg.norm(X,axis=1)[:,None]
        energy,_=energy_grad(X.ravel())
        records.append(dict(trial=trial,energy=energy,gap=energy-reference,
                            gradient_norm=float(np.linalg.norm(energy_grad(X.ravel())[1])),
                            optimizer_success=bool(result.success),message=str(result.message),
                            coordinates=X.tolist()))
    best=min(records,key=lambda r:r['energy'])
    check(f"N{n} seeded multistart found no energy below candidate",best['energy']>=reference-2e-9,
          dict(starts=starts,seed=seed,best_energy=best['energy'],candidate=reference,
               worst_energy=max(r['energy'] for r in records),
               near_candidate=sum(abs(r['gap'])<1e-8 for r in records)),
          "heuristic numerical experiment, not a proof of global optimality")
    DETAILS[f'n{n}_multistart']=records


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=Path(__file__).resolve().parent/'Thomson-26-10.30.tex')
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'verification-results.json')
    parser.add_argument('--starts',type=int,default=64)
    parser.add_argument('--seed',type=int,default=20261102)
    args=parser.parse_args()
    if args.starts<1:
        parser.error('--starts must be at least 1')
    mp.mp.dps=90
    started=time.time()
    source=args.source.read_bytes() if args.source.exists() else None
    status='failed'
    try:
        if source is not None:
            check('Source manuscript matches the audited version',
                  hashlib.sha256(source).hexdigest()==EXPECTED_SOURCE_SHA256,
                  EXPECTED_SOURCE_SHA256,'SHA-256 version guard')
        else:
            print('[INFO] No manuscript source supplied: checking transcribed formulas; source-version guard skipped.',flush=True)
        e7=n7_exact()
        y8,e8=n8_checks()
        y7=np.array([[math.cos(2*math.pi*k/5),math.sin(2*math.pi*k/5),0] for k in range(5)]+[[0,0,1],[0,0,-1]],float)
        moment_checks(y7,'N7 candidate')
        rng=np.random.default_rng(args.seed)
        random=rng.normal(size=(7,3)); random/=np.linalg.norm(random,axis=1)[:,None]
        moment_checks(random,'N7 random configuration')
        moment_checks(y8,'N8 candidate')
        numerical_search(7,e7,args.starts,args.seed)
        numerical_search(8,e8,max(16,args.starts//2),args.seed+1)
        status='passed'
    finally:
        data=dict(status=status,source=args.source.name if source else None,
                  source_sha256=hashlib.sha256(source).hexdigest() if source else None,
                  verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__,mpmath=mp.__version__),
                  precision_digits=mp.mp.dps,elapsed_seconds=time.time()-started,
                  checks=CHECKS,details=DETAILS,
                  limitations=['No explicit SDP dual certificate (lambda, mu_1, mu_3, Z_j) is present in the supplied manuscript or this program.',
                               'The conditional global bound in the manuscript is not independently certified here.',
                               'Floating-point PSD tests, high-precision grids and random optimization are numerical evidence, not exact or interval certificates.',
                               'M,D,T are checked on selected compressed monomial bases; F_ell uses all stated basis monomials.',
                               'External references, Lean proofs and bibliography claims are not audited.'])
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
        print(f"Saved {args.output.name}; status={status}; checks={len(CHECKS)}",flush=True)
    return 0


if __name__=='__main__':
    main()
