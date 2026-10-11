"""Study 2 only: paired covariance restrictions; existing project code unchanged."""
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from laplace_learning.direct import analytic_matrix, fd_matrix

ROOT = 2026101102
METHODS = ['D', 'F', 'O', 'L', 'LP', 'H1', 'signed', 'bounded']
FIXED = ['D_at_F', 'F_at_F', 'O_at_F', 'L_at_L', 'LP_at_L']
LAMBDAS = np.logspace(-8, 1, 20)
SCENARIOS = [
    ('primary', 'half', .01, 2.5, .25),
    ('full', 'full', .01, 2.5, .25),
    ('even32', 'even32', .01, 2.5, .25),
    ('quarter', 'quarter', .01, 2.5, .25),
    ('even16', 'even16', .01, 2.5, .25),
    ('noise_low', 'half', .003, 2.5, .25),
    ('noise_high', 'half', .03, 2.5, .25),
    ('rough', 'half', .01, 1.5, .25),
    ('height5', 'half', .01, 2.5, .2625),
    ('height10', 'half', .01, 2.5, .275),
]


def rng(partition, game, obj, channel=0, repetition=0):
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(
        [ROOT, partition, game, obj, channel, repetition])))


def profiles(partition, game, count, p=2.5):
    weights = (1 + np.arange(128, dtype=float)**2)**(-p/2)
    weights /= np.linalg.norm(weights)
    return np.array([rng(partition, game, i).normal(size=128)*weights for i in range(count)])


def noises(partition, game, count, repetitions):
    return np.array([[rng(partition, game, i, 1, j).normal(size=64)
                      for j in range(repetitions)] for i in range(count)])


def mask(name):
    return {'full': np.arange(64), 'half': np.arange(32),
            'quarter': np.arange(16), 'even32': np.arange(0,64,2),
            'even16': np.arange(0,64,4)}[name]


def observations(c, eps, indices, sigma, height=.25):
    signal = c @ analytic_matrix(128, H=height).T
    return (signal[:, None, indices] + sigma*eps[:, :, indices])/8


def solve(M, B):
    return cho_solve(cho_factor(M, lower=True), B)


def covariance_family(A, c, y, rank=12):
    """No validation/test arguments. c is N x 64; y is N x m."""
    if c.shape[0] < 2 or c.shape[0] != y.shape[0] or c.shape[1] != A.shape[1] or y.shape[1] != A.shape[0]:
        raise ValueError('incompatible training shapes')
    if not all(np.isfinite(x).all() for x in (A,c,y)):
        raise ValueError('nonfinite training input')
    u,s,vh = np.linalg.svd(A,full_matrices=False)
    if rank > min(A.shape): raise ValueError('rank exceeds dimensions')
    u,v = u[:,:rank],vh[:rank].T
    mc,my = c.mean(0),y.mean(0)
    cc,yy = c-mc,y-my
    b,z = cc@v,yy@u
    den = len(c)-1
    Cbz,Czz,Ccz = b.T@z/den,z.T@z/den,cc.T@z/den
    Ccy,Cyy = cc.T@yy/den,yy.T@yy/den
    t,ty = np.trace(Czz)/rank,np.trace(Cyy)/len(A)
    if t<=0 or ty<=0: raise ValueError('zero covariance scale')
    P=v@v.T
    candidates=[]
    def add(method,index,param,make):
        try:
            R=make()
            if not np.isfinite(R).all(): raise FloatingPointError('nonfinite R')
            norm=float(np.linalg.norm(R,2))
            valid=np.isfinite(norm)
            failure=None if valid else 'nonfinite norm'
        except (np.linalg.LinAlgError,ValueError,FloatingPointError) as exc:
            R=np.full((A.shape[1],len(A)),np.nan);norm=float('inf');valid=False;failure=str(exc)
        candidates.append(dict(method=method,index=index,parameter=float(param),R=R,
            norm=norm,valid=bool(valid),failure=failure))
    for j,l in enumerate(LAMBDAS):
        # Both covariance blocks are restricted, with the SAME ridge scale.
        add('D',j,l,lambda: (v*(np.diag(Cbz)/(np.diag(Czz)+l*t)))@u.T)
        add('F',j,l,lambda: v@solve(Czz+l*t*np.eye(rank),Cbz.T).T@u.T)
        add('O',j,l,lambda: solve(Czz+l*t*np.eye(rank),Ccz.T).T@u.T)
        add('L',j,l,lambda: solve(Cyy+l*ty*np.eye(len(A)),Ccy.T).T)
        add('LP',j,l,lambda: P@solve(Cyy+l*ty*np.eye(len(A)),Ccy.T).T)
    penalty=1+(np.pi*np.arange(A.shape[1]))**2;penalty=penalty/penalty.max()
    for j,a in enumerate(np.logspace(-10,1,20)):
        add('H1',j,a,lambda: solve(A.T@A+a*np.diag(penalty),A.T))
    gains=np.divide(np.diag(Cbz),np.diag(Czz),out=np.zeros(rank),where=np.diag(Czz)>0)
    for j,cap in enumerate(np.logspace(0,4,20)):
        add('signed',j,cap,lambda: (v*np.clip(gains,-cap,cap))@u.T)
        add('bounded',j,cap,lambda: (v*np.clip(gains,0,np.minimum(cap,1/s[:rank])))@u.T)
    return dict(candidates=candidates,mc=mc,my=my,u=u,v=v,P=P,s=s,Cbz=Cbz,Czz=Czz,t=t)


def predict(family,candidate,y):
    return family['mc']+(y-family['my'])@candidate['R'].T


def losses(prediction,truth,P):
    """Full 128-mode L2 squared loss; last modes are zero in reconstruction."""
    e=prediction-truth[...,:64]
    visible=e@P
    projected=np.sum(visible**2,axis=-1)
    complementary=np.sum((e-visible)**2,axis=-1)+np.sum(truth[...,64:]**2,axis=-1)
    return np.stack([projected+complementary,projected,complementary],axis=-1)


def select(family,yval,cval):
    """Only explicitly supplied validation truths; no test truth is read here."""
    scores=[]
    for candidate in family['candidates']:
        score=float('inf')
        if candidate['valid']:
            p=predict(family,candidate,yval)
            score=float(np.mean(np.sum((p-cval[:,:64])**2,axis=1)+np.sum(cval[:,64:]**2,axis=1)))
        candidate=dict(candidate,score=score)
        if not np.isfinite(score): candidate.update(valid=False,failure=candidate['failure'] or 'nonfinite validation')
        scores.append(candidate)
    selected={}
    for method in METHODS:
        group=[a for a in scores if a['method']==method and a['valid']]
        if not group:
            selected[method]=dict(next(a for a in scores if a['method']==method),valid=False)
            continue
        minimum=min(a['score'] for a in group);tol=1e-12*max(1,abs(minimum))
        selected[method]=min((a for a in group if a['score']-minimum<=tol),key=lambda a:(a['norm'],a['index']))
    fixed=[]
    for method,parent in [('D','F'),('F','F'),('O','F'),('L','L'),('LP','L')]:
        fixed.append(next(a for a in scores if a['method']==method and a['index']==selected[parent]['index']))
    return selected,fixed,scores
