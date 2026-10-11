import hashlib,json,unittest,time
from pathlib import Path
import numpy as np
from studies.study2.src.model import *
from studies.study2.src.analysis import paired_bootstrap

class Audit(unittest.TestCase):
    def setUp(self):
        z=np.random.default_rng(2196)
        self.A=z.normal(size=(5,7));self.c=z.normal(size=(80,7));self.y=self.c@self.A.T+z.normal(size=(80,5))*.2
        self.f=covariance_family(self.A,self.c,self.y,rank=3)
    def get(self,m,j=8): return next(a for a in self.f['candidates'] if a['method']==m and a['index']==j)
    def test_ridge_matches_augmented_least_squares(self):
        f=self.f;z=(self.y-f['my'])@f['u'];b=(self.c-f['mc'])@f['v'];l=LAMBDAS[8]
        X=np.vstack([z,np.sqrt((len(z)-1)*l*f['t'])*np.eye(3)])
        Y=np.vstack([b,np.zeros((3,3))])
        oracle=np.linalg.lstsq(X,Y,rcond=None)[0].T
        np.testing.assert_allclose(f['v'].T@self.get('F')['R']@f['u'],oracle,atol=1e-12)
    def test_diagonal_covariance_full_equals_diagonal(self):
        # Exact orthogonal columns with zero means, not a near-diagonal random sample.
        Q=np.linalg.qr(np.column_stack([np.ones(20),np.random.default_rng(3).normal(size=(20,4))]))[0][:,1:]
        c=Q*np.array([1,2,3,4]);A=np.diag([4.,3.,2.,1.]);y=c@A.T
        f=covariance_family(A,c,y,rank=4)
        for j in range(20):
            a={k:next(v for v in f['candidates'] if v['method']==k and v['index']==j)['R'] for k in ['D','F']}
            np.testing.assert_allclose(a['D'],a['F'],atol=1e-12)
    def test_projected_operator_identities(self):
        for j in [0,8,19]:
            P=self.f['P'];np.testing.assert_allclose(P@self.get('O',j)['R'],self.get('F',j)['R'],atol=1e-12)
            np.testing.assert_allclose(P@self.get('L',j)['R'],self.get('LP',j)['R'],atol=1e-12)
    def test_diagonal_is_diagonal_quadratic_optimum(self):
        f=self.f;g=np.diag(f['v'].T@self.get('D')['R']@f['u'])
        derivative=g*(np.diag(f['Czz'])+LAMBDAS[8]*f['t'])-np.diag(f['Cbz'])
        np.testing.assert_allclose(derivative,0,atol=1e-12)
    def test_sign_change_does_not_change_diagonal_operator(self):
        f=self.f;sign=np.array([-1,1,-1]);u=f['u']*sign;v=f['v']*sign
        bz=(sign[:,None]*f['Cbz'])*sign;zz=(sign[:,None]*f['Czz'])*sign
        R=(v*(np.diag(bz)/(np.diag(zz)+LAMBDAS[8]*f['t'])))@u.T
        np.testing.assert_allclose(R,self.get('D')['R'],atol=1e-12)
    def test_parseval_decomposition_and_shared_projection(self):
        z=np.random.default_rng(1);truth=z.normal(size=(5,128));pred=z.normal(size=(5,64));V=np.linalg.qr(z.normal(size=(64,12)))[0];P=V@V.T
        loss=losses(pred,truth,P)
        padded=np.pad(pred,((0,0),(0,64)))
        np.testing.assert_allclose(loss[:,0],np.sum((padded-truth)**2,axis=1),atol=1e-11)
        np.testing.assert_allclose(loss[:,0],loss[:,1]+loss[:,2])
    def test_partition_disjoint_and_nested(self):
        a=profiles(10,0,5);np.testing.assert_array_equal(a,profiles(10,0,6)[:5])
        for part in [20,30,31,110,120,130]: self.assertFalse(np.array_equal(a,profiles(part,0,5)))
        self.assertFalse(np.array_equal(a,profiles(10,1,5)))
        eps=noises(30,0,5,4);np.testing.assert_array_equal(eps,noises(30,0,5,4))
    def test_noise_pairs_and_forward_scaling(self):
        c=profiles(30,0,5);e=noises(30,0,5,4)
        full=observations(c,e,mask('full'),.01);half=observations(c,e,mask('half'),.01)
        np.testing.assert_array_equal(full[:,:,:32],half)
        d=observations(c,np.zeros_like(e),mask('full'),.01)
        np.testing.assert_allclose(full-d,.01*e/8,atol=1e-16)
    def test_bootstrap_identity_and_failure(self):
        a=np.arange(1,21,dtype=float).reshape(5,4)
        result,b=paired_bootstrap(a,a,0,400);np.testing.assert_allclose(b,1)
        result,b=paired_bootstrap(.9*a,a,0,400);self.assertTrue(result['criteria_met']);np.testing.assert_allclose(b,.9)
        a[0,0]=np.nan;result,_=paired_bootstrap(a,np.ones_like(a),0,400);self.assertFalse(result['valid'])
    def test_invalid_training_explicit(self):
        self.c[0,0]=np.nan
        with self.assertRaises(ValueError):covariance_family(self.A,self.c,self.y,rank=3)
    def test_selection_budget_and_exact_shared_indices(self):
        # All 8 families have exactly 20 candidates; oracle validation uses only these arguments.
        f=self.f
        for m in METHODS:self.assertEqual(sum(a['method']==m for a in f['candidates']),20)
        # Selection supports n=64 physical setting.
        c=profiles(10,0,64);A=fd_matrix(64)[mask('half')]/8;y=observations(c,noises(10,0,64,1),mask('half'),.01)[:,0]
        f=covariance_family(A,c[:,:64],y);v=profiles(20,0,7);yv=observations(v,noises(20,0,7,2),mask('half'),.01).reshape(-1,32)
        selected,fixed,scores=select(f,yv,np.repeat(v,2,axis=0))
        self.assertEqual(len(scores),160)
        self.assertEqual([a['index'] for a in fixed[:3]],[selected['F']['index']]*3)
        self.assertEqual([a['index'] for a in fixed[3:]],[selected['L']['index']]*2)
        for m in METHODS:self.assertLessEqual(selected[m]['score'],min(a['score'] for a in scores if a['method']==m)+1e-12)

if __name__=='__main__':
    start=time.process_time();suite=unittest.defaultTestLoader.loadTestsFromTestCase(Audit);result=unittest.TextTestRunner(verbosity=2).run(suite)
    base=Path(__file__).resolve().parents[1]
    report=dict(success=result.wasSuccessful(),tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),cpu_seconds=time.process_time()-start,
        source_hashes={str(p.relative_to(base)):hashlib.sha256(p.read_bytes()).hexdigest() for folder in ['src','tests'] for p in (base/folder).glob('*.py')})
    (base/'results/tests.json').write_text(json.dumps(report,indent=2)+'\n')
    raise SystemExit(not result.wasSuccessful())
