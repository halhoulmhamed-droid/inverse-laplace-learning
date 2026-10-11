import unittest
import numpy as np
from studies.study3.src.data import profiles,noises,observations,mask,losses,ROOT
from studies.study3.src.analysis import bootstrap
from studies.study2.src.model import ROOT as ROOT2
from laplace_learning.data import ROOT_SEED as ROOT1

class DataAnalysis(unittest.TestCase):
 def test_namespace_disjoint_and_nested(self):
  self.assertEqual(len({ROOT,ROOT1,ROOT2}),3)
  x=profiles(10,0,6);np.testing.assert_array_equal(x[:4],profiles(10,0,4))
  for part in [20,30,31,110,120,130]:self.assertFalse(np.array_equal(x,profiles(part,0,6)))
  self.assertFalse(np.array_equal(x,profiles(10,1,6)))
 def test_noise_pairing_and_scale(self):
  c=profiles(30,0,3);e=noises(30,0,3,4)
  full=observations(c,e,mask('full'),.01);part=observations(c,e,mask('half'),.01)
  np.testing.assert_array_equal(full[:,:,:32],part)
  clean=observations(c,np.zeros_like(e),mask('full'),.01)
  np.testing.assert_allclose(full-clean,e*.01/8,atol=1e-16)
 def test_parseval_and_high_modes(self):
  random=np.random.default_rng(97);c=random.normal(size=(5,128));p=random.normal(size=(5,64));A=random.normal(size=(9,64));y=random.normal(size=(5,9))
  actual=losses(p,c,y,A);padded=np.pad(p,((0,0),(0,64)))
  np.testing.assert_allclose(actual[:,0],((padded-c)**2).sum(1))
  np.testing.assert_allclose(actual[:,1],((padded-c)[:,8:]**2).sum(1))
  np.testing.assert_allclose(actual[:,2],((p@A.T-y)**2).sum(1))
 def test_bootstrap_shared_objects_preserved(self):
  base=np.arange(1,25,dtype=float).reshape(6,4)
  a=np.stack([.9*base,base],axis=-1);b=np.stack([1.01*base,base],axis=-1)
  r,v=bootstrap(a,b,0,400);self.assertTrue(r['criteria_met'])
  np.testing.assert_allclose(v,np.tile([.9,1.01],(400,1)))
  r,v=bootstrap(np.stack([base,base],axis=-1),b,0,400);self.assertFalse(r['criteria_met'])
 def test_bootstrap_invalid_explicit(self):
  a=np.ones((5,4,2));a[0,0,0]=np.nan
  r,v=bootstrap(a,np.ones_like(a),0,400);self.assertFalse(r['valid']);self.assertTrue(np.isnan(v).all())
  with self.assertRaises(ValueError):bootstrap(np.ones((5,4,3)),np.ones((5,4,2)),0)
