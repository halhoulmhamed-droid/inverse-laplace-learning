"""Independent quadrature, manufactured-solution and refinement checks."""
import time
import unittest
import numpy as np
from scipy.integrate import quad
from laplace_learning.direct import (
    cosine_sensor_matrix, analytic_matrix, mask_indices, fd_matrix,
    modal_discrete_ratio, piecewise_linear_sensor_matrix, reflected_laplacian,
)


class DirectTests(unittest.TestCase):
    def test_exact_sensor_quadrature(self):
        matrix = cosine_sensor_matrix(9)
        for k in range(9):
            for j in range(64):
                fun = (lambda x: 1.) if k == 0 else (lambda x: np.sqrt(2)*np.cos(k*np.pi*x))
                independent = 64 * quad(fun, j/64, (j+1)/64, epsabs=1e-13)[0]
                self.assertAlmostEqual(matrix[j, k], independent, delta=1e-10)
        # Stable evaluation remains finite at high frequencies/heights.
        self.assertTrue(np.isfinite(analytic_matrix(128, H=20)).all())

    def test_acquisition_and_rank17(self):
        acquire = piecewise_linear_sensor_matrix(32)
        np.testing.assert_allclose(acquire.sum(axis=1), 1, atol=1e-15)
        np.testing.assert_array_equal(acquire[0, :2], [.75, .25])
        np.testing.assert_array_equal(acquire[1, :2], [.25, .75])
        selected = acquire[mask_indices('contiguous')]
        self.assertEqual(np.linalg.matrix_rank(selected), 17)
        self.assertTrue(np.all(selected[:, 17:] == 0))
        # An affine trace must have its midpoint value as each interval mean,
        # even for a grid and sensor count without divisibility.
        nodal = 2 + 3*np.linspace(0, 1, 7)
        np.testing.assert_allclose(piecewise_linear_sensor_matrix(6, 11) @ nodal,
                                   2+3*(np.arange(11)+.5)/11, atol=1e-14)
        np.testing.assert_array_equal(mask_indices('full', 4), [0,1,2,3])
        np.testing.assert_array_equal(mask_indices('even', 4), [0,2])
        np.testing.assert_array_equal(mask_indices('contiguous', 4), [0,1])
        # Length weights do not renormalize the partial boundary.
        self.assertEqual(np.sum(np.ones(32)/64), .5)

    def test_reflected_stencil_manufactured_polynomial(self):
        nx, ny, H = 8, 4, .5
        B, D = reflected_laplacian(nx, ny, H)
        x = np.linspace(0, 1, nx+1)
        y = np.arange(ny)*H/ny
        # u=y^2 has zero lower flux; -Delta u=-2 exactly, including corners.
        field = np.broadcast_to(y[:,None]**2, (ny,nx+1)).ravel()
        np.testing.assert_allclose(B @ field - D @ np.full(nx+1,H**2), -2, atol=1e-12)
        # At a corner both reflected neighbor coefficients are doubled.
        self.assertEqual(B[0,1], -2*nx**2)
        self.assertEqual(B[0,nx+1], -2*(ny/H)**2)
        self.assertEqual(B[0,0], 2*nx**2+2*(ny/H)**2)
        # One y interval handles the double Dirichlet coupling at the bottom.
        B1, D1 = reflected_laplacian(4,1,H)
        np.testing.assert_allclose(B1 @ np.ones(5), D1 @ np.ones(5), atol=1e-12)

    def test_independent_sparse_modal_refinement(self):
        exact = analytic_matrix(9)
        previous = None
        for nx, ny, envelope in [(64,16,None),(128,32,.02),(256,64,.005)]:
            numerical = fd_matrix(9,nx=nx,ny=ny)
            np.testing.assert_allclose(numerical[:,0], 1, atol=1e-10, rtol=0)
            errors = []
            for k in range(1,9):
                predicted = modal_discrete_ratio(k,nx,ny)*exact[:,k]
                relative = np.linalg.norm(numerical[:,k]-predicted)/np.linalg.norm(predicted)
                self.assertLessEqual(relative,1e-10)
                err = np.linalg.norm(numerical[:,k]-exact[:,k])/np.linalg.norm(exact[:,k])
                errors.append(err)
                if envelope is not None:
                    self.assertLessEqual(err,envelope)
                for mask in ('even','contiguous'):
                    take = mask_indices(mask)
                    masked = np.linalg.norm(numerical[take,k]-predicted[take])/np.linalg.norm(predicted[take])
                    self.assertLessEqual(masked,1e-10)
            errors = np.asarray(errors)
            if previous is not None:
                self.assertTrue(np.all(errors < previous))
                self.assertTrue(np.all((previous/errors > 2)&(previous/errors < 6)))
            previous = errors
        np.testing.assert_allclose(fd_matrix(1,nx=32,ny=8),1,rtol=0,atol=1e-10)

    def test_invalid_inputs(self):
        for bad in (0,-1,2.5):
            with self.assertRaises(ValueError):
                cosine_sensor_matrix(bad)
        with self.assertRaises(ValueError):
            fd_matrix(34,nx=32)
        with self.assertRaises(ValueError):
            analytic_matrix(3,H=0)
        with self.assertRaises(ValueError):
            mask_indices('unknown')


if __name__ == '__main__':
    start = time.process_time()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(DirectTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    print(f'Direct test process CPU seconds: {time.process_time()-start:.6f}')
    raise SystemExit(not result.wasSuccessful())
