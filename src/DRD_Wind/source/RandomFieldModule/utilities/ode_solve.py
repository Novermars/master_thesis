import sys
import numpy as np
from scipy.linalg import solve_banded, solveh_banded, cholesky_banded
from scipy.sparse import spdiags
from collections.abc import Iterable, Callable

from math import *
import scipy
import numba
# import numba_scipy
from time import time

# import sys, petsc4py
# petsc4py.init(sys.argv)
# from petsc4py import PETSc


############################################################################

class Grid1D:

    def __init__(self, size, lower_bound=0, upper_bound=1, grid=None):
        if grid is not None:
            self.grid = np.array(grid)
            self.size = self.grid.size
            self.lower_bound = np.min(self.grid)
            self.upper_bound = np.max(self.grid)
            self.h = np.diff(self.grid)
        else:
            self.size = size
            self.lower_bound = lower_bound
            self.upper_bound = upper_bound
            self.h = (self.upper_bound-self.lower_bound) / (self.size-1)
            self.grid = None


    def __getitem__(self, index):
        if self.grid is not None:
            return self.grid[index]
        else:
            if isinstance(index, slice):
                start = 0 if index.start is None else index.start
                stop = self.size if index.stop is None else index.stop
                step = 1 if index.step is None else index.step
                return self.lower_bound + self.h*np.arange(start,stop,step)
            else: 
                return self.lower_bound + self.h*index

    def __len__(self):
        return self.size


############################################################################

class FEM_coefficient_matrix_generator:

    def quadrature_points(self):
        # return (self.grid[1:]-self.grid.h/2)
        q = np.zeros([2, len(self.grid)-1])
        q[0,:] = (self.grid[1:] - self.grid.h/2) - 1/sqrt(3) * self.grid.h/2
        q[1,:] = (self.grid[1:] - self.grid.h/2) + 1/sqrt(3) * self.grid.h/2
        return q

    def __init__(self, grid, coef):
        if isinstance(coef, Iterable):
            if len(coef)==1:
                self.coef = [coef[0]] * 3
            elif len(coef)==2:
                self.coef = [ coef[0], coef[0], coef[1] ]
            else:
                self.coef = coef
        else:
            self.coef = [coef] * 3
        self.grid=grid
        self.qp = self.quadrature_points()
        # self.t=t

        # create functions which return the coefficients in the operator
        # term_M1 = self._term_M1(self.coef)
        # term_M2 = self._term_M2(self.coef)
        # term_D1 = self._term_D1(self.coef)
        # term_D2 = self._term_D2(self.coef)
        # term_D0 = self._term_D0(self.coef)
        # create the matrices from the coefficient functions
        self.Mass  = self.mass_matrix(self.func_const1())
        self.Mass1 = self.mass_matrix(self.coef[0])
        self.Mass2 = self.mass_matrix(self.coef[1])
        self.Mass3 = self.mass_matrix(self.coef[2])
        self.Cross = self.crossterm_diffusion_matrix(self.coef[2])
        self.Diff  = self.diffusion_matrix(self.coef[2])
        det = ( self.coef[0](self.grid[:]) * self.coef[1](self.grid[:]) * self.coef[2](self.grid[:]) )**(1/3)
        self.L2_vec = det**(17/12)
        self.sqrtMass = cholesky_banded(self.Mass[1:,:],lower=True)
        # self.sqrtM = self.mass_sqrt_matrix(term_M1)

    def assemble(self, d, k1, k2, t_func):
        if isinstance(t_func, Callable):
            t = t_func(k1, k2, 0)
        else:
            t = t_func
        M = self.Mass
        M1 = self.Mass1
        M2 = self.Mass2
        M3 = self.Mass3
        C = self.Cross
        D = self.Diff
        return (d+1)*M + (k1**2)*M1 + (k2**2)*M2 + ((t*k1)**2)*M3 + (t*k1)*C + D

    def __call__(self, d, k1, k2, t=0.0):
        '''create FEM matrix in so-called "matrix diagonal ordered form" '''
        # create functions which return the coefficients in the operator
        reaction_function = self._reaction_function(self.coef, d, k1, k2, t)
        diffusion_function = self._diffusion_function(self.coef, t)
        crossterm_diffusion_function = self._crossterm_diffusion_function(self.coef, k1, t)
        # create the matrices from the coefficient functions
        M = self.mass_matrix(reaction_function)
        D1 = self.diffusion_matrix(diffusion_function)
        D2 = self.crossterm_diffusion_matrix(crossterm_diffusion_function)
        return M+D1+D2

    # @staticmethod
    # def _reaction_function(kappa, d, k1, k2, t):
    #     return (lambda z : d + 1.0 + (1.0 + t**2)*kappa(z)*k1**2 + kappa(z)*k2**2)

    # @staticmethod
    # def _diffusion_function(kappa, t):
    #     return (lambda z : kappa(z))

    # @staticmethod
    # def _crossterm_diffusion_function(kappa, k1, t):
    #     return (lambda z : t*kappa(z)*k1)

    ###
    @staticmethod
    def func_const1():
        return (lambda z : 1)

    # @staticmethod
    # def _term_M1(kappa):
    #     return (lambda z : 1)

    # @staticmethod
    # def _term_M2(kappa):
    #     return (lambda z : kappa(z))

    # @staticmethod
    # def _term_D1(kappa):
    #     return (lambda z : 1)

    # @staticmethod
    # def _term_D2(kappa):
    #     return (lambda z : kappa(z))

    # @staticmethod
    # def _term_D0(kappa):
    #     return (lambda z : kappa(z))

    def mass_matrix(self, reaction_function):
        '''creates the mass matrix in matrix diagonal ordered form'''
        M = np.zeros((3,len(self.grid)))
        fun_qp = reaction_function(self.qp)
        if np.isscalar(fun_qp): fun_qp *= np.ones_like(self.qp)
        w = 1/2
        x1 = (1-1/np.sqrt(3))/2
        x2 = (1+1/np.sqrt(3))/2
        coeff_element_0 = (fun_qp[0,:]*(x1**2) + fun_qp[1,:]*(x2**2)) * self.grid.h * w
        coeff_element_1 = (fun_qp[0,:]*(x1*(1-x1)) + fun_qp[1,:]*(x2*(1-x2))) * self.grid.h * w
        M[0,1:]  = coeff_element_1
        M[1,:-1] = coeff_element_0
        M[1,1:] += coeff_element_0
        M[2,:-1] = M[0,1:]
        return M
    
    def diffusion_matrix(self, diffusion_function):
        '''creates the diffusion matrix in matrix diagonal ordered form'''
        D = np.zeros((3,len(self.grid)))
        fun_qp = diffusion_function(self.qp)
        if np.isscalar(fun_qp): fun_qp *= np.ones_like(self.qp)
        w = 1/2
        coeff_element = fun_qp.sum(axis=0)*w/self.grid.h
        # coeff_element_centers = diffusion_function(self.quadrature_points)
        # coeff_element_centers /= self.grid.h
        D[0,1:]  = -coeff_element
        D[1,:-1] = coeff_element
        D[1,1:] += coeff_element
        D[2,:-1] = D[0,1:]
        return D

    def crossterm_diffusion_matrix(self, crossterm_diffusion_function):
        '''creates the cross-term contributions to the diffusion matrix in matrix diagonal ordered form
           NOTE that this contribution is ONLY active when t > 0 '''
        D = np.zeros((3,len(self.grid)))
        fun_qp = crossterm_diffusion_function(self.qp)
        if np.isscalar(fun_qp): fun_qp *= np.ones_like(self.qp)
        w = 1/2
        x1 = (1-1/np.sqrt(3))/2
        x2 = (1+1/np.sqrt(3))/2
        coeff_element = 2*(fun_qp[0,:]*x1 + fun_qp[1,:]*x2) * w
        # coeff_element_centers = crossterm_diffusion_function(self.quadrature_points)
        D[0,1:]  = -coeff_element
        D[2,:-1] = -D[0,1:]
        return 1j*D

    # def mass_sqrt_matrix(self, reaction_function):
    #     '''creates sqrt of the mass matrix in matrix diagonal ordered form'''
    #     sqrtM = np.zeros((3,len(self.grid)))
    #     coeff_element_centers = reaction_function(self.quadrature_points)
    #     coeff_element_centers *= self.grid.h
    #     coeff_element_centers = np.sqrt(coeff_element_centers)
    #     # sqrtM[0,1:]  = coeff_element_centers
    #     sqrtM[1,:-1] = coeff_element_centers
    #     # sqrtM[1,1:] += coeff_element_centers
    #     sqrtM[2,:-1] = coeff_element_centers
    #     return 1/2*sqrtM

############################################################################

class FEM_load_vector_generator:

    def quadrature_points(self, order):
        # return (self.grid[1:]-self.grid.h/2)
        x = quadrature_points_local(order)
        q = np.zeros([len(x), len(self.grid)-1])
        for i in range(len(x)):
            q[i,:] = (self.grid[1:] - self.grid.h) + x[i]*self.grid.h
        return q

    def __init__(self, grid):
        self.grid=grid

    def __call__(self, fun, order=5):
        if callable(fun):
            fun = np.vectorize(fun)
            fun_qp = fun(self.quadrature_points(order))
        else:
            print('ERROR: incompatible vector size.\n')
            print('Expected vector of length ' + str(len(self.grid)) + '\t received vector of length ' + str(len(f)))
            sys.exit()
        x = quadrature_points_local(order)
        w = quadrature_weights_local(order)
        coeff_element1 = np.zeros(len(self.grid)-1, dtype=np.complex)
        coeff_element2 = np.zeros(len(self.grid)-1, dtype=np.complex)
        for i in range(order): 
            coeff_element1 += fun_qp[i,:] * (1-x[i]) * w[i] * self.grid.h
            coeff_element2 += fun_qp[i,:] * x[i] * w[i] * self.grid.h
        rhs = np.append(coeff_element1,[0])
        rhs[1:] += coeff_element2
        return rhs

    def Reynolds_stress_load(self, k1, k2, **kwargs):
        stress_comp = kwargs.get('stress_comp')
        psi_comp = kwargs.get('component')
        z = kwargs.get('z')
        self.delta_load = self._delta_load(z)
        self.diff_delta_load = self._diff_delta_load(z)
        if stress_comp == 0:
            if psi_comp == 0:
                return None
            elif psi_comp == 1:
                return -self.diff_delta_load
            elif psi_comp == 2:
                return 1j*k2*self.delta_load
            else:
                raise Exception("incompatible psi_comp: {}".format(psi_comp))
        elif stress_comp == 1:
            if psi_comp == 0:
                return self.diff_delta_load
            elif psi_comp == 1:
                return None
            elif psi_comp == 2:
                return -1j*k1*self.delta_load
            else:
                raise Exception("incompatible psi_comp: {}".format(psi_comp))
        elif stress_comp == 2:
            if psi_comp == 0:
                return -1j*k2*self.delta_load
            elif psi_comp == 1:
                return 1j*k1*self.delta_load
            elif psi_comp == 2:
                return None
            else:
                raise Exception("incompatible psi_comp: {}".format(psi_comp))
        else:
            raise Exception("incompatible stress_comp: {}".format(stress_comp))

    def _delta_load(self, z):
        # find index of last node before z
        ind = np.sum(1*(self.grid[:]<=z))-1

        # check if z is a node
        if not isclose(self.grid[ind], z):
            raise Exception("z is not a grid point. z = {} ".format(z))

        h = self.grid[ind+1]-self.grid[ind]
        w = 1 - (z-self.grid[ind])/h
        # construct load vector
        rhs = np.zeros(len(self.grid), dtype=np.complex)
        # rhs[ind] = w
        # rhs[ind+1] = (1-w)
        if ind == 0:
            rhs[ind] = 2/3
            rhs[ind+1] = 1/3
        else:
            h0 = self.grid[ind]-self.grid[ind-1]
            h1 = self.grid[ind+1]-self.grid[ind]
            rhs[ind-1] = 1/3*h0/(h0+h1)
            rhs[ind] = 2/3
            rhs[ind+1] = 1/3*h1/(h0+h1)
        return rhs
    
    def _diff_delta_load(self, z):
        # find index of last node before z
        ind = np.sum(1*(self.grid[:]<=z))-1

        # check if z is a node
        if not isclose(self.grid[ind], z):
            raise Exception("z is not a grid point. z = {} ".format(z))

        # h = np.mean(np.diff(self.grid[:]))
        # construct load vector
        rhs = np.zeros(len(self.grid), dtype=np.complex)
        # if ind is not 0:
        #     h = 2*(self.grid[ind]-self.grid[ind-1])
        #     rhs[ind-1] = -1/h
        #     rhs[ind+1] = 1/h
        # else:
        # rhs[ind] = -1/h
        # rhs[ind+1] = 1/h
        if ind == 0:
            h = self.grid[ind+1]-self.grid[ind]
            rhs[ind] = -1/h
            rhs[ind+1] = 1/h
        else:
            h0 = self.grid[ind]-self.grid[ind-1]
            h1 = self.grid[ind+1]-self.grid[ind]
            rhs[ind-1] = -1/(h0+h1)
            rhs[ind+1] = 1/(h0+h1)

        return rhs
    
    def _delta_load_smooth(self, z, epsilon=1e-2, order=5):

        # find index of last node before z
        ind = np.sum(1*(self.grid[:]<=z))-1

        # check if z is a node
        if not isclose(self.grid[ind], z):
            raise Exception("z is not a grid point. z = {} ".format(z))
        
        # two edge cases that have to be delt with
        if epsilon >= h[ind]:
            fun = np.vectorize(self.mollifier(epsilon, z0=z))
            return self.__call__(fun,order=order)
        elif ind > 1:
            if epsilon >= h[ind-1]:
                fun = np.vectorize(self.mollifier(epsilon, z0=z))
                return self.__call__(fun,order=order)

        # the normal situation
        fun = np.vectorize(self.mollifier(epsilon))
        x = quadrature_points_local(order)
        w = quadrature_weights_local(order)
        if callable(fun):
            qp = np.append(x*epsilon)
            fun_qp = fun(qp)
            x_qp_right = x*epsilon/h[ind]
            if ind > 1:
                x_qp_left = x*epsilon/h[ind-1]
        else:
            print('ERROR: incompatible vector size.\n')
            print('Expected vector of length ' + str(len(self.grid)) + '\t received vector of length ' + str(len(f)))
            sys.exit()
        const = np.dot(fun_qp,w)*epsilon
        coeff_shape_function_left = np.zeros(len(self.grid)-1, dtype=np.complex)
        coeff_shape_function_right = np.zeros(len(self.grid)-1, dtype=np.complex)
        for i in range(order): 
            coeff_shape_function_left[ind] += fun_qp[i] * (1-x_qp_right[i]) * w[i] * epsilon
            coeff_shape_function_right[ind] += fun_qp[i] * x_qp_right[i] * w[i] * epsilon
            if ind > 1:
                coeff_shape_function_left[ind-1] += fun_qp[i] * x_qp_left[i] * w[i] * epsilon
                coeff_shape_function_right[ind-1] += fun_qp[i] * (1-x_qp_left[i]) * w[i] * epsilon
                
        rhs = np.append(coeff_shape_function_left,[0])
        rhs[1:] += coeff_shape_function_right

        if ind > 1:
            rhs /= 2*C
        else:
            rhs /= C
        return rhs

    @staticmethod
    def mollifier(epsilon, z0=0.0):
        return lambda z: 0.0 if (z-z0)**2 >= epsilon**2 else np.exp(-epsilon**2/(epsilon**2-(z-z0)**2))

############################################################################

class ode_solve:

    def __init__(self, dof, coef, domain_height=1, grid=None):
        self.grid = Grid1D(dof, upper_bound=domain_height, grid=grid)
        self.FEM_coefficient_matrix_generator = FEM_coefficient_matrix_generator(self.grid, coef)
        self.FEM_load_vector_generator = FEM_load_vector_generator(self.grid)
    
    def __call__(self, d, f, k1, k2, t=0.0, Robin_const=None, adjoint=False, **kwargs):
        test = kwargs.get('test', False)

        # M1 = self.FEM_coefficient_matrix_generator.M1
        # M2 = self.FEM_coefficient_matrix_generator.M2
        # D1 = self.FEM_coefficient_matrix_generator.D1
        # D2 = self.FEM_coefficient_matrix_generator.D2
        # D0 = self.FEM_coefficient_matrix_generator.D0
        # L2 = self.FEM_coefficient_matrix_generator.L2_vec
        sqrtM = self.FEM_coefficient_matrix_generator.sqrtMass
        # A = assemble(d, (1+t**2)*k1**2 + k2**2, 0, t*k1, M1, M2, D1, D2, D0)
        A = self.FEM_coefficient_matrix_generator.assemble(d, k1, k2, t)


        if test:
            rhs = self.FEM_load_vector_generator(f)
        else:
            if callable(f):
                rhs = f(self.grid[:])
            else:
                rhs = 1*f  ### variant of copying vector (important fot not modifying the input)
            # if not adjoint:
            #     rhs = 1*f
                # rhs = mult(sqrtM, rhs) ### this is applied in Sampling_method now


        if Robin_const is not None:
            if Robin_const == 'infty' or np.isinf(Robin_const):
                # Dirichlet BCs at z=0
                A[0,1] = 0.0
                A[1,0] = 1.0
                rhs[0] = 0
            elif Robin_const < 0:
                
                B = np.array(spdiags(A, [1, 0, - 1], len(rhs), len(rhs)).toarray())
                B[0,0] += B[-1,-1]
                B[0,-2] = B[-1,-2]
                B[-2,0] = B[-2,-1]
                B = B[:-1,:-1]

                sol = np.linalg.solve(B, rhs[:-1])
                return np.append(sol,[sol[0]])

                # sol = np.linalg.solve(B, rhs[:-1])
                # sol = np.append(sol,[sol[0]])
                # sol -= np.flip(sol)
                # return 0.5*sol

            else:
                # Robin BCs at z=0
                A[1,0] += Robin_const

        # NOTE: ONLY FOR TESTING
        if test:
            ### Dirichlet BCs at z=H
            A[1,-1] = 1.0
            A[2,-2] = 0.0
            rhs[-1] = 1

        sol = solve_banded((1, 1), A, rhs)

        return sol



    # def apply_matvec(self, d, f, k1, k2, t=0.0, Robin_const=None):
        
    #     # t = self.t
    #     M1 = self.FEM_coefficient_matrix_generator.M1
    #     M2 = self.FEM_coefficient_matrix_generator.M2
    #     D1 = self.FEM_coefficient_matrix_generator.D1
    #     D2 = self.FEM_coefficient_matrix_generator.D2
    #     D0 = self.FEM_coefficient_matrix_generator.D0
    #     L2 = self.FEM_coefficient_matrix_generator.L2_vec
    #     A = assemble(d, (1+t**2)*k1**2 + k2**2, 0, t*k1, M1, M2, D1, D2, D0)

    #     # if callable(f):
    #     #     rhs = self.FEM_load_vector_generator(f)
    #     # else:
    #     #     rhs = f*self.grid.h
    #     rhs = f / self.grid.h

    #     return mult(A, rhs)


    def apply_Mass(self, x):    
        M = self.FEM_coefficient_matrix_generator.Mass
        return mult(M, x)

    def apply_sqrtMass(self, x):    
        sqrtM = self.FEM_coefficient_matrix_generator.sqrtMass
        return mult_transpose(sqrtM, x)


# @numba.jit(nogil=True, cache=True)
# def solve(A, rhs):
# #     return scipy.linalg.solve_banded((1, 1), A, rhs)

# @numba.jit(nopython=True, nogil=True, cache=True)
# def load(f):
#     coeff_element_centers = (f[1:]+f[:-1])/2
#     rhs = np.append(coeff_element_centers,[0])
#     rhs[1:] += coeff_element_centers
#     return rhs/2

### (Legacy)
# @numba.jit(nopython=True, nogil=True, cache=True) 
# def assemble(d, k, ikt1, ikt2, M1, M2, D1, D2, D0):
#     # return d*M1 + k*M2 + ikt1*D1 + ikt2*D2  + D0
#     return (d+1)*M1 + k*M2 + ikt1*D1 + ikt2*D2  + D0


@numba.jit(nogil=True, cache=True, forceobj=True)
def mult(A, x):
    y = np.zeros_like(x, dtype=np.complex)
    if np.shape(A)[0] == 3:
        y[:-1]+= A[0,1:]*x[1:]
        y[:]  += A[1,:]*x
        y[1:] += A[2,:-1]*x[:-1]
    else:
        y[:]  += A[0,:]*x
        y[1:] += A[1,:-1]*x[:-1]
    return y

@numba.jit(nogil=True, cache=True, forceobj=True)
def mult_transpose(A, x):
    y = np.zeros_like(x, dtype=np.complex)
    if np.shape(A)[0] == 3:
        y[:-1]+= A[2,:-1]*x[1:]
        y[:]  += A[1,:]*x
        y[1:] += A[0,1:]*x[:-1]
    else:
        y[:-1]+= A[1,:-1]*x[1:]
        y[:]  += A[0,:]*x
    return y



def quadrature_points_local(order):
    if order == 1:
        x = np.array([0])
    if order == 2:
        x = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
    if order == 3:
        x = np.array([-np.sqrt(3/5), 0, np.sqrt(3/5)])
    if order == 4:
        x = np.array([-np.sqrt(3/7+2/7*np.sqrt(6/5)), -np.sqrt(3/7-2/7*np.sqrt(6/5)), np.sqrt(3/7-2/7*np.sqrt(6/5)), np.sqrt(3/7+2/7*np.sqrt(6/5))])
    if order == 5:
        x = np.array([-np.sqrt(5+2*np.sqrt(10/7)), -np.sqrt(5-2*np.sqrt(10/7)), 0, np.sqrt(5-2*np.sqrt(10/7)), np.sqrt(5+2*np.sqrt(10/7))])/3
    return (x+1)/2


def quadrature_weights_local(order):
    if order == 1:
        w = np.array([2])
    if order == 2:
        w = np.array([1, 1])
    if order == 3:
        w = np.array([5, 8, 5])/9
    if order == 4:
        w = np.array([18-np.sqrt(30), 18+np.sqrt(30), 18+np.sqrt(30), 18-np.sqrt(30)])/36
    if order == 5:
        w = np.array([(322-13*np.sqrt(70))/900, (322+13*np.sqrt(70))/900, 128/225, (322+13*np.sqrt(70))/900, (322-13*np.sqrt(70))/900])
    return w/2









############################################################################
############################################################################

if __name__ == "__main__":

    import matplotlib.pyplot as plt

    '''TEST by method of manufactured solutions'''

    # parameters
    L1 = lambda z : z
    L2 = lambda z : z**2
    L3 = lambda z : z**3
    dL3dz = lambda z : 3*z**2
    t = 100
    dof = 2**10
    domain_height = 0.9
    grid = (np.linspace(0,domain_height,dof)/domain_height)**2 *domain_height

    d = 1
    k1 = 1
    k2 = 1

    ## manufactured solution and derivatives (u=0 at z=0  dudz=0 at z=domain_height)
    n = 2
    const = (n + 1/2)*np.pi/domain_height
    psi = lambda z : np.sin(const*z)
    dpsidz = lambda z : const*np.cos(const*z)
    d2psidz2 = lambda z : -const**2*np.sin(const*z)

    # const = 1/domain_height**n
    # psi = lambda z : const*z**n
    # dpsidz = lambda z : const*n*z**(n-1)
    # d2psidz2 = lambda z : const*n*(n-1)*z**(n-2)*0

    # f = lambda z : ((d + 1.0 + L(z)**2*(1.0+t**2)*k1**2 + L(z)**2*k2**2 - 2j*t*L(z)*dLdz(z)*k1)*psi(z)
    #                 - 2.0*(1j*t*L(z)**2*k1 + L(z)*dLdz(z))*dpsidz(z) - L(z)**2*d2psidz2(z))
    f = lambda z : ( (d + 1.0 + (L1(z)**2 + t**2*L3(z)**2)*k1**2 + L2(z)**2*k2**2 - (1j*t*k1)*(2*L3(z)*dL3dz(z)))*psi(z)
                    - 2.0*(1j*t*k1*L3(z)**2 + L3(z)*dL3dz(z))*dpsidz(z) - L3(z)**2*d2psidz2(z) )

    kappa = [ np.vectorize(lambda z: L1(z)**2), np.vectorize(lambda z: L2(z)**2), np.vectorize(lambda z: L3(z)**2) ]
    ode_solve_inst = ode_solve(dof, kappa, domain_height=domain_height, grid=grid)
    psi_approx = ode_solve_inst(d, f, k1, k2, t=t, test=True, Robin_const=np.infty)
    z = np.linspace(0,domain_height,len(psi_approx))
    z = grid

    plt.plot(z,psi_approx.real)
    # plt.plot(z,psi_approx.real, 'o-')
    plt.plot(z,psi(z).real)
    plt.show()

    print(np.sum(psi_approx.imag**2)/np.sum(psi_approx.real**2))
    plt.plot(z,psi_approx.imag)
    plt.show()

