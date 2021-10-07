import sympy as sp
import pystencils as ps

from lbmpy import LBMConfig, LBMOptimisation, LBStencil, Method, Stencil

from lbmpy.creationfunctions import create_lb_update_rule
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
from lbmpy.boundaries import NoSlip, SimpleExtrapolationOutflow, UBB

from pystencils_walberla import CodeGeneration, generate_sweep, generate_pack_info_from_kernel
from lbmpy_walberla import generate_boundary


#   ========================
#      General Parameters
#   ========================

stencil = 'D3Q27'
omega = sp.Symbol('omega')
ubb = sp.symbols('ubbX, ubbY, ubbZ')
layout = 'fzyx'

#   PDF Fields
pdfs, pdfs_tmp = ps.fields('pdfs(27), pdfs_tmp(27): [3D]', layout=layout)

#   Velocity Output Field
velocity = ps.fields("velocity(3): [3D]", layout=layout)
output = {'velocity': velocity}

#   Optimization
optimization = {'cse_global': True,
                'symbolic_field': pdfs,
                'symbolic_temporary_field': pdfs_tmp,
                'field_layout': layout}


#   ==================
#      Method Setup
#   ==================

lbm_params = {'stencil': stencil,
              'method': 'cumulant',
              'relaxation_rate': omega,
              'compressible': False}

lbm_update_rule = create_lb_update_rule(optimization=optimization,
                                        output=output,
                                        **lbm_params)

lbm_method = lbm_update_rule.method

#   ========================
#      PDF Initialization
#   ========================

initial_rho = sp.Symbol('rho_0')

pdfs_setter = macroscopic_values_setter(lbm_method,
                                        initial_rho,
                                        velocity.center_vector,
                                        pdfs.center_vector)

#   =====================
#      Code Generation
#   =====================

with CodeGeneration() as ctx:
    target = ps.Target.CPU

    #   LBM Sweep
    generate_sweep(ctx, "CumulantMRTSweep", lbm_update_rule, field_swaps=[(pdfs, pdfs_tmp)], target=target)

    #   Pack Info
    generate_pack_info_from_kernel(ctx, "CumulantMRTPackInfo", lbm_update_rule, target=target)

    #   Macroscopic Values Setter
    generate_sweep(ctx, "InitialPDFsSetter", pdfs_setter, target=target)

    #   NoSlip Boundary
    generate_boundary(ctx, "CumulantMRTNoSlip", NoSlip(), lbm_method, target=target)

    # UBB with constant velocity
    generate_boundary(ctx, "CumulantMRTSimpleUBB", UBB(ubb), lbm_method, target=target)

    stencil2 = LBStencil(Stencil.D3Q27)
    outflow = SimpleExtrapolationOutflow(stencil2[4], stencil)

    # Outflow
    generate_boundary(ctx, "CumulantMRTOutflow", outflow, lbm_method, target=target)