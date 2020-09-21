#!/usr/bin/env python3
#
# Authors: Daniel Richtmann 2020
#          Christoph Lehner 2020
#
# Desc.: Test multigrid for clover
#
import gpt as g
import numpy as np

# setup rng, mute
g.default.set_verbose("random", False)
rng = g.random("test_mg")

# just run with larger volume
L = [8, 8, 16, 16]

# setup gauge field
U = g.qcd.gauge.random(g.grid(L, g.single), rng)
g.message("Plaquette:", g.qcd.gauge.plaquette(U))

# quark
w = g.qcd.fermion.wilson_clover(
    U,
    {
        "kappa": 0.137,
        "csw_r": 0,
        "csw_t": 0,
        "xi_0": 1,
        "nu": 1,
        "isAnisotropic": False,
        "boundary_phases": [1.0, 1.0, 1.0, 1.0],
    },
)

# default grid
grid = U[0].grid

# create source
src = g.vspincolor(grid)
src[:] = g.vspincolor([[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]])

# abbreviations
i = g.algorithms.inverter
mg = i.multi_grid
p = g.qcd.fermion.preconditioner

# NOTE: mg params structure
# - list -> configure each level by itself explicitly (correct lengths asserted inside)
# - scalar value (= not a list) -> broadcast parameter to every level


# mg setups
# mg_setup_2lvl = mg.setup(
#     w,
#     {
#         "block": [[2, 2, 2, 2]],
#         "northo": 1,
#         "nbasis": 30,
#         "make_hermitian": False,
#         "save_links": True,
#         "vecstype": "null",
#         "preortho": False,
#         "postortho": False,
#         "solver": i.fgmres(
#             {"eps": 1e-3, "maxiter": 50, "restartlen": 25, "checkres": False}
#         ),
#         "distribution": rng.cnormal,
#     },
# )
mg_setup_3lvl = mg.setup(
    w,
    {
        "block": [[2, 2, 2, 2], [1, 2, 2, 2]],
        "northo": 1,
        "nbasis": 30,
        "make_hermitian": False,
        "save_links": True,
        "vecstype": "null",
        "preortho": False,
        "postortho": False,
        "solver": i.fgmres(
            {"eps": 1e-3, "maxiter": 50, "restartlen": 25, "checkres": False}
        ),
        "distribution": rng.cnormal,
    },
)
# mg_setup_4lvl = mg.setup(
#     w,
#     {
#         "block": [[2, 2, 2, 2], [1, 2, 1, 1], [1, 1, 2, 2]],
#         "northo": 1,
#         "nbasis": 30,
#         "make_hermitian": False,
#         "save_links": True,
#         "vecstype": "null",
#         "preortho": False,
#         "postortho": False,
#         "solver": i.fgmres(
#             {"eps": 1e-3, "maxiter": 50, "restartlen": 25, "checkres": False}
#         ),
#         "distribution": rng.cnormal,
#     },
# )

# g.message(f"mg_setup_2lvl = {mg_setup_2lvl.params}")
g.message(f"mg_setup_3lvl = {mg_setup_3lvl.params}")
# g.message(f"mg_setup_4lvl = {mg_setup_4lvl.params}")

# mg inner solvers
wrapper_solver = i.fgmres(
    {"eps": 1e-1, "maxiter": 10, "restartlen": 5, "checkres": False}
)
smooth_solver = i.fgmres(
    {"eps": 1e-14, "maxiter": 8, "restartlen": 4, "checkres": False}
)
coarsest_solver = i.fgmres(
    {"eps": 5e-2, "maxiter": 50, "restartlen": 25, "checkres": False}
)

# mg solver/preconditioner objects
# mg_2lvl_vcycle = mg.inverter(
#     mg_setup_2lvl,
#     {
#         "coarsest_solver": coarsest_solver,
#         "smooth_solver": smooth_solver,
#         "wrapper_solver": None,
#     },
# )
# mg_2lvl_kcycle = mg.inverter(
#     mg_setup_2lvl,
#     {
#         "coarsest_solver": coarsest_solver,
#         "smooth_solver": smooth_solver,
#         "wrapper_solver": wrapper_solver,
#     },
# )
# mg_3lvl_vcycle = mg.inverter(
#     mg_setup_3lvl,
#     {
#         "coarsest_solver": coarsest_solver,
#         "smooth_solver": smooth_solver,
#         "wrapper_solver": None,
#     },
# )
mg_3lvl_kcycle = mg.inverter(
    mg_setup_3lvl,
    {
        "coarsest_solver": coarsest_solver,
        "smooth_solver": smooth_solver,
        "wrapper_solver": wrapper_solver,
    },
)
# mg_4lvl_vcycle = mg.inverter(
#     mg_setup_4lvl,
#     {
#         "coarsest_solver": coarsest_solver,
#         "smooth_solver": smooth_solver,
#         "wrapper_solver": None,
#     },
# )
# mg_4lvl_kcycle = mg.inverter(
#     mg_setup_4lvl,
#     {
#         "coarsest_solver": coarsest_solver,
#         "smooth_solver": smooth_solver,
#         "wrapper_solver": wrapper_solver,
#     },
# )

# preconditioners
smoother_prec = mg_3lvl_kcycle.smooth_solver[0]

# outer solver
fgmres_params = {"eps": 1e-6, "maxiter": 1000, "restartlen": 20}

# preconditioned inversion (using only smoother, w/o coarse grid correction)
fgmres_outer = i.fgmres(fgmres_params, prec=smoother_prec)
sol_smooth = g.eval(fgmres_outer(w) * src)
eps2 = g.norm2(w * sol_smooth - src) / g.norm2(src)
niter_prec_smooth = len(fgmres_outer.history)
g.message("Test resid/iter smoother prec fgmres:", eps2, niter_prec_smooth)
assert eps2 < 1e-10

# preconditioned inversion (2lvl mg -- vcycle)
# fgmres_outer = i.fgmres(fgmres_params, prec=mg_2lvl_vcycle)
# sol_prec_2lvl_mg_vcycle = g.eval(fgmres_outer(w) * src)
# eps2 = g.norm2(w * sol_prec_2lvl_mg_vcycle - src) / g.norm2(src)
# niter_prec_2lvl_mg_vcycle = len(fgmres_outer.history)
# g.message(
#     "Test resid/iter 2lvl vcycle mg prec fgmres:", eps2, niter_prec_2lvl_mg_vcycle
# )
# assert eps2 < 1e-10
# assert niter_prec_2lvl_mg_vcycle < niter_prec_smooth

# # preconditioned inversion (2lvl mg -- kcycle)
# fgmres_outer = i.fgmres(fgmres_params, prec=mg_2lvl_kcycle)
# sol_prec_2lvl_mg_kcycle = g.eval(fgmres_outer(w) * src)
# eps2 = g.norm2(w * sol_prec_2lvl_mg_kcycle - src) / g.norm2(src)
# niter_prec_2lvl_mg_kcycle = len(fgmres_outer.history)
# g.message(
#     "Test resid/iter 2lvl kcycle mg prec fgmres:", eps2, niter_prec_2lvl_mg_kcycle
# )
# assert eps2 < 1e-10
# assert niter_prec_2lvl_mg_kcycle == niter_prec_2lvl_mg_vcycle  # equivalent for 2 lvls

# # preconditioned inversion (3lvl mg -- vcycle)
# fgmres_outer = i.fgmres(fgmres_params, prec=mg_3lvl_vcycle)
# sol_prec_3lvl_mg_vcycle = g.eval(fgmres_outer(w) * src)
# eps2 = g.norm2(w * sol_prec_3lvl_mg_vcycle - src) / g.norm2(src)
# niter_prec_3lvl_mg_vcycle = len(fgmres_outer.history)
# g.message(
#     "Test resid/iter 3lvl vcycle mg prec fgmres:", eps2, niter_prec_3lvl_mg_vcycle
# )
# assert eps2 < 1e-10
# assert niter_prec_3lvl_mg_vcycle < niter_prec_smooth

# preconditioned inversion (3lvl mg -- kcycle)
fgmres_outer = i.fgmres(fgmres_params, prec=mg_3lvl_kcycle)
sol_prec_3lvl_mg_kcycle = g.eval(fgmres_outer(w) * src)
eps2 = g.norm2(w * sol_prec_3lvl_mg_kcycle - src) / g.norm2(src)
niter_prec_3lvl_mg_kcycle = len(fgmres_outer.history)
g.message(
    "Test resid/iter 3lvl kcycle mg prec fgmres:", eps2, niter_prec_3lvl_mg_kcycle
)
assert eps2 < 1e-10
assert niter_prec_3lvl_mg_kcycle <= niter_prec_smooth

# # preconditioned inversion (4lvl mg -- vcycle)
# fgmres_outer = i.fgmres(fgmres_params, prec=mg_4lvl_vcycle)
# sol_prec_4lvl_mg_vcycle = g.eval(fgmres_outer(w) * src)
# eps2 = g.norm2(w * sol_prec_4lvl_mg_vcycle - src) / g.norm2(src)
# niter_prec_4lvl_mg_vcycle = len(fgmres_outer.history)
# g.message(
#     "Test resid/iter 4lvl vcycle mg prec fgmres:", eps2, niter_prec_4lvl_mg_vcycle
# )
# assert eps2 < 1e-10
# assert niter_prec_4lvl_mg_vcycle <= niter_prec_3lvl_mg_vcycle

# # preconditioned inversion (4lvl mg -- kcycle)
# fgmres_outer = i.fgmres(fgmres_params, prec=mg_4lvl_kcycle)
# sol_prec_4lvl_mg_kcycle = g.eval(fgmres_outer(w) * src)
# eps2 = g.norm2(w * sol_prec_4lvl_mg_kcycle - src) / g.norm2(src)
# niter_prec_4lvl_mg_kcycle = len(fgmres_outer.history)
# g.message(
#     "Test resid/iter 4lvl kcycle mg prec fgmres:", eps2, niter_prec_4lvl_mg_kcycle
# )
# assert eps2 < 1e-10
# assert niter_prec_4lvl_mg_kcycle <= niter_prec_4lvl_mg_vcycle

# print contributions to mg setup runtime
g.message("Contributions to time spent in MG setups")
for name, t in [
    #    ("2lvl", mg_setup_2lvl.t),
    ("3lvl", mg_setup_3lvl.t),
    #    ("4lvl", mg_setup_4lvl.t),
]:
    g.message(name + ":")
    for lvl in reversed(range(len(t))):
        g.message(t[lvl])

# print contributions to mg solve runtime
g.message("Contributions to time spent in MG preconditioners")
for name, t in [
    #    ("2lvl_vcycle", mg_2lvl_vcycle.t),
    #    ("2lvl_kcycle", mg_2lvl_kcycle.t),
    #    ("3lvl_vcycle", mg_3lvl_vcycle.t),
    ("3lvl_kcycle", mg_3lvl_kcycle.t),
    #    ("4lvl_vcycle", mg_4lvl_vcycle.t),
    #    ("4lvl_kcycle", mg_4lvl_kcycle.t),
]:
    g.message(name + ":")
    for lvl in reversed(range(len(t))):
        g.message(t[lvl])

# print average iteration counts / time per level
g.message("Average iteration counts of inner solvers")
for name, h in [
    #    ("2lvl_vcycle", mg_2lvl_vcycle.history),
    #    ("2lvl_kcycle", mg_2lvl_kcycle.history),
    #    ("3lvl_vcycle", mg_3lvl_vcycle.history),
    ("3lvl_kcycle", mg_3lvl_kcycle.history),
    #    ("4lvl_vcycle", mg_4lvl_vcycle.history),
    #    ("4lvl_kcycle", mg_4lvl_kcycle.history),
]:
    for lvl in reversed(range(len(h))):
        for k, v in h[lvl].items():
            stats = list(map(lambda l: sum(l) / len(l), zip(*v)))
            if stats:
                g.message(f"{name}: lvl {lvl}: {k:10s} = {int(stats[0])}")
