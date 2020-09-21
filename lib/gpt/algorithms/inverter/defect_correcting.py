#
#    GPT - Grid Python Toolkit
#    Copyright (C) 2020  Christoph Lehner (christoph.lehner@ur.de, https://github.com/lehner/gpt)
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
import gpt as g
import sys


class defect_correcting:

    #
    # Less numerically stable (leads to suppression of critical low-mode space before inner_mat^{-1}):
    #
    # outer_mat = inner_mat (1 + eps defect)     (*)
    #
    # outer_mat^{-1} = (1 - eps defect + eps^2 defect^2 - eps^3 defect^3 + ...) inner_mat^{-1}   (geometric series)
    #
    # lhs = outer_mat^{-1} rhs
    #     = (1 - eps defect + eps^2 defect^2 - eps^3 defect^3 + ...) inner_mat^{-1} rhs
    #
    # Now defining lhs^{(0)} = inner_mat^{-1} rhs, we have
    #
    #     lhs = lhs^{(0)} - eps defect lhs^{(0)} + eps^2 defect^2 lhs^{(0)} + ...
    #
    # We can rearrange (*) to find
    #
    #     eps defect = inner_mat^{-1} outer_mat - 1
    #

    #
    # More numerically stable (implemented here):
    #
    # outer_mat =  (1 + eps defect) inner_mat     (*)
    #
    # outer_mat^{-1} = inner_mat^{-1} (1 - eps defect + eps^2 defect^2 - eps^3 defect^3 + ...)   (geometric series)
    #
    # lhs = outer_mat^{-1} rhs
    #     = inner_mat^{-1} (1 - eps defect + eps^2 defect^2 - eps^3 defect^3 + ...) rhs
    #
    # Now defining lhs^{(0)} = inner_mat^{-1} rhs, we have
    #
    #     lhs = lhs^{(0)} - inner_mat^{-1} eps defect rhs + ...
    #
    # We can rearrange (*) to find
    #
    #     eps defect = outer_mat inner_mat^{-1} - 1
    #

    @g.params_convention(eps=1e-15, maxiter=1000000)
    def __init__(self, inner_inverter, params):
        self.params = params
        self.eps = params["eps"]
        self.maxiter = params["maxiter"]
        self.history = None
        self.inner_inverter = inner_inverter

    def __call__(self, outer_mat):

        inner_inv_mat = self.inner_inverter(outer_mat)

        def inv(psi, src):

            # verbosity
            verbose = g.default.is_verbose("dci")
            t_start = g.time()

            # leading order
            n = len(src)
            _s = [g.copy(x) for x in src]
            _d = [g.copy(x) for x in psi]

            self.history = []
            for i in range(self.maxiter):

                # correction step
                t0 = g.time()
                for j in range(n):
                    _s[j] -= outer_mat * _d[j]
                t1 = g.time()
                _d = g.eval(inner_inv_mat * _s)
                t2 = g.time()
                for j in range(n):
                    psi[j] += _d[j]

                # true resid
                eps = max(
                    [g.norm2(outer_mat * psi[j] - src[j]) ** 0.5 for j in range(n)]
                )
                self.history.append(eps)

                if verbose:
                    g.message(
                        "Defect-correcting inverter:  eps[",
                        i,
                        "] =",
                        eps,
                        ".  Timing:",
                        t2 - t1,
                        "s (innver_inv), ",
                        t1 - t0,
                        "s (outer_mat)",
                    )

                if eps < self.eps:
                    if verbose:
                        g.message(
                            "Defect-correcting inverter: converged at iteration",
                            i,
                            "after",
                            g.time() - t_start,
                            "s",
                        )
                    break

        otype, grid, cb = None, None, None
        if type(outer_mat) == g.matrix_operator:
            otype, grid, cb = outer_mat.otype, outer_mat.grid, outer_mat.cb

        return g.matrix_operator(
            mat=inv,
            inv_mat=outer_mat,
            otype=otype,
            accept_guess=(True, False),
            grid=grid,
            cb=cb,
            accept_list=True,
        )
