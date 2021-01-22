#
#    GPT - Grid Python Toolkit
#    Copyright (C) 2020  Christoph Lehner (christoph.lehner@ur.de, https://github.com/lehner/gpt)
#                  2020  Daniel Richtmann (daniel.richtmann@ur.de)
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


class cg:
    @g.params_convention(eps=1e-15, maxiter=1000000)
    def __init__(self, params):
        self.params = params
        self.eps = params["eps"]
        self.maxiter = params["maxiter"]
        self.history = None

    def __call__(self, mat):

        otype, grid, cb = None, None, None
        if type(mat) == g.matrix_operator:
            otype, grid, cb = mat.otype, mat.grid, mat.cb
            mat = mat.mat
            # remove wrapper for performance benefits

        def inv(psi, src):
            assert src != psi
            self.history = []
            verbose = g.default.is_verbose("cg")
            t = g.timer("cg")
            t("setup")
            p, mmp, r = g.copy(src), g.copy(src), g.copy(src)
            mat(mmp, psi)  # in, out
            d = g.inner_product(psi, mmp).real
            b = g.norm2(mmp)
            r @= src - mmp
            p @= r
            a = g.norm2(p)
            cp = a
            ssq = g.norm2(src)
            if ssq == 0.0:
                assert a != 0.0  # need either source or psi to not be zero
                ssq = a
            rsq = self.eps ** 2.0 * ssq
            for k in range(1, self.maxiter + 1):
                c = cp
                t("mat")
                mat(mmp, p)
                t("inner")
                dc = g.inner_product(p, mmp)
                d = dc.real
                a = c / d
                t("axpy_norm")
                cp = g.axpy_norm2(r, -a, mmp, r)
                t("linearcomb")
                b = cp / c
                psi += a * p
                p @= b * p + r
                t("other")
                self.history.append(cp)
                if verbose:
                    g.message("cg: res^2[ %d ] = %g, target = %g" % (k, cp, rsq))
                if cp <= rsq:
                    if verbose:
                        t()
                        g.message(
                            "cg: converged in %d iterations, took %g s" % (k, t.total)
                        )
                        g.message(t)
                    break

        return g.matrix_operator(
            mat=inv,
            inv_mat=mat,
            otype=otype,
            accept_guess=(True, False),
            grid=grid,
            cb=cb,
        )
