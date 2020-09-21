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
import numpy as np


class fgmres:
    @g.params_convention(
        eps=1e-15, maxiter=1000000, restartlen=20, checkres=True, prec=None
    )
    def __init__(self, params):
        self.params = params
        self.eps = params["eps"]
        self.maxiter = params["maxiter"]
        self.restartlen = params["restartlen"]
        self.checkres = params["checkres"]
        self.prec = params["prec"]
        self.history = None

    @g.params_convention()
    def modified(self, params):
        return fgmres({**self.params, **params})

    def qr_update(self, s, c, H, gamma, i):
        # apply previous givens to matrix
        for j in range(i):
            tmp = -s[j] * H[j, i] + c[j] * H[j + 1, i]
            H[j, i] = np.conjugate(c[j]) * H[j, i] + np.conjugate(s[j]) * H[j + 1, i]
            H[j + 1, i] = tmp

        # compute new rotation matrix
        den = (np.absolute(H[i, i]) ** 2 + np.absolute(H[i + 1, i]) ** 2) ** 0.5
        c[i] = H[i, i] / den
        s[i] = H[i + 1, i] / den

        # apply new givens to matrix
        H[i, i] = den
        H[i + 1, i] = 0.0

        # apply new givens to vector
        gamma[i + 1] = -s[i] * gamma[i]
        gamma[i] *= np.conjugate(c[i])

    def update_psi(self, psi, gamma, H, y, V, i):
        # backward substitution
        for j in reversed(range(i + 1)):
            y[j] = (gamma[j] - np.dot(H[j, j + 1 : i + 1], y[j + 1 : i + 1])) / H[j, j]

        for j in range(i + 1):
            psi += y[j] * V[j]

    def restart(self, mat, psi, mmpsi, src, r, V, Z, gamma):
        r2 = self.calc_res(mat, psi, mmpsi, src, r)
        gamma[0] = r2 ** 0.5
        V[0] @= r / gamma[0]
        if Z is not None:
            for z in Z:
                z[:] = 0.0
        return r2

    def calc_res(self, mat, psi, mmpsi, src, r):
        mat(mmpsi, psi)
        return g.axpy_norm2(r, -1.0, mmpsi, src)

    def __call__(self, mat):

        otype, grid, cb = None, None, None
        if type(mat) == g.matrix_operator:
            otype, grid, cb = mat.otype, mat.grid, mat.cb
            mat = mat.mat
            # remove wrapper for performance benefits

        prec = self.prec(mat) if self.prec is not None else None

        def inv(psi, src):
            self.history = []
            # verbosity
            verbose = g.default.is_verbose("fgmres")

            # timing
            t = g.timer("fgmres")
            t("setup")

            # parameters
            rlen = self.restartlen

            # tensors
            dtype = g.double.complex_dtype
            H = np.zeros((rlen + 1, rlen), dtype)
            c = np.zeros((rlen + 1), dtype)
            s = np.zeros((rlen + 1), dtype)
            y = np.zeros((rlen + 1), dtype)
            gamma = np.zeros((rlen + 1), dtype)

            # fields
            mmpsi, r = (
                g.copy(src),
                g.copy(src),
            )
            V = [g.lattice(src) for i in range(rlen + 1)]
            Z = (
                [g.lattice(src) for i in range(rlen + 1)] if prec is not None else None
            )  # save vectors if unpreconditioned

            # initial residual
            r2 = self.restart(mat, psi, mmpsi, src, r, V, Z, gamma)

            # source
            ssq = g.norm2(src)
            if ssq == 0.0:
                assert r2 != 0.0  # need either source or psi to not be zero
                ssq = r2

            # target residual
            rsq = self.eps ** 2.0 * ssq

            for k in range(self.maxiter):
                # iteration within current krylov space
                i = k % rlen

                # iteration criteria
                reached_maxiter = k + 1 == self.maxiter
                need_restart = i + 1 == rlen

                t("prec")
                if prec is not None:
                    prec(Z[i], V[i])

                t("mat")
                if prec is not None:
                    mat(V[i + 1], Z[i])
                else:
                    mat(V[i + 1], V[i])

                t("ortho")
                g.default.push_verbose("orthogonalize", False)
                g.orthogonalize(V[i + 1], V[0 : i + 1], H[:, i])
                g.default.pop_verbose()

                t("linalg")
                H[i + 1, i] = g.norm2(V[i + 1]) ** 0.5
                if H[i + 1, i] == 0.0:
                    g.message("fgmres: breakdown, H[%d, %d] = 0" % (i + 1, i))
                    break
                V[i + 1] /= H[i + 1, i]

                t("qr")
                self.qr_update(s, c, H, gamma, i)

                t("other")
                r2 = np.absolute(gamma[i + 1]) ** 2
                self.history.append(r2)

                if verbose:
                    g.message(
                        "fgmres: res^2[ %d, %d ] = %g, target = %g" % (k, i, r2, rsq)
                    )

                if r2 <= rsq or need_restart or reached_maxiter:
                    t("update_psi")
                    if prec is not None:
                        self.update_psi(psi, gamma, H, y, Z, i)
                    else:
                        self.update_psi(psi, gamma, H, y, V, i)
                    comp_res = r2 / ssq

                    if r2 <= rsq:
                        if verbose:
                            t()
                            g.message(
                                "fgmres: converged in %d iterations, took %g s"
                                % (k + 1, t.dt["total"])
                            )
                            g.message(t)
                            if self.checkres:
                                res = self.calc_res(mat, psi, mmpsi, src, r) / ssq
                                g.message(
                                    "fgmres: computed res = %g, true res = %g, target = %g"
                                    % (comp_res ** 0.5, res ** 0.5, self.eps)
                                )
                            else:
                                g.message(
                                    "fgmres: computed res = %g, target = %g"
                                    % (comp_res ** 0.5, self.eps)
                                )
                        break

                    if reached_maxiter:
                        if verbose:
                            t()
                            g.message(
                                "fgmres: did NOT converge in %d iterations, took %g s"
                                % (k + 1, t.dt["total"])
                            )
                            g.message(t)
                            if self.checkres:
                                res = self.calc_res(mat, psi, mmpsi, src, r) / ssq
                                g.message(
                                    "fgmres: computed res = %g, true res = %g, target = %g"
                                    % (comp_res ** 0.5, res ** 0.5, self.eps)
                                )
                            else:
                                g.message(
                                    "fgmres: computed res = %g, target = %g"
                                    % (comp_res ** 0.5, self.eps)
                                )
                        break

                    if need_restart:
                        t("restart")
                        r2 = self.restart(mat, psi, mmpsi, src, r, V, Z, gamma)
                        if verbose:
                            g.message("fgmres: performed restart")

        return g.matrix_operator(
            mat=inv,
            inv_mat=mat,
            otype=otype,
            accept_guess=(True, False),
            grid=grid,
            cb=cb,
        )
