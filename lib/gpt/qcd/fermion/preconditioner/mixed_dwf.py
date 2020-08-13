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
#
# This file implements the MADWF algorithm of
#
#  https://arxiv.org/abs/1111.5059
#
# to create an approximate solution of an outer DWF inverse
# using a different inner DWF inverse.
#
import gpt as g


class mixed_dwf:
    def __init__(self, solver, solver_pv, dwf_inner):
        self.solver = solver
        self.solver_pv = solver_pv
        self.dwf_inner = dwf_inner
        self.dwf_inner_pv = dwf_inner.modified(mass=1.0)

    def __call__(self, dwf_outer):
        dwf_inner = self.dwf_inner
        dwf_inner_pv = self.dwf_inner_pv
        dwf_outer_pv = dwf_outer.modified(mass=1.0)

        inv_dwf_outer_pv = self.solver_pv(dwf_outer_pv)
        inv_dwf_inner = self.solver(dwf_inner)

        def sep(x):
            return g.separate(x, dimension=0)

        def mrg(x):
            return g.merge(x, dimension=0)

        def _P(dst, src, offset):
            src_s = sep(src)
            Ls = len(src_s)
            Pplus = 0.5 * (g.gamma["I"] + g.gamma[5])
            Pminus = 0.5 * (g.gamma["I"] - g.gamma[5])
            dst @= mrg(
                [
                    g(Pminus * src_s[s] + Pplus * src_s[(s + Ls + offset) % Ls])
                    for s in range(Ls)
                ]
            )

        def _P_mat(dst, src):
            _P(dst, src, 1)

        def _P_inv_mat(dst, src):
            _P(dst, src, -1)

        P = g.matrix_operator(
            mat=_P_mat, adj_mat=_P_inv_mat, adj_inv_mat=_P_mat, inv_mat=_P_inv_mat
        )

        def inv(dst_outer, src_outer):
            Ls_inner = dwf_inner.F_grid.fdimensions[0]
            zero4d = g.lattice(dwf_outer.U_grid, src_outer.otype)
            zero4d[:] = 0
            c_s = sep(g.adj(P) * inv_dwf_outer_pv * src_outer)
            y0prime = sep(
                g.adj(P)
                * inv_dwf_inner
                * dwf_inner_pv
                * P
                * mrg([c_s[0]] + [zero4d] * (Ls_inner - 1))
            )[0]
            dst_outer @= P * mrg(
                [y0prime]
                + sep(
                    g.adj(P)
                    * inv_dwf_outer_pv
                    * dwf_outer
                    * P
                    * mrg([g(-y0prime)] + c_s[1:])
                )[1:]
            )

        return g.matrix_operator(
            mat=inv,
            inv_mat=dwf_outer,
            otype=dwf_outer.otype,
            zero=(True, False),
            grid=dwf_outer.F_grid,
            cb=None,
        )
