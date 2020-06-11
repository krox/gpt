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
import gpt

# First EO preconditioning (upper triangular left)
#
#      ( EE EO )   ( 1  EO OO^-1 ) ( N 0 ) ( EE 0  )
#  M = ( OE OO ) = ( 0    1      ) ( 0 1 ) ( OE OO )
#
#  N = 1 - EO OO^-1 OE EE^-1
#
#  Verify:  ( N   EO OO^-1  ) ( EE 0  )   ( N EE + EO OO^-1 OE    EO )   ; N EE + EO OO^-1 OE = EE - EO OO^-1 OE + EO OO^-1 OE = EE
#           ( 0   1         ) ( OE OO ) = ( OE                    OO )
#
# Then:
#
#        ( EE 0  )^-1 ( N^dag^-1 (N^dag N)  0 )^-1 ( 1     EO OO^-1  )^-1
# M^-1 = ( OE OO )    ( 0                   1 )    ( 0         1     )
#
#        ( EE^-1              0     )  ( (N^dag N)^-1 N^dag   0 )  ( 1   - EO OO^-1 )
#      = ( -OO^-1 OE EE^-1    OO^-1 )  ( 0                    1 )  ( 0       1      )
#
# M^-1 = L (N^dag N)^-1 R + S
#
# R = N^dag ( 1   - EO OO^-1 )
#
#     ( EE^-1           ) 
# L = ( -OO^-1 OE EE^-1 )
#
#     ( 0   0     )
# S = ( 0   OO^-1 )
#

class eo2:
    def __init__(self, op):
        self.op = op
        self.otype = op.otype
        self.F_grid_eo = op.F_grid_eo
        self.F_grid = op.F_grid
        self.tmp = gpt.lattice(self.F_grid_eo,self.otype)
        self.tmp2 = gpt.lattice(self.F_grid_eo,self.otype)
        self.ImportPhysicalFermionSource = self.op.ImportPhysicalFermionSource
        self.ExportPhysicalFermionSolution = self.op.ExportPhysicalFermionSolution

        def _N(oe, ie):
            self.op.Mooee.inv()(self.tmp2,ie)
            self.op.Meooe(oe,self.tmp2)
            self.op.Mooee.inv()(self.tmp2,oe)
            self.op.Meooe(oe,self.tmp2)
            oe @= ie - oe

        def _NDag(oe, ie):
            self.op.Meooe.adj()(self.tmp2,ie)
            self.op.Mooee.inv().adj()(oe,self.tmp2)
            self.op.Meooe.adj()(self.tmp2,oe)
            self.op.Mooee.inv().adj()(oe,self.tmp2)
            oe @= ie - oe

        def _NDagN(oe, ie):
            _N(self.tmp,ie)
            _NDag(oe,self.tmp)

        self.N = gpt.matrix_operator(mat = _N, adj_mat = _NDag, otype = op.otype)
        self.NDagN = gpt.matrix_operator(mat = _NDagN, adj_mat = _NDagN, otype = op.otype)

    def R(self, oe, ie, io):
        self.op.Mooee.inv()(self.tmp,io)
        self.op.Meooe(oe,self.tmp)
        self.tmp @= ie - oe
        self.N.adj()(oe,self.tmp)

    def L(self, oe, oo, ie):
        self.op.Mooee.inv()(oe,ie)
        self.op.Meooe(self.tmp,oe)
        self.op.Mooee.inv()(oo,self.tmp)
        oo @= - oo

    def S(self, oe, oo, ie, io):
        self.op.Mooee.inv()(oo,io)
        oe[:]=0

