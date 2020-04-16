#!/usr/bin/env python3
#
# Authors: Christoph Lehner 2020
#
# Desc.: Illustrate core concepts and features
#
import gpt as g
import sys

# load configuration
fine_grid=g.grid([8,8,8,16],g.single)
coarse_grid=g.grid([2,2,2,4],g.single)

# basis
n=30
basis=[ g.vcolor(fine_grid) for i in range(n) ]
g.random.seed("block_seed_string_13")
g.random.cnormal(basis)
for i in range(2):
    g.message("Ortho step %d" % i)
    g.block.orthogonalize(coarse_grid,basis)

# test coarse vector
lcoarse=g.vcomplex(coarse_grid,n)
g.random.cnormal(lcoarse)

# temporary fine and coarse vectors
tmpf=g.lattice(basis[0])
lcoarse2=g.lattice(lcoarse)

# coarse-to-fine-to-coarse
g.block.promote(lcoarse,tmpf,basis)
g.block.project(lcoarse2,tmpf,basis)

# report error
g.message(g.norm2(lcoarse-lcoarse2))
