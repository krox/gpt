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
import gpt, numpy

# matrix exponential
def exp(i, order=19):
    if i.grid.precision != gpt.double:
        x = gpt.convert(i, gpt.double)
    else:
        x = gpt.copy(i)
    n = gpt.norm2(x) ** 0.5 / x.grid.gsites
    maxn = 0.05
    ns = 0
    if n > maxn:
        ns = int(numpy.log2(n / maxn))
        x /= 2 ** ns
    o = gpt.lattice(x)
    o[:] = 0
    nfac = 1.0
    xn = gpt.copy(x)
    o[:] = numpy.identity(o.otype.shape[0], o.grid.precision.complex_dtype)
    o += xn
    for j in range(2, order + 1):
        nfac /= j
        xn @= xn * x
        o += xn * nfac
    for j in range(ns):
        o @= o * o
    if i.grid.precision != gpt.double:
        r = gpt.lattice(i)
        gpt.convert(r, o)
        o = r
    return o


def log(i, convergence_threshold=0.5):
    # i = n*(1 + x), log(i) = log(n) + log(1+x)
    # x = i/n - 1, |x|^2 = <i/n - 1, i/n - 1> = |i|^2/n^2 + |1|^2 - (<i,1> + <1,i>)/n
    # d/dn |x|^2 = -2 |i|^2/n^3 + (<i,1> + <1,i>)/n^2 = 0 -> 2|i|^2 == n (<i,1> + <1,i>)
    if i.grid.precision != gpt.double:
        x = gpt.convert(i, gpt.double)
    else:
        x = gpt.copy(i)
    I = numpy.identity(x.otype.shape[0], x.grid.precision.complex_dtype)
    lI = gpt.lattice(x)
    lI[:] = I
    n = gpt.norm2(x) / gpt.inner_product(x, lI).real
    x /= n
    x -= lI
    n2 = gpt.norm2(x) ** 0.5 / x.grid.gsites
    order = 8 * int(16 / (-numpy.log10(n2)))
    assert n2 < convergence_threshold
    o = gpt.copy(x)
    xn = gpt.copy(x)
    for j in range(2, order + 1):
        xn @= xn * x
        o -= xn * ((-1.0) ** j / j)
    o += lI * numpy.log(n)
    if i.grid.precision != gpt.double:
        r = gpt.lattice(i)
        gpt.convert(r, o)
        o = r
    return o

def inv(i):
    """ site-local matrix inversion for color-matrix/spin-matrix/spin-color-matrix """
    tmp = i[:]  # get numpy array of local sites

    # simple: just a color- or spin-matrix
    if type(i.otype) == gpt.ot_matrix_su3_fundamental or type(i.otype) == gpt.ot_matrix_color or type(i.otype) == gpt.ot_matrix_spin:
        tmp = numpy.linalg.inv(tmp)
    # little more complex: combined color-spin-matrix
    elif type(i.otype) == gpt.ot_matrix_spin_color:
        # GPT uses shape = (Ns, Ns, Nc, Nc)
        # for numpy inversion we need shape = (Ns*Nc, Ns*Nc)
        ns = i.otype.shape[0]  # spin dimensions (i.e. 4)
        nc = i.otype.shape[2]  # color dimensions (i.e. 3)
        sites = tmp.shape[0]   # number of sites on this rank
        tmp = numpy.transpose(tmp, (0, 1, 3, 2, 4))
        tmp = tmp.reshape((sites, ns * nc, ns * nc))
        tmp = numpy.linalg.inv(tmp)
        tmp = tmp.reshape((sites, ns, nc, ns, nc))
        tmp = numpy.transpose(tmp, (0, 1, 3, 2, 4))
        tmp = numpy.copy(tmp, order='C')  # GPT expects this order
    else:
        assert False

    o = gpt.lattice(i)
    o[:] = tmp
    return o

def det(i):
    """ site-local matrix determinant for color-matrix """
    tmp = i[:]
    assert len(i.otype.shape) == 2
    o = gpt.complex(i.grid)
    o[:] = numpy.linalg.det(tmp)
    return o
