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
import cgpt
import gpt


def ferm_to_prop(p, f, s, c):
    assert len(f.v_obj) == 1 and len(p.v_obj) == 1
    return cgpt.util_ferm2prop(f.v_obj[0], p.v_obj[0], s, c, True)

def prop_to_ferm(f, p, s, c):
    assert len(f.v_obj) == 1 and len(p.v_obj) == 1
    return cgpt.util_ferm2prop(f.v_obj[0], p.v_obj[0], s, c, False)

def reunitize(U):
    if type(U) == list:
        for a in U:
            reunitize(a)
    elif type(U) == gpt.lattice:
        assert len(U.v_obj) == 1
        return cgpt.util_reunitize(U.v_obj[0], U.v_obj[0])
    else:
        assert False
