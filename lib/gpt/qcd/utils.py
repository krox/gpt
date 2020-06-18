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
    assert(len(f.v_obj) == 1 and len(p.v_obj) == 1)
    return cgpt.util_ferm2prop(f.v_obj[0],p.v_obj[0],s,c,True)

def prop_to_ferm(f, p, s, c):
    assert(len(f.v_obj) == 1 and len(p.v_obj) == 1)
    return cgpt.util_ferm2prop(f.v_obj[0],p.v_obj[0],s,c,False)

def project_on_group(dst, src=None):
    if src is None:
        return project_on_group(dst, dst)
    elif type(dst) == list and type(src) == list:
        assert len(dst) == len(src)
        for a,b in zip(dst,src):
            project_on_group(a,b)
    elif type(dst) == gpt.lattice and type(src) == gpt.lattice:
        assert len(dst.v_obj) == 1 and len(src.v_obj) == 1
        return cgpt.util_project_on_group(dst.v_obj[0], src.v_obj[0])
    else:
        assert False

def project_on_algebra(dst, src=None):
    if src is None:
        return project_on_algebra(dst, dst)
    elif type(dst) == list and type(src) == list:
        assert len(dst) == len(src)
        for a,b in zip(dst,src):
            project_on_algebra(a,b)
    elif type(dst) == gpt.lattice and type(src) == gpt.lattice:
        assert len(dst.v_obj) == 1 and len(src.v_obj) == 1
        return cgpt.util_project_on_algebra(dst.v_obj[0], src.v_obj[0])
    else:
        assert False
