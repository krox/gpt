/*
    GPT - Grid Python Toolkit
    Copyright (C) 2020  Christoph Lehner (christoph.lehner@ur.de, https://github.com/lehner/gpt)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/
template<typename T>
class cgpt_fermion_operator : public cgpt_fermion_operator_base {
public:
  T* op;

  cgpt_fermion_operator(T* _op) : op(_op) {
  }

  virtual ~cgpt_fermion_operator() {
    delete op;
  }

  virtual RealD unary(int opcode, cgpt_Lattice_base* in, cgpt_Lattice_base* out) {
    return cgpt_fermion_operator_unary<T>(*op,opcode,in,out);
  }

  virtual void deriv(std::array<cgpt_Lattice_base*,4> force, cgpt_Lattice_base* X, cgpt_Lattice_base* Y, int dag) {
    typedef typename T::FermionField::vector_object vobj_ferm; // (vectorized) spin-vector x color-vector
    typedef typename T::GaugeLinkField::vector_object vobj_gauge; // (vectorized) color-matrix

    // Ugly workaround: Grid's FermOp.MDeriv() interface expects the gauge-config as a single
    // "lorentz-vector x color-matrix" field. But GPT uses 4 separate "color-matrix" fields instead.
    // So we need a temporary here, which is kinda ugly
    typename T::GaugeField tmp(X->get_grid());
    op->MDeriv(tmp, compatible<vobj_ferm>(X)->l, compatible<vobj_ferm>(Y)->l, dag);
    for(int mu = 0; mu < 4; ++mu)
        compatible<vobj_gauge>(force[mu])->l = peekLorentz(tmp, mu);
  }

};
