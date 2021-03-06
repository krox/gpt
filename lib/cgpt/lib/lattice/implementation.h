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
template<class T>
class cgpt_Lattice : public cgpt_Lattice_base {
public:
  Lattice<T> l;

  typedef typename Lattice<T>::vector_object vobj;
  typedef typename vobj::scalar_object sobj;
  typedef typename Lattice<T>::vector_type vCoeff_t;
  typedef typename Lattice<T>::scalar_type Coeff_t;

  cgpt_Lattice(GridBase* grid) : l(grid) {
  }

  virtual ~cgpt_Lattice() {
    //std::cout << "Deallocate" << std::endl;
  }

  cgpt_Lattice_base* create_lattice_of_same_type() {
    return new cgpt_Lattice<T>(l.Grid());
  }

  virtual std::string type() {
    return typeid(T).name();
  }

  virtual int singlet_rank() {
    return ::singlet_rank(l);
  }

  virtual PyObject* to_decl() {   
    return PyTuple_Pack(3,PyLong_FromVoidPtr(this),
			PyUnicode_FromString(get_otype(l).c_str()),
			PyUnicode_FromString(get_prec(l).c_str()));
  }

  void set_to_zero() {
    l = Zero();
  }

  // use norm2 convention for squared norm, talked to Peter, Grid may also change to this cleaner notation
  virtual RealD axpy_norm2(ComplexD a, cgpt_Lattice_base* x, cgpt_Lattice_base* y) {
    return ::axpy_norm(l,(Coeff_t)a,compatible<T>(x)->l,compatible<T>(y)->l);
  }

  virtual void axpy(ComplexD a, cgpt_Lattice_base* x, cgpt_Lattice_base* y) {
    return ::axpy(l,(Coeff_t)a,compatible<T>(x)->l,compatible<T>(y)->l);
  }

  virtual RealD norm2() {
    return ::norm2(l);
  }

  virtual void rank_inner_product(ComplexD* res, std::vector<cgpt_Lattice_base*> & left, std::vector<cgpt_Lattice_base*> & right, long n_virtual, bool use_accelerator) {
    PVector<Lattice<T>> _left;
    PVector<Lattice<T>> _right;
    cgpt_basis_fill(_left,left);
    cgpt_basis_fill(_right,right);
    if (use_accelerator) {
      ::rankInnerProduct(res,_left,_right,n_virtual);
    } else {
      ::rankInnerProductCpu(res,_left,_right,n_virtual);
    }
  }

  virtual void inner_product_norm2(ComplexD& ip, RealD& a2, cgpt_Lattice_base* other) {
    ::innerProductNorm(ip,a2,l,compatible<T>(other)->l);
  }

  // ac == { true : add result to dst, false : replace dst }
  virtual cgpt_Lattice_base* mul(cgpt_Lattice_base* dst, bool ac, cgpt_Lattice_base* b, int unary_a, int unary_b, int unary_expr) {
    if (typeid(T) == typeid(iSinglet<vCoeff_t>) && b->type() != type()) {
      // singlet multiplication always commutes, can save half cost of instantiation
      return b->mul(dst,ac,this,unary_b,unary_a,unary_expr);
    }
    return cgpt_lattice_mul(dst,ac,unary_a,l,unary_b,b,unary_expr);
  }

  virtual cgpt_Lattice_base* compatible_linear_combination(cgpt_Lattice_base* dst,bool ac, std::vector<cgpt_lattice_term>& f, int unary_factor, int unary_expr) {
    return cgpt_compatible_linear_combination(l,dst,ac,f,unary_factor,unary_expr);
  }

  virtual cgpt_Lattice_base* matmul(cgpt_Lattice_base* dst, bool ac, PyArrayObject* b, std::string& bot, int unary_b, int unary_a, int unary_expr, bool reverse) {
    return cgpt_lattice_matmul(dst,ac,unary_a,l,b,bot,unary_b,unary_expr,reverse);
  }

  virtual cgpt_Lattice_base* gammamul(cgpt_Lattice_base* dst, bool ac, Gamma::Algebra gamma, int unary_a, int unary_expr, bool reverse) {
    return cgpt_lattice_gammamul(dst,ac,unary_a,l,gamma,unary_expr,reverse);
  }

  virtual void copy_from(cgpt_Lattice_base* _src) {
    cgpt_Lattice<T>* src = compatible<T>(_src);
    l = src->l;
  }

  virtual void fft_from(cgpt_Lattice_base* src, const std::vector<int> & dims, int sign) {
    FFT fft((GridCartesian*)l.Grid());
    Lattice<T> tmp = compatible<T>(src)->l;
    for (long i=0;i<dims.size();i++) {
      fft.FFT_dim(l,tmp,dims[i],sign);
      if (i != dims.size()-1)
	tmp = l;
    }
  }

  virtual void unary_from(cgpt_Lattice_base* src, PyObject* params) {
    cgpt_unary_from(l,compatible<T>(src)->l,params);
  }

  virtual void cshift_from(cgpt_Lattice_base* _src, int dir, int off) {
    cgpt_Lattice<T>* src = compatible<T>(_src);
    l = Cshift(src->l, dir, off);
  }

  virtual PyObject* sum() {
    return cgpt_numpy_export( ::sum(l) );
  }

  virtual PyObject* to_str() {
    return PyUnicode_FromString(cgpt_lattice_to_str(l).c_str());
  }

  virtual void convert_from(cgpt_Lattice_base* src) {
    cgpt_lattice_convert_from(l,src);
  }

  virtual PyObject* slice(int dim) {
    return cgpt_lattice_slice(l,dim);
  }

  virtual void ferm_to_prop(cgpt_Lattice_base* prop, int spin, int color, bool f2p) {
    cgpt_ferm_to_prop(l,prop,spin,color,f2p);
  }

  virtual void pick_checkerboard_from(int cb, cgpt_Lattice_base* src) {
    pickCheckerboard(cb, l, compatible<T>(src)->l);
  }

  virtual void set_checkerboard_from(cgpt_Lattice_base* src) {
    setCheckerboard(l, compatible<T>(src)->l);
  }

  virtual void change_checkerboard(int cb) {
    l.Checkerboard() = cb;
  }

  virtual int get_checkerboard() {
    return l.Checkerboard();
  }

  virtual void basis_rotate(std::vector<cgpt_Lattice_base*> &_basis,RealD* Qt,int j0, int j1, int k0,int k1,int Nm) {
    PVector<Lattice<T>> basis;
    cgpt_basis_fill(basis,_basis);
    cgpt_basis_rotate(basis,Qt,j0,j1,k0,k1,Nm);
  }

  virtual void basis_rotate(std::vector<cgpt_Lattice_base*> &_basis,ComplexD* Qt,int j0, int j1, int k0,int k1,int Nm) {
    PVector<Lattice<T>> basis;
    cgpt_basis_fill(basis,_basis);
    cgpt_basis_rotate(basis,Qt,j0,j1,k0,k1,Nm);
  }

  virtual void linear_combination(std::vector<cgpt_Lattice_base*> &_basis,ComplexD* Qt) {
    PVector<Lattice<T>> basis;
    cgpt_basis_fill(basis,_basis);
    cgpt_linear_combination(l,basis,Qt);
  }

  virtual PyObject* memory_view() {
    auto v = l.View(CpuWrite);
    size_t sz = v.size() * sizeof(v[0]);
    char* ptr = (char*)&v[0];
    v.ViewClose();

    // this marks Cpu as dirty, so data will be copied to Gpu; this is not fully safe
    // and the ViewClose should be moved to the destructor of the PyMemoryView object.
    // Do this in the same way as currently done in mview() in the future.
    return PyMemoryView_FromMemory(ptr,sz,PyBUF_WRITE);
  }

  virtual PyObject* memory_view_coordinates() {
    return cgpt_memory_view_coordinates(l.Grid(),l.Checkerboard());
  }

  virtual void describe_data_layout(long & Nsimd, long & word, long & simd_word, std::vector<long> & ishape) {
    GridBase* grid = l.Grid();
    Nsimd = grid->Nsimd();
    word = sizeof(sobj);
    simd_word = sizeof(Coeff_t);
    ishape.resize(0);
    cgpt_numpy_data_layout(sobj(),ishape);
    if (ishape.size() == 0) // treat complex numbers as 1d array with one element
      ishape.push_back(1);
  }
  
  virtual int get_numpy_dtype() {
    return infer_numpy_type(Coeff_t());
  }

  virtual cgpt_block_map_base* block_map(GridBase* coarse, 
					 std::vector<cgpt_Lattice_base*>& basis,
					 long basis_n_virtual, long basis_virtual_size, long basis_n_block,
					 cgpt_Lattice_base* mask) {

    ASSERT(basis.size() > 0 && basis.size() % basis_n_virtual == 0);

#define BASIS_SIZE(n) if (n == basis_virtual_size) { return new cgpt_block_map<T, iVSinglet ## n<vCoeff_t> >(coarse,basis,basis_n_virtual,basis_n_block,mask); }
#include "../basis_size.h"
#undef BASIS_SIZE
    
    { ERR("Unknown basis size %d",(int)basis_virtual_size); }

  }

  virtual GridBase* get_grid() {
    return l.Grid();
  }

  virtual PyObject* advise(std::string type) {
    if (type == "infrequent_use") {
      l.Advise() = AdviseInfrequentUse;
    } else {
      ERR("Unknown advise %s",type.c_str());
    }    
    return PyLong_FromLong(0);
  }

  virtual PyObject* prefetch(std::string type) {
    if (type == "accelerator") {
      //l.AcceleratorPrefetch();
    } else if (type == "host") {
      //l.HostPrefetch();
    } else {
      ERR("Unknown prefetch %s",type.c_str());
    }
    return PyLong_FromLong(0);
  }

};
