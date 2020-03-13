/*
  CGPT

  Authors: Christoph Lehner 2020
*/

template<class T> class cgpt_Lattice;

template<typename T>
cgpt_Lattice<T>* compatible(cgpt_Lattice_base* other) {
  ASSERT(typeid(T).name() == other->type());
  return (cgpt_Lattice<T>*)other;
}

template<class T>
class cgpt_Lattice : public cgpt_Lattice_base {
public:
  Lattice<T> l;

  typedef typename Lattice<T>::vector_object vobj;
  typedef typename vobj::scalar_object sobj;
  typedef typename Lattice<T>::vector_type vCoeff_t;
  typedef typename Lattice<T>::scalar_type Coeff_t;

  cgpt_Lattice(GridCartesian* grid) : l(grid) {
  }

  virtual ~cgpt_Lattice() {
    //std::cout << "Deallocate" << std::endl;
  }

  cgpt_Lattice_base* create_lattice_of_same_type() {
    return new cgpt_Lattice<T>((GridCartesian*)l.Grid());
  }

  virtual std::string type() {
    return typeid(T).name();
  }

  virtual PyObject* to_decl() {   
    return PyTuple_Pack(3,PyLong_FromVoidPtr(this),
			PyUnicode_FromString(::get_otype(l)),
			PyUnicode_FromString(::get_prec(l)));
  }

  virtual RealD axpy_norm(ComplexD a, cgpt_Lattice_base* x, cgpt_Lattice_base* y) {
    return ::axpy_norm(l,(Coeff_t)a,compatible<T>(x)->l,compatible<T>(y)->l);
  }

  virtual RealD norm2() {
    return ::norm2(l);
  }

  virtual ComplexD innerProduct(cgpt_Lattice_base* other) {
    return ::innerProduct(l,compatible<T>(other)->l);
  }

  // ac == { true : add result to dst, false : replace dst }
  virtual cgpt_Lattice_base* mul(cgpt_Lattice_base* dst, bool ac, cgpt_Lattice_base* b, int unary_a, int unary_b, int unary_expr) {
    return cgpt_lattice_mul(dst,ac,unary_a,l,unary_b,b,unary_expr);
  }

  virtual cgpt_Lattice_base* compatible_linear_combination(cgpt_Lattice_base* dst,bool ac, std::vector<cgpt_lattice_term>& f, int unary_factor, int unary_expr) {
    return cgpt_compatible_linear_combination(l,dst,ac,f,unary_factor,unary_expr);
  }

  virtual void copy_from(cgpt_Lattice_base* _src) {
    cgpt_Lattice<T>* src = compatible<T>(_src);
    l = src->l;
  }

  virtual void cshift_from(cgpt_Lattice_base* _src, int dir, int off) {
    cgpt_Lattice<T>* src = compatible<T>(_src);
    l = Cshift(src->l, dir, off);
  }

  virtual void set_val(std::vector<int>& coor, ComplexD& val) {
    int nc = (int)coor.size();
    if (!nc && abs(val) == 0.0) {
      l = Zero();
    } else {
      cgpt_lattice_poke_value(l,coor,val);
    }
  }

  virtual PyObject* to_str() {
    std::stringstream st;
    st << l;
    return PyUnicode_FromString(st.str().c_str());
  }
  
};