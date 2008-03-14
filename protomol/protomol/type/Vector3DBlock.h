/*  -*- c++ -*-  */
#ifndef VECTOR3DBLOCK_H
#define VECTOR3DBLOCK_H

#include <vector>

#include <protomol/type/Vector3D.h>
#include <protomol/base/Proxy.h>


namespace ProtoMol {
  //_____________________________________________________________ Vector3DBlock
  /**
   * Container holding a vector (array) of 3D coordinates/vectors
   */
  class Vector3DBlock : public Proxy {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types & enum's
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    typedef std::vector<Vector3D>::value_type value_type;
    typedef std::vector<Vector3D>::pointer pointer;
    typedef std::vector<Vector3D>::reference reference;
    typedef std::vector<Vector3D>::const_reference const_reference;
    typedef std::vector<Vector3D>::size_type size_type;
    typedef std::vector<Vector3D>::difference_type difference_type;
    typedef std::vector<Vector3D>::iterator iterator;
    typedef std::vector<Vector3D>::const_iterator const_iterator;
    typedef std::vector<Vector3D>::reverse_iterator reverse_iterator;
    typedef std::vector<Vector3D>::const_reverse_iterator
    const_reverse_iterator;
  private:
    enum {LIMIT = 30};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Vector3DBlock()  : Proxy(), vec() {}
    explicit Vector3DBlock(size_type n) : Proxy(), vec(n) {}
    Vector3DBlock(size_type n, const Vector3D &t) : Proxy(), vec(n, t) {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  New methods of class Vector3DBlock from std::vector
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    iterator begin() {return vec.begin();}
    iterator end() {return vec.end();}
    const_iterator begin() const {return vec.begin();}
    const_iterator end() const {return vec.end();}
    reverse_iterator rbegin() {return vec.rbegin();}
    reverse_iterator rend() {return vec.rend();}
    const_reverse_iterator rbegin() const {return vec.rbegin();}
    const_reverse_iterator rend() const {return vec.rend();}
    size_type size() const {return vec.size();}
    size_type max_size() const {return vec.max_size();}
    size_type capacity() const {return vec.capacity();}
    bool empty() const {return vec.empty();}
    reference operator[](size_type n) {return vec[n];}
    const_reference operator[](size_type n) const {return vec[n];}
    void reserve(size_t n) {vec.reserve(n);}
    reference front() {return vec.front();}
    const_reference front() const {return vec.front();}
    reference back() {return vec.back();}
    const_reference back() const {return vec.back();}
    void push_back(const Vector3D &t) {return vec.push_back(t);}
    void pop_back() {return vec.pop_back();}
    void swap(Vector3DBlock &x) {vec.swap(x.vec);}
    iterator insert(iterator pos, const Vector3D &x) {
      return vec.insert(pos, x);
    }
    void insert(iterator pos, size_type n, const Vector3D &x) {
      vec.insert(pos, n, x);
    }
    iterator erase(iterator pos) {return vec.erase(pos);}
    iterator erase(iterator first, iterator last) {
      return vec.erase(first, last);
    }
    void clear() {vec.clear();}
    void resize(size_type n, Vector3D t = Vector3D(0.0, 0.0, 0.0)) {
      vec.resize(n, t);
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Vector3DBlock
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /// Clear (set to zero) each element of the vector.
    void zero(int n = -1) {
      if (n >= 0)
        vec.resize(n);
      const unsigned int count = vec.size();
      for (unsigned int i = 0; i < count; ++i)
        vec[i] = Vector3D(0, 0, 0);
    }

    /// copy one Vector3DBlock to another one.
    Vector3DBlock &intoAssign(const Vector3DBlock &x) {
      vec = x.vec;
      return *this;
    }

    /// Add another Vector3DBlock into this one.
    Vector3DBlock &intoAdd(const Vector3DBlock &x) {
      const unsigned int count = vec.size();
      for (unsigned int i = 0; i < count; ++i)
        vec[i] += x.vec[i];

      return *this;
    }

    /// Subtract another Vector3DBlock from this one.
    Vector3DBlock &intoSubtract(const Vector3DBlock &x) {
      const unsigned int count = vec.size();
      for (unsigned int i = 0; i < count; ++i)
        vec[i] -= x.vec[i];

      return *this;
    }

    /// Add a scalar multiple of another Vector3DBlock into this one.
    Vector3DBlock &intoWeightedAdd(Real weight, const Vector3DBlock &x) {
      const unsigned int count = vec.size();
      for (unsigned int i = 0; i < count; ++i)
        vec[i] += x.vec[i] * weight;

      return *this;
    }

    /// Assign a scalar multiple of another Vector3DBlock into this one.
    Vector3DBlock &intoWeighted(Real weight, const Vector3DBlock &x) {
      const unsigned int count = vec.size();
      for (unsigned int i = 0; i < count; ++i)
        vec[i] = x.vec[i] * weight;

      return *this;
    }

    /// Subtract a scalar multiple of another Vector3DBlock from this one.
    Vector3DBlock &intoWeightedSubtract(Real weight, const Vector3DBlock &x) {
      const unsigned int count = vec.size();
      for (unsigned int i = 0; i < count; ++i)
        vec[i] -= x.vec[i] * weight;

      return *this;
    }

    /// Compute the boinding box over all elements
    void boundingbox(Vector3D &minbb, Vector3D &maxbb) const;

    /// Compute the sum over all elements
    Vector3D sum() const;

    /// Compute regression plane by SVD
    bool fitplane(Vector3D &normal, Real &d, Real &err, int limit =
                    LIMIT) const;

  public:
    friend bool operator==(const Vector3DBlock &a, const Vector3DBlock &b);
    friend bool operator<(const Vector3DBlock &a, const Vector3DBlock &b);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    std::vector<Vector3D> vec;
  };

  inline bool operator==(const Vector3DBlock &a, const Vector3DBlock &b) {
    return a.vec == b.vec;
  }
  inline bool operator<(const Vector3DBlock &a, const Vector3DBlock &b) {
    return a.vec == b.vec;
  }
}

inline void swap(ProtoMol::Vector3DBlock &a, ProtoMol::Vector3DBlock &b) {
  a.swap(b);
}

#endif /* VECTOR3DBLOCK_H */
