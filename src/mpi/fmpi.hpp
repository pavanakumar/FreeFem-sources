#ifndef FMPI_HPP

#define FMPI_HPP

// Remark: on mipich MPI_Comm, MPI_Resquest, MPI_Group, MPI_Op are int
// => encapsulation

//static long verbosity = 1000;
template<class MPI_type, int DIFF>
struct fMPI {
  MPI_type v;
  operator MPI_type &() { return v; }
  operator MPI_type *() { return &v; }
  operator MPI_type () const { return v; }

  // MPI_type * operator &() { return &v; }
  void operator = (const MPI_type vv) { v=vv; }
  fMPI(const MPI_type vv=0) : v(vv) {}
  bool operator != (MPI_type vv) const { return vv != v; }
  bool operator == (MPI_type vv) const { return vv == v; }
};

// the encapsulation for the for MPI type (int on mpich )
typedef fMPI<MPI_Comm, 1> fMPI_Comm;
typedef fMPI<MPI_Group, 2> fMPI_Group;
typedef fMPI<MPI_Request, 3> fMPI_Request;
typedef fMPI<MPI_Op, 4> fMPI_Op;
// end of encapsulation ..

#endif

