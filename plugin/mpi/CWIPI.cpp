#include <algorithm>
//ff-c++-LIBRARY-dep: cxx11 [petsc] cwipi mpi
//ff-c++-cpp-dep:
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <ff++.hpp>
#include <cwipi.h>
#if defined(WITH_petsc)
#include <petsc.h>
#endif

#include "AFunction_ext.hpp"
#include "AFunction.hpp"
extern KN<String>* pkarg;
MPI_Comm cwipi_comm;

using namespace Fem2D;
using namespace std;
std::vector<int> global_connectivity;
std::vector<int> global_connectivity_offset;

long cwipiComm(pcommworld const &comm) {
  MPI_Comm* mpi_comm = (MPI_Comm*)comm;
  if(mpi_comm && *mpi_comm != MPI_COMM_NULL)
      MPI_Comm_free(mpi_comm);
  MPI_Comm_dup(cwipi_comm, mpi_comm);
  return 0L;
}

long CwipiWaitIsSendFfpp(string *const &coupling_name, long const &request) {
  cwipi_wait_issend(coupling_name->c_str(), request);
  return 0L;
}

long CwipiWaitIrecvFfpp(string *const &coupling_name, long const &request) {
  cwipi_wait_irecv(coupling_name->c_str(), request);
  return 0L;
}
long CwipiLocateFfpp(string *const &coupling_id) {
  cwipi_locate(coupling_id->c_str());
  return 0L;
}
long CwipiGetNLocatedPointsFfpp(string *const &coupling_id) {
  cwipi_get_n_located_points(coupling_id->c_str());
  return 0L;
}
long CwipiGetNNotLocatedPointsFfpp(string *const &coupling_id) {
  cwipi_get_n_not_located_points(coupling_id->c_str());
  return 0L;
}

long CwipiDefineMeshFfpp(string *const &coupling_id, long const &n_vertex,
                         long const &n_element, KN<double> *const &coordinates,
                         KN<long> *const &connectivity_offset,
                         KN<long> *const &connectivity) {
  global_connectivity_offset.resize(connectivity_offset->size());
  for (int i = 0; i < connectivity_offset->size(); ++i)
    global_connectivity_offset[i] = (*connectivity_offset)[i];
  global_connectivity.resize(connectivity->size());
  for (int j = 0; j < connectivity->size(); ++j)
    global_connectivity[j] = (*connectivity)[j];
  cwipi_define_mesh(coupling_id->c_str(), n_vertex, n_element,
                    &((*coordinates)[0]), global_connectivity_offset.data(), global_connectivity.data());
  return 0L;
}
long CwipiHoDefineMeshFfpp(string *const &coupling_id, long const &n_vertex,
                           long const &n_element, long const &order,
                           KN<double> *const &coordinates,
                           KN<long> *const &connectivity_offset,
                           KN<long> *const &connectivity) {
  global_connectivity_offset.resize(connectivity_offset->size());
  for (int i = 0; i < connectivity_offset->size(); ++i)
    global_connectivity_offset[i] = (*connectivity_offset)[i];
  global_connectivity.resize(connectivity->size());
  for (int j = 0; j < connectivity->size(); ++j)
    global_connectivity[j] = (*connectivity)[j];
  cwipi_ho_define_mesh(coupling_id->c_str(), n_vertex, n_element, order,
                       &((*coordinates)[0]), global_connectivity_offset.data(), global_connectivity.data());
  return 0L;
}

long CwipiDeleteCouplingFfpp(string *const &coupling_id) {
  cwipi_delete_coupling(coupling_id->c_str());
  return 0L;
}

long CwipiIsSendFfpp(string *const &coupling_name, string *const &exchange_name,
                     long const &tag, long const &stride, long const &time_step,
                     double const &time_value,
                     string *const &sending_field_name,
                     KN<double> *const &sending_field, long *const &request) {
  int c_request;
  cwipi_issend(coupling_name->c_str(), exchange_name->c_str(), tag, stride,
               time_step, time_value, sending_field_name->c_str(),
               &((*sending_field)[0]), &c_request);
  *request = c_request;
  return 0L;
}
long CwipiIrecvFfpp(string *const &coupling_name, string *const &exchange_name,
                    long const &tag, long const &stride, long const &time_step,
                    double const &time_value,
                    string *const &receiving_field_name,
                    KN<double> *const &receiving_field, long *const &request) {
  int r_request;
  cwipi_irecv(coupling_name->c_str(), exchange_name->c_str(), tag, stride,
              time_step, time_value, receiving_field_name->c_str(),
              &((*receiving_field)[0]), &r_request);
  *request = r_request;
  return 0L;
}
long CwipiCreateCouplingFfpp(string *const &coupling_name,
                             string *const &coupled_application,
                             long const &entitiesDim, double const &tolerance,
                             long const &output_frequency,
                             string *const &output_format,
                             string *const &output_format_option) {
  cwipi_create_coupling(
      coupling_name->c_str(), CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
      coupled_application->c_str(), entitiesDim, tolerance, CWIPI_STATIC_MESH,
      CWIPI_SOLVER_CELL_VERTEX, output_frequency, output_format->c_str(),
      output_format_option->c_str());
  return 0L;
}

long CwipiInitFfpp(pcommworld const &cwipi_comm,
                   string *const &application_name,
                   pcommworld const &application_comm) {
  cwipi_init(*(MPI_Comm*)cwipi_comm, application_name->c_str(), (MPI_Comm*)application_comm);
  return 0L;
}

long GetUniqueConnectivity(KN<long> *const &arr, KN<long> *const &conn_index,
                           KN<long> *const &conn) {
  long ne = conn_index->size() - 1;
  set<long> s;
  for (int i = 0; i < ne; ++i) {
    for (int j = (*conn_index)[i]; j < (*conn_index)[i + 1]; ++j) {
      s.insert((*conn)[j]);
    }
  }
  long num_unique = s.size();
  std::cout << "Found totally " << num_unique << " vertices\n";
  arr->resize(num_unique);
  long counter = 0;
  for (const auto &item : s ) { 
     (*arr)[counter] = item;
     counter++;
  }

  // Create the hash table of GID to LID
  std::map<long, long> gid_to_lid;
  counter = 0;
  for (int i = 0; i < num_unique; ++i)
    gid_to_lid[(*arr)[i]] = ++counter;
  for (int i = 0; i < ne; ++i) {
    for (int j = (*conn_index)[i]; j < (*conn_index)[i + 1]; ++j) {
      (*conn)[j] = gid_to_lid[(*conn)[j]];
    }
  }
#if 0 // For printing purpose
        for( int i = 0; i < ne; ++i ) {
          for( int j = (*conn_index)[i]; j < (*conn_index)[i+1]; ++j ) {
            std::cout << (*conn)[j] << " ";
          }
          std::cout << "\n";
        }
#endif
  return 0L;
}

long CwipiDumpApplicationPropertiesFfpp () {
     cwipi_dump_application_properties();
     return 0L;
     }

/*long CwipiExchange(string *const &coupling_id,string *const &exchange_name,long const &stride,long const &time_step,
     long const &time_value,string *const &sending_field_name,double *const &sending_field,string *const &receiving_field_name,
     double *const &receiving_field,long *const &n_not_located_points) {
    
     cwipi_exchange(coupling_id->c_str(),exchange_name->c_str(),stride,time_step,
     time_value,sending_field_name->c_str(),sending_field,reinterpret_cast<char*>(receiving_field_name),receiving_field,
     reinterpret_cast<int *>(n_not_located_points));
     return 0L;
}*/

void finalize( ) {
  global_connectivity.clear();
  global_connectivity.shrink_to_fit();
  global_connectivity_offset.clear();
  global_connectivity_offset.shrink_to_fit();
  cwipi_finalize();
#if defined(PETSC_HAVE_ELEMENTAL)
  PetscElementalFinalizePackage();
#endif
}

void init() {
  const std::string str_cwipi_prefix("-cwipi-prefix");
  std::stringstream cat;
  int argc = pkarg->n;
  char** argv = new char*[argc];
  for (int i = 0; i < argc; ++i) argv[i] = const_cast< char* >((*(*pkarg)[i].getap( ))->c_str( ));
  cat << "FreeFem++";
  for(int i=0; i<argc-1; ++i) {
    if( str_cwipi_prefix.compare(argv[i]) == 0 ) {
      cat << "-" << argv[i+1];
      break;
    }
  }
  delete[] argv;
  cwipi_init(MPI_COMM_WORLD, cat.str().c_str(), &cwipi_comm);
#if defined(WITH_petsc)
#if defined(PETSC_HAVE_ELEMENTAL)
  PetscElementalInitializePackage();
#endif
  PETSC_COMM_WORLD = cwipi_comm;
#endif
  ff_atend(finalize);
  Global.Add("cwipiComm", "(",
             new OneOperator1_<long, pcommworld>(cwipiComm));
  Global.Add("CwipiWaitIsSendFfpp", "(",
             new OneOperator2_<long, string *, long>(CwipiWaitIsSendFfpp));
  Global.Add("CwipiWaitIrecvFfpp", "(",
             new OneOperator2_<long, string *, long>(CwipiWaitIrecvFfpp));
  Global.Add("CwipiLocateFfpp", "(",
             new OneOperator1_<long, string *>(CwipiLocateFfpp));
  Global.Add("CwipiGetNLocatedPointsFfpp", "(",
             new OneOperator1_<long, string *>(CwipiGetNLocatedPointsFfpp));
  Global.Add("CwipiGetNNotLocatedPointsFfpp", "(",
             new OneOperator1_<long, string *>(CwipiGetNNotLocatedPointsFfpp));
  Global.Add("CwipiDefineMeshFfpp", "(",
             new OneOperator6_<long, string *, long, long, KN<double> *,
                               KN<long> *, KN<long> *>(CwipiDefineMeshFfpp));
  Global.Add("CwipiHoDefineMeshFfpp", "(",
             new OneOperator7_<long, string *, long, long, long, KN<double> *,
                               KN<long> *, KN<long> *>(CwipiHoDefineMeshFfpp));
  Global.Add("CwipiDeleteCouplingFfpp", "(",
             new OneOperator1_<long, string *>(CwipiDeleteCouplingFfpp));
  
  Global.Add(
      "CwipiIsSendFfpp", "(",
      new OneOperator9_<long, string *, string *, long, long, long, double,
                        string *, KN<double> *, long *>(CwipiIsSendFfpp));
  Global.Add(
      "CwipiIrecvFfpp", "(",
      new OneOperator9_<long, string *, string *, long, long, long, double,
                        string *, KN<double> *, long *>(CwipiIrecvFfpp));
                        
  Global.Add("CwipiCreateCouplingFfpp", "(",
             new OneOperator7_<long, string *, string *, long, double, long,
                               string *, string *>(CwipiCreateCouplingFfpp));
  Global.Add(
      "CwipiInitFfpp", "(",
      new OneOperator3_<long, pcommworld, string *, pcommworld>(
          CwipiInitFfpp));
  Global.Add("GetUniqueConnectivity", "(",
             new OneOperator3_<long, KN<long> *, KN<long> *, KN<long> *>(
                 GetUniqueConnectivity));
  Global.Add("CwipiDumpApplicationPropertiesFfpp", "(",
             new OneOperator0<long>(CwipiDumpApplicationPropertiesFfpp));
  //Global.Add("CwipiExchange", "(", new OneOperator10_<long, string *, string*,
  //long,long,long,string *, KN<long> *, string * ,KN<long> *,long>(CwipiExchange));
}
LOADFUNC(init);
