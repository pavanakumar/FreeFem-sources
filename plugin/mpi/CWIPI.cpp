#include "cwipi.h"
#include "AFunction.hpp"
#include "AFunction_ext.hpp"
#include "algorithm"
#include <bits/stdc++.h>
#include <ff++.hpp>
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include "../../src/mpi/fmpi.hpp"

// cwipi_solver_type_t cellvertex =CWIPI_SOLVER_CELL_VERTEX;
// int solvert = static_cast<int>(cellvertex);

extern MPI_Comm ff_global_comm_world;
extern KN<String>* pkarg;
extern long mpirank;
extern long mpisize;

using namespace Fem2D;
using namespace std;

#ifdef ENABLE_CWIPI
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
                         KN<long> *const &connectivity_index,
                         KN<long> *const &connectivity) {
  std::vector<int> temp(connectivity_index->size());
  for (int i = 0; i < connectivity_index->size(); ++i)
    temp[i] = (*connectivity_index)[i];
  //  std::copy( *(*connectivity_index), *(*connectivity_index) +
  //  connectivity_index->size(), temp.begin());
  std::vector<int> temp1(connectivity->size());
  for (int j = 0; j < connectivity->size(); ++j)
    temp1[j] = (*connectivity)[j];
  //  std::vector<double> temp2(coordinates->size());
  //  for( int k=0; k < coordinates->size(); ++k)
  //    temp2[k] = (*coordinates)[k];
  //  std::copy( *(*connectivity), *(*connectivity) + connectivity->size(),
  //  temp1.begin());
  cwipi_define_mesh(coupling_id->c_str(), n_vertex, n_element,
                    &((*coordinates)[0]), temp.data(), temp1.data());
  return 0L;
}
long CwipiHoDefineMeshFfpp(string *const &coupling_id, long const &n_vertex,
                           long const &n_element, long const &order,
                           KN<double> *const &coordinates,
                           KN<long> *const &connectivity_index,
                           KN<long> *const &connectivity) {
  std::vector<int> temp2(connectivity_index->size());
  for (int i = 0; i < connectivity_index->size(); ++i)
    temp2[i] = (*connectivity_index)[i];
  //  std::copy( *(*connectivity_index), *(*connectivity_index) +
  //  connectivity_index->size(), temp.begin());
  std::vector<int> temp3(connectivity->size());
  for (int j = 0; j < connectivity->size(); ++j)
    temp3[j] = (*connectivity)[j];
  //  std::vector<double> temp2(coordinates->size());
  //  for( int k=0; k < coordinates->size(); ++k)
  //    temp2[k] = (*coordinates)[k];
  //  std::copy( *(*connectivity), *(*connectivity) + connectivity->size(),
  //  temp1.begin());
  cwipi_ho_define_mesh(coupling_id->c_str(), n_vertex, n_element, order,
                       &((*coordinates)[0]), temp2.data(), temp3.data());
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
                     double *const &sending_field, long *const &request) {
  // int tag = dtag;
  // int stride = dstride;
  // int time_step = dtime_step;
  cwipi_issend(coupling_name->c_str(), exchange_name->c_str(), tag, stride,
               time_step, time_value, sending_field_name->c_str(),
               sending_field, reinterpret_cast<int *>(request));
  return 0L;
}
long CwipiIrecvFfpp(string *const &coupling_name, string *const &exchange_name,
                    long const &tag, long const &stride, long const &time_step,
                    double const &time_value,
                    string *const &receiving_field_name,
                    double *const &receiving_field, long *const &request) {
  cwipi_irecv(coupling_name->c_str(), exchange_name->c_str(), tag, stride,
              time_step, time_value, receiving_field_name->c_str(),
              receiving_field, reinterpret_cast<int *>(request));
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
// long CwipiInitFfpp(string *const &application_name){ // MPI_Comm *const
// &application_comm){ MPI_Comm const &cwipi_comm, string *const
// &application_name, MPI_Comm *const &application_comm)
//   MPI_Comm temp;
//   MPI_Comm_dup(MPI_COMM_WORLD, &temp);
//   cwipi_init( MPI_COMM_WORLD, (const char *) application_name, (&temp) );
//   return 0L;
//}

void CwipiFinalizeFfpp() {
  cwipi_finalize();
}

void CwipiInitFfpp() {
  int argc = pkarg->n;
  char** argv = new char*[argc];
  for (int i = 0; i < argc; ++i) argv[i] = const_cast< char* >((*(*pkarg)[i].getap( ))->c_str( ));  
  const std::string str_cwipi_prefix("--cwipi-prefix");
  std::stringstream cat;
  cat << "FreeFem++";
  for(int i=0; i<argc; ++i) {
      if( str_cwipi_prefix.compare(argv[i]) == 0 ) {
          cat << "-" << argv[i+1];
          break;
      }
    }
  std::cout << "Application name = " << cat.str().c_str() << "\n";
  cwipi_init(MPI_COMM_WORLD, cat.str().c_str(), &ff_global_comm_world);
  int rank, size;
  MPI_Comm_rank(ff_global_comm_world, &rank);
  MPI_Comm_size(ff_global_comm_world, &size);
  mpisize = size;
  mpirank = rank;
  ff_atend(CwipiFinalizeFfpp);
}

bool LinkTest(KN<double> *const &testInput1, KN<long> *const &testInput2) {
  cout << "This is a shared library test..." << endl;
  cout << *testInput1 << endl;
  cout << *testInput2 << endl;
  // cwipi_init_ffpp();
  return 0;
}

long GetUniqueConnectivity(KN<long> *const &arr, KN<long> *const &conn_index,
                           KN<long> *const &conn) {
  long ne = conn_index->size() - 1;
  set<long> s;
  for (int i = 0; i < ne; ++i) {
    // std::cout << "offset " << (*conn_index)[i] << " to " <<
    // (*conn_index)[i+1] << "\n";
    for (int j = (*conn_index)[i]; j < (*conn_index)[i + 1]; ++j) {
      s.insert((*conn)[j]);
    }
  }
  long num_unique = s.size();
  std::cout << "Found totally " << num_unique << " vertices\n";
  long *temp = (long *)malloc(sizeof(long) * num_unique);
  std::copy(s.begin(), s.end(), temp);
  arr->unset();
  arr->set(temp, num_unique);
  // Create the hash table of GID to LID
  std::map<long, long> gid_to_lid;
  long counter = 0;
  for (int i = 0; i < num_unique; ++i)
    gid_to_lid[temp[i]] = ++counter;
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

static void init() {
  // Init the CWIPI plugin
  CwipiInitFfpp();
  static fMPI_Comm fcwipiWorld=ff_global_comm_world;
  //Global.New("mpiCommWorld",CConstant<fMPI_Comm*>(&fcwipiWorld));
  Global.New("cwipiCommWorld",CConstant<fMPI_Comm*>(&fcwipiWorld));
  Global.Add("LinkTest", "(",
             new OneOperator2_<bool, KN<double> *, KN<long> *>(LinkTest));
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
                        string *, double *, long *>(CwipiIsSendFfpp));
  Global.Add(
      "CwipiIrecvFfpp", "(",
      new OneOperator9_<long, string *, string *, long, long, long, double,
                        string *, double *, long *>(CwipiIrecvFfpp));
  Global.Add("CwipiCreateCouplingFfpp", "(",
             new OneOperator7_<long, string *, string *, long, double, long,
                               string *, string *>(CwipiCreateCouplingFfpp));
  Global.Add("GetUniqueConnectivity", "(",
             new OneOperator3_<long, KN<long> *, KN<long> *, KN<long> *>(
                 GetUniqueConnectivity));
}
LOADFUNC(init);
#endif
