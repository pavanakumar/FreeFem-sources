// ORIG-DATE:     Fev 2010
// -*- Mode : c++ -*-
//
// SUMMARY  : liaison medit freefem++ : adaptmesh in 3d 
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
// AUTHOR   : Jacques Morice
// E-MAIL   : jacques.morice@ann.jussieu.fr
//
//   for automatic  compilation with ff-c++
//ff-c++-LIBRARY-dep:   mmg3d
//ff-c++-cpp-dep: 
//  

/* 
 This file is part of Freefem++
 
 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.
 
 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

 Thank to the ARN ()  FF2A3 grant
 ref:ANR-07-CIS7-002-01 
 */


//  ./ff-c++ mmg3d.cpp -I../src/libMesh/ -I../download/include/mmg3d/ -L../download/lib/mmg3d/ -lmmg3d -L../src/libMesh/ -lMesh

#include "ff++.hpp" 
#include "msh3.hpp"
//#define ADAPTLIBRARY
#include "libmmg3d.h"

using namespace  Fem2D;
using namespace  mmg3d;

Mesh3 * MMG_pMesh_to_msh3(MMG_pMesh meshMMG){
  int i;
   
  cout << "transformation maillage --> msh3 " << endl;
  cout << "meshMMG->np =" << meshMMG->np << endl;
  cout << "meshMMG->ne =" << meshMMG->ne << endl;
  cout << "meshMMG->nt =" << meshMMG->nt << endl;
  Vertex3 *v = new Vertex3[meshMMG->np];
  Tet *t  = new Tet[meshMMG->ne];
  Tet *tt = t;
  Triangle3 *b  = new Triangle3[meshMMG->nt];
  Triangle3 *bb = b;
     
  int k;
  MMG_pPoint ppt;
  for (k=1; k<=meshMMG->np; k++) {
    ppt = &meshMMG->point[k];
    v[k-1].x   = ppt->c[0] ;
    v[k-1].y   = ppt->c[1] ;
    v[k-1].z   = ppt->c[2] ;
    v[k-1].lab = ppt->ref  ;
  }

  MMG_pTetra ptetra;
  for (k=1; k<=meshMMG->ne; k++) {
    int iv[4],lab;
    ptetra = &meshMMG->tetra[k];
    iv[0]=ptetra->v[0]-1;
    iv[1]=ptetra->v[1]-1;
    iv[2]=ptetra->v[2]-1;
    iv[3]=ptetra->v[3]-1;
    lab  =ptetra->ref;
    (*tt++).set( v, iv, lab);   
  }
 
  MMG_pTria ptriangle;
  for (k=1; k<=meshMMG->nt; k++) {
    int iv[3],lab;
    ptriangle = &meshMMG->tria[k];
    iv[0]=ptriangle->v[0]-1;
    iv[1]=ptriangle->v[1]-1;
    iv[2]=ptriangle->v[2]-1;
    lab=ptriangle->ref;
    (*bb++).set( v, iv, lab);
  }


 Mesh3 *T_TH3 = new Mesh3(meshMMG->np, meshMMG->ne, meshMMG->nt, v, t, b);
 cout << "transformation maillage --> msh3 " << endl;
 cout << "meshMMG->np =" << meshMMG->np << endl;
 cout << "meshMMG->ne =" << meshMMG->ne << endl;
 cout << "meshMMG->nt =" << meshMMG->nt << endl;
 cout << "T_TH3" << T_TH3->nv  << " " << T_TH3->nt  << " " << T_TH3->nbe <<endl;
 return T_TH3;
}

MMG_pMesh mesh3_to_MMG_pMesh(const Mesh3 &Th3, const int & nvmax, const int &ntrimax, const int & ntetmax , const bool & boolMoving , KN<double> &Moving){
  MMG_pMesh meshMMG;
  meshMMG = (MMG_pMesh)calloc(1,sizeof(MMG_Mesh)) ;

  meshMMG->np = Th3.nv;
  meshMMG->nt = Th3.nbe;
  meshMMG->ne = Th3.nt;

  meshMMG->npmax = nvmax;
  meshMMG->ntmax = ntrimax;
  meshMMG->nemax = ntetmax;

  meshMMG->point = (MMG_pPoint)calloc(meshMMG->npmax+1,sizeof(MMG_Point));
  meshMMG->tetra = (MMG_pTetra)calloc(meshMMG->nemax+1,sizeof(MMG_Tetra));
  meshMMG->tria = (MMG_pTria) calloc(meshMMG->ntmax+1,sizeof(MMG_Tria));
  //meshMMG->disp = NULL;
 
  if( boolMoving ){
    MMG_pDispl ppd;
    meshMMG->disp = (MMG_pDispl)calloc(1,sizeof(MMG_Displ));
    ppd = meshMMG->disp;
    ppd->np = meshMMG->np;
    ppd->mv = (double *) calloc(3*(meshMMG->np+1),sizeof(double));
    assert(mesh->disp->mv);
    ppd->alpha = (short *) calloc( meshMMG->np+1,sizeof(short));
    assert(mesh->disp->alpha);
    for(int ii=0; ii < meshMMG->np; ii++){      
      ppd->mv[3*ii+1] = Moving[3*ii];
      ppd->mv[3*ii+2] = Moving[3*ii+1];
      ppd->mv[3*ii+3] = Moving[3*ii+2];
    }
  }
  else{
    meshMMG->disp = NULL;
  }
  meshMMG->adja = (int*)calloc(4*meshMMG->nemax+5,sizeof(int));
  
  int k;
  MMG_pPoint ppt;
  for (k=1; k<=meshMMG->np; k++) {
    ppt = &meshMMG->point[k];
    ppt->c[0] = Th3.vertices[k-1].x;
    ppt->c[1] = Th3.vertices[k-1].y;
    ppt->c[2] = Th3.vertices[k-1].z;
    ppt->ref  = Th3.vertices[k-1].lab;
  }

  
  MMG_pTetra ptetra;
  for (k=1; k<=meshMMG->ne; k++) {
    const Tet & K(Th3.elements[k-1]);
    ptetra = &meshMMG->tetra[k];
    ptetra->v[0] = Th3.operator()(K[0])+1;
    ptetra->v[1] = Th3.operator()(K[1])+1;
    ptetra->v[2] = Th3.operator()(K[2])+1;
    ptetra->v[3] = Th3.operator()(K[3])+1;
    ptetra->ref = K.lab;
  }

  MMG_pTria ptriangle;
  for (k=1; k<=meshMMG->nt; k++) {
    const Triangle3 & K(Th3.be(k-1));
    ptriangle = &meshMMG->tria[k];
    ptriangle->v[0] = Th3.operator()(K[0])+1;
    ptriangle->v[1] = Th3.operator()(K[1])+1;
    ptriangle->v[2] = Th3.operator()(K[2])+1;
    ptriangle->ref = K.lab;
  }
  
  /*  // pour le deplacement des corps rigides
  MMG_pDispl pd;
  for (k=1; k<mesh->np; k++) {
    pd = &mesh->disp[k];
    pd->mv[0] = depx;
    pd->mv[1] = depy;
    pd->mv[2] = depz;
  }
  */
  return meshMMG;
}

MMG_pSol metric_mmg3d(const int & nv, const int & nvmax, const KN<double> &metric){
  static const int wrapperMetric[6]={0,1,2,3,4,5};
  MMG_pSol sol;
  
  sol= (MMG_pSol)calloc(1,sizeof(MMG_Sol)) ;

  cout << metric.N() << " taille de la metrique "<< endl;
  if( metric.N() > 0){
    sol->np = nv;
    sol->npmax=nvmax;
    
    if(metric.N() == nv){ 
      const int ic=1;
      char newvalue[sizeof(int)];
      sprintf(newvalue,"%s", (char*)&ic);
      sol->offset = *newvalue;
    }
    else{
      const int ic=6;
      char newvalue[sizeof(int)];
      sprintf(newvalue,"%s", (char*)&ic);
      sol->offset = *newvalue;
    }
    
    sol->met = (double*)calloc(sol->npmax+1,sol->offset*sizeof(double));
    int k,isol,i;
    double tmp;
    for (k=1; k<=sol->np; k++) {
      isol = (k-1)*sol->offset + 1;
      for (i=0; i< sol->offset; i++){
	sol->met[isol + i] = metric[(isol-1)+i];
      }
      // MMG_swap data
      tmp                = sol->met[isol + 2];
      sol->met[isol + 2] = sol->met[isol + 3];
      sol->met[isol + 3] = tmp;
    }
  }
  else{
    sol->np=0;
    sol->offset=1;
  }
  return sol;
}

void metric_mmg3d_to_ff_metric(MMG_pSol sol, KN<double> &metric){
  static const int wrapperMetric[6]={0,1,2,3,4,5};

  int k,isol,i;
  for (k=1; k<=sol->np; k++) {
    isol = (k-1)*sol->offset + 1;
    for (i=0; i< sol->offset; i++)
      metric[(isol-1)+i]= sol->met[isol + i];
  }
}



class mmg3d_Op: public E_F0mps 
{
public:
  Expression eTh,xx,yy,zz;
  static const int n_name_param = 5; // 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  KN_<double>  arg(int i,Stack stack,KN_<double> a ) const
  { return nargs[i] ? GetAny<KN_<double> >( (*nargs[i])(stack) ): a;}
  double  arg(int i,Stack stack,double a ) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  int  arg(int i,Stack stack, int a ) const{ return nargs[i] ? GetAny< int >( (*nargs[i])(stack) ): a;}
  
  
public:
  mmg3d_Op(const basicAC_F0 &  args ,Expression tth) 
    : eTh(tth)
  {
    //if(verbosity >1) 
    cout << "mmg3d"<< endl;
    args.SetNameParam(n_name_param,name_param,nargs);
    
    const E_Array * a1=0 ;
    if(nargs[2])  a1  = dynamic_cast<const E_Array *>(nargs[2]);
  
    if(a1) {
      if(a1->size() !=3) 
	CompileError("mmg3d(Th,displacement=[X,Y,Z],) ");
      xx=to<double>( (*a1)[0]); 
      yy=to<double>( (*a1)[1]);
      zz=to<double>( (*a1)[2]);
    }    
  } 
  
  AnyType operator()(Stack stack)  const ;
};


basicAC_F0::name_and_type  mmg3d_Op::name_param[]= {
  /*
  {  "nvmax", &typeid(long)},     // 0
  {  "ntrimax", &typeid(long)},   // 1
  {  "ntetmax", &typeid(long)},   // 2
  */
  {  "options", &typeid(KN_<long>)},    // 3 -> 0
  {  "metric", &typeid(KN_<double>)},   // 4 -> 1
  {  "displacement", &typeid(E_Array)},  // 5 -> 2
  {  "displVect", &typeid(KN_<double>)},   // 6 ->3 
  {  "memory", &typeid(long)} // 7->4
};

class mmg3d_ff : public OneOperator { public:  
     mmg3d_ff() : OneOperator( atype<pmesh3>(), atype<pmesh3>() ) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
	return  new  mmg3d_Op( args,t[0]->CastTo(args[0]) ); 
  }
};

AnyType mmg3d_Op::operator()(Stack stack)  const 
{
  // initialisation
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh3 * pTh= GetAny<Mesh3 *>((*eTh)(stack));
  ffassert( pTh );
  Mesh3 &Th3=*pTh;
  int nv=Th3.nv;
  int nt=Th3.nt;
  int nbe=Th3.nbe;

  // default value for max nv,ntri,ntet  max
  // division of tetrahedrons in the middle of edges 
  int defaultnvmax=    500000;       
  int defaultntrimax= 1000000;   
  int defaultntetmax= 3000000;     
  
  KN<long> defaultopt(6);
  defaultopt(0)= 1;
  defaultopt(1)= 0;
  defaultopt(2)= 64;
  defaultopt(3)= 0;
  defaultopt(4)= 0;
  defaultopt(5)= 3;

  int  nvmax;   //(arg(0,stack,-1));  
  int  ntrimax; //(arg(1,stack,-1));
  int  ntetmax; //(arg(2,stack,-1));
  KN<int> opt(arg(0,stack,defaultopt));

  int memory(arg(4,stack,-1));
  
  if( memory < 0 ){
    nvmax   = max( (int) 1.5*nv,defaultnvmax);
    ntrimax = max( (int) 1.5*nbe,defaultntrimax);
    ntetmax = max( (int) 1.5*nt,defaultntetmax);
  }     
  else{
    int million = 1048576L;
    int bytes = sizeof(MMG_Point)   + 0.2*sizeof(MMG_Tria)	
      + 6*sizeof(MMG_Tetra) + 4*sizeof(int) 
      + sizeof(MMG_Sol) + sizeof(MMG_Displ) 
      + sizeof(int) + 5*sizeof(int);
    int npask = (double)memory / bytes * million;
    nvmax   = max( (int) 1.5*nv,npask);
    ntetmax = max( (int) 1.5*nt,6*npask);
    ntrimax = max( (int) 1.5*nbe,(int)(0.3*npask));

    cout << " nvmax " << nvmax << endl;
    cout << " ntrimax " << ntrimax << endl;
    cout << " ntetmax " << ntetmax << endl;
  }										   
  
  KN<double> metric;
  if(metric.N() != 0) exit(1);
  // definiton d'une metric par default
  if( nargs[1]  ){ 
    metric = GetAny<KN_<double> >( (*nargs[1])(stack) ); 
    assert(metric.N()==Th3.nv || metric.N()==6*Th3.nv);
  }

 
  bool BoolMoving=0;
  KN<double> Moving(0);
  
  if( nargs[2] || nargs[3] ){
    BoolMoving=1;
    if( nargs[3] ){
      Moving = GetAny<double>( (*nargs[3])(stack) );
      assert( Moving.N() == 3*Th3.nv );
      if( Moving.N() != 3*Th3.nv ){ cerr << " Displacement vector is of size 3*Th.nv" << endl; exit(1);} 
    }
    else{ 
      MeshPoint *mp3(MeshPointStack(stack));
      Moving.resize(3*Th3.nv);
      for( int i=0; i<Th3.nv; ++i){
	mp3->set( Th3.vertices[i].x, Th3.vertices[i].y, Th3.vertices[i].z );
	if(xx) Moving[3*i]   = GetAny<double>((*xx)(stack)); 
	if(yy) Moving[3*i+1] = GetAny<double>((*yy)(stack));
	if(zz) Moving[3*i+2] = GetAny<double>((*zz)(stack));  
      }
    }
    //if(verbosity > 2) 
    if(verbosity >2) cout << "displacement vector is loading" << endl;
  }
  

  MMG_pMesh MMG_Th3=mesh3_to_MMG_pMesh(Th3, nvmax, ntrimax, ntetmax,BoolMoving, Moving);

  /*
    
    MMG_pMesh MMG_Th3=mesh3_to_MMG_pMesh(Th3, nvmax, ntrimax, ntetmax);
  MMG_pMesh meshMMG;
  meshMMG = (MMG_pMesh)calloc(1,sizeof(MMG_Mesh)) ;
  
  meshMMG->np = Th3.nv;
  meshMMG->nt = Th3.nbe;
  meshMMG->ne = Th3.nt;

  meshMMG->npmax = nvmax;
  meshMMG->ntmax = ntrimax;
  meshMMG->nemax = ntetmax;

  meshMMG->point = (MMG_pPoint)calloc(meshMMG->npmax+1,sizeof(MMG_Point));
  meshMMG->tetra = (MMG_pTetra)calloc(meshMMG->nemax+1,sizeof(MMG_Tetra));
  meshMMG->tria = (MMG_pTria) calloc(meshMMG->ntmax+1,sizeof(MMG_Tria));
  //meshMMG->disp = (MMG_pDispl)calloc(meshMMG->npmax+1,sizeof(MMG_Displ));
  meshMMG->disp = NULL;
  meshMMG->adja = (int*)calloc(4*meshMMG->nemax+5,sizeof(int));
  
  int k;
  MMG_pPoint ppt;
  for (k=1; k<=meshMMG->np; k++) {
    ppt = &meshMMG->point[k];
    ppt->c[0] = Th3.vertices[k-1].x;
    ppt->c[1] = Th3.vertices[k-1].y;
    ppt->c[2] = Th3.vertices[k-1].z;
    ppt->ref  = Th3.vertices[k-1].lab;
  }

  
  MMG_pTetra ptetra;
  for (k=1; k<=meshMMG->ne; k++) {
    const Tet & K(Th3.elements[k-1]);
    ptetra = &meshMMG->tetra[k];
    ptetra->v[0] = Th3.operator()(K[0])+1;
    ptetra->v[1] = Th3.operator()(K[1])+1;
    ptetra->v[2] = Th3.operator()(K[2])+1;
    ptetra->v[3] = Th3.operator()(K[3])+1;
    ptetra->ref = K.lab;
  }

  MMG_pTria ptriangle;
  for (k=1; k<=meshMMG->nt; k++) {
    const Triangle3 & K(Th3.be(k-1));
    ptriangle = &meshMMG->tria[k];
    ptriangle->v[0] = Th3.operator()(K[0])+1;
    ptriangle->v[1] = Th3.operator()(K[1])+1;
    ptriangle->v[2] = Th3.operator()(K[2])+1;
    ptriangle->ref = K.lab;
  }
  
 
  MMG_Mesh *MMG_Th3= meshMMG;
  */

  MMG_pSol sol=metric_mmg3d(nv, nvmax, metric);
  
  int res=MMG_mmg3dlib( opt, MMG_Th3, sol);

  if( res > 0){
    cout << " problem of remeshing with mmg3d :: error" <<  res << endl; 
    exit(1);
  }
  
  Mesh3 *Th3_T = MMG_pMesh_to_msh3( MMG_Th3 );
  if(verbosity > 10) cout << "buildGtree" << endl;
  Th3_T->BuildGTree();

  /* free mem */
  if(verbosity > 10) cout << "mesh" << endl;
  free( MMG_Th3->point );
  free( MMG_Th3->tria  );
  free( MMG_Th3->tetra );
  /*la desallocation de ce pointeur plante dans certains cas...*/
  if(verbosity > 10) cout << "mesh: adja" << endl;
  free( MMG_Th3->adja);
  if(verbosity > 10) cout << "mesh: disp" << endl;
  if( BoolMoving ){
    free( MMG_Th3->disp->alpha );
    free( MMG_Th3->disp->mv );
  }
  free( MMG_Th3->disp );
  free( MMG_Th3 );

  if(verbosity > 10) cout << "sol" << endl;
  free(sol->met);
  free(sol);
  
  *mp=mps;
  Add2StackOfPtr2FreeRC(stack,Th3_T);
  return Th3_T;
}


class Init1 { public:
  Init1();
};

static Init1 init1;  //  une variable globale qui serat construite  au chargement dynamique 

Init1::Init1(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 
  
  //if (verbosity)
  if(verbosity) cout << " load: mmg3d  " << endl;
  
  Global.Add("mmg3d","(",new mmg3d_ff);
  
}


#define  WITH_NO_INIT
#include "msh3.hpp" 
