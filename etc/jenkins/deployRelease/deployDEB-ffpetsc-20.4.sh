#!/usr/bin/env bash

set -x
set -u
set -e

## Parameters
TOKEN=$1
ORGANIZATION="FreeFem"
REPOSITORY="FreeFem-sources"
VERSION=`grep AC_INIT configure.ac | cut -d"," -f2`
RELEASE_TAG_NAME="v$VERSION"
distrib=`uname -s`-`uname -r`
OSRELEASE=`lsb_release -r|awk '{print $2}'`
DISTRIB="Ubuntu"

DEB_NAME="freefem-${VERSION}-amd64-${OSRELEASE}"
GH_DEB_NAME="FreeFEM-${VERSION}-amd64-${OSRELEASE}.deb"
if $1 = "no" ; then 
## DEB build
autoreconf -i
./configure --enable-download --enable-optim --enable-generic
./3rdparty/getall -a -o PETSc,Ipopt,NLopt,freeYams,FFTW,Gmm++,MMG3D,mshmet,MUMPS,htool
## compile and install f, libhdf5-dev (>=1.10.4)f-petsc
cd 3rdparty/ff-petsc && make petsc-slepc && cd -
./reconfigure

make -j4
make -j4 install
fi
#create FreeFEM Debian package  with ff-petsc

mkdir $DEB_NAME
mkdir $DEB_NAME/DEBIAN
touch $DEB_NAME/DEBIAN/control
echo "Package: freefem" >> $DEB_NAME/DEBIAN/control
echo "Version: "$VERSION >> $DEB_NAME/DEBIAN/control
echo "Section: custom" >> $DEB_NAME/DEBIAN/control
echo "Architecture: amd64" >> $DEB_NAME/DEBIAN/control
echo "Depends: libc6 (>= 2.31), g++ (>= 9.3), gcc (>= 9.3), gfortran (>= 9.3), libgsl-dev (>=2.5), libhdf5-dev (>=1.10.4), liblapack-dev (>= 3.9), libopenmpi-dev (>=4.0.3) ,freeglut3-dev (>= 2.8.1) ">> $DEB_NAME/DEBIAN/control
echo "Maintainer: FreeFEM, Frédéric Hecht <frederic.hecht@sorbonne-universite.fr> " >> $DEB_NAME/DEBIAN/control
echo "Description: FreeFEM, Finite Element Language software" >> $DEB_NAME/DEBIAN/control
echo "Homepage: https://freefem.org" >> $DEB_NAME/DEBIAN/control
mkdir -p $DEB_NAME/usr/local
mkdir -p $DEB_NAME/usr/local/share/FreeFEM
mkdir -p $DEB_NAME/usr/local/bin
mkdir -p $DEB_NAME/usr/local/lib/ff++
mkdir -p $DEB_NAME/usr/share/doc/freefem

cp -r /usr/local/ff-petsc/ $DEB_NAME/usr/local/ff-petsc
cp -r /usr/local/lib/ff++/$VERSION $DEB_NAME/usr/local/lib/ff++/$VERSION
cp -r /usr/local/share/FreeFEM/$VERSION $DEB_NAME/usr/local/share/FreeFEM/$VERSION
cp -r /usr/local/bin/FreeFem++ $DEB_NAME/usr/local/bin/FreeFem++
cp -r /usr/local/bin/FreeFem++-mpi $DEB_NAME/usr/local/bin/FreeFem++-mpi
cp -r /usr/local/bin/FreeFem++-nw $DEB_NAME/usr/local/bin/FreeFem++-nw
cp -r /usr/local/bin/bamg $DEB_NAME/usr/local/bin/bamg
cp -r /usr/local/bin/cvmsh2 $DEB_NAME/usr/local/bin/cvmsh2
cp -r /usr/local/bin/ff-c++ $DEB_NAME/usr/local/bin/ff-c++
cp -r /usr/local/bin/ff-get-dep $DEB_NAME/usr/local/bin/ff-get-dep
cp -r /usr/local/bin/ff-mpirun $DEB_NAME/usr/local/bin/ff-mpirun
cp -r /usr/local/bin/ff-pkg-download $DEB_NAME/usr/local/bin/ff-pkg-download
cp -r /usr/local/bin/ffglut $DEB_NAME/usr/local/bin/ffglut
cp -r /usr/local/bin/ffmaster $DEB_NAME/usr/local/bin/ffmaster
cp -r /usr/local/bin/ffmedit $DEB_NAME/usr/local/bin/ffmedit
cp  AUTHORS $DEB_NAME/usr/share/doc/freefem/AUTHOR
cp  README.md $DEB_NAME/usr/share/doc/freefem/README.md

dpkg-deb --build $DEB_NAME/
mv $DEB_NAME.deb $GH_DEB_NAME
exit 0
## Deploy in GitHub release
RELEASE=`curl 'https://api.github.com/repos/'$ORGANIZATION'/'$REPOSITORY'/releases/tags/'$RELEASE_TAG_NAME`
UPLOAD_URL=`printf "%s" "$RELEASE" | jq -r '.upload_url'`

if [ -x $UPLOAD_URL ]
then
    echo "Release does not exists"
    exit 1
else
    RESPONSE=`curl --data-binary "@$GH_DEB_NAME" -H "Authorization: token $TOKEN" -H "Content-Type: application/octet-stream" "$UPLOAD_URL=$GH_DEB_NAME"`
fi

# clean the VM
rm -rf $DEB_NAME
rm $GH_DEB_NAME

. ./bin/uninstall-ff++
