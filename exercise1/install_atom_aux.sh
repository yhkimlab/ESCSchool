wdir=$(pwd)

wget https://launchpad.net/xmlf90/trunk/1.5/+download/xmlf90-1.5.0.tgz // 1.5.0버전
wget https://launchpad.net/libgridxc/trunk/0.8/+download/libgridxc-0.8.5.tgz // 0.8.5버전
tar -xvzf xmlf90-1.5.0.tgz
tar -xvzf libgridxc-0.8.5.tgz
cd xmlf90-1.5.0
mkdir Gfortran
cd Sys
cp gfortran.make ../Gfortran/fortran.mk
cd ../Gfortran
sh ../config.sh
make

cd $wdir
cd libgridxc-0.8.5
mkdir Gfortran
cd Gfortran
cp ../extra/fortran.mk .
sh ../src/config.sh
make clean
make

cd $wdir
wget https://siesta-project.org/SIESTA_MATERIAL/Pseudos/Code/atom-4.2.7-100.tgz // 4.2.7버전
tar xvzf atom-4.2.7-100.tgz
cd atom-4.2.7-100
cp arch.make.sample arch.make
#vi arch.make

