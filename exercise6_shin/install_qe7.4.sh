mkdir install
cd install
wget https://gitlab.com/QEF/q-e/-/archive/qe-7.4/q-e-qe-7.4.tar
tar -xvf q-e-qe-7.4.tar
cd q-e-qe-7.4
./configure
make pw pp