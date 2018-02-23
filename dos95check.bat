echo off
nuweb binit.w
echo id31sum
g95 id31sum.f90 -o test/id31sum.exe 
echo id31offsets
g95 id31offsets.f90 -o test/id31offsets.exe 
echo id31sumall
g95 id31sumall.f90 -o test/id31sumall.exe 
rem df %OPTS% /exe:sifit.exe       sifit.f90        