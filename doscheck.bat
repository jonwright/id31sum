echo off
set OPTS=/warn:all /assume:accuracy_sensitive /fast  
nuweb binit.w
latex binit
nuweb binit.w
latex binit
dvipdfm -p a4 binit
dvips binit
df %OPTS% /exe:id31sum.exe     id31sum.f90
df %OPTS% /exe:id31sumall.exe  id31sumall.f90
df %OPTS% /exe:id31offsets.exe id31offsets.f90      
df %OPTS% /exe:sifit.exe       sifit.f90        