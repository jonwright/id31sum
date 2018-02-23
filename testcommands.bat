testspec test.dat > testspec.out
testbins1         > testbins1.out
testbins2         > testbins2.out
echo bindumpedfile.out > in
bindump test.dat 0.01 1  < in > bindump.out
id31sum test.dat 0.01 1 1 
id31sumall test.dat 0.01 1 1
                                                                             