--I want to test another example on Jason's code, so I stole that for here
restart
kerDD = method()
kerDD(Module) := M -> (
     DM = coker dual presentation M;
     Ext^1(DM, ring M)
     )

cokerDD = method()
cokerDD(Module) := M -> (
     DM = coker dual presentation M;
     Ext^2(DM, ring M)
     )

doubleDualMap = method()
doubleDualMap(Module):= M -> (
     map(dual (dual M), M, id_(ambient M))
--     map(dual (dual M), M, id_(cover M))
     )


--testing another example to try to solve the problem
--this example has doubledual module a submodule of 0
--but M itself is a submodule of ZZ^2.
--This (presumedly) is why it breaks the doubleDualMap code
--for either choice of ambient or cover above
R=ZZ
L=coker matrix{{2}}
N=coker matrix{{3}}
M=L++N
doubleDualMap(M)

--learning sufficiently much about the truncate command to implement it
restart
R=QQ[x,y,z]
M=coker matrix{{x,y,z}}
C=resolution(M,LengthLimit=>10)
loadPackage "SpectralSequences" --don't have this...
help truncate

--creating method to truncate, dualize, and shift a complex
truncateDualShift = method()
truncateDualShift(ZZ,ChainComplex) := (g,P) -> (
     if g > max C-1

