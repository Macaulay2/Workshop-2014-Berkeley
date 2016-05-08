--*************************************************
--*************************************************
--This file contains the documentation for the 
--Fsing package.
--*************************************************
--*************************************************

beginDocumentation()

doc ///
   Key
      Fsing
   Headline
      A package for calculations of singularities in positive characteristic 
   Description
      Text    
         This will do a lot of cool stuff someday. 
///

--***********************************************
--***********************************************
--Documentation for IntegerComps.m2
--***********************************************
--***********************************************

doc ///
     Key
     	floorlog
     Headline
        Computes the floor of the log base b of x
     Usage
     	 floorlog(b,x)
     Inputs 
     		b:ZZ
		x:ZZ		
     Outputs
         :ZZ
     Description
	Text
	    This differs from floor(log_b(x)) in that it corrects problems due to rounding.
///

doc ///
     Key
     	multOrder
     	(multOrder, ZZ, ZZ)
     Headline
        Computes the multiplicative order of a modulo b
     Usage
     	 multOrder(a,b)
     Inputs 
     		a:ZZ
		b:ZZ		
     Outputs
         :ZZ
     Description
	Text
	    This computes the multiplicative order of a modulo b.  If a and b are not relatively prime, it returns an error.
///

