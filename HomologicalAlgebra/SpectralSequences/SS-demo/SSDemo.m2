-- Spectral Sequence Demo
--restart
--uninstallPackage"SpectralSequences"
--installPackage"SpectralSequences"
--installPackage("SpectralSequences", RemakeAllDocumentation => true)
--check "SpectralSequences"


restart
needsPackage"SpectralSequences"
viewHelp"SpectralSequences"
-- Example 1 --
-- We compute the Serre Spectral Sequence arising from the Hopf Fibration
-- SS^1 --> SS^3 --> SS^2.  This example is made possible by 
-- the minimal triangualtion of this fibration given 
-- by K.V. Madahar and K.S Sarkaria. Geom Dedicata, 2000.
restart
needsPackage "SpectralSequences"
needs "~/Desktop/Workshop-2014-Berkeley/HomologicalAlgebra/SpectralSequences/SS-demo/Hopf-Preload.m2";
K = filteredComplex({S3,F1S3,F0S3}, ReducedHomology => false);
-- the filtered complex K above arises from a simplicial 
-- realization of the hopf fibration SS^1 --> SS^3 --> SS^2
E = prune spectralSequence K		    
E^0		    
E^0 .dd
E^1
E^1 .dd
E^2
-- can check that the E^2 page has been computed correctly
E^2 .dd
E^3
prune HH K_infinity
--
-- Example 2 --
-----------------------------------------------------------------
-- Examples related to g^1_3's on a general curve of genus 4. 
-----------------------------------------------------------------
-- An example related to g^1_3s on a canoncially embedded genus 4 curve 
-- C \subseteq X := PP^1 x PP^1
-- The following uses the hypercohomology spectral sequence to
-- compute the cohomology groups of a g^1_3 on C
-- Of course the hypercohomology spectral sequence is not
-- very interesting but it is still a proof of concept.
restart
needsPackage"SpectralSequences"
	         -- C \subseteq PP^1 x PP^1 type (3,3)
		 -- Use hypercohomology to compute HH OO_C(1,0) 
		 R = ZZ/101[a_0..b_1, Degrees=>{2:{1,0},2:{0,1}}]; -- PP^1 x PP^1
		 B = intersect(ideal(a_0,a_1),ideal(b_0,b_1)) ; -- irrelevant ideal
		 B = B_*/(x -> x^5)//ideal ; -- Sufficentily high Frobenius power 
		 G = res image gens B ;
	  	 I = ideal random(R^1, R^{{-3,-3}}) ; -- ideal of C
		 F = res comodule I ;
		 K = Hom(G , filteredComplex (F ** R^{{1,0}})) ; -- Twist F by a line of ruling and make filtered complex whose ss abuts to HH OO_C(1,0) 
		 E = prune spectralSequence K ; --the spectral sequence degenerates on the second page 
		 E^1 
		 E^2 ; -- output is a mess
--	       The cohomology groups we want are obtained as follows.
		 basis({0,0}, E^2_{0,0}) --  == HH^0 OO_C(1,0)
		 basis({0,0}, E^2_{1,-2}) --  == HH^1 OO_C(1,0)	 
		 basis({0,0}, E^2)
-- Example 3 --
restart
needsPackage"SpectralSequences"
--     	  I-adic filtrations of chain complexes and their spectral sequences
--	  By multiplying a chain complex by sucessive powers of an ideal we obtain a 
--        filtered complex.       
	      B = QQ[a..d]
	      J = ideal vars B
	      C = complete res monomialCurveIdeal(B,{1,3,4})
	      K = filteredComplex(J,C,4)
--      Here are higher some pages of the associated spectral sequence:
	       E = prune spectralSequence K
	       E^2
	       E^3
	       E^3 .dd
	       E^3 .dd_{-1,2}
    	       source E^3 .dd_{-1,2}	
    	       target E^3 .dd_{-1,2}
	       E^4
	       E^4 .dd
    	       E^infinity
    	       hilbertPolynomial(E^4)    	    	
               prune HH K_infinity
    	       hilbertPolynomial(HH_0 K_infinity)
    	       E^infinity
    	       prune associatedGradedHomologyObject(-4,0,K)
	       E^infinity _{-4,4}
    	       prune associatedGradedHomologyObject(-3,0,K)
    	       E^infinity _{-3,3}
    	       --etc
	       
--- scratch related to change of rings for Tor

	       