newPackage("NCAlgebraV2",
     Headline => "Additional functions for NCAlgebra",
     Version => "0.1",
     Date => "January 6th, 2014",
     Authors => {
	  {Name => "Frank Moore",
	   HomePage => "http://www.math.wfu.edu/Faculty/Moore.html",
	   Email => "moorewf@wfu.edu"},
	  {Name => "Andrew Conner",
	   HomePage => "http://www.math.wfu.edu/Faculty/Conner.html",
	   Email => "connerab@wfu.edu"},
          {Name => "Courtney Gibbons",
	   HomePage => "",
	   Email => "crgibbon@hamilton.edu"}},
     AuxiliaryFiles => true,
     DebuggingMode => true
     )

export {}

needsPackage "NCAlgebra"

--- bug fix/performance/interface improvements
------------------------------------
--- Testing!

--- additions in the near future
------------------------------------
--- skewPolynomialRing with NCRing bases
--- Finish right mingens
--- Implement left kernels and mingens etc (opposite ring now done)
--- Make sure that trivial ideals are handled correctly
--- Make sure constructions over the tensor algebra are handled correctly
--- isFiniteDimensional?
--- 'basis' for f.d. algebras?
--- Generating set for algebras not over a field
--- fix coefficients to return a pair again.

--- other things to add or work on in due time
-----------------------------------
--- notation to refer to a certain graded piece of an algebra e.g. A_3
--- a dimension function for graded pieces dim(A,17)
--- Dare I say it, Diamond Lemma?
--- Make Quotients of Quotients work.
--- NCRingMap kernels (to a certain degree)  -- Not sure I can do this with
---   Bergman, can't use block orders in bergman.
--- anick resolution
--- NCModules (?) (including module gb (via simple), hilbert series, modulebettinumbers)
--- Work on reduction code a bit?
--- Hilbert series for modules (and sided ideals)
--- enveloping algebras, tensor products of algebra, Hochschild (co)homology?

--- symbols in hash tables of exported types
----
restart
uninstallPackage "NCAlgebraV2"
installPackage "NCAlgebraV2"
needsPackage "NCAlgebraV2"
viewHelp "NCAlgebraV2"

