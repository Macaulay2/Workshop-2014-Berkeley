--------------------------------------------------------------------------
-- PURPOSE : Visualize package for Macaulay2 provides the ability to 
-- visualize various algebraic objects in java script using a 
-- modern browser.
--
-- Copyright (C) 2013 Branden Stone and Jim Vallandingham
--
-- This program is free software; you can redistribute it and/or
-- modify it under the terms of the GNU General Public License version 2
-- as published by the Free Software Foundation.
--
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU General Public License for more details.
--------------------------------------------------------------------------


newPackage(
	"Visualize",
    	Version => "0.2", 
    	Date => "October 2, 2013",
    	Authors => {       
	     {Name => "Elliot Korte", Email => "ek2872@bard.edu"},	     
	     {Name => "Will Smith", Email => "smithw12321@gmail.com"},		
	     {Name => "Branden Stone", Email => "bstone@bard.edu", HomePage => "http://www.bard.edu/~bstone/"},	     
	     {Name => "Jim Vallandingham", Email => "vlandham@gmail.com", HomePage => "http://vallandingham.me/"}
	     },
    	Headline => "Visualize",
    	DebuggingMode => true
    	)

export {
    
    -- Options
     "Path",
     "visTemplate",
    
    -- Methods
     "visIntegralClosure",
     "visIdeal",
     "visGraph",
     "runServer", --helper
     "toArray", --helper
     "getCurrPath" --helper
}

needsPackage"Graphs"


------------------------------------------------------------
-- METHODS
------------------------------------------------------------

-- Input: None.
-- Output: String containing current path.

getCurrPath = method()
installMethod(getCurrPath, () -> (local currPath; currPath = get "!pwd"; substring(currPath,0,(length currPath)-1)|"/"))


--input: A list of lists
--output: an array of arrays
--
-- would be nice if we could use this on any nesting of lists/seq
--
toArray = method() 
toArray(List) := L -> (
     return new Array from apply(L, i -> new Array from i);
     )
    


--input: A path
--output: runs a server for displaying objects
--
runServer = method(Options => {Path => getCurrPath()})
runServer(String) := opts -> (visPath) -> (
    return run visPath;
    )


--input: Three Stings. The first is a key word to look for.  The second
--    	 is what to replace the key word with. The third is the path 
--    	 where template file is located.
--output: A file with visKey replaced with visString.
--
visOutput = method(Options => {Path => getCurrPath()})
visOutput(String,String,String) := opts -> (visKey,visString,visTemplate) -> (
    local fileName; local openFile; local PATH;
    
    fileName = (toString currentTime() )|".html";
    PATH = opts.Path|fileName;
    openOut PATH << 
    	replace(visKey, visString , get visTemplate) << 
	close;
                  
    return show new URL from { "file://"|PATH };
    )



--input: A monomial ideal of a polynomial ring in 2 or 3 variables.
--output: The newton polytope of the of the ideal.
--
visIdeal = method(Options => {Path => getCurrPath()|"/temp-files/", visTemplate => getCurrPath() |"/templates/visIdeal/visIdeal"})
visIdeal(Ideal) := opts -> J -> (
    local R; local arrayList; local arrayString; local numVar; local visTemp;
    
    R = ring J;
    numVar = rank source vars R;
    
    if ((numVar != 2) and (numVar != 3)) then (error "Ring needs to have either 2 or 3 variables.");
    
    if numVar == 2 
    then (
    	visTemp = opts.visTemplate|"2D.html";
	arrayList = apply( flatten entries gens J, m -> flatten exponents m);	
	arrayList = toArray arrayList;
	arrayString = toString arrayList;
    )
    else (
	visTemp = opts.visTemplate|"3D.html";
    	arrayList = apply(flatten entries basis(0,infinity, R/J), m -> flatten exponents m );
    	arrayList = toArray arrayList;
    	arrayString = toString arrayList;
    );
	
    return visOutput( "visArray", arrayString, visTemp, Path => opts.Path );
    )

--input: A graph
--output: the graph in the browswer
--
visGraph = method(Options => {Path => getCurrPath()|"/temp-files/", visTemplate => getCurrPath() | "/templates/visGraph/visGraph-template.html"})
visGraph(Graph) := opts -> G -> (
    local A; local arrayList; local arrayString;
    
    A = adjacencyMatrix G;
    arrayList = toArray entries A;
    arrayString = toString arrayList;
    
    return visOutput( "visArray", arrayString, opts.visTemplate, Path => opts.Path );
    )


--------------------------------------------------
-- DOCUMENTATION
--------------------------------------------------

-- use simple doc
beginDocumentation()

document {
     Key => Visualize,
     Headline => "A package to help visualize algebraic objects in the browser using javascript.",
     
     "Lots of cool stuff happens here.",
     
     PARA{}, "For the mathematical background see ",

     
     UL {
	  {"Winfried Bruns and Jürgen Herzog.", EM " Cohen-Macaulay Rings."},
	},
     
     }

-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------

end

-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------

-----------------------------
-----------------------------
-- Stable Tests
-----------------------------
-----------------------------

restart
loadPackage"Graphs"
loadPackage"Visualize"

G = graph({{x_0,x_1},{x_0,x_3},{x_0,x_4},{x_1,x_3},{x_2,x_3}},Singletons => {x_5})
visGraph G

R = QQ[x,y,z]
I = ideal"x4,xyz3,yz,xz,z6,y5"
visIdeal I
visIdeal( I, Path => "/Users/bstone/Desktop/Test/")

S = QQ[x,y]
I = ideal"x4,xy3,y50"
visIdeal I
visIdeal( I, Path => "/Users/bstone/Desktop/Test/")

-----------------------------
-----------------------------
-- Demo
-----------------------------
-----------------------------


restart
loadPackage"Graphs"
loadPackage"Visualize"

G = graph({{x_0,x_1},{x_0,x_3},{x_0,x_4},{x_1,x_3},{x_2,x_3}},Singletons => {x_5})
visGraph G

R = QQ[x,y,z]
I = ideal"x4,xyz3,yz,xz,z6,y5"
visIdeal I
visIdeal( I, Path => "/Users/bstone/Desktop/Test/")

S = QQ[x,y]
I = ideal"x4,xy3,y50"
visIdeal I
visIdeal( I, Path => "/Users/bstone/Desktop/Test/")


