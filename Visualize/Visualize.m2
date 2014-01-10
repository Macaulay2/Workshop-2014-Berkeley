---------------------------------------------------------------------------
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
     	     {Name => "Brett Barwick", Email => "Brett@barwick.edu", HomePage => "http://www.bard.edu/~bstone/"},	     
	     {Name => "Elliot Korte", Email => "ek2872@bard.edu"},	     
	     {Name => "Will Smith", Email => "smithw12321@gmail.com"},		
	     {Name => "Branden Stone", Email => "bstone@bard.edu", HomePage => "http://www.bard.edu/~bstone/"},	     
	     {Name => "Jim Vallandingham", Email => "vlandham@gmail.com", HomePage => "http://vallandingham.me/"}
	     },
    	Headline => "Visualize",
    	DebuggingMode => true,
	AuxiliaryFiles => true,
	Configuration => {"DefaultPath" => concatenate(currentDirectory(),"temp-files/") } 
    	)

export {
    
    -- Options
     "VisPath",
     "VisTemplate",
    
    -- Methods
     "visIntegralClosure",
     "visIdeal",
     "visGraph",
     "copyJS",
     
    -- Helpers 
     "runServer",
     "toArray", 
     "getCurrPath", 
     "copyTemplate"     

}

needsPackage"Graphs"


defaultPath = (options Visualize).Configuration#"DefaultPath"

-- (options Visualize).Configuration

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
runServer = method(Options => {VisPath => currentDirectory()})
runServer(String) := opts -> (visPath) -> (
    return run visPath;
    )

--- add methods for output here:
--







--input: Three Stings. The first is a key word to look for.  The second
--    	 is what to replace the key word with. The third is the path 
--    	 where template file is located.
--output: A file with visKey replaced with visString.
--
visOutput = method(Options => {VisPath => currentDirectory()})
visOutput(String,String,String) := opts -> (visKey,visString,visTemplate) -> (
    local fileName; local openFile; local PATH;
    
    fileName = (toString currentTime() )|".html";
    PATH = opts.VisPath|fileName;
    openOut PATH << 
    	replace(visKey, visString , get visTemplate) << 
	close;
                  
    return (show new URL from { "file://"|PATH }, fileName);
    )

-- input: path to ah html file
-- output: a copy of the input file in a temporary folder
--
copyTemplate = method()
copyTemplate String := src -> (
    local fileName; local dirPath;
    
    fileName = (toString currentTime() )|".html";
    
    dirPath = temporaryFileName();
    makeDirectory dirPath;
    dirPath = concatenate(dirPath,"/",fileName);
    
    copyFile( src, dirPath);
    
    return dirPath;
)


-- input:
-- output:
searchReplace = method(Options => {VisPath => currentDirectory()})
searchReplace(String,String,String) := opts -> (oldString,newString,visSrc) -> (
    local visFilePathTemp;
    
    visFilePathTemp = temporaryFileName();
    copyFile(visSrc,visFilePathTemp);
    openOut visSrc << 
    	replace(oldString, newString , get visFilePathTemp) << 
	close;
	
    return visSrc;
    )




--input: A monomial ideal of a polynomial ring in 2 or 3 variables.
--output: The newton polytope of the of the ideal.
--
visIdeal = method(Options => {VisPath => defaultPath, VisTemplate => currentDirectory() |"Visualize/templates/visIdeal/visIdeal"})
visIdeal(Ideal) := opts -> J -> (
    local R; local arrayList; local arrayString; local numVar; local visTemp;
    local varList;
    -- local A;
    
    R = ring J;
    numVar = rank source vars R;
    varList = flatten entries vars R;
        
    if ((numVar != 2) and (numVar != 3)) then (error "Ring needs to have either 2 or 3 variables.");
    
    if numVar == 2 
    then (
	visTemp = copyTemplate(opts.VisTemplate|"2D.html");
	copyJS(replace(baseFilename visTemp, "", visTemp));

	arrayList = apply( flatten entries gens J, m -> flatten exponents m);	
	arrayList = toArray arrayList;
	arrayString = toString arrayList;
	
	searchReplace("visArray",arrayString, visTemp);
--	searchReplace("XXX",toString(varList_0), visTemp);
--	searchReplace("YYY",toString(varList_1), visTemp);
--	searchReplace("ZZZ",toString(varList_2), visTemp)
    )
    else (
	visTemp = copyTemplate(opts.VisTemplate|"3D.html");
	copyJS(replace(baseFilename visTemp, "", visTemp));
	
	
    	arrayList = apply(flatten entries basis(0,infinity, R/J), m -> flatten exponents m );
    	arrayList = toArray arrayList;
    	arrayString = toString arrayList;
	
	searchReplace("visArray",arrayString, visTemp);
	searchReplace("XXX",toString(varList_0), visTemp);
	searchReplace("YYY",toString(varList_1), visTemp);
	searchReplace("ZZZ",toString(varList_2), visTemp)
    );
    
    show new URL from { "file://"|visTemp };
--    A = visOutput( "visArray", arrayString, visTemp, VisPath => opts.VisPath );
    
    return visTemp;--opts.VisPath|A_1;
    )

--input: A graph
--output: the graph in the browswer
--
visGraph = method(Options => {VisPath => defaultPath, VisTemplate => currentDirectory() | "Visualize/templates/visGraph/visGraph-template.html"})
visGraph(Graph) := opts -> G -> (
    local A; local arrayString; local vertexString; local visTemp;
    local tempG; local vertexSet;
    
    A = adjacencyMatrix G;
    arrayString = toString toArray entries A; -- Turn the adjacency matrix into a nested array (as a string) to copy to the template html file.
    
    -- Add this back in when we figure out how to deal with the old
    -- Graphs package not knowing what G.vertexSet means.
    
    if value((options Graphs).Version) == 0.1 then (
	 vertexString = toString new Array from apply(keys(G#graph), i -> "\""|toString(i)|"\""); -- Create a string containing an ordered list of the vertices in the older Graphs package.
	 tempG = new MutableHashTable from G;
	 tempG#vertexSet = {};
	 G = new HashTable from tempG;
    ) else (
         vertexString = toString new Array from apply(G.vertexSet, i -> "\""|toString(i)|"\""); -- Create a string containing an ordered list of the vertices in the newer Graphs package.
    );
    
    visTemp = copyTemplate(currentDirectory()|"Visualize/templates/visGraph/visGraph-template.html"); -- Copy the visGraph template to a temporary directory.
    
    searchReplace("visArray",arrayString, visTemp); -- Replace visArray in the visGraph html file by the adjacency matrix.
    searchReplace("visLabels",vertexString, visTemp); -- Replace visLabels in the visGraph html file by the ordered list of vertices.

    copyJS(replace(baseFilename visTemp, "", visTemp)); -- Copy the javascript libraries to the temp folder.
    
    show new URL from { "file://"|visTemp };
    
    return visTemp;
)


--input: a String of a path to a directory
--output: Copies the js library to path
--
--caveat: Checks to see if files exist. If they do exist, the user
--        must give permission to continue. Continuing will overwrite
--        current files and cannont be undone.
copyJS = method()
copyJS(String) := dst -> (
    local jsdir; local ans; local quest;
    
    dst = dst|"js/";    
    
    -- get list of filenames in js/
    jsdir = delete("..",delete(".",
	    readDirectory(currentDirectory()|"Visualize/js/")
	    ));
    
    -- test to see if files exist in target
    if (scan(jsdir, j -> if fileExists(concatenate(dst,j)) then break true) === true)
    then (
    	   quest = concatenate(" -- Some files in ",dst," will be overwritten.\n -- This action cannot be undone.");
	   print quest;
	   ans = read "Would you like to continue? (yes or no):  ";
	   while (ans != "yes" and ans != "no") do (
	       ans = read "Would you like to continue? (yes or no):  ";
	       );  
	   if ans == "no" then (
	       error "Process was aborted."
	       );
    	);
    
    copyDirectory(currentDirectory()|"Visualize/js/",dst);
    
    return "Created directory "|dst;
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
-- branden
restart 
loadPackage"Graphs"
loadPackage"Visualize"
(options Visualize).Configuration

searchReplace("visArray","kickass string", testFile)
searchReplace("XXX","kickass string", testFile)
searchReplace("YYY","kickass string", testFile)
searchReplace("ZZZ","kickass string", testFile)

-- Old Graphs
G = graph({{x_0,x_1},{x_0,x_3},{x_0,x_4},{x_1,x_3},{x_2,x_3}},Singletons => {x_5})
visGraph G
H = graph({{x_1, x_0}, {x_3, x_0}, {x_3, x_1}, {x_4, x_0}}, Singletons => {x_2, x_5, 6, cat_sandwich})
visGraph H

-- New Graphs
G = graph(toList(0..5),{{0,1},{0,3},{0,4},{1,3},{2,3}},Singletons => {5},EntryMode => "edges")
G = graph(toList(0..5),{0,{1,2,3,4}},Singletons => {5})--,EntryMode => "edges")
visGraph G
visGraph( G, VisPath => "/Users/bstone/Desktop/Test/")
S = G.vertexSet
toString S

R = QQ[a,b,c]
I = ideal"a2,ab,b2c,c5,b4"
I = ideal"x4,xyz3,yz,xz,z6,y5"
visIdeal I
visIdeal( I, VisPath => "/Users/bstone/Desktop/Test/")

S = QQ[x,y]
I = ideal"x4,xy3,y5"
visIdeal I
visIdeal( I, VisPath => "/Users/bstone/Desktop/Test/")


copyJS "/Users/bstone/Desktop/Test/"
yes
copyJS ( currentDirectory()|"temp


-----------------------------
-- Julio's tests
-----------------------------



-----------------------------
-- end Julio's Test
-----------------------------



-----------------------------
-----------------------------
-- Demo
-----------------------------
-----------------------------


restart
loadPackage"Visualize"

-- Creates staircase diagram 
-- 2 variables
S = QQ[x,y]
I = ideal"x4,xy3,y5"
visIdeal I

-- User can choose where to place files
visIdeal( I, VisPath => "/Users/bstone/Desktop/Test/")

-- 3 variables
R = QQ[x,y,z]
J = ideal"x4,xyz3,yz2,xz3,z6,y5"
visIdeal J
visIdeal( J, VisPath => "/Users/bstone/Desktop/Test/")

restart
needsPackage"Graphs"
loadPackage"Visualize"

-- we are also focusing on graphs
G = graph({{x_0,x_1},{x_0,x_3},{x_0,x_4},{x_1,x_3},{x_2,x_3}},Singletons => {x_5})
-- displayGraph A
visGraph G

M = 
A = graph M
visGraph A
