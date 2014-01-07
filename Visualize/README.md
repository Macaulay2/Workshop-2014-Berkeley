Visualize.m2
========

A [Macaulay2](https://github.com/Macaulay2/M2) package to help visualize algebraic objects in the browser using javascript.

Current Methods
----

This package is still in development.

* `visIdeal` This allows the user to view an interactive image of the lower boundary of the exponent set of a monomial ideal in three variables.




Usage
=====

Assuming [Macaulay2](https://github.com/Macaulay2/M2) is installed on your machine, the following directions will help in the downloading and running of `Visualize.m2`.

First Clone the Repository
------

The easiest way to download the needed files is to clone the entire repository. You will need to install [Git](https://help.github.com/articles/set-up-git) for this.

```
git clone https://github.com/b-stone/Visualize-M2.git
cd Visualize-M2
```

You only need the `Visualize.m2` file and the folders `js` and `templates`.

Running in M2
----

First make sure that the file `Visualize.m2` is on the load [path](http://www.math.uiuc.edu/Macaulay2/doc/Macaulay2-1.6/share/doc/Macaulay2/Macaulay2Doc/html/_path.html). To run, execute the following. (This is assuming that `visualize.m2` is on the path `./` )

```
loadPackage "Visualize"
R = QQ[x,y,z]
I = ideal"x4,xy,yz,xz,z3,y3"
visIdeal( I,  Path => "./temp-files/" )
```

At this point your browser should open and you should have an interactive image of the lower boundary of the ideals exponent set (as can bee seen [here](http://math.bard.edu/bstone/visideal/)). As is, the `html` file that is created is saved in the `./Visualize-M2/temp-files/` directory. If you wish to create the file elsewhere, change the path, but make sure the JavaScript files are moved to your target directory.

For visualization of graphs, use the following command. (An example is found [here])(http://math.bard.edu/~bstone/visgraph/))

```
loadPackage"Graphs"
G = graph({{x_0,x_1},{x_0,x_3},{x_0,x_4},{x_1,x_3},{x_2,x_3}},Singletons => {x_5})
visGraph( G, Path => "./temp-files/" )
```



