
function generate2dGraph(rawData) {

    var w = window.innerWidth;
    var h = window.innerHeight - 125;
    var color = d3.scale.category10();


    var data = rawData;

    //convert adjacency matrix into edge pairs for d3 to process
    var nodeList = [];
    var edgeList = [];
    for (var i = 0; i<data.length; i++) {

        nodeList.push( {name: i } );
        //console.log("I " + i);

    }
    for (var i = 0; i<data.length; i++) {
        for (var j = 0; j< data[i].length; j++) {
            if (data[i][j] == 1) {
                edgeList.push( { source: i, target: j} );
                //console.log("I " + i + "  J " + j);
            }    
        }
    }
    
    //begin d3
    
    var svg = d3.select("body")
    .append("svg")
    .attr("height", h)
    .attr("width", w)
    .attr("id", "canvasElement");

    var dataset = {
        nodes: nodeList,
        edges: edgeList
    };

    var force = d3.layout.force()
    .nodes(dataset.nodes)
    .links(dataset.edges)
    .size([w, h])
    .gravity(.1)
    .linkDistance([150])
    .charge([-1000])
    .start();

    // https://github.com/mbostock/d3/wiki/Force-Layout#wiki-drag
    // implement 'dragstart' function that fixes nodes so they are
    // not repositioned by force layout.
    var drag = force.drag()
    .on("dragstart", dragstart);

    function dragstart(d) {
        d.fixed = true;
        d3.select(this).classed("fixed", true);
    }


    var edges = svg.selectAll("line")
    .data(dataset.edges)
    .enter()
    .append("line")
    .style("stroke", "#ccc")
    .style("stroke-width", 2);

    var nodes = svg.selectAll("g")
    .data(dataset.nodes)
    .enter()
    .append("g")
    .call(force.drag);

    nodes.append("circle")
    .attr("r", 10)
    .style("fill", function(d, i) {
        return color(i);
    });

    nodes.append("text")
    .attr("font-family", "sans-serif")
    .attr("font-size", "12px")
    .attr("fill", "white")
    .attr("text-anchor", "middle")
    .attr("dy", "4px")
    .text(function(d) { return d.name; });


    force.on("tick", function() {

        edges.attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });

        nodes.attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; });                 
        
    });

}

function updateWindowSize2d() {
    if (curView == "2d") {
        console.log("uhhhh"); 
        var svg = document.getElementById("canvasElement");
        svg.style.width = window.innerWidth;
        svg.style.height = window.innerHeight - 125;
        svg.width = window.innerWidth;
        svg.height = window.innerHeight - 125;
    }

}
