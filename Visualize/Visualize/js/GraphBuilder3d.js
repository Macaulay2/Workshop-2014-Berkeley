var renderer = null,
scene = null,
camera = null,
container = null,
isMouseDown = false,
bigObject = null,
windowWidth = 0,
windowHeight = 0,
clock = null,
nodes3d = [],
mouse = new THREE.Vector2(),
objects = [],
offset = new THREE.Vector3(),
projector = new THREE.Projector(),
objectSelected, objectIntersected,
plane = null,
currentObject = null;

function getRandomColor() {
    var letters = '0123456789ABCDEF'.split('');
    var color = '#';
    for (var i = 0; i < 6; i++ ) {
        color += letters[Math.round(Math.random() * 15)];
    }
    return color;
}

function initializeBuilder3d() {

    container = document.getElementById("canvasElement3d");
    windowWidth = window.innerWidth;
    windowHeight = window.innerHeight -  125;

    renderer = new THREE.WebGLRenderer( {antialias:true} );


    renderer.setSize(windowWidth, windowHeight);
    renderer.setClearColor(0xffffff, 1);

    container.appendChild( renderer.domElement );

    // Create a new Three.js scene
    scene = new THREE.Scene();

    // Put in a camera
    camera = new THREE.PerspectiveCamera( 60, windowWidth/ windowHeight, 1, 20000 );
    camera.position.set( 0, 0, 600 );

    // EVENTS
    THREEx.WindowResize(renderer, camera);
    THREEx.FullScreen.bindKey({ charCode : 'f'.charCodeAt(0) });
    // CONTROLS
    controls = new THREE.OrbitControls( camera, renderer.domElement );

    bigObject = new THREE.Object3D();
    
    
    var material = new THREE.LineBasicMaterial({ color: 0xff0000, linewidth: 1 });
    var geometry = new THREE.Geometry();
    geometry.vertices.push(new THREE.Vector3(120,0,0));
    geometry.vertices.push(new THREE.Vector3(0,0,0));
    var xAxis = new THREE.Line(geometry, material);
    bigObject.add(xAxis);
    
    material = new THREE.LineBasicMaterial({ color: 0x00ff00, linewidth: 1 });
    geometry = new THREE.Geometry();
    geometry.vertices.push(new THREE.Vector3(0,120,0));
    geometry.vertices.push(new THREE.Vector3(0,0,0));
    var yAxis = new THREE.Line(geometry, material);
    bigObject.add(yAxis);
    
    
    material = new THREE.LineBasicMaterial({ color: 0x0000ff, linewidth: 1 });
    geometry = new THREE.Geometry();
    geometry.vertices.push(new THREE.Vector3(0,0,120));
    geometry.vertices.push(new THREE.Vector3(0,0,0));
    var zAxis = new THREE.Line(geometry, material);
    bigObject.add(zAxis);
    

    addMouseHandler(); 
    //var light = new THREE.DirectionalLight( 0x00ffff, 1.5);
    //light.position.set(0, 0, 1);
    //scene.add( light );

    scene.add(bigObject);

    plane = new THREE.Mesh( new THREE.PlaneGeometry( 2000, 2000, 8, 8 ),
        new THREE.MeshBasicMaterial( { color: 0x000000, opacity: 0.25, transparent: true, wireframe: true } ) );
    plane.visible = false;
    scene.add( plane );


}

function generateGraph() {
    //takes an ajacency matrix representation of a graph g and draws it in Three dimensional space
    renderer.setClearColor(0xffffff, 1);
    var size = dataData[0].length;
    var numRows = Math.floor(size / 4);
    var curNode = null;
    for (var i = 0; i<size; i++) {
        curNode = new node(labelData[i], getRandomColor());
        curNode.initialize();
        curNode.setCoords((i*40) - 280, (i % numRows)*80 - 160, Math.floor(Math.random()*600) - 300 );
        //console.log("x:" + (-240 + i*80)+ "   y:" + (((i % numRows) * 160)  + "   z:" + (Math.floor(Math.random()*480) - 240) );
        nodes3d.push(curNode);
    }
    for (var i = 0; i<size; i++) {
        for (var j = i; j < size; j++) {
            if (dataData[i][j]==1) {
                addEdge(nodes3d[i],nodes3d[j]);
            }
        }
    }



    }

    function node(n, c) {
        this.color = c;
        this.name = n;
        this.x = 0;
        this.y = 0;
        this.z = 0;
        this.mesh = null;
        this.label = null;
        this.edgeVectors = [];

        this.initialize = function() {
            var material = new THREE.MeshBasicMaterial({color: c});
            var geometry = new THREE.SphereGeometry(16, 32, 32);
            this.mesh = new THREE.Mesh(geometry, material);
            this.mesh["nodey"] = this;
            this.label = makeTextSprite( this.name, { fontsize: 24, backgroundColor: {r:255, g:100, b:100, a:1} } );
            bigObject.add( this.label);
            bigObject.add(this.mesh);
            objects.push(this.mesh);
        };

        this.setX = function(xx) {
            this.x = xx;
            this.mesh.position.x = xx;
        };
        this.setY = function(yy) {
            this.y = yy;
            this.mesh.position.y = yy;
        };
        this.setZ = function(zz) {
            this.z = zz;
            this.mesh.position.z = zz;
        };
        this.setCoords = function(xx,yy,zz) {
            this.x = xx;
            this.mesh.position.x = xx;
            this.y = yy;
            this.mesh.position.y = yy;
            this.z = zz;
            this.mesh.position.z = zz;
            var zAdd = (Math.abs(zz) * (-15 / 240) + 15);
                    //if (zz<0) 
                        //zAdd = -zAdd;
                    //this.label.position = (new THREE.Vector3(xx + (Math.abs(xx) * (-5 / 240) + 5), yy + (Math.abs(yy) * (-5 / 240) + 5), zz + zAdd)).multiplyScalar(1.1);
            this.label.position = (new THREE.Vector3());
            updateLabels();
        };
        this.updateLocation = function() {
            this.x = this.mesh.position.x;
            this.y = this.mesh.position.y;
            this.z = this.mesh.position.z;
            for (var i = 0; i<this.edgeVectors.length; i++) {
                console.log(this.edgeVectors[i]);
                this.edgeVectors[i].set(this.mesh.position.x, this.mesh.position.y,this.mesh.position.z);
            }
            this.label.position = (new THREE.Vector3());
            updateLabels();
        };

            }

    function updateLabels() {
        var cameraDistance = Math.sqrt(Math.pow(camera.position.x,2) + Math.pow(camera.position.y,2) + Math.pow(camera.position.z,2));

        var scaledDownCameraX = (camera.position.x);
        var scaledDownCameraY = (camera.position.y);
        var scaledDownCameraZ = (camera.position.z);


        //console.log("scaledx: " + scaledDownCameraX + "  scaledy: " + scaledDownCameraY + "  scaledz: " + scaledDownCameraZ);
        
        var curNode = null;
        for (var i = 0; i < nodes3d.length; i++) {
            curNode = nodes3d[i];
            var scalar = 20 / Math.sqrt(Math.pow(scaledDownCameraX - curNode.x, 2) + Math.pow(scaledDownCameraY - curNode.y, 2) +  Math.pow(scaledDownCameraZ - curNode.z, 2) );
            
            console.log("huur\n");
            console.log(curNode);
            curNode.label.position.x = curNode.x - (scalar * (-scaledDownCameraX + curNode.x));
            curNode.label.position.y = curNode.y - (scalar * (-scaledDownCameraY + curNode.y));
            curNode.label.position.z = curNode.z - (scalar * (-scaledDownCameraZ + curNode.z));
            //console.log("Node["+i+"]  x: " + curNode.x + "   labelX: " + curNode.label.position.x);
            //console.log("Node["+i+"]  y: " + curNode.y + "   labelY: " + curNode.label.position.y);
            //console.log("Node["+i+"]  z: " + curNode.z + "   labelZ: " + curNode.label.position.z);

        }
    }

function addEdge(node1, node2) {

    var material = new THREE.LineBasicMaterial({color:0x000000, linewidth: 2 });

    var geometry = new THREE.Geometry();
    var firstVector = new THREE.Vector3(node1.x,node1.y,node1.z);
    node1.edgeVectors.push(firstVector);
    geometry.vertices.push(firstVector);
    var secondVector = new THREE.Vector3(node2.x,node2.y,node2.z);
    node1.edgeVectors.push(secondVector);
    geometry.vertices.push(secondVector);

    var line = new THREE.Line(geometry, material);
    bigObject.add(line);
    
}

function makeTextSprite( message, parameters )
{
    if ( parameters === undefined ) parameters = {};
    
    var fontface = parameters.hasOwnProperty("fontface") ? 
    parameters["fontface"] : "Arial";
    
    var fontsize = parameters.hasOwnProperty("fontsize") ? 
    parameters["fontsize"] : 18;
    
    var borderThickness = parameters.hasOwnProperty("borderThickness") ? 
    parameters["borderThickness"] : 4;
    
    var borderColor = parameters.hasOwnProperty("borderColor") ?
    parameters["borderColor"] : { r:0, g:0, b:0, a:1.0 };
    
    var backgroundColor = parameters.hasOwnProperty("backgroundColor") ?
    parameters["backgroundColor"] : { r:255, g:255, b:255, a:1.0 };

    //var spriteAlignment = parameters.hasOwnProperty("alignment") ?
    //  parameters["alignment"] : THREE.SpriteAlignment.topLeft;

    var spriteAlignment = THREE.SpriteAlignment.topLeft;
    

    var canvas = document.createElement('canvas');
    var context = canvas.getContext('2d');
    context.font = "Bold " + fontsize + "px " + fontface;
    
    // get size data (height depends only on font size)
    var metrics = context.measureText( message );
    var textWidth = metrics.width;
    
    // background color
    context.fillStyle   = "rgba(" + backgroundColor.r + "," + backgroundColor.g + ","
      + backgroundColor.b + "," + backgroundColor.a + ")";
    // border color
    context.strokeStyle = "rgba(" + borderColor.r + "," + borderColor.g + ","
      + borderColor.b + "," + borderColor.a + ")";

context.lineWidth = borderThickness;
roundRect(context, borderThickness/2, borderThickness/2, textWidth + borderThickness, fontsize * 1.4 + borderThickness, 6);
    // 1.4 is extra height factor for text below baseline: g,j,p,q.
    
    // text color
    context.fillStyle = "rgba(0, 0, 0, 1.0)";

    context.fillText( message, borderThickness, fontsize + borderThickness);
    
    // canvas contents will be used for a texture
    var texture = new THREE.Texture(canvas) 
    texture.needsUpdate = true;

    var spriteMaterial = new THREE.SpriteMaterial( 
        { map: texture, useScreenCoordinates: false, alignment: spriteAlignment } );
    var sprite = new THREE.Sprite( spriteMaterial );
    sprite.scale.set(100,50,1.0);
    return sprite;  
}

// function for drawing rounded rectangles
function roundRect(ctx, x, y, w, h, r) 
{
    ctx.beginPath();
    ctx.moveTo(x+r, y);
    ctx.lineTo(x+w-r, y);
    ctx.quadraticCurveTo(x+w, y, x+w, y+r);
    ctx.lineTo(x+w, y+h-r);
    ctx.quadraticCurveTo(x+w, y+h, x+w-r, y+h);
    ctx.lineTo(x+r, y+h);
    ctx.quadraticCurveTo(x, y+h, x, y+h-r);
    ctx.lineTo(x, y+r);
    ctx.quadraticCurveTo(x, y, x+r, y);
    ctx.closePath();
    ctx.fill();
    ctx.stroke();   
}

function addMouseHandler()
{
    var dom = renderer.domElement;

    dom.addEventListener( 'mousemove', onMouseMove, false);
    dom.addEventListener( 'mousedown', onMouseDown, false);
    dom.addEventListener( 'mouseup' , onMouseUp, false);
}

function onMouseDown(event)
{
    event.preventDefault();
    isMouseDown = true;
    var vector = new THREE.Vector3( mouse.x, mouse.y, 0.5 );
    projector.unprojectVector( vector, camera );

    var raycaster = new THREE.Raycaster( camera.position, vector.sub( camera.position ).normalize() );

    var intersects = raycaster.intersectObjects( objects );

    if ( intersects.length > 0 ) {

        controls.enabled = false;

        objectSelected = intersects[ 0 ].object;

        var intersects = raycaster.intersectObject( plane );
        offset.copy( intersects[ 0 ].point ).sub( plane.position );
        console.log(intersects[0].object);
        for (var i = 0; i<nodes3d.length;i++) {
            if (intersects[0]["object"] == nodes3d[i].mesh) {

                currentObject = nodes3d[i];
                
            }
        }

        container.style.cursor = 'move';

    }
    
}   
function onMouseUp (event) {
    event.preventDefault();

    controls.enabled = true;

    if ( objectIntersected ) {

        plane.position.copy( objectIntersected.position );

        objectSelected = null;

    }

    container.style.cursor = 'auto';
    isMouseDown = false;

}

function onMouseMove( event ) {
    event.preventDefault();

    if (isMouseDown) {
        updateLabels();
    }

    mouse.x = ( event.clientX / window.innerWidth  ) * 2 - 1;
    mouse.y = - ( event.clientY / (window.innerHeight - 125) ) * 2 + 1;

    //

    var vector = new THREE.Vector3( mouse.x, mouse.y, 0.5 );
    projector.unprojectVector( vector, camera );

    var raycaster = new THREE.Raycaster( camera.position, vector.sub( camera.position ).normalize() );

    if ( objectSelected ) {

        var intersects = raycaster.intersectObject( plane );
        objectSelected.position.copy( intersects[ 0 ].point.sub( offset ) );
        objectSelected["nodey"].updateLocation();
        return;

    }

    var intersects = raycaster.intersectObjects( objects );

    if ( intersects.length > 0 ) {

        if ( objectIntersected != intersects[ 0 ].object ) {

            if ( objectIntersected ) objectIntersected.material.color.setHex( objectIntersected.currentHex );

            objectIntersected = intersects[ 0 ].object;
            objectIntersected.currentHex = objectIntersected.material.color.getHex();

            plane.position.copy( objectIntersected.position );
            

            plane.lookAt( camera.position );

        }

        container.style.cursor = 'pointer';

    } 
    else {

        if ( objectIntersected ) objectIntersected.material.color.setHex( objectIntersected.currentHex );

        objectIntersected = null;

        container.style.cursor = 'auto';

    }
}


function runIt() {
   // Render the scene
   
   requestAnimationFrame(runIt);
   
   
   renderer.render( scene, camera );

   controls.update();
   
}







