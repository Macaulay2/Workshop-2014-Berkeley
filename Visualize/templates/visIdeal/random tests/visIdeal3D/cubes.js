var dataset = [[5,0,0], [0,4,0], [0,0,6], [1, 2, 4], [3,3,2]]

var cubeArray = []

//get size of cube field
xMax = 0;
yMax = 0;
zMax = 0;
for ( i = 0; i < dataset.length; i++ ) {
    if ( dataset[i][0] > xMax ) { xMax = dataset[i][0]; }
    if ( dataset[i][1] > yMax ) { yMax = dataset[i][1]; }
    if ( dataset[i][2] > zMax ) { zMax = dataset[i][2]; }
}

//figure out where to put cubes
for ( i = 0; i < xMax; i++ ) {
    console.log(i);
    for ( j = 0; j < yMax; j++ ) {
	for ( k = 0; k < zMax; k++ ) {
	    check = true;
	    console.log(check);
	    for ( m = 0; m < dataset.length; m++ ) {
		if ( i >= dataset[m][0] && j >= dataset[m][1] && k >= dataset[m][2] ) { check = false; }
	    }
	    if( check ) {
		cubeArray.push([i, j, k]);
	    }
	}
    }
}


// make 3js model, etc.

var camera, scene, renderer;
var geometry, material, mesh;

// make a cube with given x, y, and z coordinates

var makeCube = function(xPos, yPos, zPos) {
    geometry = new THREE.CubeGeometry( 1, 1, 1 );
    if ( (xPos + yPos + zPos) % 2 === 0 ) {
	material = new THREE.MeshBasicMaterial( { color: 0x82A2E8, wireframe: false, vertexColors: true } );
    } else {
	material = new THREE.MeshBasicMaterial( { color: 0xE8C882, wireframe: false, vertexColors: true } );
    }
    mesh = new THREE.Mesh( geometry, material );
    mesh.position.x = xPos;
    mesh.position.y = yPos;
    mesh.position.z = zPos;
    return mesh;
}

var distance = Math.max(Math.max(xMax, yMax), zMax) + 5;

function init() {

    camera = new THREE.PerspectiveCamera( 75, window.innerWidth / window.innerHeight, 0.01, 10000 );
    camera.position.z = distance;
    
    scene = new THREE.Scene();
    
    cubes = []
    for ( i = 0; i < cubeArray.length; i++ ) {
	cubes.push(makeCube(cubeArray[i][0], cubeArray[i][1], cubeArray[i][2]));
    }


    for( i = 0; i < cubes.length; i++ ) {
	scene.add(cubes[i]);
    }

    var light = new THREE.PointLight( 0xffffff, 1, 100 );
    light.position.set( 6, 6, 6 );
    scene.add( light );
    
    renderer = new THREE.CanvasRenderer();
    renderer.setSize( window.innerWidth, window.innerHeight );
    
    document.body.appendChild( renderer.domElement );
    controls = new THREE.OrbitControls(camera, renderer.domElement);
    
}

function animate() {

    requestAnimationFrame( animate );
    renderer.render( scene, camera );
    controls.update();
}


init();
animate();