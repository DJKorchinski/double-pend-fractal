const  L = 404;
var grid_buffer1 = new ArrayBuffer(L*L);
var grid = new Uint8ClampedArray(grid_buffer1); // identity of the zero found. 
var grid_buffer2 = new ArrayBuffer(L*L);
var grid_iter = new Uint8ClampedArray(grid_buffer2); // number of (scaled) iterations to find a zero.
const scale_factor = 1;
const Npix = L*L*scale_factor*scale_factor;
var imgData = new ImageData(L*scale_factor,L*scale_factor); 
var imgDataMap_buffer = new ArrayBuffer(Npix*4); //4 bytes per pixel, for 32 bit integer addressing. 
var imgDataMap = new Int32Array(imgDataMap_buffer);

GRID_TO_MS_FACTOR = 2/100; // View grid  width / height in units of math


const xcenter_grid = Math.floor(L/2), ycenter_grid = Math.floor(L/2);
function ind_to_grid_x(ind) { return ind%L; }
function ind_to_grid_y(ind) { return Math.floor(ind/L); }
function ind_to_x(ind) { return (ind_to_grid_x(ind)-xcenter_grid) * GRID_TO_MS_FACTOR; }
function ind_to_y(ind) { return (ind_to_grid_y(ind)-ycenter_grid) * GRID_TO_MS_FACTOR; }
function grid_xy_to_ind(x,y) { return x+y*L; }

//our equation is: f(z) = z^n - 1 = 0
// df/dz = n*z^(n-1)
// -f(z) / (df/dz) = dz = dx + i dy
// let's specialize to n = 3 for now.

// (f0r + f0i) / ( f1r + f1i) = (f0r + f0i) * (f1r- f1i) / (( f1r + f1i) * (f1r- f1i))
// =  (f0r * f1r + f0i*f1i) / (f1r^2 + f1i^2) + i (f0i*f1r - f0r*f1i) / (f1r^2 + f1i^2)

MAX_ITER = 100;
var ind_to_compute = 0; 
function tick(){
    var x = ind_to_x(ind_to_compute), y = ind_to_y(ind_to_compute);
    for(var iter = 0; iter < MAX_ITER; iter ++) {
        f0_real = x*x*x - 3*x*y*y - 1;
        f0_imag = 3*x*x*y - y*y*y;
        f1_real = 3*x*x - 3*y*y;
        f1_imag = 6*x*y;
        // -f(z) / (df/dz) = dz = dx + i dy
        dx = -(f0_real*f1_real + f0_imag*f1_imag) / (f1_real*f1_real + f1_imag*f1_imag);
        dy = -(f0_imag*f1_real - f0_real*f1_imag) / (f1_real*f1_real + f1_imag*f1_imag);
        x += dx;
        y += dy;
        if(dx*dx + dy*dy < 1e-20){ break; }
    }
    grid_iter[ind_to_compute] = 255*Math.pow(Math.log(iter+1)/Math.log(MAX_ITER+1),2); 
    //identifying which root we found with a number.
    if(x > 0){
        grid[ind_to_compute] = 0;
    } else if (y > 0){
        grid[ind_to_compute] = 1;
    } else {
        grid[ind_to_compute] = 2;
    }
    //increment the index to compute.
    ind_to_compute = (ind_to_compute +1 ) % (L*L);
}

function draw() {
    //blitting every site!
    //color each pixel based on the source in the grid. 
    for(var ind = 0; ind < Npix; ind++){
        // *4 b/c each pixel is rgba. +1 to get green component. 
        // imgData.data[ind * 4 + 1] = 32;
        // imgData.data[ind * 4 + 0] =  ((grid[imgDataMap[ind]] == 0)? 0 : 255);
        // imgData.data[ind * 4 + 1] =  ((grid[imgDataMap[ind]] == 1)? 0 : 255);
        // imgData.data[ind * 4 + 2] =  ((grid[imgDataMap[ind]] == 2)? 0 : 255);

        imgData.data[ind * 4 + 0] =  ((grid[imgDataMap[ind]] == 0)? 0 : 255-(grid_iter[imgDataMap[ind]]));
        imgData.data[ind * 4 + 1] =  ((grid[imgDataMap[ind]] == 1)? 0 : 255-(grid_iter[imgDataMap[ind]]));
        imgData.data[ind * 4 + 2] =  ((grid[imgDataMap[ind]] == 2)? 0 : 255-(grid_iter[imgDataMap[ind]]));

    }
    ctx.putImageData(imgData,0,0)
}

var lastframems = 0;
var dtSinceLastFrame = 0;
var frameRate = 1;
var MAX_GRID_PER_FRAME = L*L/100;
function main_loop(timestamp){
    var dt = timestamp - lastframems;
    dtSinceLastFrame += dt; 
    // console.log(dt)
    if(dtSinceLastFrame > 0.){ //can reduce frame rate using this. 
        for (var i = 0; i < MAX_GRID_PER_FRAME; i++){ tick(); if(window.performance.now() - timestamp > 1000/61) break; }
        draw();
        dtSinceLastFrame = 0;
    }
    lastframems = timestamp;
    requestAnimationFrame(main_loop);
}



function init(){
    //fixing the alpha values in the img_buffer.
    for(var i =3 ; i < Npix*4; i+=4){ imgData.data[i] = 255; }
    //making it green!
    for(var i =1 ; i < Npix*4; i+=4){ imgData.data[i] = 128; }
    //building the image data map: 
    var ind = 0;
    for(var y=0; y < L*scale_factor; y++){
        for(var x = 0; x < L * scale_factor; x++){
            ind = y*L*scale_factor + x;
            imgDataMap[ind] = grid_xy_to_ind(Math.floor(x / scale_factor),Math.floor(y/scale_factor));
        }
    }

    requestAnimationFrame(main_loop);
}


const can = document.getElementById('can');
const ctx = can.getContext("2d");
document.addEventListener('DOMContentLoaded', init, false);
