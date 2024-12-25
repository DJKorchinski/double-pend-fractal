const  L = 400;
const Npix = L*L;
var imgData = new ImageData(L,L); 


var theta = Math.PI * 0.5;
var thetadot = 0;
var thetaddot = 0;
const g = 9.8; R = 50; T_FACTOR = 3;
var theta2 = Math.PI * 0.5;
var theta2dot = 0; 
var theta2ddot = 0;
function tick(dttot) {
    dt = T_FACTOR * dttot / 1000;
    theta = theta + thetadot * dt;
    thetadot = thetadot + thetaddot * dt;
    thetaddot = -Math.sin(theta) * g / R;
    first=false; 
}

const FADE_INTERVAL = 4;
var frameCount = 0;
var lastFrameImage;
function draw() {
    if(frameCount % FADE_INTERVAL == 0){
        for(var ind = 0; ind < Npix; ind++){
            imgData.data[ind * 4 + 0] =  lastFrameImage.data[ind * 4 + 0] - 1
            imgData.data[ind * 4 + 1] =  lastFrameImage.data[ind * 4 + 1] - 1
            imgData.data[ind * 4 + 2] =  lastFrameImage.data[ind * 4 + 2] - 1
        }
    }
    x1 = Math.sin(theta) * R + L/2;
    y1 = Math.cos(theta) * R + L/2;
    imgData.data[(Math.floor(x1) + Math.floor(y1)*L) * 4 + 0] = 255;
    x2 = Math.sin(theta2) * R + x1;
    y2 = Math.cos(theta2) * R + y1;
    imgData.data[(Math.floor(x2) + Math.floor(y2)*L) * 4 + 0] = 255;
    ctx.putImageData(imgData,0,0)
    //capture all of the fade-able pixels in the last frame.
    lastFrameImage = ctx.getImageData(0,0,L,L);
    //Drawing the pendulum bob
    drawBob(L/2, L/2, x1, y1, "rgba(150,30,20,1.0)");
    drawBob(x1,y1,x2,y2,"rgba(150,30,20,1.0)");
    frameCount++;
}

function drawBob(x0,y0,x1,y1,style){
    ctx.fillStyle = style;
    ctx.strokeStyle = style;
    ctx.beginPath();
    ctx.moveTo(x0, y0);
    ctx.lineTo(x1, y1);
    ctx.stroke();
    ctx.beginPath();
    ctx.arc(x1,y1,5,0,Math.PI*2);
    ctx.fill();
}

var lastframems = 0;
var dtSinceLastFrame = 0;
var frameRate = 1;
var MAX_GRID_PER_FRAME = L*L/100;
function main_loop(timestamp){
    var dt = timestamp - lastframems;
    tick(dt);
    dtSinceLastFrame += dt; 
    // console.log(dt)
    if(dtSinceLastFrame > 0.){ //can reduce frame rate using this. 
        draw(dtSinceLastFrame);
        dtSinceLastFrame = 0;
    }
    lastframems = timestamp;
    requestAnimationFrame(main_loop);
}



function init(){
    //fixing the alpha values in the img_buffer.
    //making it green!
    initFrame = ctx.getImageData(0,0,L,L);
    for(var i =3 ; i < Npix*4; i+=4){ initFrame.data[i] = 255; imgData.data[i] = 255; } //maxing the alpha 
    for(var i =1 ; i < Npix*4; i+=4){ initFrame.data[i] = 50; } //setting it green
    ctx.putImageData(initFrame,0,0);
    lastFrameImage = ctx.getImageData(0,0,L,L);

    requestAnimationFrame(main_loop);
}


const can = document.getElementById('can');
const ctx = can.getContext("2d");
document.addEventListener('DOMContentLoaded', init, false);
