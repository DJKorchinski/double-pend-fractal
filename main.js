const  L = 400;
const Npix = L*L;
var imgData = new ImageData(L,L); 

const g = 9.8; const R = 50; const T_FACTOR = 3; const RSQ = R*R; 
var theta1 = Math.PI * 0.5;
var theta1dot = 0;
var theta2 = Math.PI * 0.5;
var theta2dot = 0;
var t = 0;
function tick(dttot) {
    dttick = T_FACTOR * dttot / 1000;
    dt = Math.min(0.0025,dttick)
    dttot = 0.0;
    while(dttot < dttick){
        dt = Math.min(0.0025,dttick - dttot);
        dttot += dt;
    // const NSTEP = 10; 
    // for (var i = 0; i < NSTEP; i++){
        //Runge-Kutta integration
        //BUTCHER TABLEAU
        // 0 | 0
        // 1/3 | 1/3
        // 2/3 | -1/3 1
        // 1 | 1 -1 1
        // ----------------
        //   1/8 3/8 3/8 1/8
        const [k11d, k11dd, k12d, k12dd] = derivs(theta1,theta1dot,theta2,theta2dot,R);
        const [k21d, k21dd, k22d, k22dd] = derivs(theta1 + 1./3 * dt * k11d, theta1dot + 1./3 * dt * k11dd, theta2 + 1./3 * dt * k12d, theta2dot + 1./3 * dt * k12dd, R);
        const [k31d, k31dd, k32d, k32dd] = derivs(theta1 - 1./3 * dt * k11d + dt * k21d, theta1dot - 1./3 * dt * k11dd + dt * k21dd, theta2 - 1./3 * dt * k12d + dt * k22d, theta2dot - 1./3 * dt * k12dd + dt * k22dd, R);
        const [k41d, k41dd, k42d, k42dd] = derivs(theta1 + dt * ( k11d - k21d + k31d), theta1dot + dt * ( k11dd - k21dd + k31dd), theta2 + dt * ( k12d - k22d + k32d), theta2dot + dt * ( k12dd - k22dd + k32dd), R);
        
        theta1 = theta1 + dt * (1./8 * k11d + 3./8 * k21d + 3./8 * k31d + 1./8 * k41d);
        theta1dot = theta1dot + dt * (1./8 * k11dd + 3./8 * k21dd + 3./8 * k31dd + 1./8 * k41dd);
        theta2 = theta2 + dt * (1./8 * k12d + 3./8 * k22d + 3./8 * k32d + 1./8 * k42d);
        theta2dot = theta2dot + dt * (1./8 * k12dd + 3./8 * k22dd + 3./8 * k32dd + 1./8 * k42dd);
    }
    t+=dttick; 
    console.log(t,KE(theta1,theta1dot,theta2,theta2dot) + PE(theta1,theta1dot,theta2,theta2dot));
}

function derivs(theta1,thetadot1,theta2,thetadot2,r){
    var thetadotdot1 = -((-2*g*Math.sin(theta1) + 
        g*Math.cos(theta1 - theta2) * Math.sin(theta2) - 
        r*Math.cos(theta1 - theta2)*Math.sin(theta1 - theta2) * thetadot1*thetadot1 
        - r*Math.sin(theta1 - theta2)*thetadot2*thetadot2 )/(
       r*(-2 + Math.cos(theta1 - theta2)*Math.cos(theta1 - theta2) )));
    var thetadotdot2 = -((2*g*Math.cos(theta1 - theta2)*Math.sin(theta1) - 
        2 *g* Math.sin(theta2)+ 
        2 *r* Math.sin(theta1 - theta2) * thetadot1*thetadot1 +
         r*Math.cos(theta1 - theta2) * Math.sin(theta1 - theta2)* thetadot2*thetadot2)/(
       r * (-2 + Math.cos(theta1 - theta2)*Math.cos(theta1 - theta2))));
    return [thetadot1,thetadotdot1,thetadot2,thetadotdot2] //returning the derivatives
}

function KE(theta1,thetadot1,theta2,thetadot2){
    return 0.5 * R*R* (2*thetadot1*thetadot1 + thetadot2*thetadot2 + 2*thetadot1*thetadot2*Math.cos(theta1 - theta2));
}

function PE(theta1,thetadot1,theta2,thetadot2){
    return -g * R * (2*Math.cos(theta1) + Math.cos(theta2));
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
    x1 = Math.sin(theta1) * R + L/2;
    y1 = Math.cos(theta1) * R + L/2;
    imgData.data[(Math.floor(x1) + Math.floor(y1)*L) * 4 + 0] = 255;
    x2 = Math.sin(theta2) * R + x1;
    y2 = Math.cos(theta2) * R + y1;
    imgData.data[(Math.floor(x2) + Math.floor(y2)*L) * 4 + 0] = 255;
    ctx.putImageData(imgData,0,0)
    //capture all of the fade-able pixels in the last frame.
    lastFrameImage = ctx.getImageData(0,0,L,L);
    //Drawing the pendulum bobs
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
