const  L = 512;
const Npix = L*L;
var imgData = new ImageData(L,L); 
var fractalData = new ImageData(L,L);

var NUM_PENDULA_SUBDIVISIONS = 3;
const g = 9.8; const R = 50; const T_FACTOR = 100; const RSQ = R*R; 
const TMAX = 60; 
var theta1s = [Math.PI * 0.5,Math.PI * 0.5 ]; var theta2s = [Math.PI * 0.5,Math.PI * 0.495 ]; var theta1dots = [0.0, 0.]; var theta2dots = [0.,0.];
var totalFlips = [0,0]; var timeToFlip = [0,0]; 
var t = 0;
var firstFlipIndsThisTick = [];
var indsToRedraw = [];
function tick(dttot) {
    dttot = Math.min(dttot, 1000); //cap the dt to no lower than 1 fps
    dttick = T_FACTOR * dttot / 1000;
    dt = Math.min(0.0025,dttick)
    firstFlipIndsThisTick = [];
    for (var i = 0; i < theta1s.length; i++){
        theta1 = theta1s[i]; theta2=theta2s[i]; theta1dot = theta1dots[i]; theta2dot = theta2dots[i];
        dttot = 0.0;
        while(dttot < dttick){
            dt = Math.min(0.0025,dttick - dttot);
            dttot += dt;
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
        //checking for flips of either bob
        if(Math.abs(theta1) > 2*Math.PI){
            if(totalFlips[i] == 0) {timeToFlip[i] = t; firstFlipIndsThisTick.push(i);}
            totalFlips[i] += 1;
            theta1 -= 2*Math.PI*Math.sign(theta1);
        }
        if(Math.abs(theta2) > 2*Math.PI){
            if(totalFlips[i] == 0) {timeToFlip[i] = t; firstFlipIndsThisTick.push(i);}
            totalFlips[i] += 1;
            theta2 -= 2*Math.PI*Math.sign(theta2);
        }
        theta1s[i] = theta1; theta2s[i] = theta2; theta1dots[i] = theta1dot; theta2dots[i] = theta2dot;
    }
    t+=dttick; 
    // console.log(t,KE(theta1,theta1dot,theta2,theta2dot) + PE(theta1,theta1dot,theta2,theta2dot));
    //Now, we need to process the flipped pendula, and update the grid.
    for (var i = 0; i < firstFlipIndsThisTick.length; i++){
        pendula_ind = firstFlipIndsThisTick[i];
        grid_ind = grid_inds_being_computed[pendula_ind];
        current_grid.mspoints[grid_ind].computed = true;
        current_grid.mspoints[grid_ind].tfirst = timeToFlip[pendula_ind];
        indsToRedraw.push(grid_ind);
    }

    //Now, check to see  if we have exceeded our simulation time!
    if(t > TMAX){
        for(var i = 0; i < inds_to_compute.length; i++){
            grid_ind = inds_to_compute[i];
            if(current_grid.mspoints[grid_ind].computed == false){
                current_grid.mspoints[grid_ind].computed = true;
                indsToRedraw.push(grid_ind);
            }
        }
        should_reinit = true;
    }
}

//Computing the derivatives for a double pendulum: 
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
//Computing the kinetic energy
function KE(theta1,thetadot1,theta2,thetadot2){
    return 0.5 * R*R* (2*thetadot1*thetadot1 + thetadot2*thetadot2 + 2*thetadot1*thetadot2*Math.cos(theta1 - theta2));
}
//Computing the potential energy
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
    for (var i = 0; i < theta1s.length; i++){
        theta1 = theta1s[i]; theta2 = theta2s[i];
        x1 = Math.sin(theta1) * R + L/2;
        y1 = Math.cos(theta1) * R + L/2;
        imgData.data[(Math.floor(x1) + Math.floor(y1)*L) * 4 + 0] = 255;
        x2 = Math.sin(theta2) * R + x1;
        y2 = Math.cos(theta2) * R + y1;
        imgData.data[(Math.floor(x2) + Math.floor(y2)*L) * 4 + 0] = 255;
    }
    ctx.putImageData(imgData,0,0)
    //capture all of the fade-able pixels in the last frame.
    lastFrameImage = ctx.getImageData(0,0,L,L);
    //Drawing the pendulum bobs
    for (var i = 0; i < theta1s.length; i++){
        theta1 = theta1s[i]; theta2 = theta2s[i];
        x1 = Math.sin(theta1) * R + L/2;
        y1 = Math.cos(theta1) * R + L/2;
        x2 = Math.sin(theta2) * R + x1;
        y2 = Math.cos(theta2) * R + y1;
        drawBob(L/2, L/2, x1, y1, "rgba(150,30,20,1.0)");
        drawBob(x1,y1,x2,y2,"rgba(150,30,20,1.0)");
    }
    //Now, we need to draw the fractal! 
    ctx.putImageData(fractalData,0,L);
    for (var i = 0; i < indsToRedraw.length; i++){
        let ind = indsToRedraw[i];
        let x = current_grid.mspoints[ind].pixel_extent.sx;
        let y = current_grid.mspoints[ind].pixel_extent.sy + L;
        let w = current_grid.mspoints[ind].pixel_extent.w;
        let h = current_grid.mspoints[ind].pixel_extent.h;
        let tFlip = current_grid.mspoints[ind].tfirst;
        ctx.fillStyle = `rgba(0,${Math.min(Math.floor(255*tFlip/TMAX),255)},0,1.0)`;
        ctx.fillRect(x,y,w,h);
    }
    if(indsToRedraw.length > 0){
        fractalData = ctx.getImageData(0,L,L,L);
        indsToRedraw = [];
    }
    //Let's also draw the location of the pendula being computed!
    for (var i = 0; i < grid_inds_being_computed.length; i ++){
        let ind = grid_inds_being_computed[i];
        let x = current_grid.mspoints[ind].pixel_extent.sx;
        let y = current_grid.mspoints[ind].pixel_extent.sy + L;
        let w = current_grid.mspoints[ind].pixel_extent.w;
        let h = current_grid.mspoints[ind].pixel_extent.h;
        x+=w/2; y+=h/2;
        ctx.fillStyle = `rgba(0,0,255,1.0)`;
        ctx.beginPath();
        ctx.arc(x,y,5,0,Math.PI*2);
        ctx.fill();
        ctx.closePath();
    }
    frameCount++;
}

//A function to draw a pendulum bob.
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
    ctx.closePath();
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
    if(dtSinceLastFrame > 0.){ //can reduce frame rate using this, by selecting a value higher than 0.0 (in ms).
        draw(dtSinceLastFrame);
        dtSinceLastFrame = 0;
    }
    if(should_reinit && indsToRedraw.length == 0){
        should_reinit = false;
        needNewGrid = !initPendula(NUM_PENDULA_SUBDIVISIONS);
        if(needNewGrid){
            current_grid = setupGrid(current_grid.subdivisions + 1);
            initPendula(NUM_PENDULA_SUBDIVISIONS);
        }
    }

    lastframems = timestamp;
    requestAnimationFrame(main_loop);
}

//fractal sampling grid set up
const MIN_THETA = -Math.PI; const MAX_THETA = Math.PI;
function setupGrid(numDivisions) {
    let pixels_per_division = L/2**numDivisions;
    var grid = {subdivisions: numDivisions, axissamples:Math.pow(2,numDivisions), mspoints: []}
    for (var i = 0; i < grid.axissamples; i++){
        ms_x = (i+0.5)/(grid.axissamples) * (MAX_THETA - MIN_THETA) + MIN_THETA;
        for (var j = 0; j < grid.axissamples; j++){
            ms_y = (j+0.5)/(grid.axissamples) * (MAX_THETA - MIN_THETA) + MIN_THETA;
            //okay, let's map each ms_grid to the pixel grid.
            extent = { sx : pixels_per_division*i, sy : pixels_per_division*j, w : pixels_per_division, h : pixels_per_division };
            grid.mspoints.push({x: ms_x, y: ms_y,computed: false, tfirst: -1, nflips: -1, pixel_extent: extent});
        }
    }
    return grid;
}

var grid_inds_being_computed = [];
var should_reinit = false;
function initPendula(numDivisions){
    theta1s = []; theta2s = []; theta1dots = []; theta2dots = [];
    
    //okay, let's iterate along until we find a cell that has not been computed.
    inds_to_compute = [];
    found_uncomputed = false;
    for (var i = 0; i < current_grid.mspoints.length; i++){
        if(current_grid.mspoints[i].computed == false){
            for (var j = 0; j < 2**numDivisions; j++){
                for(var k = 0; k < 2**numDivisions; k++){
                    inds_to_compute.push(i+j+current_grid.axissamples*k);
                }
            }
            found_uncomputed = true;
            break;
        }
    }
    if(!found_uncomputed){
        return false;
    }
    grid_inds_being_computed = inds_to_compute;
    console.log(inds_to_compute);
    
    totalFlips = []; timeToFlip = [];
    for (var i = 0; i < inds_to_compute.length; i++){
        ind = inds_to_compute[i];
        theta1s.push(current_grid.mspoints[ind].x); theta2s.push(current_grid.mspoints[ind].y);
        theta1dots.push(0.0); theta2dots.push(0.0);
        totalFlips.push(0); timeToFlip.push(0.);
    }

    t = 0;
    return true;
}

var current_grid = setupGrid(Math.max(2,NUM_PENDULA_SUBDIVISIONS));
function init(){
    //fixing the alpha values in the img_buffer.
    //making it green!
    initFrame = ctx.getImageData(0,0,L,L);
    for(var i =3 ; i < Npix*4; i+=4){ initFrame.data[i] = 255; imgData.data[i] = 255; } //maxing the alpha 
    for(var i =1 ; i < Npix*4; i+=4){ initFrame.data[i] = 50; } //setting it green
    ctx.putImageData(initFrame,0,0);
    lastFrameImage = ctx.getImageData(0,0,L,L);
    fractalData = ctx.getImageData(0,L,L,L);
    for(var i =3 ; i < Npix*4; i+=4){ fractalData.data[i] = 255; } //maxing the alpha 
    initPendula(NUM_PENDULA_SUBDIVISIONS);
    requestAnimationFrame(main_loop);
}


const can = document.getElementById('can');
const ctx = can.getContext("2d");
document.addEventListener('DOMContentLoaded', init, false);
