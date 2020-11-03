/** 
Project logic credits to:
Jos Stam:  http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf (if link is broken search for his academic works, they are brilliant)
Oliver Hunt: https://nerget.com/fluidSim/ (incredible exemplar for the implementation of Jos Stam's work)
Evgeny Demidov: https://www.ibiblio.org/e-notes/webgl/gpu/fluid.htm (exemplar using webgl)
Karl Sims: https://www.karlsims.com/fluid-flow.html (interpretation of Jos Stam's works)
*/

function FluidField(canvas) {
    function setFlag(flagName){
        if(flagName === "free"){
            resetFlags();
            flags[0] = 1;
        } else if(flagName === "heart"){
            resetFlags();
            flags[1] = 1;
        } else if(flagName === "doubleHeart"){
            resetFlags();
            flags[2] = 1;
        } else if (flagName === "opposingSources"){
            resetFlags();
            flags[3] = 1;
        } else if (flagName === "fourWaySources"){
            resetFlags();
            flags[4] = 1;
        }

        function resetFlags(){
            flags = [0, 0, 0, 0, 0];
        }
    }

    function addFields(x, s, dt)
    {
        for (var i=0; i<size ; i++ ) x[i] += dt*s[i];
    }

    function set_bnd(b, x)
    {
        //core function for blank slate and counter negation
        if (flags[0] === 1){
            if (b===1) { //horizontal 
                for (var i = 1; i <= width; i++) {
                    x[i] =  x[i + rowSize]; //top
                    x[i + (height+1) *rowSize] = x[i + height * rowSize]; //bottom
                }
    
                for (var j = 1; i <= height; i++) {
                    x[j * rowSize] = -x[1+j * rowSize]; //left 
                    x[(width + 1) + j * rowSize] = -x[width + j * rowSize]; //right
                }
            } else if (b === 2) { //vertical
                for (var i = 1; i <= width; i++) {
                    x[i] = -x[i + rowSize]; // top
                    x[i + (height + 1) * rowSize] = -x[i + height * rowSize]; //bottom 
                }
    
                for (var j = 1; j <= height; j++) {
                   x[j * rowSize] =  x[1 + j * rowSize]; //left
                    x[(width + 1) + j * rowSize] =  x[width + j * rowSize]; //right
                }
            } else { // called for density and some of project, advect 
                for (var i = 1; i <= width; i++) {
                    x[i] =  x[i + rowSize]; // top (becomes one row below)
                    x[i + (height + 1) * rowSize] = x[i + height * rowSize]; // bottom? 
                } 
                // height is excluding the 2 pixel padding, L R; 
                //rowSize is including 2 pixel padding L R
                //width is excluding 2 pixel padding LR 
    
                for (var j = 1; j <= height; j++) {
                    x[j * rowSize] =  x[1 + j * rowSize]; // left side
                    x[(width + 1) + j * rowSize] =  x[width + j * rowSize]; // right side
                }
            }
        } else if (flags[1] === 1){

            let localX = 20;
            for (var j = 1; j <= height; ++j){
                x[j * rowSize ] = 0; //left
                if ( (j == height/2 - 1)){
                    x[j*rowSize + localX] = 0.5*(x[(j-1)*rowSize + localX] + x[j*rowSize + localX - 1]);

                }
                if((j >= height/2 + 2) && (j <= height/2 + 4)){
                    if (b == 1) {
                        x[j*rowSize + localX] = -1* x[j*rowSize + localX - 1]; // left end of block
                        x[j*rowSize + localX + 1] = x[j*rowSize + localX];
                        x[j*rowSize + localX + 2] = -1 * x[j*rowSize + localX + 3] ;// right end of the block
                    }
                }
                
                x[width + j * rowSize ] = 0; //right
            }
            for (var i = 1; i <= width; ++i){
                x[i ] = 0

                x[i + (height + 1)*rowSize ] = 0
            }
        } else if (flags[2] === 1){ // opposing bifurcation sources 
            let localX = 20;

            //neutralizing velocities
            for (var j = 1; j <= height; ++j){
                x[j * rowSize ] = 0; //left
                if ( (j == height/2 - 1)){
                    x[j*rowSize + localX] = 0.5*(x[(j-1)*rowSize + localX] + x[j*rowSize + localX -1]); //leftSource
                    x[j*rowSize + width - localX] = 0.5*(x[(j-1)*rowSize + width - localX] + x[j*rowSize + width - localX + 1]); //rightSource
                }
                if((j >= height/2 + 2) && (j <= height/2 + 4)){
                    if (b == 1) {
                        //leftsource
                        x[j*rowSize + localX] = -1* x[j*rowSize + localX - 1]; // left end of block
                        x[j*rowSize + localX + 1] = x[j*rowSize + localX];
                        x[j*rowSize + localX + 2] = -1 * x[j*rowSize + localX + 3] ;// right end of the block
                    

                        //rightsource
                        x[j*rowSize + width - localX] = -1* x[j*rowSize + width - localX + 1]; // left end of block
                        x[j*rowSize + width - localX - 1] = x[j*rowSize + width - localX];
                        x[j*rowSize + width - localX - 2] = -1 * x[j*rowSize + width - localX - 3] ;// right end of the block
                    
                    
                    }

                }
                x[width + j * rowSize ] = 0; //right
            }
            for (var i = 1; i <= width; ++i){
                x[i ] = 0
                x[i + (height + 1)*rowSize ] = 0
            }
        } else if (flags[3] === 1 || flags[4] === 1){//opposing dual sources (3) or 4 way source (4)
            
            for (var j = 1; j <= height; ++j){
                x[j * rowSize ] = 0; //left
                x[width + j * rowSize ] = 0; //right
            }
            for (var i = 1; i <= width; ++i){
                x[i ] = 0
                x[i + (height + 1)*rowSize ] = 0
            }
        }

        //corners
        var maxEdge = (height + 1) * rowSize;
        x[0]                 = 0.5 * (x[1] + x[rowSize]);
        x[maxEdge]           = 0.5 * (x[1 + maxEdge] + x[height * rowSize]);
        x[(width+1)]         = 0.5 * (x[width] + x[(width + 1) + rowSize]);
        x[(width+1)+maxEdge] = 0.5 * (x[width + maxEdge] + x[(width + 1) + height * rowSize]);
        
    }

    function lin_solve(b, x, x0, a, c)
    {
        if (a === 0 && c === 1) {
            for (var j=1 ; j<=height; j++) {
                var currentRow = j * rowSize;
                ++currentRow;
                for (var i = 0; i < width; i++) {
                    x[currentRow] = x0[currentRow];
                    ++currentRow;
                }
            }
            set_bnd(b, x);
        } else {
            var invC = 1 / c;
            for (var k=0 ; k<iterations; k++) {
                for (var j=1 ; j<=height; j++) {
                    var lastRow = (j - 1) * rowSize;
                    var currentRow = j * rowSize;
                    var nextRow = (j + 1) * rowSize;
                    var lastX = x[currentRow];
                    ++currentRow;
                    for (var i=1; i<=width; i++)
                        lastX = x[currentRow] = (x0[currentRow] + a*(lastX+x[++currentRow]+x[++lastRow]+x[++nextRow])) * invC;
                }
                set_bnd(b, x);
            }
        }
    }
    
    function diffuse(b, x, x0, dt)
    {
        var a = 0;
        lin_solve(b, x, x0, a, 1 + 4*a);
    }
    
    function lin_solve2(x, x0, y, y0, a, c)
    {
        if (a === 0 && c === 1) {
            for (var j=1 ; j <= height; j++) {
                var currentRow = j * rowSize;
                ++currentRow;
                for (var i = 0; i < width; i++) {
                    x[currentRow] = x0[currentRow];
                    y[currentRow] = y0[currentRow];
                    ++currentRow;
                }
            }
            set_bnd(1, x);
            set_bnd(2, y);
        } else {
            var invC = 1/c;
            for (var k=0 ; k<iterations; k++) {
                for (var j=1 ; j <= height; j++) {
                    var lastRow = (j - 1) * rowSize;
                    var currentRow = j * rowSize;
                    var nextRow = (j + 1) * rowSize;
                    var lastX = x[currentRow];
                    var lastY = y[currentRow];
                    ++currentRow;
                    for (var i = 1; i <= width; i++) {
                        lastX = x[currentRow] = (x0[currentRow] + a * (lastX + x[currentRow] + x[lastRow] + x[nextRow])) * invC;
                        lastY = y[currentRow] = (y0[currentRow] + a * (lastY + y[++currentRow] + y[++lastRow] + y[++nextRow])) * invC;
                    }
                }
                set_bnd(1, x);
                set_bnd(2, y);
            }
        }
    }
    
    function diffuse2(x, x0, y, y0, dt) // only called by velocity components
    {
        var a = 0;
        lin_solve2(x, x0, y, y0, a, 1 + 4 * a);
    }
    
    function advect(b, d, d0, u, v, dt)
    {
        var Wdt0 = dt * width;
        var Hdt0 = dt * height;
        var Wp5 = width + 0.5;
        var Hp5 = height + 0.5;
        for (var j = 1; j<= height; j++) {
            var pos = j * rowSize;
            for (var i = 1; i <= width; i++) {
                var x = i - Wdt0 * u[++pos]; 
                var y = j - Hdt0 * v[pos];
                if (x < 0.5)
                    x = 0.5;
                else if (x > Wp5)
                    x = Wp5;
                var i0 = x | 0;
                var i1 = i0 + 1;
                if (y < 0.5)
                    y = 0.5;
                else if (y > Hp5)
                    y = Hp5;
                var j0 = y | 0;
                var j1 = j0 + 1;
                var s1 = x - i0;
                var s0 = 1 - s1;
                var t1 = y - j0;
                var t0 = 1 - t1;
                var row1 = j0 * rowSize;
                var row2 = j1 * rowSize;
                d[pos] = s0 * (t0 * d0[i0 + row1] + t1 * d0[i0 + row2]) + s1 * (t0 * d0[i1 + row1] + t1 * d0[i1 + row2]);
            }
        }
        set_bnd(b, d);
    }
    
    function project(u, v, p, div)
    {
        var h = -0.5 / Math.sqrt(width * height);
        for (var j = 1 ; j <= height; j++ ) {
            var row = j * rowSize;
            var previousRow = (j - 1) * rowSize;
            var prevValue = row - 1;
            var currentRow = row;
            var nextValue = row + 1;
            var nextRow = (j + 1) * rowSize;
            for (var i = 1; i <= width; i++ ) {
                div[++currentRow] = h * (u[++nextValue] - u[++prevValue] + v[++nextRow] - v[++previousRow]);
                p[currentRow] = 0;
            }
        }
        set_bnd(0, div);
        set_bnd(0, p); 
        
        lin_solve(0, p, div, 1, 4 );
        var wScale = 0.5 * width;
        var hScale = 0.5 * height;
        for (var j = 1; j<= height; j++ ) {
            var prevPos = j * rowSize - 1;
            var currentPos = j * rowSize;
            var nextPos = j * rowSize + 1;
            var prevRow = (j - 1) * rowSize;
            var currentRow = j * rowSize;
            var nextRow = (j + 1) * rowSize;

            for (var i = 1; i<= width; i++) {
                u[++currentPos] -= wScale * (p[++nextPos] - p[++prevPos]);
                v[currentPos]   -= hScale * (p[++nextRow] - p[++prevRow]);
            }
        }
        set_bnd(1, u);
        set_bnd(2, v);
    }
    
    function dens_step(x, x0, u, v, dt)
    {
        addFields(x, x0, dt);
        diffuse(0, x0, x, dt ); // else statement called in set_bnd (b = 0)
        advect(0, x, x0, u, v, dt );
    }
    
    function vel_step(u, v, u0, v0, dt)
    {
        addFields(u, u0, dt );
        addFields(v, v0, dt );
        var temp = u0; u0 = u; u = temp;
        var temp = v0; v0 = v; v = temp;
        diffuse2(u,u0,v,v0, dt);
        project(u, v, u0, v0);
        var temp = u0; u0 = u; u = temp; 
        var temp = v0; v0 = v; v = temp;
        advect(1, u, u0, u0, v0, dt);
        advect(2, v, v0, u0, v0, dt);
        project(u, v, u0, v0 );
    }
    var uiCallback = function(d,u,v) {};

    function Field(dens, u, v) {
        this.setDensity = function(x, y, d) {
             dens[(x + 1) + (y + 1) * rowSize] = d;
        }
        this.getDensity = function(x, y) {
             return dens[(x + 1) + (y + 1) * rowSize];
        }
        this.setVelocity = function(x, y, xv, yv) {
             u[(x + 1) + (y + 1) * rowSize] = xv;
             v[(x + 1) + (y + 1) * rowSize] = yv;
        }
        this.getXVelocity = function(x, y) {
             return u[(x + 1) + (y + 1) * rowSize];
        }
        this.getYVelocity = function(x, y) {
             return v[(x + 1) + (y + 1) * rowSize];
        }
        this.width = function() { return width; }
        this.height = function() { return height; }
        this.getFlagState = function() {
            if(flags[0] === 1){
                return "free";
            } else if (flags[1] === 1) {
                return "heart";
            } else if(flags[2] === 1) {
                return "doubleHeart";
            } else if (flags[3] === 1) {
                return "opposingSources";
            } else if (flags[4] === 1) {
                return "fourWaySources";
            } else {
                console.log("getFlagState failed");
                return false;
            }
        }
    }
    function queryUI(d, u, v)
    {
        for (var i = 0; i < size; i++)
            u[i] = v[i] = d[i] = 0.0;
        uiCallback(new Field(d, u, v));
    }

    this.update = function () {
        queryUI(dens_prev, u_prev, v_prev);
        vel_step(u, v, u_prev, v_prev, dt);
        dens_step(dens, dens_prev, u, v, dt);
        displayFunc(new Field(dens, u, v));
    }
    this.setDisplayFunction = function(func) {
        displayFunc = func;
    }
    
    this.iterations = function() { return iterations; }
    this.setIterations = function(iters) {
        if (iters > 0 && iters <= 100)
           iterations = iters;
    }
    this.setUICallback = function(callback) {
        uiCallback = callback;
    }
    var iterations = 25;
    var dt = 0.1;
    var dens;
    var dens_prev;
    var u;
    var u_prev;
    var v;
    var v_prev;
    var width;
    var height;
    var rowSize;
    var size;
    var displayFunc;
    var flags = [1, 0, 0, 0, 0]; // flag 0 = default free for all; 1= heart 


    function reset(configuration = "undefined")
    {
        if (configuration !== "undefined"){

            if(configuration == "heart"){
                console.log("changed to heart");
                setFlag("heart");
                document.getElementById("stateHeader").innerHTML = "Bifurcating object";
            } else if (configuration === "free"){
                console.log("changed to free");
                setFlag('free');
                document.getElementById("stateHeader").innerHTML = "Free draw canvas";

            } else if (configuration === "doubleHeart"){
                console.log("changed to doubleHeart");
                setFlag('doubleHeart');
                document.getElementById("stateHeader").innerHTML = "Opposing bifurcations";

            } else if (configuration === "opposingSources"){
                console.log("changed to opposingSources");
                setFlag('opposingSources');
                document.getElementById("stateHeader").innerHTML = "Opposing sources";

            } else if (configuration === "fourWaySources"){
                console.log("changed to fourWaySources");
                setFlag("fourWaySources");
                document.getElementById("stateHeader").innerHTML = "Four way sources";

            }

        } else {
            console.log("no change to config");
        }   
            rowSize = width + 2;
            size = (width+2)*(height+2);
            dens = new Array(size);
            dens_prev = new Array(size);
            u = new Array(size);
            u_prev = new Array(size);
            v = new Array(size);
            v_prev = new Array(size);
            for (var i = 0; i < size; i++)
                dens_prev[i] = u_prev[i] = v_prev[i] = dens[i] = u[i] = v[i] = 0;
        
    }
    this.reset = reset;
    this.setResolution = function (hRes, wRes)
    {
        var res = wRes * hRes;
        if (res > 0 && res < 1000000 && (wRes != width || hRes != height)) {
            width = wRes;
            height = hRes;
            reset();
            return true;
        }
        return false;
    }
}

