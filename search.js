//var text1, text2;

var EPSILON = 0.01;

function assert(condition, errorMsg) {
    if (! condition) throw new Error(errorMsg);
}

function zeroedArray(n) {
    var a = new Array(n);
    for (var i = 0; i < n; ++i) a[i] = 0;
    return a;
}

function factorial(n) {
    var r = 1;
    for (var i = n; i > 1; --i)
        r *= i;
    return r;
}

function gaussCumulative(parms, x) {
    return jstat.pnorm(x, parms.mean, parms.sigma);
}
function gaussCumulativeBetween(parms, x, y) {
    assert(x < y, "Bad arguments to 'gaussCumulativeBetween'");
    var r = gaussCumulative(parms, y) - gaussCumulative(parms, x);
    assert(r >= 0 && r <= 1, "Bad probability as return value of 'gaussCumulativeBetween': " + r);
    return r;
}

function expCumulative(lambda, x) {
    if (x < 0)
        return 0;
    else
        return 1-Math.exp(-1*lambda*x);
}

function getWhiteNoiseParms() {
    var parms = { };
    parms.f = 1;
    parms.m = 20;
    // Fig 17b:
//    parms.phi = [0.19, 0.91, 0.37, 0.82, 0.42, 0.31, 0.54, 0.53, 0.27, 0.73, 0.27, 0.10, 0.59, 0.75, 0.65, 0.50, 0.78, 0.61, 0.46, 0.55];
    // Figs 15 etc.
//    parms.phi = [0.46, 0.52, 0.03, 0.56, 0.19, 0.03, 0.40, 0.02, 0.53, 0.73, 0.22, 0.39, 0.38, 0.14, 0.86, 0.17, 0.25, 0.27, 0.66, 0.51];
    parms.phi = new Array(parms.m);
    for (var i = 0; i < parms.m; ++i) {
        parms.phi[i] = Math.random(); // Between 0 and 1.
        assert(parms.phi[i] >= 0 && parms.phi[i] <= 1, "Bad random number");
    }
    return parms;
}

function whiteNoise(parms, p) {
    var sum = 0;
    for (var i = 1; i <= parms.m; ++i) {
        sum += Math.sin((2*Math.PI*i*parms.f*p)/500 + parms.phi[i-1]);
    }
    var r = (1/parms.m) * sum;
    assert(r >= -1 && r <= 1, "Bad return value for 'whiteNoise': " + r);
    return r;
}

function getSearchParms() {
    var numRadii = 10, radiusSize = 40;
    return {
        maxSimulationSteps: 1000,
        gaussParms: { mean: 0, sigma: Math.sqrt(numRadii) },
        whiteNoiseParms: getWhiteNoiseParms(),
        pDiscover: 0.7,
        radiusSize: radiusSize,
        stepSize: radiusSize / 4,
        stepTime: 20, // Milliseconds
        numRadii: numRadii,
        scalingParm: 10,
        defaultTanRat: 0.9,
        defaultRadRat: 0.1
    }
}

function getInitialSearchState() {
    var parms = getSearchParms();
    return {
        parms: parms,
        step: 0,
        timeInRadii: zeroedArray(parms.numRadii),
        radiiMinXs: zeroedArray(parms.numRadii),
        radiiMaxXs: zeroedArray(parms.numRadii),
        radiiMinYs: zeroedArray(parms.numRadii),
        radiiMaxYs: zeroedArray(parms.numRadii),
        inOrOut: 1, // 1 for out, 0 for stay in the same radius, -1 for in.
        posX: 0, // From origin (not the same as processing coord).
        posY: 0, // From origin (not the same as processing coord).
        oldPosX: 0,
        oldPosY: 0,
        currentRadius: 0,
        previousRadius: null,
        excursionNumber: 1,
        ticksSinceLastIncrementOfExcursionNumber: 0
    };
}

var pInRadii;
(function () {
    var memoized;
    var memoizedForStep = null;
    pInRadii = function (state) {
        if (state.step === memoizedForStep) {
            return memoized;
        }
        else {
            var fracs = new Array(state.parms.numRadii);
            var withDiscounts = new Array(state.parms.numRadii);
            var fracsTotal = 0;
            var discountTotal = 0;
            var priorTotal = 0;
            for (var i = 0; i < state.parms.numRadii; ++i) {
                var circ = (i+1) * state.parms.radiusSize * 2 * Math.PI;
                var fractionInspectedEstimate = (state.timeInRadii[i] * state.parms.stepSize) / circ;
                fractionInspectedEstimate = expCumulative(-Math.log(1-state.parms.pDiscover), fractionInspectedEstimate);
                fracs[i] = 1-fractionInspectedEstimate;
                fracsTotal += fracs[i];
                
                // Multiplying by two because we're only using one tail of the distribution.
                var prior = gaussCumulativeBetween(state.parms.gaussParms, i, i+1) * 2;
                priorTotal += prior;
                withDiscounts[i] = prior * (1 - fractionInspectedEstimate);
                discountTotal += prior - withDiscounts[i];
            }
            assert(Math.abs(priorTotal-1) < EPSILON, "Bad total probability (X)");
            
            if (fracsTotal != 0) {
                for (var i = 0; i < withDiscounts.length; ++i) {
                    var extra = (fracs[i] / fracsTotal) * discountTotal;
                    withDiscounts[i] += extra;
                }
            }

            memoizedForStep = state.step;
            memoized = withDiscounts;
            return withDiscounts;
        }
    }
})();

function pInRadius(state, radius) { return pInRadii(state)[radius]; }
function mostProbableRadius(state) {
    var ps = pInRadii(state);
    var max = 0;
    var maxr = 0;
    for (var i = 0; i < ps.length; ++i) {
        if (ps[i] > max) {
            max = ps[i];
            maxr = i;
        }
    }
    return maxr;
}

// Given the ant's current distance from the origin,
// this function decides whether the next step in the search
// should move toward or away from the origin. It returns 1
// to indicate 'away', 0 to indicate 'stay at current distance',
// and -1 to indicate 'toward'.
function inOrOut(state) {
    var maxr = mostProbableRadius(state);

    if (state.currentRadius > maxr)
        return -1;
    else if (state.currentRadius < maxr)
        return 1;
    else
        return 0;
}

function assertConsistentRadiusPs(state) {
    var tot = 0;
    for (var i = 0; i < state.parms.numRadii; ++i) {
        tot += pInRadius(state, i);
    }
    assert(Math.abs(tot-1.0) < EPSILON, "Bad total probability: " + tot);
}

function updateSearchState(state) {
    if (state.step > state.parms.maxSimulationSteps)
        return false;
    state.step++;

    assertConsistentRadiusPs(state);

    var tanRadRat = whiteNoise(state.parms.whiteNoiseParms, state.step) * (state.parms.scalingParm/state.excursionNumber) * (state.currentRadius+1);
    var trrs = tanRadRat * tanRadRat;
    var rad = 1/Math.sqrt(1+trrs);
    var tan = tanRadRat/Math.sqrt(1+trrs);

    var vecFromRadX = state.posX == 0 ? 1 : state.posX;
    var vecFromRadY = state.posY;
    var l = Math.sqrt((vecFromRadX*vecFromRadX) + (vecFromRadY*vecFromRadY));
    state.currentRadius = Math.floor(l / state.parms.radiusSize);
    if (state.currentRadius >= state.parms.numRadii)
        return false;

    state.inOrOut = inOrOut(state);
    if (state.inOrOut == 0) {
        rad = 0;//state.parms.defaultRadRat;
        tan = 1;//state.parms.defaultTanRat;
//        var d = Math.sqrt(state.posX*state.posX + state.posY*state.posY);
//        var middle = state.currentRadius * state.parms.radiusSize + state.parms.radiusSize*0.5;
//        if (d <= middle)
//            state.inOrOut = -1;
//        else
//            state.inOrOut = 1;
    }

    state.oldPosX = state.posX;
    state.oldPosY = state.posY;

    var sin = (state.inOrOut == -1 ? vecFromRadY * -1 : vecFromRadY) / l;
    var cos = (state.inOrOut == -1 ? vecFromRadX * -1 : vecFromRadX) / l;
    state.posX += rad * state.parms.stepSize * cos;
    state.posY += rad * state.parms.stepSize * sin;
        
    var perpX = -vecFromRadY;
    var perpY = vecFromRadX;
        
    var sin = perpY / l;
    var cos = perpX / l;
    state.posX += tan * cos * state.parms.stepSize;
    state.posY += tan * sin * state.parms.stepSize;

    // Check that the ant hasn't gone outside the area,
    // and move it back inside if it has.
    var dsqd = state.posX*state.posX + state.posY*state.posY;
    var radsqd = state.parms.numRadii * state.parms.radiusSize * state.parms.numRadii * state.parms.radiusSize;
    if (dsqd > radsqd) {
        var d = Math.sqrt(dsqd);
        var sin = state.posY / d;
        var cos = state.posX / d;
        state.posX = state.parms.numRadii * state.parms.radiusSize * cos;
        state.posY = state.parms.numRadii * state.parms.radiusSize * sin;
    }

    if (state.posX < state.radiiMinXs[state.currentRadius])
        state.radiiMinXs[state.currentRadius] = state.posX;
    else if (state.posX > state.radiiMaxXs[state.currentRadius])
        state.radiiMaxXs[state.currentRadius] = state.posX;

    if (state.posY < state.radiiMinYs[state.currentRadius])
        state.radiiMinYs[state.currentRadius] = state.posY;
    else if (state.posY > state.radiiMaxYs[state.currentRadius])
        state.radiiMaxYs[state.currentRadius] = state.posY;

    if (state.previousRadius === null) {
        state.timeInRadii[state.currentRadius] *= tan;
    }
    else {
        for (var i = state.previousRadius; i <= state.currentRadius; ++i) {
            state.timeInRadii[i] += 1/(Math.abs(state.currentRadius-state.previousRadius)+1);
        }
    }

    if (state.previousRadius > 0 && state.currentRadius == 0) {
        if (state.ticksSinceLastIncrementOfExcursionNumber > 10) {
            state.excursionNumber++;
            state.ticksSinceLastIncrementOfExcursionNumber = -1;
        }
    }
    state.ticksSinceLastIncrementOfExcursionNumber++;

    state.previousRadius = state.currentRadius;

    if (((state.oldPosX < 0 && state.posX > 0) || (state.posX < 0 && state.oldPosX > 0) || (state.posX == 0 && state.oldPosX == 0)) &&
        ((state.oldPosY < 0 && state.posY > 0) || (state.posY < 0 && state.oldPosY > 0) || (state.posY == 0 && state.oldPosY == 0))) {
        state.inOrOut *= -1;
    }

    return true;
}

var SINCOS45 = 0.707106781467586;
var SIN65 = Math.sin((65/360)*2*Math.PI);
var COS65 = Math.cos((65/360)*2*Math.PI);
var SIN35 = Math.sin((25/360)*2*Math.PI);
var COS35 = Math.cos((25/360)*2*Math.PI);

function sketchProc(p) {
    var state = getInitialSearchState();
    var positions = new Array(state.parms.maxSimulationSteps * 2);
    positions[0] = 0;
    positions[1] = 0;

    var WIDTH = 820, HEIGHT = 820;

    p.size(WIDTH, HEIGHT);

    var intervalId = null;
    var intermediateX = state.posX;
    var intermediateY = state.posY;
    var INTERVAL = 10;
    assert(INTERVAL <= state.parms.stepTime, "Bad INTERVAL");
    var ticks = 0;
    function draw(antX, antY, targetX, targetY) {
        // Draw black square as background.
        p.fill(0, 0, 0);
        p.stroke(0, 0, 0);
        p.rect(0, 0, WIDTH, HEIGHT);

        // Draw probability circles.
        var psInRadii = new Array(state.parms.numRadii);
        for (var i = 0; i < state.parms.numRadii; ++i) {
            psInRadii[i] = pInRadius(state, i);            
        }
        var maxr = mostProbableRadius(state);

        for (var i = state.parms.numRadii-1; i >= 0; --i) {
            var prob = psInRadii[i];
            var v = expCumulative(5, prob)*255;
            p.stroke(255, 255, 255);
            if (i == maxr)
                p.fill(218, 165, 32);
            else
                p.fill(v,v,v);
            p.ellipse(WIDTH/2, HEIGHT/2, (i+1) * state.parms.radiusSize * 2, (i+1) * state.parms.radiusSize * 2);
        }

        for (var i = 1; i < state.parms.numRadii; ++i) {
            // Arrow to indicate whether this is more or less probable than inner circle.
            var prob = psInRadii[i];
            var innerProb = psInRadii[i-1];
            p.stroke(0,255,0);
            p.fill(0,255,0);
            var ms = [ [1, 1, 1, 1,     -1, -1, 1, 1, 1, 1,   1, 1, -1, -1, -1, -1 ],
                       [1, 1, -1, 1,    -1, 1, 1, -1, 1, -1,  1, -1, -1, 1, -1, 1  ],
                       [-1, 1, -1, 1,   1, 1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1   ],
                       [-1, 1, 1, 1,   1, -1, -1, 1, -1, 1,  -1, 1, 1, -1, 1, -1   ] ];
            for (var j = 0; j < ms.length; ++j) {
                var m = ms[j];

                var arrowPosX = (i * state.parms.radiusSize * SINCOS45 * m[0]) + (WIDTH/2)*m[1];
                var arrowPosY = (i * state.parms.radiusSize * SINCOS45 * m[2]) + (HEIGHT/2)*m[3];

                if (innerProb < prob) {
                    // Arrow pointing inward.
                    var pointPosX = arrowPosX + 5*m[4];
                    var pointPosY = arrowPosY + 5*m[5];
                    var p1X = pointPosX + COS65*10*m[6];
                    var p1Y = pointPosY + SIN65*10*m[7];
                    var p2X = pointPosX + COS35*10*m[8];
                    var p2Y = pointPosY + SIN35*10*m[9];
                    p.triangle(pointPosX, pointPosY, p1X, p1Y, p2X, p2Y);
                }
                else if (innerProb > prob) {
                    // Arrow pointing outward.
                    var pointPosX = arrowPosX + 5*m[10];
                    var pointPosY = arrowPosY + 5*m[11];
                    var p1X = pointPosX + COS65*10*m[12];
                    var p1Y = pointPosY + SIN65*10*m[13];
                    var p2X = pointPosX + COS35*10*m[14];
                    var p2Y = pointPosY + SIN35*10*m[15];
                    p.triangle(pointPosX, pointPosY, p1X, p1Y, p2X, p2Y);
                }
                else {
                    // Circle to indicate that they're equal.
                    p.ellipse(arrowPosX, arrowPosY, 10, 10);
                }
            }
        }

        var PRECISION = 2;
        for (var i = 1; i < state.parms.numRadii; ++i) {
            var sprob = psInRadii[i].toFixed(PRECISION); // Not perfect rounding (uses binary rep), but good enough.
            p.fill(0, 255, 0);
            // TODO: Dimensions should be computed properly from paramaters, this is very
            // fragile at present.
            p.text(sprob, WIDTH/2 - 10, HEIGHT/2 - (i * state.parms.radiusSize) - state.parms.radiusSize*0.4);
            p.text(sprob, WIDTH/2 - 10, HEIGHT/2 + (i * state.parms.radiusSize) + state.parms.radiusSize*0.4);
            p.text(sprob, WIDTH/2 - 30 - (i * state.parms.radiusSize), HEIGHT/2);
            p.text(sprob, WIDTH/2 + 10 + (i * state.parms.radiusSize), HEIGHT/2);
        }
        p.text(psInRadii[0].toFixed(PRECISION), WIDTH/2 - 10, HEIGHT/2);

        p.stroke(255, 0, 0);
        p.fill(255, 0, 0);
        p.ellipse(antX + WIDTH/2, HEIGHT/2-antY, 20, 20);
    }
    draw(0, 0);

    if (! updateSearchState(state))
        return;
    function update() {
        if (ticks >= state.parms.stepTime / INTERVAL) {
            positions[state.step * 2] = state.posX;
            positions[state.step * 2 + 1] = state.posY;
            draw(state.posX, state.posY);
            intermediateX = state.posX;
            intermediateY = state.posY;
            if (! updateSearchState(state)) {
                clearInterval(intervalId);
                intervalId = "STOPPED";
                // Draw the entire path.
                p.stroke(255, 0, 0);
                p.noFill();
                p.strokeWeight(4);
                p.beginShape();
                for (var i = 0; i < positions.length; i += 2) {
                    var x = positions[i]; var y = positions[i+1];
                    p.vertex(x+WIDTH/2, HEIGHT/2-y);
                }
                p.endShape();
            }
            ticks = 0;
        }
        else {
            var frac = INTERVAL / state.parms.stepTime;
            var xd = state.posX - intermediateX;
            var yd = state.posY - intermediateY;
            assert(!(xd == 0 && yd == 0), "Unexpected zero value (1)");
            var l = Math.sqrt(xd*xd + yd*yd);
            assert(l > 0, "Unexpected zero value (2)");
            var sin = yd/l;
            var cos = xd/l;
            var distanceMoving = frac * state.parms.stepSize;
            intermediateX += cos * distanceMoving;
            intermediateY += sin * distanceMoving;
            draw(intermediateX, intermediateY, state.posX, state.posY);
            ticks++;
        }

        $("#inorout").text(state.inOrOut);
        $("#xpos").text(state.posX.toFixed(2));
        $("#ypos").text(state.posY.toFixed(2));
    }
    intervalId = setInterval(update, INTERVAL);

    $("#pause").click(function (e) {
        e.preventDefault();
        if (intervalId != "STOPPED") {
            if (intervalId == "PAUSED") {
                intervalId = setInterval(update, INTERVAL);
                $("#pause").val('Pause');
            }
            else {
                clearInterval(intervalId);
                intervalId = "PAUSED";
                $("#pause").val('Go');
            }
        }
    });
}

$(document).ready(function () {
    var canvas = document.getElementById("canvas1");
//    text1 = document.getElementById("text1");
//    text2 = document.getElementById("text2");

    var processingInstance = new Processing(canvas, sketchProc);
});
