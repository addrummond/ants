//var text1, text2;

var EPSILON = 0.01;

function assert(condition, errorMsg) {
    if (! condition) throw new Error(errorMsg);
}

function mysum(fromInclusive, toInclusive, f) {
    var total = 0;
    for (var i = fromInclusive; i <= toInclusive; ++i) {
        total += f(i);
    }
    return total;
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
    parms.phi = new Array(parms.m);
    for (var i = 0; i < parms.m; ++i) {
        parms.phi[i] = Math.random(); // Between 0 and 1.
    }
    return parms;
}

function whiteNoise(parms, p) {
    var r = (1/parms.m) *
            mysum(1, parms.m, function (i) {
                return Math.sin(2*Math.PI*((i*parms.f*p)+parms.phi[i-1]));
            });
    assert(r >= -1 && r <= 1, "Bad return value for 'whiteNoise': " + r);
    return r;
}

function getSearchParms() {
    var numRadii = 10, radiusSize = 40;
    return {
        maxSimulationSteps: 1000,
        gaussParms: { mean: 0, sigma: Math.sqrt(numRadii) },
        whiteNoiseParms: getWhiteNoiseParms(),
        pDiscover: 0.9,
        radiusSize: radiusSize,
        stepSize: radiusSize/4,
        numRadii: numRadii,
        scalingParm: 10//numRadii * radiusSize * 5
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
        inOrOut: 1, // 1 for out, -1 for in.
        posX: 0, // From origin (not the same as processing coord).
        posY: 0, // From origin (not the same as processing coord).
        oldPosX: 0,
        oldPosY: 0
    };
}

function pInRadius(state, radius) {
    var fracs = new Array(state.parms.numRadii);
    var withDiscounts = new Array(state.parms.numRadii);
    var fracsTotal = 0;
    var discountTotal = 0;
    var priorTotal = 0;
    for (var i = 0; i < state.parms.numRadii; ++i) {
        var circ = (i+1) * state.parms.radiusSize * 2 * Math.PI;
        var fractionInspectedEstimate = (state.timeInRadii[i] * state.parms.stepSize * 0.5) / circ;
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

    if (fracsTotal == 0) {
        return withDiscounts[radius];
    }
    else {
        var extra = (fracs[radius] / fracsTotal) * discountTotal;
        return withDiscounts[radius] + extra;
    }
}

// Given the ant's current distance from the origin,
// this function decides whether the next step in the search
// should move toward or away from the origin. It returns 1
// to indicate 'away', -1 to indicate 'toward', and 0 to
// indicate that the ant should stay in the current radius.
function inOrOut(state, currentRadius) {
    if (currentRadius == 0)
        return 1; // Out
    if (currentRadius == state.parms.numRadii-1)
        return -1; // In

//    var inCurrentR = pInRadius(state, currentRadius);
    var inOuterR = pInRadius(state, currentRadius+1);
    var inInnerR = pInRadius(state, currentRadius-1);
    if (inOuterR > inInnerR)
        return 1;
    else
        return -1
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

    var tanRadRat = whiteNoise(getWhiteNoiseParms(), state.step) * state.parms.scalingParm;
    var trrs = tanRadRat * tanRadRat;
    var rad = 1/Math.sqrt(1+trrs);
    var tan = tanRadRat/Math.sqrt(1+trrs);

    var vecFromRadX = state.posX == 0 ? 1 : state.posX;
    var vecFromRadY = state.posY;
    var l = Math.sqrt((vecFromRadX*vecFromRadX) + (vecFromRadY*vecFromRadY));
    var currentRadius = state.posX == 0 && state.posY == 0 ? 0 : Math.floor(l / state.parms.radiusSize);
    if (currentRadius >= state.parms.numRadii)
        return false;

    state.inOrOut = inOrOut(state, currentRadius);

    state.oldPosX = state.posX;
    state.oldPosY = state.posY;

    if (state.inOrOut == -1) {
        vecFromRadX *= -1;
        vecFromRadY *= -1;
    }

    var sin = vecFromRadY / l;
    var cos = vecFromRadX / l;
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
    if (dsqd > radsqd/2) {
        var d = Math.sqrt(dsqd);
        var sin = state.posY / d;
        var cos = state.posX / d;
        state.posX = rad * cos;
        state.posY = rad * sin;
    }

    if (state.posX < state.radiiMinXs[currentRadius])
        state.radiiMinXs[currentRadius] = state.posX;
    else if (state.posX > state.radiiMaxXs[currentRadius])
        state.radiiMaxXs[currentRadius] = state.posX;

    if (state.posY < state.radiiMinYs[currentRadius])
        state.radiiMinYs[currentRadius] = state.posY;
    else if (state.posY > state.radiiMaxYs[currentRadius])
        state.radiiMaxYs[currentRadius] = state.posY;

    state.timeInRadii[currentRadius]++;

    if (((state.oldPosX < 0 && state.posX > 0) || (state.posX < 0 && state.oldPosX > 0) || (state.posX == 0 && state.oldPosX == 0)) &&
        ((state.oldPosY < 0 && state.posY > 0) || (state.posY < 0 && state.oldPosY > 0) || (state.posY == 0 && state.oldPosY == 0))) {
        state.inOrOut *= -1;
    }

    return true;
}

function sketchProc(p) {
    var state = getInitialSearchState();

    var WIDTH = 820, HEIGHT = 820;

    p.size(WIDTH, HEIGHT);

    var MAXP = 1.0; //maxPInRadius(state);

    var antImage = p.loadImage("ant.png");
    var ANT_IMAGE_WIDTH = 20;
    var ANT_IMAGE_HEIGHT = 20;

    var id;
    function update() {
        // Draw probability circles.
        var psInRadii = new Array(state.parms.numRadii);
        for (var i = state.parms.numRadii-1; i >= 0; --i) {
            var prob = pInRadius(state, i);
            psInRadii[i] = prob;
            var v = (prob*255)/MAXP;
            p.stroke(255, 255, 255);
            p.fill(v,v,v);
            p.ellipse(WIDTH/2, HEIGHT/2, (i+1) * state.parms.radiusSize * 2, (i+1) * state.parms.radiusSize * 2);
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

        p.image(antImage, state.posX + WIDTH/2 - ANT_IMAGE_WIDTH/2, HEIGHT-(state.posY+HEIGHT/2)+ANT_IMAGE_HEIGHT/2, ANT_IMAGE_WIDTH, ANT_IMAGE_HEIGHT);
        if (! updateSearchState(state)) {
            clearInterval(id);
        }
    }
    update();
    id = setInterval(update, 100);
}

window.onload = function () {
    var canvas = document.getElementById("canvas1");
//    text1 = document.getElementById("text1");
//    text2 = document.getElementById("text2");

    var processingInstance = new Processing(canvas, sketchProc);
};
