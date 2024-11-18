// ***** JavaScript D'ni Clock *****
// Author: Jehon the Scribe (mail: Jehon[at]gmx.net)
// Site: http://jehon.mystfans.com (long dead)
// Last Update: 13.05.2002
// If you use this code somewhere, please give proper attribution
var vileeFullName = true; // show the vilee as digit (false) or as name like "lEfo" (true)
var usePartavotee = true; // show pahrtahvo (0-24) and tahvo5 (0-4) instead of gahrtahvo (0-4) and tahvo (0-24)
var special0 = false;     // use "cyclical" 0/25 hybrid digit ]/[ instead of plain zero ].[ for time?
var sep1 = '&bull;';      // separator for D'ni date digits (only if vileeFullName==false)
var sep2 = ' - ';         // separator between date and time
var sep3 = '&bull;';      // separator for D'ni time digits
var sep1r = '.';          // separator for romanized date digits (only if vileeFullName==false)
var sep2r = ' - ';        // separator between romanized date and time
var sep3r = ':';          // separator for romanized time digits
// DniNum[25] is the single 25 digit ]X[, and DniNum[26] the special 25/0 hybrid
var DniNum = new Array('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ')', '!', '@', '#', '$', '%', '^', '&', '*', '(', '[', ']', '{', '}', '\\', '|', ':');
var vileeScriptPrefix = 'lE';
var vileeScript = new Array('fo', 'bro', 'san', 'tar', 'vot', 'vofo', 'vobro', 'vosan', 'votar', 'novU');
var vileeNamesPrefix = 'Lee';
var vileeNames = new Array('fo', 'bro', 'sahn', 'tahr', 'vot', 'vofo', 'vobro', 'vosahn', 'votahr', 'novoo');

// reference date: 9656 leefo 1 0.0.0.0 = 2000-04-20 22:12:48 UTC
//var RefMillisec = Date.parse('2000-04-20') + ((22 * 60 + 12) * 60 + 48) * 1000;
//var RefHar = 9656;
// reference date: 9647 leefo 1 0.0.0.0 = 1991-04-21 17:54:00 UTC
// (time stamp on the original HyperCard Stack file for MYST) - source: RAWA
var RefMillisec = Date.parse('1991-04-21') + ((17 * 60 + 54) * 60 + 00) * 1000;
var RefHar = 9647;
var MillisecPerHar = 31556925216; // Mean Solar Tropical Year in 1995, source: RAWA
//var MillisecPerHar = 31556924886.4;
//var MillisecPerHar = 365.2421875 * 24 * 60 * 60 * 1000; // = 31556925000
var ProranteePerHar = 10 * 29 * 5 * 25 * 25 * 25;       // = 22656250
var PrevTimeString = '';
var har, vilee, yar, gartavo, partavo, tavoMod5, tavo, goran, proran;

function Earth2Dni() {
    var D = new Date();
    var delta = D.getTime() - RefMillisec;

    har = Math.floor(delta / MillisecPerHar);
    delta -= har * MillisecPerHar;
    delta *= ProranteePerHar / MillisecPerHar;
    vilee = Math.floor(delta / (29 * 5 * 25 * 25 * 25));
    delta -= vilee * (29 * 5 * 25 * 25 * 25);
    yar = Math.floor(delta / (5 * 25 * 25 * 25));
    delta -= yar * (5 * 25 * 25 * 25);
    gartavo = Math.floor(delta / (25 * 25 * 25));
    delta -= gartavo * (25 * 25 * 25);
    tavo = Math.floor(delta / (25 * 25));
    delta -= tavo * (25 * 25);
    goran = Math.floor(delta / 25);
    delta -= goran * 25;
    proran = Math.floor(delta);
    // correct potential value underflow
    if (proran < 0) { proran += 25; goran--; }
    if (goran < 0) { goran += 25; tavo--; }
    if (tavo < 0) { tavo += 25; gartavo--; }
    if (gartavo < 0) { gartavo += 5; yar--; }
    if (yar < 0) { yar += 29; vilee--; }
    if (vilee < 0) { vilee += 10; har--; }
    // calculate pahrtahvo from gahrtahvo and tahvo
    partavo = gartavo * 5 + Math.floor(tavo / 5);
    tavoMod5 = tavo % 5;
    // add reference D'ni hahr (year) and make vilee (month) & yahr (day) 1-based instead of 0-based
    har += RefHar;
    vilee++;
    yar++;
}

function ToBase25(myint, useSpecialZero25Hybrid) {
    if (useSpecialZero25Hybrid === undefined) useSpecialZero25Hybrid = false;
    if (myint == 0) return useSpecialZero25Hybrid ? DniNum[26] : DniNum[0];
    var tempStr = '';
    while (myint >= 1) {
        tempStr = DniNum[myint % 25] + tempStr;
        myint = Math.floor(myint / 25);
    }
    return tempStr;
}

function ToDniScript() {
    var temp = ToBase25(har);
    if (vileeFullName) temp += '&nbsp; ' + vileeScriptPrefix + vileeScript[vilee - 1] + ' '
    else temp += sep1 + ToBase25(vilee) + sep1;
    temp += ToBase25(yar) + sep2;
    if (usePartavotee) temp += ToBase25(partavo, special0) + sep3 + ToBase25(tavoMod5)
    else temp += ToBase25(gartavo) + sep3 + ToBase25(tavo, special0);
    temp += sep3 + ToBase25(goran, special0) + sep3 + ToBase25(proran, special0);
    return temp;
}

function ToRoman() {
    var temp = har;
    if (vileeFullName) temp += ' ' + vileeNamesPrefix + vileeNames[vilee - 1] + ' '
    else temp += sep1r + vilee + sep1r;
    temp += yar + sep2r;
    if (usePartavotee) temp += partavo + sep3r + tavoMod5
    else temp += gartavo + sep3r + tavo;
    temp += sep3r + goran + sep3r + proran;
    return temp;
}

function JehonDniClock() {
    Earth2Dni();
    var TimeString = '&nbsp;' + ToDniScript() + '&nbsp;';
    if (TimeString != PrevTimeString) {
        var TimeStringRoman = ToRoman();
        if (document.getElementById) {
            document.getElementById('DniClock').innerHTML = TimeString;
            document.getElementById('DniClock').title = "D'ni Date and Time:\n" + TimeStringRoman;
            document.getElementById('DniClockRoman').innerHTML = TimeStringRoman;
            vileeFullName = document.getElementById("chkVileeNames").checked;
            usePartavotee = document.getElementById("chkPahrtahvotee").checked;
            special0 = document.getElementById("chkSpecial0").checked;
        } else if (document.all) {
            DniClock.innerHTML = TimeString;
            DniClock.title = TimeStringRoman;
            DniClockRoman.innerHTML = TimeStringRoman
            vileeFullName = chkVileeNames.checked;
            usePartavotee = chkPahrtahvotee.checked;
            special0 = chkSpecial0.checked;
        }
        PrevTimeString = TimeString;
    }
    // update at least 10 times per second (every 100 millisec)
    // because a prorahn is 1.39 seconds - otherwise the clock doesn't tick smoothly
    setTimeout("JehonDniClock()", 100);
}

function Init() {
    if (document.getElementById) {
        document.getElementById("chkVileeNames").checked = vileeFullName;
        document.getElementById("chkPahrtahvotee").checked = usePartavotee;
        document.getElementById("chkSpecial0").checked = special0;
    }
}

window.onload = function() {
    Init();
    if (document.getElementById || document.all) JehonDniClock();
}
