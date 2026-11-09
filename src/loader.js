// Add cache buster to eliminate possibility of browser fetching a cached version of a local file instead of the actual one.
// Important for when you're making iterative edits and all of a sudden your changes aren't taking effect.
var cacheBuster = Math.floor(Math.random() * 99999999);

// var styleElem = document.createElement("link");
// styleElem.rel = "stylesheet";
// styleElem.type = "text/css";
// styleElem.href = "./src/app.css?v=" + cacheBuster;
// document.getElementsByTagName('head')[0].appendChild(styleElem);

var srcFiles = [
    "config",
    "model",
    "api",
    "view",
    "main"
];

for ( var i = 0; i < srcFiles.length; i++ ) {
    var ithSrcFile = srcFiles[i];
    var ithScriptElem = document.createElement("script");
    ithScriptElem.defer = true;
    ithScriptElem.async = false;
    ithScriptElem.src = "./src/"+ithSrcFile+".js?v=" + cacheBuster;
    document.body.appendChild(ithScriptElem);

    // setTimeout(function(scriptElem) {
        
    // }, i*100, ithScriptElem);
}