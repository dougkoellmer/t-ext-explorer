var RESPONSE_MESSAGE_PART__NETWORK_ERROR = "Network error";
var RESPONSE_MESSAGE_PART__CANCELLED = "cancelled";


function isResponseCancellation(responseMessage) {
    return responseMessage.toLocaleLowerCase().includes(RESPONSE_MESSAGE_PART__CANCELLED.toLocaleLowerCase());
}

function isResponseNetworkError(responseMessage) {
    return responseMessage.toLocaleLowerCase().includes(RESPONSE_MESSAGE_PART__NETWORK_ERROR.toLocaleLowerCase());
}

function isRunningOnLocalHost() {
  if (location.hostname === "localhost" || location.hostname === "127.0.0.1" || location.hostname === "10.0.2.2") {
    return true;
  } else {
    return false;
  }
}

function resizeInputEquationCells(inputElem, elem) {
    var equationWrapper = inputElem.parentElement;
    var equationTableCol = inputElem.parentElement.parentElement;
    equationTableCol.style.width = "100%";

    var renderElem = getMathFieldRenderElem(elem);

    if ( renderElem ) {
        // all good
    } else {
        elem.render(); // try to render.

        renderElem = getMathFieldRenderElem(elem);
    }

    if ( renderElem ) {
        elem.parentElement.style.width = renderElem.clientWidth + "px";
    }

    var equationTableColWidth = equationTableCol.clientWidth;

    if ( equationTableColWidth < window.innerWidth/3 ) {
        equationTableCol.style.width = window.innerWidth/3 + "px";
    }

    if ( equationWrapper.clientHeight < equationTableCol.clientHeight ) {
        equationWrapper.style.height = equationTableCol.clientHeight + "px";
    } else {
        equationWrapper.style.height = "100%;"
    }
}

function setCondensedLatex(inputElem, elem, latex, wasAbleToCondense, hadInputErrors, showTippy = true) {
    var previousLatex = elem.getAttribute("last_condensed_latex");
    elem.setAttribute("last_condensed_latex", latex);

    if ( latex && elem ) {
        setBackEndResultsField(elem, "= " + latex);

        resizeInputEquationCells(inputElem, elem);
    }

    var isSimplifyAnimationTarget = true;

    var condensedLatexChanged = previousLatex != latex;

    if ( condensedLatexChanged ) {
        if ( wasAbleToCondense ) {
            if ( hadInputErrors ) {
                elem.style.color = "#CC5500";
            } else {
                elem.style.color = "darkgreen";
            }
            
            showSimplifyOrExpansionStatus(inputElem, elem, true, isSimplifyAnimationTarget, hadInputErrors, showTippy);
        } else {
            elem.style.color = "red";
            showSimplifyOrExpansionStatus(inputElem, elem, false, isSimplifyAnimationTarget, hadInputErrors, showTippy);
        }
    }
}

function handlePhiPsiInputValidationResponse(endpoint, responseJson, onInputPassed, onInputFailed) {
    
    var inputErrorElem = appState.results0;

    clearInputErrorContainer(0);

    var hadInputErrors = false;
    if ( responseJson.input_errors && responseJson.input_errors.length > 0 ) {
        
        setVisibleStatus(VISIBLE_STATUS__INPUT_ERRORS);

        responseJson.input_errors.forEach(ith => {
            var errorHtml = "<div class='output_error'>"+ith+"</div>";
            addToOutput(inputErrorElem, errorHtml);
        });

        hadInputErrors = true;
    } else {
        var message = "";
        if ( (endpoint == "module_cruncher" && responseJson.success === true) ) {
            message = "No input errors. Click <b>Run</b> in upper-right to re-evaluate.";
        } else if (endpoint == "validate_phi_psi_input" && appState.hasValidOutput ) {
            message = "No input errors. Results below may be stale. Click <b>Run</b> in upper-right to re-evaluate.";
        } else {
            message = "No input errors. Click <b>Run</b> in upper-right to evaluate.";
        }
        
        var successHtml = "<div class='output_success'>"+message+"</div>";
        addToOutput(inputErrorElem, successHtml);

        hadInputErrors = false;
    }

    var showTippy = endpoint != "module_cruncher" && appState.inputValidationRequestId > 1;
    setCondensedLatex(appState.equation0, appState.resultsPhiModuleCondensed, responseJson.phi_module_condensed_latex, responseJson.was_able_to_condense_phi, responseJson.was_able_to_validate_phi == false, showTippy);
    setCondensedLatex(appState.equation1, appState.resultsPsiModuleCondensed, responseJson.psi_module_condensed_latex, responseJson.was_able_to_condense_psi, responseJson.was_able_to_validate_psi == false, showTippy);

    if ( hadInputErrors ) {
        onInputFailed();
    } else {
        onInputPassed();
    }
}

function onPrimaryMathFieldOrLatexChange(equationOrdinal) {

    var mathField = equationOrdinal === 0 ? appState.equation0 : appState.equation1;
    var latexElem = equationOrdinal === 0 ? appState.latex0 : appState.latex1;

    var expandedLatex = mathField.value;
    var lastExpandedLatex = mathField.getAttribute("last_expanded_latex");

    if ( lastExpandedLatex != null && expandedLatex == lastExpandedLatex ) {
        // Likely getting here from a click-off after the input has been changed.
        // We already know it changed, so ignore.
        return;
    }

    // console.log("Doing another change event: simplifiedLatex="+simplifiedLatex+", lastSimplifiedLatex="+lastSimplifiedLatex )

    mathField.setAttribute("last_expanded_latex", expandedLatex);

    onInputChanged();
}

function markCurrentRequestId() {
    appState.inputValidationRequestId += 1;

    var currentRequestId = appState.inputValidationRequestId;

    return currentRequestId;
}

function shouldProcessResponse(currentRequestId) {
    // Ensure that an old request that took longer doesn't stomp out a newer request that returned sooner.
    // This also makes it so that if a few requests are sent right in a row, practically speaking only
    // the last response is handled.
    if ( currentRequestId /*still*/ == appState.inputValidationRequestId ) {
        return true;
    } else {
        return false;
    }
}

const onPrimaryMathFieldOrLatexChange_debounced = debounce(onPrimaryMathFieldOrLatexChange, DEBOUNCED_IN_MILLIS__EXPAND_OR_SIMPLIFY);

function validatePhiPsiInput(phiLatex, psiLatex) {

    setVisibleStatus(VISIBLE_STATUS__CHECKING_INPUT);

    var simulatedDelay = 0;
    var currentRequestId = markCurrentRequestId();
    
    setTimeout(function() {
        var requestJson = {
        };

        buildCommonRequestJson(requestJson, phiLatex, psiLatex);

        var endpoint = "validate_phi_psi_input";

        return callApiEndpoint(endpoint, requestJson, function(responseJson) {

            if ( shouldProcessResponse(currentRequestId) ) {
                setVisibleStatus(VISIBLE_STATUS__RESULTS_STALE);

                handlePhiPsiInputValidationResponse(endpoint, responseJson, function() {
                }, function() {
                });
            }
        });
    }, simulatedDelay)
}

function buildCommonRequestJson(requestJson, phiLatex, psiLatex) {
    requestJson["phi_module"] = phiLatex;
    requestJson["psi_module"] = psiLatex;
    requestJson["show_reduction_steps"] = false; //getCheckboxOptionValue("show_reduction_steps");
    requestJson["replace_c_with_tau"] = true; //getCheckboxOptionValue("replace_c_with_tau");
    requestJson["replace_q_with_integer"] = true; //getCheckboxOptionValue("replace_q_with_integer");
    requestJson["group_coefficients"] = false; // getCheckboxOptionValue("group_coefficients");
    requestJson["strip_integer_exponents_for_symbolic_execution"] = true; //getCheckboxOptionValue("strip_integer_exponents_for_symbolic_execution");
    requestJson["execute_symbolically"] = isSymbolicExecutionEnabled();
    requestJson["prime_characteristic"] = appState.primeCharacteristicInput.value;
    requestJson["extension_degree"] = appState.extensionDegreeInput.value;
    // requestJson["only_generate_ext"] = true;
}

function runCalculationOnServer(phiLatex, psiLatex, callback) {

    var variant = getVariant();

    var requestJson = {
        variant: variant
    };

    if ( variant == VARIANT__PHI_PSI_EMBEDDING ) {
        buildCommonRequestJson(requestJson, phiLatex, psiLatex);
    } else if ( variant == VARIANT__DIAGONAL ) {
        requestJson["upsilon_module"] = flippedPhiModule;
    } else if ( variant == VARIANT__RIGHTMOST_BASIS ) {
        requestJson["phi_module"] = flippedPhiModule;
    }

    var currentRequestId = markCurrentRequestId();

    var endpoint = "module_cruncher";
    return callApiEndpoint(endpoint, requestJson, function(responseJson) {

        if ( shouldProcessResponse(currentRequestId) ) {
            handlePhiPsiInputValidationResponse(endpoint, responseJson, function() {
                callback(responseJson); // Continue on if input is OK.
            }, function() {
                showErrorOutput("Error: Could not validate input.");
            });
        }
    });
}

function callApiEndpoint(endpoint, requestJson, callback) {
    var apiRoot = "";

    if ( isRunningOnLocalHost() ) {
        apiRoot = "http://localhost:6789";
    } else {
        apiRoot = "https://t-module-cruncher-hnkpe.ondigitalocean.app";
    }

    var url = apiRoot + '/' + endpoint;

    var xhr = new XMLHttpRequest();
    xhr.open("POST", url, true);
    xhr.setRequestHeader("Content-Type", "application/json");

    xhr.onload = function() {
        // 1. HTTP-Level & Application Failures
        if (xhr.status >= 200 && xhr.status < 300) {
            // console.log("API response successful, parsing JSON...");
            var jsonResponse = {};
            try {
                
                var rawText = xhr.responseText
                var jsonStartIndex = rawText.indexOf('{');
            
                if (jsonStartIndex === -1) {
                    throw new Error("No JSON object found in response");
                }
                
                // 2. Extract and trim only the JSON portion
                var cleanJsonText = rawText.substring(jsonStartIndex).trim();

                cleanJsonText = cleanJsonText.endsWith("0") ? cleanJsonText.slice(0, -1) : cleanJsonText;
                
                jsonResponse = JSON.parse(cleanJsonText);

                // console.log("API response successful, processing...")
            } catch (error ) {
                callback("Error handling response: " + error);
                return;
            }

            // console.log("JSON parse successful, returning result...");
            callback(jsonResponse);
        } else {
            // e.g., 401 Unauthorized, 404 Not Found, 500 Server Error
            console.warn(`HTTP Error: ${xhr.status} ${xhr.statusText}`);
            callback(`HTTP Error: ${xhr.status} ${xhr.statusText}`);
        }
    };

    // 2. Network-Level Failures
    xhr.onerror = function(e) {
        console.error('Network Error (DNS failure, lost connection, CORS issue):', e);
        callback(RESPONSE_MESSAGE_PART__NETWORK_ERROR + ". Please check your Internet connection and try again.");
    };

    // 3. Timeout Failures
    xhr.ontimeout = function() {
        console.error('Request timed out.');
        callback("Request timed out.");
    };

    // 4. Abort Failures
    xhr.onabort = function() {
        // console.error('Request was aborted by the user or script.');
        callback("Calculation was "+RESPONSE_MESSAGE_PART__CANCELLED+".");
    };
    
    xhr.send(JSON.stringify(requestJson));

    return xhr;
}
