var appState = {};


appState.nInput = document.getElementById("n");
appState.mInput = document.getElementById("m");

appState.executeSymbolicallyCheckbox = document.getElementById("execute_symbolically");
appState.primeCharacteristicInput = document.getElementById("prime_characteristic");
appState.extensionDegreeInput = document.getElementById("extension_degree");

appState.dForNInput = document.getElementById("d-for-n");
appState.dForMInput = document.getElementById("d-for-m");
appState.equation0 = document.getElementById("equation-0");
appState.equation1 = document.getElementById("equation-1");
appState.addTermButton0 = document.getElementById("add-term-button-0");
appState.addTermButton1 = document.getElementById("add-term-button-1");
appState.latex0 = document.getElementById("latex-0");
appState.latex1 = document.getElementById("latex-1");
appState.resetButton = document.getElementById("reset-button");
appState.recalcButton = document.getElementById("recalc-button");
appState.results0 = document.getElementById("results-0");
appState.results1 = document.getElementById("results-1");
appState.results2 = document.getElementById("results-2");
appState.results3 = document.getElementById("results-3");
appState.results4 = document.getElementById("results-4");
appState.results5 = document.getElementById("results-5");
appState.resultsAMatrices = document.getElementById("results-a-matrices");
appState.resultsMultiplyByA = document.getElementById("results-multiply-by-a");
appState.resultsPhiModuleCondensed = document.getElementById("equation-0-simplified");
appState.resultsPsiModuleCondensed = document.getElementById("equation-1-simplified");
// appState.resultsMaxDegrees = document.getElementById("results-max-degrees");
appState.resultsEMatricesByDegrees = document.getElementById("results-e-matrices-by-degrees");
appState.unitTestResults = document.getElementById("results-unit-tests");
appState.metricsElem = document.getElementById("metrics");
appState.tmodule0 = [];
appState.tmodule1 = [];
appState.computeEngine = null;

appState.hasValidOutput = false;

appState.inputValidationRequestId = 0;

appState.reductionSteps = document.getElementById("reduction-steps");
appState.finalReduction = document.getElementById("final-reduction");

appState.analysisSteps = document.getElementById("analysis-steps");
appState.columnLists = document.getElementById("optimized-steps");
appState.sortedColumnLists = document.getElementById("sorted-column-lists");
appState.flagList = document.getElementById("flag-list");
appState.rightmostBasis = document.getElementById("rightmost-basis");


// https://auth.thirdiron.com/v2/api-tokens, query params: {"libraryId": "1931", "returnPreproxy": true, "client": "native", "success": "/success", "failure": "/failure"}



function init() {
    var variant = getVariant();

    appState.computeEngine = new ComputeEngine.ComputeEngine();

    var outputWrapper = document.getElementById("output-wrapper");
    outputWrapper.style.display = "none";

    var recalcButton = document.createElement("button");
    appState.recalcButton = recalcButton;
    recalcButton.id = "recalc-button";
    recalcButton.style.marginRight = "10px";
    recalcButton.style.display = "inline-block";
    recalcButton.innerHTML = "Run";

    var exportButton = document.getElementById("export-button");
    if ( exportButton ) {
        exportButton.addEventListener("click", function() {
            showExportDialog();
        });
    }

    var exampleButton = document.getElementById("examples-button");
    if ( exampleButton ) {
        exampleButton.addEventListener("click", function() {
            showExampleChooserDialog();
        });
    }

    var optionsButton = document.getElementById("options-button");
    if ( optionsButton ) {
        optionsButton.addEventListener("click", function() {
            showOptionsChooserDialog();
        });
    }

    appState.resetButton.parentElement.insertBefore(recalcButton, appState.resetButton);

    disableRunButton();

    recalcButton.addEventListener("click", function() {
        onInputChanged(true);
    });

    // initNInputField(appState.nInput, 0);
    // initNInputField(appState.mInput, 1);

    initSymbolicExecutionCheckbox();
    initCharacteristicInputField(appState.primeCharacteristicInput, DEFAULT_VALUE__PRIME_CHARACTERISTIC);
    initCharacteristicInputField(appState.extensionDegreeInput, DEFAULT_VALUE__EXTENSION_DEGREE);

    refreshCharacteristicFields();

    initDInputField(appState.dForNInput, 0);
    initDInputField(appState.dForMInput, 1);

    initEquationField(0);
    initEquationField(1);

    appState.resetButton.addEventListener("click", function() {
        showDialog("Are you sure you want to reset?<br>Equations will be cleared of current values.", "cancel", "yes reset", function() {
            resetEquations();
        });
    });

    loadMathFieldValuesFromLocalStorage();

    // var forceServerCalculation = false;
    // var firstStartUp = true;
    // onInputChanged(forceServerCalculation, firstStartUp);
}

function initIntegerInputOnlyField(inputElem) {
    inputElem.addEventListener("keydown", function(event) {
        if( ['-', '.', 'e', 'E'].includes(event.key)) {
            event.preventDefault();
        }
    });
}

function resizeMatricesInInputElem(equationOrdinal, newD) {
    var equationElem = equationOrdinal == 0 ? appState.equation0 : appState.equation1;
    var latexInputElem = equationOrdinal == 0 ? appState.latex0 : appState.latex1;

    var oldLatexString = equationElem.value;

    var matrixSize = getLargestDimensionOfLatexMatrix(oldLatexString)

    if ( matrixSize < 0 ) {
        oldLatexString = wrapInLaTeXMatrix([oldLatexString]);
    }

    var newLatexString = resizeLatexMatrices(oldLatexString, newD, newD, "0");

    setEquationInputProgramatically(equationOrdinal, newLatexString, false);
}

function enableCharacteristicInputElement(inputElem, isSymbolicExecutionEnabledTrueOrFalse) {
    if ( isSymbolicExecutionEnabledTrueOrFalse === false ) {
        inputElem.disabled = false;
        inputElem.parentElement.style.opacity = 1;
    } else {
        inputElem.disabled = true;
        inputElem.parentElement.style.opacity = .5;
    }
}

function isSymbolicExecutionEnabled() {
    var enabled = appState.executeSymbolicallyCheckbox ? appState.executeSymbolicallyCheckbox.checked : false;

    return enabled;
}

function refreshCharacteristicFields() {
    var enabled = isSymbolicExecutionEnabled();

    if ( appState.primeCharacteristicInput ) {
        enableCharacteristicInputElement(appState.primeCharacteristicInput, enabled);
    }

    if ( appState.extensionDegreeInput ) {
        enableCharacteristicInputElement(appState.extensionDegreeInput, enabled);
    }
}

function initSymbolicExecutionCheckbox() {
    if ( appState.executeSymbolicallyCheckbox ) {
        var defautToChecked = true;
        appState.executeSymbolicallyCheckbox.checked = getCheckboxOptionValue("execute_symbolically", defautToChecked);

        appState.executeSymbolicallyCheckbox.addEventListener("change", function() {
            setCheckboxOptionValue("execute_symbolically", appState.executeSymbolicallyCheckbox.checked);
            setVisibleStatus(VISIBLE_STATUS__RESULTS_STALE);
            refreshCharacteristicFields();

            clearCondensedLatexAttributes();

            doContinuationAfterInputChange(function() {
                onInputChanged();
            });
        });
    }
}

function initCharacteristicInputField(inputElem, defaultValue) {
    if ( inputElem ) {
        initIntegerInputOnlyField(inputElem);

        var valueFromLocalStorage = window.localStorage.getItem(inputElem.id);

        if ( valueFromLocalStorage ) {
            defaultValue = parseInt(valueFromLocalStorage);
        }
        
        setInputValueProgrammaticallyForIntInput(inputElem, defaultValue);

        inputElem.addEventListener('change', function() {
            var newValue = parseInt(this.value);
            newValue = newValue >= 1 ? newValue : 1;

            if ( inputElem == appState.primeCharacteristicInput && isPrime(newValue) == false ) {
                setInputValueProgrammaticallyForIntInput(inputElem, inputElem.getAttribute("old_value"));
                showIntInputError(inputElem, "Must be prime.");
            } else {
                window.localStorage.setItem(inputElem.id, newValue);
                setInputValueProgrammaticallyForIntInput(inputElem, newValue);
                setVisibleStatus(VISIBLE_STATUS__RESULTS_STALE);
            }
        });
    }
}

function initDInputField(inputElem, equationOrdinal) {
    if ( inputElem ) {
        initIntegerInputOnlyField(inputElem);
        
        setInputValueProgrammaticallyForIntInput(inputElem, DEFAULT_VALUE__D);

        inputElem.addEventListener('change', function() {
            var newD = parseInt(this.value);
            newD = newD >= 1 ? newD : 1;

            if ( newD > MAX_INPUT_MATRIX_DIMENSION ) {
                setInputValueProgrammaticallyForIntInput(inputElem, inputElem.getAttribute("old_value"));

                var equationElem = equationOrdinal == 0 ? appState.equation0 : appState.equation1;
                showMatrixExceededMaxDimensionsError(equationElem);
            } else {
                setInputValueProgrammaticallyForIntInput(inputElem, newD);

                resizeMatricesInInputElem(equationOrdinal, newD);

                doContinuationAfterInputChange(function() {
                    onInputChanged();
                });
            }
        });
    }
}

function doContinuationAfterInputChange(continuation) {
    setTimeout(() => {
        if (document.activeElement == appState.recalcButton) {
            // Skip input validation since we're running the whole calculation anyway.
            console.log("Skipping input change continuation since we're already running calculation.");
        } else {
            console.log("Calling input change continuation.");
            continuation();
        }
    }, 0);
}

function clearCondensedLatexAttributes() {
    if ( appState.resultsPhiModuleCondensed ) {
        appState.resultsPhiModuleCondensed.setAttribute("last_condensed_latex", "");
    }

    if ( appState.resultsPsiModuleCondensed ) {
        appState.resultsPsiModuleCondensed.setAttribute("last_condensed_latex", "");
    }
}
    

function initEquationField(equationOrdinal) {
    if ( equationOrdinal >= 0 ) {

        var mathField = equationOrdinal === 0 ? appState.equation0 : appState.equation1;
        var addTermButton = equationOrdinal === 0 ? appState.addTermButton0 : appState.addTermButton1;

        mathField.mathVirtualKeyboardPolicy = "manual";
        // mathField.addEventListener("focusin", () =>  mathVirtualKeyboard.show());
        // mathField.addEventListener("focusout", () =>  mathVirtualKeyboard.hide());

        var latexElem = equationOrdinal === 0 ? appState.latex0 : appState.latex1;

        function onMathFieldChange_continuation() {
            doContinuationAfterInputChange(function() {
                onPrimaryMathFieldOrLatexChange_debounced(equationOrdinal);
            });
        }

        function onMathFieldChange(inputElem) {
            latexElem.value = inputElem.value;
            
            onMathFieldChange_continuation();
        }

        function onLatexElemChange(inputElem) {
            mathField.value = inputElem.value;

            onMathFieldChange_continuation();
        }

        // mathField.addEventListener("input", function() {
        //     onMathFieldChange(this);
        // });

        // latexElem.addEventListener("input", function() {
        //     onLatexElemChange(this);
        // });

        mathField.addEventListener("change", function() {
            onMathFieldChange(this);
        });

        latexElem.addEventListener("change", function() {
            onLatexElemChange(this);
        });

        addTermButton.addEventListener("click", function() {
            var existingLatex = mathField.value;
            var matrixSize = getLargestDimensionOfLatexMatrix(existingLatex);

            var constrainedMatrixSize = matrixSize;
            constrainedMatrixSize = constrainedMatrixSize >= 1 ? constrainedMatrixSize : 1;
            constrainedMatrixSize = constrainedMatrixSize <= MAX_INPUT_MATRIX_DIMENSION ? constrainedMatrixSize : MAX_INPUT_MATRIX_DIMENSION;

            var nextTerm = "";
            var startingPlusSign = "";

            if ( existingLatex ) {
                var largestDegreeOfTau = getLargestExponentInLatex(existingLatex, TAU_LATEX);
                var largestDegreeOfTheta = getLargestExponentInLatex(existingLatex, THETA_ON_DIAGONAL);
                var containsTheta = largestDegreeOfTheta >= 0; // Assume that any kind of theta means that it's "correctly" there.
                var isUsingMatrices = matrixSize > 1;

                if ( containsTheta === false ) {
                    if ( isUsingMatrices === true ) {
                        nextTerm = createLatexMatrix(constrainedMatrixSize, constrainedMatrixSize);
                        // "resize" to the same size just to re-use a method that already puts thetas on the diagonal for us.
                        nextTerm = resizeLatexMatrices(nextTerm, constrainedMatrixSize, constrainedMatrixSize, DEFAULT_LATEX_FOR_NEW_MATRIX_ELEMENTS, true)
                    } else {
                        nextTerm = THETA_ON_DIAGONAL;
                    }
                } else {
                    if ( largestDegreeOfTau >= 1 ) {
                        largestDegreeOfTau = largestDegreeOfTau + 1;
                    } else {
                        largestDegreeOfTau = 1;
                    }

                    var latexExponent = largestDegreeOfTau > 9 ? "{" + largestDegreeOfTau + "}" : largestDegreeOfTau;
                    var tauToThePowLatex = "\\tau^" + latexExponent;

                    if ( isUsingMatrices === true ) {
                        nextTerm = createLatexMatrix(constrainedMatrixSize, constrainedMatrixSize) + " \\cdot " + tauToThePowLatex;
                    } else {
                        nextTerm = tauToThePowLatex;
                    }
                }

                startingPlusSign = " + ";
                
            } else {
                nextTerm = THETA_ON_DIAGONAL;
            }

            var newTermLatex = existingLatex + startingPlusSign + nextTerm;
            
            setEquationInputProgramatically(equationOrdinal, newTermLatex);

            scrollMathFieldToEnd(mathField);
        });

        setUpLatexCopyButton(latexElem);
        
        function onLatexCleared(latexElemThatWasCleared) {
            // Clear them first in addition to setEquationInputProgramatically() so focus() is called on empty field.
            latexElem.value = '';
            mathField.value = '';

            latexElemThatWasCleared.focus();

            setEquationInputProgramatically(equationOrdinal, "", true);
        }

        var latexClearButton = document.getElementById(latexElem.id + "-clear");

        latexClearButton.addEventListener("click", function() {
            onLatexCleared(latexElem);
        });
    }
}

window.addEventListener('resize', (event) => {
  onWindowResize();
});

function doResizeRoutine() {
    if ( appState.equation0 ) {
        resizeInputEquationCells(appState.equation0, appState.resultsPhiModuleCondensed);
    }

    if ( appState.equation1 ) {
        resizeInputEquationCells(appState.equation1, appState.resultsPsiModuleCondensed);
    }
}

const onWindowResize = debounce(function() {
    doResizeRoutine();
}, 100);

function getParamInputValue(inputElem) {
    var rawValue = inputElem && inputElem.value ? inputElem.value : DEFAULT_VALUE__D;
    var parsedValue = parseInt(rawValue);

    parsedValue = parsedValue >= 1 ? parsedValue : DEFAULT_VALUE__D;

    return parsedValue;
}

function setEquationInputProgramatically(equationOrdinal, latex, callOnChange = true) {
    var equationInput = equationOrdinal === 0 ? appState.equation0 : appState.equation1;
    var latexInput = equationOrdinal === 0 ? appState.latex0 : appState.latex1;
    equationInput.value = latex;
    latexInput.value = latex;

    if ( callOnChange ) {
        onPrimaryMathFieldOrLatexChange(equationOrdinal);
    } else {
        equationInput.setAttribute("last_expanded_latex", latex);
    }
}

init();