var DEFAULT_VALUE__D = 2;
var DEFAULT_VALUE__N = 2;

var DEFAULT_VALUE__PRIME_CHARACTERISTIC = 3;
var DEFAULT_VALUE__EXTENSION_DEGREE = 2;

var TAU_LATEX = "\\tau";
var THETA_ON_DIAGONAL = "\\theta";
var MATRIX_ELEM__ZERO = "0";

var MATRIX_NAME__W = "MATRIX_NAME__W";
var MATRIX_NAME__V = "MATRIX_NAME__V";
var MATRIX_NAME__WV = "MATRIX_NAME__WV";

// NOTE: None of the uses of this constant actually ended up working.
// I had to find/replace "maxMatrixCols:10" with "maxMatrixCols:64" in mathlive.js.
// I'm sure there's a more elegant way to configure it from the outside, but I gave up
// after trying at least 10 different things. I left some of the attempts in case
// they ever become relevant.
var MAX_INPUT_MATRIX_DIMENSION = 64;

var DEFAULT_LATEX_FOR_NEW_MATRIX_ELEMENTS = "0";
var DEFAULT_LATEX_MATRIX_ENVIRONMENT = "pmatrix";

MathfieldElement.maxMatrixCols = MAX_INPUT_MATRIX_DIMENSION;

MathfieldElement.macros = {
  ...MathfieldElement.macros,
  MaxMatrixCols: '' + MAX_INPUT_MATRIX_DIMENSION,
  maxMatrixCols: '' + MAX_INPUT_MATRIX_DIMENSION,
};

var VARIANT__PHI_PSI_EMBEDDING = "VARIANT__PHI_PSI_EMBEDDING";
var VARIANT__DIAGONAL = "VARIANT__DIAGONAL";
var VARIANT__RIGHTMOST_BASIS = "VARIANT__RIGHTMOST_BASIS";

var FEATURE_FLAG__RENDER_ALL_MATH_OUTPUT_UP_FRONT = false;
var FEATURE_FLAG__USE_LATEX_PLACEHOLDERS = false;

var R_MATRIX_CELL_SIZE_PX = 64;


function newOutputMatrixConfig(baseId, title, responseKey) {
    return {
        base_id: baseId,
        title: title,
        response_key: responseKey
    };
}

var R_AND_S_CONFIG = newOutputMatrixConfig("r-matrix", "R\\text{ and }S", "R_matrix");

var S_CONFIG = newOutputMatrixConfig("s-matrix", "S", "S_matrix");

var OUTPUT_MATRIX_CONFIGS = [
    newOutputMatrixConfig("ext-matrix", "\\text{Ext}^1", "Ext_matrix"),
    newOutputMatrixConfig("ext0-matrix", "\\text{Ext}^0", "Ext0_matrix"),
    newOutputMatrixConfig("k-matrix", "K", "K_matrix"),
    newOutputMatrixConfig("p-matrix", "P", "P_matrix"),

    R_AND_S_CONFIG,
    S_CONFIG
    
    // newOutputMatrixConfig("v-matrix", "Exact Sequence Partial Basis", "V_matrix")
];

var CHECKBOX_OPTIONS = [
    newCheckboxOption(
        "show_reduction_steps",
        "Show Reduction Steps",
        "WARNING: Enabling this may cause poor performance with large examples or slow computers.",
        false
    ),
    newCheckboxOption(
        "replace_c_with_tau",
        "Substitute <math-span class='options_chooser_math'>c^{q^k}</math-span> with <math-span>\\tau^k</math-span> in final output matrices",
        "Perform said substitution in the final Ext/Ext0 output matrices.<br>Symbol <math-span>c</math-span> by itself with no exponent becomes <math-span>1</math-span>",
        true
    ),
    newCheckboxOption(
        "replace_q_with_integer",
        "Substitute e.g. <math-span class='options_chooser_math'>\\theta^{q^3}</math-span> with <math-span class='options_chooser_math'>\\theta^{(3)}</math-span> in final output matrices",
        "NOTE: This is NOT applied to non-symbols (integers) if Option 5 is enabled and Characteristic is None (prime and n not provided).",
        true
    ),
    newCheckboxOption(
        "group_coefficients",
        "Group coefficients by <math-span class='options_chooser_math'>\\tau</math-span> between each reduction step",
        "Groups coefficients by <math-span>\\tau</math-span> between each reduction.<br>If disabled, it may lead to more reduction steps, but the final Ext/Ext0 results should be the same.",
        true
    ),
    newCheckboxOption(
        "strip_integer_exponents_for_symbolic_execution",
        "Substitute e.g. <math-span class='options_chooser_math'>5^{q^3}</math-span> or <math-span class='options_chooser_math'>5^{(3)}</math-span> with simply <math-span class='options_chooser_math'>5</math-span> for symbolic execution",
        "This setting is only applicable when Characteristic is None (prime and n not provided)",
        true
    ),
];

function getCheckboxOptionValue(localStorageKey, defaultToChecked = null) {
    var value = window.localStorage.getItem(localStorageKey);

    if ( value == "true" || value == "false" ) {
        return value == "true" ? true : false;
    } else {
        if ( defaultToChecked != null ) {
            return defaultToChecked;
        } else {
            var checkboxOption = CHECKBOX_OPTIONS.find(ith => {
                return ith.local_storage_key == localStorageKey
            });

            return checkboxOption.default_value;
        }
    }
}

function setCheckboxOptionValue(localStorageKey, value) {
    return window.localStorage.setItem(localStorageKey, value ? "true" : "false");
}

function newCheckboxOption(localStorageKey, readableLabel, tooltipTitle, defaultValue) {
    return {
        local_storage_key: localStorageKey,
        readable_label: readableLabel,
        tooltip_title: tooltipTitle,
        default_value: defaultValue
    };
}

var EXAMPLES__PHI_PSI_EMBEDDDING = [
    [
        // This one was what exposed the need for coefficient grouping a frobenius expansion.
        String.raw`\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\theta & 0\\ 0 & \theta\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^0+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}0 & 1\\ 0 & 0\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{1}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}0 & 0\\ 1 & 1\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{2}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}1 & 0\\ 0 & 0\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{3}`,
        String.raw`\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\theta & 0 \\1 & \theta \\\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot \tau^{0} + \begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\theta^{2} & 0 \\1 & \theta \\\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot \tau^{1}`
    ],
    [
        // Same as above but simplified and with subscripted a symbols.
        String.raw`\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\theta + a_{2} \tau^{3} & \tau\\a_{1} \tau^{2} & \theta + \tau^{2}\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}`,
        String.raw`\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\theta + \theta^{2} \tau & 0\\1 + \tau & \theta + \theta \tau\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}`
    ],
    [
        // This one caused the infinite loop due to misinterpretation of what "k" should be for the embedding step.
        String.raw`\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\theta & 0\\ 1 & \theta\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^0+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}0 & 0\\ 0 & 0\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{1}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}0 & 0\\ 0 & 0\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{2}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}0 & 1\\ 0 & 1\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{3}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}1 & 0\\ 0 & 0\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{4}`,
        String.raw`\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\theta & 0\\ 1 & \theta\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^0+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}0 & 0\\ 0 & 1\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{1}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}1 & 0\\ 0 & 0\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{2}`
    ],
    [
        // This one used for initial performance profiling against mathematica.
        // As an example: To match Mathematica, striving for a bit less than .1 seconds to solve for just the Ext matrix (W) on a decent machine.
        // Currently sitting at 0.12 seconds to solve this matrix. So not quite as fast as Mathematica, but at least not orders of magnitude slower.
        // Also (A) Mathematica should be about 10x faster, and (B) there are still likely some optimizations to make to the Sympy code. Optimizations
        // are dangerous though in the sense that  they may introduce bugs in certain cases, especially around equation simplification.
        String.raw`\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\theta & 0\\ 1 & \theta\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^0+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}0 & 0\\ 0 & 0\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{1}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}0 & 0\\ 0 & 0\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{2}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}0 & 1\\ 1 & 0\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{3}`,
        String.raw`\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\theta & 0\\ 1 & \theta\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^0+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}0 & 0\\ 0 & 1\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{1}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}1 & 0\\ 0 & 0\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{2}`
    ],
    [
        String.raw`\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\theta & 0\\ 1 & \theta\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^0+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}0 & 0\\ 0 & 0\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{1}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}0 & 0\\ 1 & 0\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{2}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}0 & 1\\ 0 & 0\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{3}`,
        String.raw`\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\theta & 0 \\1 & \theta \\\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^0 + \begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}0 & 0 \\0 & 1 \\\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{1}`
    ],
    [
        String.raw`\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\theta & 0 & 0\\ 0 & \theta & 0\\ 0 & 0 & \theta\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^0+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}6 & 2 & 3\\ 4 & 5 & 6\\ 2 & 9 & 8\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{1}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}0 & 0 & 3\\ 0 & 0 & 4\\ 1 & 1 & 6\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{2}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}1 & 2 & 1\\ 3 & 4 & 2\\ 1 & 2 & 1\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{3}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}1 & 4 & 2\\ 1 & 5 & 2\\ 5 & 6 & 1\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{4}`,
        String.raw`\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\theta & 0\\ 0 & \theta\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^0+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}2 & 6\\ 2 & 7\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{1}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}1 & 5\\ 4 & 9\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{2}`
    ],
    [
        String.raw`\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\theta & 0\\ 0 & \theta\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^0+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}0 & 0\\ 0 & 0\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{1}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}1 & 0\\ 1 & 0\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{2}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}a & 0\\ 0 & 1\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{3}`,
        String.raw`\theta+\tau`
    ]
];



var EXAMPLES__PHI_BASIS = [
    [
        String.raw`\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\theta & 0 & 0\\ 0 & \theta & 0\\ a & 0 & \theta\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^0+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}a & 0 & 1\\ a & 0 & 2\\ a & 0 & 3\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{1}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}g & 0 & 0\\ 0 & h & 0\\ 0 & 0 & h\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{2}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}x & l & m\\ x & l & m\\ x & l & m\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{3}`
    ],
    [
        String.raw`\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\theta & 0 & 0\\ 0 & \theta & 0\\ 0 & 0 & \theta\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^0+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}1 & 0 & 2\\ 0 & 0 & 1\\ 0 & 0 & 0\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{1}+\begin{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}0 & 1 & 0\\ 0 & 0 & 0\\ 1 & 1 & 1\end{${DEFAULT_LATEX_MATRIX_ENVIRONMENT}}\cdot\tau^{2}`
    ]
];

var DEBOUNCED_IN_MILLIS__EXPAND_OR_SIMPLIFY = 100; // Not really a concern anymore since expression isn't validated on every single change

const debounce = (callback, wait) => {
  let timeoutId = null;
  return (...args) => {
    window.clearTimeout(timeoutId);
    timeoutId = window.setTimeout(() => {
      callback.apply(null, args);
    }, wait);
  };
};

function getVariant() {

    return VARIANT__PHI_PSI_EMBEDDING;
    // let path = window.location.pathname;
    // let fileName = path.substring(path.lastIndexOf('/') + 1);

    // if ( fileName.includes("reduction_embedding") ) {
    //     return VARIANT__PHI_PSI_EMBEDDING;
    // } else if ( fileName == "index2.html" ) {
    //     return VARIANT__DIAGONAL;
    // } else if ( fileName.includes("rightmost_basis") ) {
    //     return VARIANT__RIGHTMOST_BASIS;
    // }
}

function isSingleEquationVariant() {
    var variant = getVariant();

    return variant == VARIANT__DIAGONAL || variant == VARIANT__RIGHTMOST_BASIS;
}

function nilPddMatrix() {
    return [];
}

function getDefaultForN(equationOrdinal) {
    return equationOrdinal === 0 ? DEFAULT_VALUE__N + 1 : DEFAULT_VALUE__N;
}

function getLocalStorageFieldName(equationOrdinal) {
     var fieldName = "";

     var variant = getVariant();
    
    if ( variant == VARIANT__DIAGONAL ) {
        fieldName = "equation_diagonal";
    } else if ( variant == VARIANT__PHI_PSI_EMBEDDING ) {
         fieldName = "equation_" + equationOrdinal;
    } else if ( variant == VARIANT__RIGHTMOST_BASIS ) {
         fieldName = "equation_dimensional";
    }

    return fieldName;
}

function getMathFieldValueFromLocalStorage(equationOrdinal) {
    var variant = getVariant();
    
    var fieldName = getLocalStorageFieldName(equationOrdinal);

    var latexStringFromLocalStorage = window.localStorage.getItem(fieldName);
    latexStringFromLocalStorage = latexStringFromLocalStorage ? latexStringFromLocalStorage : "";

    return latexStringFromLocalStorage;
}

function loadMathFieldValuesFromLocalStorage() {
    var latex0 = getMathFieldValueFromLocalStorage(0);
    var latex1 = getMathFieldValueFromLocalStorage(1);

    loadLatexIntoEquationInput(0, latex0, /*callOnChange=*/false);
    loadLatexIntoEquationInput(1, latex1, /*callOnChange=*/true);
}

function loadLatexIntoEquationInput(equationOrdinal, latexString, callOnChange = true) {
    if ( latexString ) {
        // All good.
    } else {
        latexString = EXAMPLES__PHI_PSI_EMBEDDDING[0][equationOrdinal];
    }

    setEquationInputProgramatically(equationOrdinal, latexString, callOnChange);
}

function saveMathFieldValuesToLocalStorage() {
    if ( appState.equation0 ) {
        var equation0LatexString = appState.equation0.value;

        var fieldName = getLocalStorageFieldName(0);
        
        localStorage.setItem(fieldName, equation0LatexString);
    }
    
    if ( appState.equation1 ) {
        var equation1LatexString = appState.equation1.value;
        var fieldName = getLocalStorageFieldName(1);
        localStorage.setItem(fieldName, equation1LatexString);
    }
}

function clearMathFieldValuesFromLocalStorage() {

    var variant = getVariant();
    var fieldName0 = getLocalStorageFieldName(0);
    var fieldName1 = getLocalStorageFieldName(0);

    localStorage.removeItem(fieldName0);

    if ( fieldName1 ) {
        localStorage.removeItem(fieldName0);
    }
}

function isPrime(num) {
  if (num <= 1) {
    return false;
  }
  
  if (num <= 3) {
    return true;
  }
  
  if (num % 2 === 0 || num % 3 === 0) {
    return false;
  }
  
  for (let i = 5; i * i <= num; i += 6) {
    if (num % i === 0 || num % (i + 2) === 0) {
      return false;
    }
  }
  
  return true;
}