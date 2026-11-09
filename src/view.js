
var VISIBLE_STATUS__AWAITING_INPUT = "VISIBLE_STATUS__AWAITING_INPUT";
var VISIBLE_STATUS__CHECKING_INPUT = "VISIBLE_STATUS__CHECKING_INPUT";
var VISIBLE_STATUS__RESULTS_STALE = "VISIBLE_STATUS__RESULTS_STALE";
var VISIBLE_STATUS__INPUT_ERRORS = "VISIBLE_STATUS__INPUT_ERRORS";
var VISIBLE_STATUS__CALCULATING = "VISIBLE_STATUS__CALCULATING";
var VISIBLE_STATUS__CALCULATION_ERRORS = "VISIBLE_STATUS__CALCULATION_ERRORS";
var VISIBLE_STATUS__CALCULATION_CANCELLED = "VISIBLE_STATUS__CALCULATION_CANCELLED";
var VISIBLE_STATUS__NETWORK_ERRORS = "VISIBLE_STATUS__NETWORK_ERRORS";
var VISIBLE_STATUS__COMPLETE = "VISIBLE_STATUS__COMPLETE";

function hideErrorOutput() {
    var errorOutputWrapper = document.getElementById("error-output-wrapper");

    if ( errorOutputWrapper ) {
        errorOutputWrapper.style.display = "none";
    }
}

function getMathFieldRenderElem(mathDivOrSpan) {
    var render = Array.from(mathDivOrSpan.shadowRoot.querySelectorAll(".ML__latex"))[0];
    
    return render;
}

function isOutputVisible() {
    var outputWrapper = document.getElementById("output-wrapper");

    return outputWrapper ? outputWrapper.checkVisibility() : false;
}

function showErrorOutput(errorMessage) {
    var errorOutputWrapper = document.getElementById("error-output-wrapper");

    if ( errorOutputWrapper ) {
        // all good
    } else {
        errorOutputWrapper = document.createElement("div");
        errorOutputWrapper.id = "error-output-wrapper";

        errorOutput = document.createElement("div");
        errorOutput.id = "error-output";

        errorOutputWrapper.appendChild(errorOutput);

        var outputWrapper = document.getElementById("output-wrapper");

        outputWrapper.parentElement.insertBefore(errorOutputWrapper, outputWrapper);
    }

    errorOutputWrapper.style.display = "block";

    var errorMessageForSimplifiedFields = "";

    if ( isResponseCancellation(errorMessage) ) {
        errorMessageForSimplifiedFields = errorHtml("Could not simplify, calculation was cancelled.");
    } else {
        errorMessageForSimplifiedFields = errorHtml("Could not simplify, see error below.");
    }
    

    // if ( appState.resultsPhiModuleCondensed ) {
    //     appState.resultsPhiModuleCondensed.innerHTML = errorMessageForSimplifiedFields;
    // }

    // if ( appState.resultsPsiModuleCondensed ) {
    //     appState.resultsPsiModuleCondensed.innerHTML = errorMessageForSimplifiedFields;
    // }

    errorOutput.innerHTML = errorMessage;
}

function removeAllEventListeners(element) {
  var deepClone = true;
  var newElement = element.cloneNode(deepClone);

  if ( element.parentNode ) {
    element.parentNode.replaceChild(newElement, element);
  }

  return newElement;
}

function showLoadingModal(xhr) {
    var loadingModal = document.getElementById("loading-modal");
    if ( loadingModal ) {

    } else {
        loadingModal = document.createElement("div");
        loadingModal.id = "loading-modal";

        var loadingModalChild = document.createElement("div");
        loadingModalChild.classList.add("loading_modal_spinner");
        
        loadingModal.appendChild(loadingModalChild);

        var loadingModalCancelButton = document.createElement("button");
        loadingModalCancelButton.innerHTML = "Cancel";
        loadingModalCancelButton.id = "loading-modal-cancel-button";
        loadingModalCancelButton.style.marginTop = "10px";
        loadingModal.appendChild(loadingModalCancelButton);

        document.body.appendChild(loadingModal);
    }

    var loadingModalCancelButton = document.getElementById("loading-modal-cancel-button");
    loadingModalCancelButton = removeAllEventListeners(loadingModalCancelButton);

    loadingModalCancelButton.addEventListener("click", function() {
        if (xhr && xhr.readyState > 0 && xhr.readyState < 4) {
            xhr.abort();
            xhr = null;
        }
    });

    loadingModal.classList.add("active");
}

// Grab your math-field element
const mf = document.querySelector('math-field');

// function getRowCountOfFirstMatrixInLatex(latexStr) {
//   // Find the first matrix environment (e.g., bmatrix, pmatrix, matrix, vmatrix, etc.)
//   const match = latexStr.match(/\\begin\{(b|p|v|V|B)?matrix\}([\s\S]*?)\\end\{(b|p|v|V|B)?matrix\}/);
//   if (!match) {
//     return -1;
//   }

//   const matrixBody = match[2].trim();
//   if (!matrixBody) {
//     return -1;
//   }

//   // Split by row separator '\\', ignoring escaped or trailing empty entries
//   const rows = matrixBody.split(/\\\\/).filter(row => row.trim() !== '');
  
//   return rows.length;
// }

function getLargestDimensionOfLatexMatrix(latexStr) {
    var biggestDimension = -1;

    traverseLatexMatrices(latexStr, (match, envType, innerContent) => {
        
        const rawRows = innerContent.split(/\\\\/).map(row => row.trim()).filter(Boolean);

        var maxCols = 0;
        for (const row of rawRows) {
            const ithColCount = row.split('&').length;
            if (ithColCount > maxCols) {
                maxCols = ithColCount;
            }
        }

        if ( rawRows.length > 0 && rawRows.length > biggestDimension ) {
            biggestDimension = rawRows.length;
        }

        if ( maxCols > 0 && maxCols > biggestDimension ) {
            biggestDimension = maxCols;
        }

        return innerContent;
    });

    return biggestDimension;
}

function wrapInLaTeXMatrix(array, type = 'pmatrix') {
    // Join elements with the LaTeX column separator
    const innerContent = array.join(' & ');
    
    // Wrap in the specified LaTeX matrix environment
    return `\\begin{${type}}${innerContent}\\end{${type}}`;
}

function traverseLatexMatrices(latexStr, traverser) {
    // Regex matches \begin{matrix} ... \end{matrix} and its bracketed variants
    const matrixRegex = /\\begin\{(matrix|bmatrix|pmatrix|vmatrix|Vmatrix|Bmatrix)\}([\s\S]*?)\\end\{\1\}/g;

    return latexStr.replace(matrixRegex, traverser);
}

function getLargestExponentInLatex(latexStr, symbol) {
  // Escape the symbol to safely use it inside a Regular Expression
  const escapedSymbol = symbol.replace(/[-\/\\^$*+?.()|[\]{}]/g, '\\$&');
  
  // Regex matches:
  // 1. The symbol itself
  // 2. Optional exponent group: ^ followed by either:
  //    - A braced integer like {-15} or {42}
  //    - A single digit (optionally negative) like ^2 or ^-3
  const pattern = new RegExp(
    `${escapedSymbol}(?:\\s*\\^\\s*(?:\\{\\s*(-?\\d+)\\s*\\}|(-?\\d)))?`, 
    'g'
  );
  
  let maxExponent = -Infinity;
  let found = false;
  let match;

  // Scan through all occurrences of the symbol in the LaTeX string
  while ((match = pattern.exec(latexStr)) !== null) {
    found = true;
    const bracedExp = match[1];
    const singleExp = match[2];
    
    let currentExponent;
    if (bracedExp !== undefined) {
      currentExponent = parseInt(bracedExp, 10);
    } else if (singleExp !== undefined) {
      currentExponent = parseInt(singleExp, 10);
    } else {
      // If the symbol exists but has no ^ indicator, its implicit exponent is 1
      currentExponent = 1;
    }
    
    if (currentExponent > maxExponent) {
      maxExponent = currentExponent;
    }
  }

  return found ? maxExponent : -1;
}

function createLatexMatrix(rows, cols, fillValue = DEFAULT_LATEX_FOR_NEW_MATRIX_ELEMENTS, env = 'pmatrix') {

  const rowContent = Array(cols).fill(fillValue).join(' & ');
  
  const matrixBody = Array(rows).fill(rowContent).join(' \\\\\n  ');
  
  return `\\begin{${env}}\n  ${matrixBody}\n\\end{${env}}`;
}

function resizeLatexMatrices(latexStr, targetRows, targetCols, padValue = DEFAULT_LATEX_FOR_NEW_MATRIX_ELEMENTS, forceThetaOnDiagonal = false) {

    var pushedThetaAlready = false;

    return traverseLatexMatrices(latexStr, (match, envType, innerContent) => {

        const rawRows = innerContent.split(/\\\\/).map(row => row.trim()).filter(Boolean);

        const currentMatrix = rawRows.map(row =>
            row.split('&').map(cell => cell.trim())
        );

        const newRows = [];

        for (let r = 0; r < targetRows; r++) {
            const newCols = [];
            for (let c = 0; c < targetCols; c++) {
                // If the cell exists in the original matrix, keep it. Otherwise, pad it.
                if (currentMatrix[r] && currentMatrix[r][c] !== undefined) {
                    if ( forceThetaOnDiagonal && r == c ) {
                        newCols.push(THETA_ON_DIAGONAL);
                    } else {
                        newCols.push(currentMatrix[r][c]);
                    }
                } else {
                    if (r == c && pushedThetaAlready == false) {
                        newCols.push(THETA_ON_DIAGONAL);
                    } else {
                        newCols.push(padValue);
                    }
                }
            }
            newRows.push(newCols.join(' & '));
        }

        pushedThetaAlready = true;

        var includeNewlines = false;
        if (includeNewlines) {
            return `\\begin{${envType}}\n  ${newRows.join(' \\\\\n  ')}\n\\end{${envType}}`;
        } else {
            return `\\begin{${envType}}  ${newRows.join(' \\\\  ')}\\end{${envType}}`;
        }
    });
}

function hideLoadingModal() {
    var loadingModal = document.getElementById("loading-modal");

    if ( loadingModal ) {
        loadingModal.classList.remove("active");
    }
}

function disableRunButton() {
    var recalcButton = document.getElementById("recalc-button");
    recalcButton.disabled = true;
    recalcButton.style.cursor = "default";
}

function enableRunButton() {
    var recalcButton = document.getElementById("recalc-button");
    recalcButton.disabled = false;
    recalcButton.style.cursor = "pointer";
}

function setVisibleStatus(visibleStatus, color) {

    var text = "";
    var color = "";

    var outputWrapper = document.getElementById("output-wrapper");

    hideErrorOutput();

    var showOutput = false;

    if ( visibleStatus == VISIBLE_STATUS__AWAITING_INPUT ) {
        // disableRunButton();
        enableRunButton();
        text = "AWAITING INPUT";
        color = "green";
    } else if ( visibleStatus == VISIBLE_STATUS__RESULTS_STALE ) {
        showOutput = true;
        enableRunButton();
        text = "RESULTS STALE";
        color = "orange";

        outputWrapper.style.opacity = ".25";
    } else if ( visibleStatus == VISIBLE_STATUS__CHECKING_INPUT ) {
        showOutput = true;
        enableRunButton();
        text = "CHECKING INPUT...";
        color = "orange";

        outputWrapper.style.opacity = ".25";
    } if ( visibleStatus == VISIBLE_STATUS__INPUT_ERRORS ) {
        showOutput = true;
        // disableRunButton();
        enableRunButton();
        hideLoadingModal();
        text = "INPUT ERROR";
        color = "red";
    } else if ( visibleStatus == VISIBLE_STATUS__CALCULATING ) {
        disableRunButton();
        text = "CALCULATING...";
        color = "orange";
    } else if ( visibleStatus == VISIBLE_STATUS__CALCULATION_CANCELLED ) {
        appState.hasValidOutput = false;
        enableRunButton();
        hideLoadingModal();
        text = "CANCELLED";
        color = "red";
    } else if ( visibleStatus == VISIBLE_STATUS__CALCULATION_ERRORS ) {
        appState.hasValidOutput = false;
        enableRunButton();
        hideLoadingModal();
        text = "OUTPUT ERROR";
        color = "red";
    } else if ( visibleStatus == VISIBLE_STATUS__NETWORK_ERRORS ) {
        appState.hasValidOutput = false;
        enableRunButton();
        hideLoadingModal();
        text = "NETWORK ERROR";
        color = "red";
    } else if ( visibleStatus == VISIBLE_STATUS__COMPLETE ) {
        appState.hasValidOutput = true;
        showOutput = true;
        enableRunButton();
        hideLoadingModal();
        text = "COMPLETE";
        color = "green";
        outputWrapper.style.opacity = "1";
    }

    if ( showOutput && appState.hasValidOutput == true ) {
        outputWrapper.style.display = "block";
    } else {
        outputWrapper.style.display = "none";
    }

    var elem = document.getElementById("staleness-indicator");

    if ( elem ) {

    } else {
        var recalcButton = document.getElementById("recalc-button");
        elem = document.createElement("div");
        elem.id = "staleness-indicator";
        elem.style.marginRight = "10px";
        elem.style.display = "inline-block";
        appState.resetButton.parentElement.insertBefore(elem, recalcButton);
    }

    elem.innerHTML = text;
    elem.style.color = color;
}

function clearInputErrorContainers() {
    clearInputErrorContainer(0);
    clearInputErrorContainer(1);
}

function isMathFieldRendered(element) {

  if ( !element ) {
    return false;
  }

  const shadowRoot = element.shadowRoot;
  if (!shadowRoot) {
    return false;
  }

  var hasInternalStructure = shadowRoot.querySelector('.ML__base') !== null;

  return hasInternalStructure;
}

function errorHtml(errorMessage) {
    return "<div class='output_error'>"+errorMessage+"</div>";
}

function noteHtml(errorMessage) {
    return "<div class='output_note'>"+errorMessage+"</div>";
}

function isMathDivOrSpan(element) {
    return element.tagName == "MATH-DIV" || element.tagName == "MATH-SPAN";
}

function getLatexFromElement(element) {
    if ( isMathDivOrSpan(element) ) {
        return element.textContent;
    } else if ( element.tagName == "TEXTAREA" || element.tagName == "MATH-FIELD" ) {
        return element.value;
    } else {
        return "";
    }
}

function renderAllMathLive(element = document) {

    function renderAllMathLive_inner(ith) {
        if ( ith.getAttribute("already_explicitly_rendered") == "true" ) {
            // skip to improve performance.
        } else {
            if ( isMathFieldRendered(ith) ) {
                // ignore.
            } else {
                ith.render();
                ith.setAttribute("already_explicitly_rendered", "true");
            }
        }
    }
    
    if ( isMathDivOrSpan(element) ) {
        renderAllMathLive_inner(element);
    } else {
        var allSpans = Array.from(element.querySelectorAll("math-span"));
        var allDivs = Array.from(element.querySelectorAll("math-div"));
        var allMathLiveElems = [...allSpans, ...allDivs];

        allMathLiveElems.forEach(ith => {
            renderAllMathLive_inner(ith);
        });
    }
}

function clearInputErrorContainer(ordinal) {
    var resultsElem = ordinal === 0 ?  appState.results0 : appState.results1;

    if ( resultsElem ) {
        var resultsElemUl = resultsElem.getElementsByTagName("ul")[0];

        resultsElemUl.innerHTML = "";
    }
}

function addToOutput(resultsElem, string) {
    var resultsElemUl = resultsElem.getElementsByTagName("ul")[0];

    if ( resultsElemUl.innerHTML ) {
        // resultsElem.innerHTML += "<br>" + string;
        resultsElemUl.innerHTML += "<li>" + string + "</li>"
    } else {
        resultsElemUl.innerHTML = "<li>" + string + "</li>"
    }
}

function getModuleNameHtml(ordinal) {
    if ( ordinal === 0 ) {
        return "<math-span>\\Phi</math-span>";
    } else {
        return "<math-span>\\Psi</math-span>";
    }
}

function fitTextSimple(elem) {
    var fontSizeString = window.getComputedStyle(elem).getPropertyValue('font-size');
    var fontSizeInPixels = parseFloat(fontSizeString);
    fitText(elem, true, fontSizeInPixels);
}

function fitText(el, largeAsPossible, fontSize = R_MATRIX_CELL_SIZE_PX) {
  const parent = el.parentElement;
  
  var baseFontSize = largeAsPossible ? fontSize : fontSize/2;
  // Reset font size to measure the native un-clamped text width
  el.style.fontSize = baseFontSize+'px'; 
  
  const maxAvailableWidth = parent.clientWidth - parseFloat(getComputedStyle(parent).paddingLeft) * 2;
  const currentTextWidth = el.offsetWidth;

  const maxAvailableHeight = parent.clientHeight - parseFloat(getComputedStyle(parent).paddingTop) * 2;
  const currentTextHeight = el.offsetHeight;
  
  // Calculate the scaling ratio
  const widthRatio = maxAvailableWidth / currentTextWidth;
  const heightRatio = maxAvailableHeight / currentTextHeight;

  const ratio = Math.min(widthRatio, heightRatio);

//   console.log("ratio: " + ratio);
  
  // Set the precise dynamic font size (bounded by 100px max)
  const idealSize = Math.min(baseFontSize, Math.floor(baseFontSize * ratio));
  el.style.fontSize = `${idealSize}px`;
}

function shortenRawSeconds(seconds) {
    return seconds.toFixed(4);
}

function applyTippy(element, content) {
    tippy(element, {
        content: content,
        allowHTML: true,
        maxWidth: 'none',
        onShow(instance) {
            // 1. Get the container element or reference
            const box = instance.popper.querySelector('.tippy-content');
            renderAllMathLive(box);
        },
        theme: 'tooltip' // Use the custom class here,
    })
}

function getMatrixNameHtml(matrixOrdinal) {
    return "<math-span>\\tau^{"+matrixOrdinal+"}</math-span>";
}

function onMathElemClicked(mathLatex) {

    showOkDialog(mathLatex);

    var dialogElem = document.getElementById("yes-no-dialog");

    if ( dialogElem ) {
        dialogElem.style.maxWidth = "90vw";
        // dialogElem.style.overflow = "scroll";
    }
}

function hideDialog(dialogElem) {
    if ( dialogElem && dialogElem.parentElement ) {
        dialogElem.parentElement.removeChild(dialogElem);
    }

    document.body.style.overflow = "scroll";
}

function showOkDialog(message) {
    showDialog(message, "OK");
}

function prepExampleLatexForDisplay(exampleLatex) {
    // exampleLatex = exampleLatex.replaceAll("\\", "");

    return exampleLatex;
}

function prepExampleLatexForChoosing(exampleLatex) {
    if ( exampleLatex.startsWith("\\Phi_\\tau = ") ) {
        exampleLatex = exampleLatex.replace("\\Phi_\\tau = ", "");
    }

    if ( exampleLatex.startsWith("\\Psi_\\tau = ") ) {
        exampleLatex = exampleLatex.replace("\\Psi_\\tau = ", "");
    }

    return exampleLatex;
}

function setUpLatexCopyButton(latexElem, addCopiedMarkToEnd = false) {
    var copyButtonElem = document.getElementById(latexElem.id + "-copy");

    var originalButtonContent = copyButtonElem.textContent;

    copyButtonElem.addEventListener("click", function() {
        
        var latex = getLatexFromElement(latexElem);

        debugger;

        if ( latexElem.id.includes(R_AND_S_CONFIG.base_id) ) {
            var textAreaIdForS = S_CONFIG.base_id + "-latex-textarea";

            var latexElemForS = document.getElementById(textAreaIdForS);

            if ( latexElemForS ) {
                latex += LATEX_EXPORT_NEWLINE + getLatexFromElement(latexElemForS);
            }
        }

        navigator.clipboard.writeText(latex);
        if ( addCopiedMarkToEnd ) {
            copyButtonElem.textContent = originalButtonContent + ' (Copied!)';
        } else {
            copyButtonElem.textContent = 'Copied!';
        }
        
        setTimeout(() => {
            copyButtonElem.textContent = originalButtonContent;
        }, 1500);
    });
}

function arrayToLatexMatrix(arr, type = DEFAULT_LATEX_MATRIX_ENVIRONMENT) {
    
  const rows = arr.map(row => row.join(' & '));
  
  const matrixBody = rows.join(' \\\\\n');
  
  return `\\begin{${type}}\n${matrixBody}\n\\end{${type}}`;
}

function arrayToLatexForBigMatrix(arr) {
  if (!arr.length || !arr[0].length) return '';
  
  const colCount = arr[0].length;
  // Creates a centering format like cccccc... for all columns
  var colFormat = "";

  for ( var i = 0; i < colCount; i++ ) {
    colFormat += "c";
  }
  
  const rows = arr.map(row => row.join(' & '));
  
  // Use \left[ and \right] to dynamically stretch square brackets around the array
  var toReturn = `\\left(\\begin{array}{${colFormat}}\n  ${rows.join(' \\\\\n  ')}\n\\end{array}\\right)`;

  toReturn = toReturn.replaceAll("\\textcolor{orange}", "");

  return toReturn;
}

var LATEX_EXPORT_NEWLINE = "\\newline\r\n\r\n\r\n";

function generateAllPhiPsiLatexForExport() {
    var phiExpandedLatexElem = document.getElementById("latex-0");
    var psiExpandedLatexElem = document.getElementById("latex-1");
    var phiSimplifiedLatexElem = document.getElementById("equation-0-simplified");
    var psiSimplifiedLatexElem = document.getElementById("equation-1-simplified");

    var toReturn = "";

    // toReturn += "\\documentclass[12pt]{article}" + newline;
    // toReturn += "\\begin{document}" + newline;

    // toReturn += "\\boldsymbol{\\Phi_{expanded}} = " + phiExpandedLatexElem.value + newline;
    // toReturn += "\\boldsymbol{\\Phi_{simple}} = " + phiSimplifiedLatexElem.value + newline;
    // toReturn += "\\boldsymbol{\\Psi_{expanded}} = " + psiExpandedLatexElem.value + newline;
    // toReturn += "\\boldsymbol{\\Psi_{simple}} = " + psiSimplifiedLatexElem.value + newline;

    toReturn += "\\Phi_{input} = " + phiExpandedLatexElem.value + LATEX_EXPORT_NEWLINE;
    toReturn += "\\Phi_{simple} " + phiSimplifiedLatexElem.textContent + LATEX_EXPORT_NEWLINE; // "=" is part of the latex
    toReturn += "\\Psi_{input} = " + psiExpandedLatexElem.value + LATEX_EXPORT_NEWLINE;
    toReturn += "\\Psi_{simple} " + psiSimplifiedLatexElem.textContent + LATEX_EXPORT_NEWLINE; // "=" is part of the latex


    if ( isOutputVisible() ) {

        uLatex = appState.results3.getAttribute("latex");
        uInverseLatex = appState.results4.getAttribute("latex");

        if ( uLatex ) {
            toReturn += "U = " + uLatex + LATEX_EXPORT_NEWLINE;
        }

        if ( uInverseLatex ) {
            toReturn += "U^{-1} = " + uInverseLatex + LATEX_EXPORT_NEWLINE;
        }

        OUTPUT_MATRIX_CONFIGS.forEach(ith => {
            var latexElem = document.getElementById(ith.base_id + "-latex-textarea");

            if ( latexElem ) {
                toReturn += ith.title + " = " + latexElem.textContent + LATEX_EXPORT_NEWLINE;
            }
        });
    }

    // toReturn += "\\end{document}" + newline;

    return toReturn;
}

function scrollMathFieldToEnd(mathField) {
    setTimeout(() => {
        const container = mathField.shadowRoot?.querySelector('.ML__content');
        if (container) {
            container.scrollLeft = container.scrollWidth;
        }
    }, 0);
}

function exportLatexToTex(latexText, filename) {
    
  const blob = new Blob([latexText], { type: 'application/x-tex;charset=utf-8' });
  
  // Create a temporary anchor element
  const link = document.createElement('a');
  link.href = URL.createObjectURL(blob);
  link.download = filename;
  
  // Trigger download and clean up
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(link.href);
}

function generateRMatrixLatexHtml(baseIdName, latex = "") {
    var toReturn = '';

    toReturn += '<div class="textarea-container">';
    // toReturn += '    <div style="position:unset; text-align:right; display:block;" class="button-group">';
    // toReturn += '        <button id="'+(baseIdName+'-textarea-copy')+'" type="button">Copy As LaTeX</button>';
    // toReturn += '    </div>';
    toReturn += '    <math-span max-matrix-cols="'+MAX_INPUT_MATRIX_DIMENSION+'" id="'+(baseIdName+'-textarea')+'" class="r_matrix_rendered">'+latex+'</math-span>';
    // toReturn += '    <math-div id="'+(baseIdName+'-textarea')+'">'+latex+'</math-div>';
    
    toReturn += '</div>';

    return toReturn;
}

function addRMatrixLatexToDom(wrapperElem, array2d ) {
    var baseIdName = wrapperElem.id;
    var textAreaId = baseIdName + "-textarea";

    // var latex = arrayToLatexForBigMatrix(arr);
    
    // latex = "blah";
    var rMatrixLatexHtml = generateRMatrixLatexHtml(baseIdName, "");

    wrapperElem.innerHTML = rMatrixLatexHtml;
    
    var textArea = document.getElementById(textAreaId);

    textArea.maxMatrixCols = MAX_INPUT_MATRIX_DIMENSION;

    var matrixLatex = "";

    if ( array2d.length > MAX_INPUT_MATRIX_DIMENSION ) {
        // Has rendering problem but will at least show something.
        matrixLatex = arrayToLatexForBigMatrix(array2d);
    } else {
        matrixLatex = arrayToLatexMatrix(array2d);
    }
    
    textArea.textContent = matrixLatex;
   

    // textArea.value = latex;

    setUpLatexCopyButton(textArea, true);
}

function showOptionsChooserDialog() {
    var messageHtml = "<div style='font-weight:bold; text-decoration:underline; margin-bottom:10px;'>Options</div>";

    var variant = getVariant();

    messageHtml += '<div style="text-align:left;">';

    for ( var i = 0; i < CHECKBOX_OPTIONS.length; i++ ) {
        let ithCheckboxOption = CHECKBOX_OPTIONS[i];

        var borderTopNone = i != 0 ? "border-top:none;" : "";
        var borderTopNone = "";

        messageHtml += '<div class="checkbox_container">';
        messageHtml += '    <span class="module_name">'+(i+1)+'.</span>';
        messageHtml += '    <input type="checkbox" id="'+ithCheckboxOption.local_storage_key+'_id" name="'+ithCheckboxOption.local_storage_key+'_name">';
        messageHtml += '    <label for="'+ithCheckboxOption.local_storage_key+'_id">'+ithCheckboxOption.readable_label+'</label>';
        messageHtml += '</div>';

        // messageHtml += "<div style='"+borderTopNone+"' class='example_row'>";
        // messageHtml += "<math-div style='display:block; pointer-events: none;' id='example-"+i+"_0'></math-div>";
        // if ( ithExample[1] ) {
        //     messageHtml += "<math-div style='display:block; pointer-events: none;' id='example-"+i+"_1'></math-div>";
        // }
        // messageHtml += "</div>";
    }

    messageHtml += '</div">';
    
    showDialog(messageHtml, "OK");

    for ( var i = 0; i < CHECKBOX_OPTIONS.length; i++ ) {
        let ithCheckboxOption = CHECKBOX_OPTIONS[i];

        let checkboxElem = document.getElementById(ithCheckboxOption.local_storage_key+"_id");

        checkboxElem.checked = getCheckboxOptionValue(ithCheckboxOption.local_storage_key);

        checkboxElem.addEventListener("change", function() {
            setVisibleStatus(VISIBLE_STATUS__RESULTS_STALE);
            setCheckboxOptionValue(ithCheckboxOption.local_storage_key, checkboxElem.checked);
        });
        
        var toolTipHtml = "<div style='text-align:center; font-size:20px;'>"+ithCheckboxOption.tooltip_title+"</div>"
        applyTippy(checkboxElem.parentElement, toolTipHtml);
    }

    var dialogElemScrollContainer = document.getElementById("dialog-element-scroll-container");

    renderAllMathLive(dialogElemScrollContainer);
}

function fadeBorder(targetElem, isGreen, hadInputErrors) {
  targetElem.classList.remove('trigger-border-fade-red');
  targetElem.classList.remove('trigger-border-fade-orange');
  targetElem.classList.remove('trigger-border-fade-green');
  
  void targetElem.offsetWidth; 
  
  if ( isGreen ) {
    if ( hadInputErrors ) {
        targetElem.classList.add('trigger-border-fade-orange');
    } else {
        targetElem.classList.add('trigger-border-fade-green');
    }
  } else {
    targetElem.classList.add('trigger-border-fade-red')
  }
}

function hideTippyOnElem(element) {
    if ( element._tippy ) {
        element._tippy.hide();
    }
}

function showToastTippyOnElem(element, messageHtml) {
    hideTippyOnElem(element);

    const tip = tippy(element, {
        trigger: 'manual',
        allowHTML: true,
        theme: 'left-arrow',
        content: 'Initial error',
        placement: 'bottom-start',
        zIndex: 100
    });

    tip.setProps({ content: messageHtml });
    tip.show();
    
    setTimeout(() => {
        tip.hide();
    }, 1500);
}

function getTippyMessageHtml(plainTextMessage, cssColor) {
    var fontSize = 25;
    var messageHtml = "<div style='color:"+cssColor+"; font-size:"+fontSize+"px;'>"+plainTextMessage+"</math-span></div>";

    return messageHtml;
}

function showSimplifyOrExpansionStatus(inputElem, targetElem, isGreen, isSimplify, hadInputErrors, showTippy = true) {

    var elementToShowTippyOn = inputElem;
    
    var symbolForFullModule = targetElem == appState.equation0 ? "\\Phi_{t}" : "\\Psi_{t}"

    var message = "";

    if ( isGreen ) {
        if ( isSimplify ) {
            if ( hadInputErrors ) {
                message = "Simplified but with errors.";
            } else {
                message = "Successfully simplified.";
            }
        } else {
            message = "Successfully expanded.";
        }
        
    } else {
        if ( isSimplify ) {
            message = "Could not simplify.";
        } else {
            message = "Could not expand.";
        }
    }

    var cssColor = isGreen ? (hadInputErrors ? "darkorange" : "lightgreen") : "lightsalmon";
    var messageHtml = getTippyMessageHtml(message, cssColor);
    if ( showTippy ) {
        showToastTippyOnElem(elementToShowTippyOn, messageHtml);
    }

    fadeBorder(targetElem.parentElement, isGreen, hadInputErrors);
}

function showMatrixExceededMaxDimensionsError(equationElem) {
    var message = "Exceeded max matrix dimension: " + MAX_INPUT_MATRIX_DIMENSION;
    var cssColor = "lightsalmon";
    var messageHtml = getTippyMessageHtml(message, cssColor);
    
    showToastTippyOnElem(equationElem, messageHtml);
}

function showIntInputError(inputElem, message) {
    var cssColor = "lightsalmon";
    var messageHtml = getTippyMessageHtml(message, cssColor);
    
    showToastTippyOnElem(inputElem, messageHtml);
}

function generateMetricsTable(responseJson) {
    
    var metricsHtml = "";
    metricsHtml += "<table id='metrics-table'>";

    metricsHtml += "    <tr>";
    metricsHtml += "        <td>";
    metricsHtml += "        </td>";
    // metricsHtml += "        <td class='metrics_table_header'>";
    // metricsHtml += "            <div>Input Validation</div>";
    // metricsHtml += "        </td>";
    // metricsHtml += "        <td class='metrics_table_header'>";
    // metricsHtml += "            <div>Preprocessing</div>";
    // metricsHtml += "        </td>";
    // metricsHtml += "        <td class='metrics_table_header'>";
    // metricsHtml += "            <div>Rightmost Basis</div>";
    // metricsHtml += "        </td>";
    metricsHtml += "        <td class='metrics_table_header'>";
    metricsHtml += "            <math-div>\\text{Ext}^1</math-div>";
    metricsHtml += "        </td>";
    metricsHtml += "        <td class='metrics_table_header'>";
    metricsHtml += "            <math-div>\\text{Ext}^0</math-div>";
    metricsHtml += "        </td>";
    metricsHtml += "        <td class='metrics_table_header'>";
    metricsHtml += "            <math-div>K</math-div>";
    metricsHtml += "        </td>";
    metricsHtml += "        <td class='metrics_table_header'>";
    metricsHtml += "            <div>Total</div>";
    metricsHtml += "        </td>";
    metricsHtml += "    </tr>";

    // metricsHtml += "    <tr>";
    // metricsHtml += "        <td class='metrics_table_value_name'><div>Execution Time</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>"+shortenRawSeconds(responseJson.seconds_to_validate_input)+" s</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>"+shortenRawSeconds(responseJson.seconds_to_generate_preprocessing)+" s</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>"+shortenRawSeconds(responseJson.seconds_to_generate_rightmost_basis)+" s</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>"+shortenRawSeconds(responseJson.Ext_matrix.total_seconds_to_generate)+" s</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>"+shortenRawSeconds(responseJson.Ext0_matrix.total_seconds_to_generate)+" s</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>"+shortenRawSeconds(responseJson.K_matrix.total_seconds_to_generate)+" s</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>"+shortenRawSeconds(responseJson.total_seconds_to_generate)+" s</div></td>";
    // metricsHtml += "    </tr>";

    // metricsHtml += "    <tr>";
    // metricsHtml += "        <td class='metrics_table_value_name'><div>Reduction Count</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>N/A</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>N/A</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>N/A</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>"+responseJson.Ext_matrix.total_reduction_count+"</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>"+responseJson.Ext0_matrix.total_reduction_count+"</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>"+responseJson.K_matrix.total_reduction_count+"</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>"+responseJson.total_reduction_count+"</div></td>";
    // metricsHtml += "    </tr>";


    metricsHtml += "    <tr>";
    metricsHtml += "        <td class='metrics_table_value_name'><div>Leaf Count</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>N/A</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>N/A</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>N/A</div></td>";
    metricsHtml += "        <td class='metrics_table_value'><div>"+responseJson.Ext_matrix.leaf_count+"</div></td>";
    metricsHtml += "        <td class='metrics_table_value'><div>"+responseJson.Ext0_matrix.leaf_count+"</div></td>";
    metricsHtml += "        <td class='metrics_table_value'><div>"+responseJson.K_matrix.leaf_count+"</div></td>";
    metricsHtml += "        <td class='metrics_table_value'><div>"+responseJson.total_leaf_count+"</div></td>";
    metricsHtml += "    </tr>";

    metricsHtml += "    <tr>";
    metricsHtml += "        <td class='metrics_table_value_name'><div>Byte Count</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>N/A</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>N/A</div></td>";
    // metricsHtml += "        <td class='metrics_table_value'><div>N/A</div></td>";
    metricsHtml += "        <td class='metrics_table_value'><div>"+responseJson.Ext_matrix.byte_count+"</div></td>";
    metricsHtml += "        <td class='metrics_table_value'><div>"+responseJson.Ext0_matrix.byte_count+"</div></td>";
    metricsHtml += "        <td class='metrics_table_value'><div>"+responseJson.K_matrix.byte_count+"</div></td>";
    metricsHtml += "        <td class='metrics_table_value'><div>"+responseJson.total_byte_count+"</div></td>";
    metricsHtml += "    </tr>";

    metricsHtml += "</table>";


    // metricsHtml += "<tr><td style='text-align:right'>Is Cached Response = </td><td>"+cachedMessage+"</td></tr><tr><td style=''>Reduction Count</td><td>"+responseJson.total_reduction_count+", </td></tr><tr><td style='text-align:right'>Total Execution Time = </td><td>"+seconds+" seconds</td></tr></table>"
    return metricsHtml;
}
function showExportDialog() {
    var messageHtml = "<div style='font-weight:bold; margin-bottom:10px;'>Choose an Export Option</div>";

    var variant = getVariant();
    var examples = [];

    var baseCopyMessage = "Copy all LaTeX to Clipboard";
    
    if ( variant == VARIANT__PHI_PSI_EMBEDDING ) {
        messageHtml += "<div class='export_option_row_wrapper'>";
            messageHtml += "<div class='export_option_row'>";
            messageHtml += '<div style="display:inline-block; color:blue; cursor:pointer;" id="copy-all-button">'+baseCopyMessage+'</div>'
            messageHtml += "</div>";
            messageHtml += "<div class='export_option_row'>";
            messageHtml += '<div style="display:inline-block; color:blue; cursor:pointer;" id="export-all-button">Export all LaTeX to .tex file</div>'
            messageHtml += "</div>";
        messageHtml += "</div>";
    } else {
        messageHtml += "Not implemented";
    }

    showDialog(messageHtml, "Close");

    var copyAllButton = document.getElementById("copy-all-button");
    var exportAllButton = document.getElementById("export-all-button");

    copyAllButton.addEventListener("click", function() {
        var allLatex = generateAllPhiPsiLatexForExport();

        navigator.clipboard.writeText(allLatex);

        copyAllButton.textContent = baseCopyMessage + ' (Copied!)';
        setTimeout(() => {
            copyAllButton.textContent = baseCopyMessage;
        }, 1500);
    });

    exportAllButton.addEventListener("click", function() {
        var allLatex = generateAllPhiPsiLatexForExport();

        exportLatexToTex(allLatex, 'phi_psi_reduction_embedding.tex');
    });
}

function loadPhiPsiExample(exampleLatex0, exampleLatex1) {
    var dialogElemScrollContainer = document.getElementById("dialog-element-scroll-container");

    var callOnChange = false;

    loadLatexIntoEquationInput(0, prepExampleLatexForChoosing(exampleLatex0), callOnChange);
    
    if ( exampleLatex1 ) {
        loadLatexIntoEquationInput(1, prepExampleLatexForChoosing(exampleLatex1), callOnChange)
    }

    if ( dialogElemScrollContainer ) {
        hideDialog(dialogElemScrollContainer);
    }

    onInputChanged();
}

function showExampleChooserDialog() {
    var messageHtml = "<div style='font-weight:bold; margin-bottom:10px;'>Click an example to choose it as input.</div>";

    messageHtml += "<div style='margin-bottom:10px;'>(examples show how t-modules can be expanded or simplified)</div>"

    var variant = getVariant();
    var examples = [];
    
    if ( variant == VARIANT__PHI_PSI_EMBEDDING ) {
        examples = EXAMPLES__PHI_PSI_EMBEDDDING;
    } else if ( variant == VARIANT__RIGHTMOST_BASIS ) {
        examples = EXAMPLES__PHI_BASIS;
    }

    for ( var i = 0; i < examples.length; i++ ) {
        let ithExample = examples[i];

        var borderTopNone = i != 0 ? "border-top:none;" : "";
        var borderTopNone = "";

        messageHtml += "<div style='"+borderTopNone+"' class='example_row'>";
        messageHtml += "<math-div style='display:block; pointer-events: none;' id='example-"+i+"_0'></math-div>";
        if ( ithExample[1] ) {
            messageHtml += "<math-div style='display:block; pointer-events: none;' id='example-"+i+"_1'></math-div>";
        }
        messageHtml += "</div>";
    }
    
    showDialog(messageHtml, "Cancel");

    var exampleRows = Array.from(document.getElementsByClassName("example_row"));

    for ( var i = 0; i < examples.length; i++ ) {
        var ithExample = examples[i];
        var ithExampleRow = exampleRows[i];

        var exampleElem0 = document.getElementById("example-"+i+"_0");
        var exampleElem1 = document.getElementById("example-"+i+"_1");

        let exampleLatex0 = prepExampleLatexForDisplay(ithExample[0]);
        let exampleLatex1 = prepExampleLatexForDisplay(ithExample[1]);

        if ( variant == VARIANT__PHI_PSI_EMBEDDING ) {
            exampleLatex0 = "\\Phi_\\tau = " + exampleLatex0;
            exampleLatex1 = "\\Psi_\\tau = " + exampleLatex1;
        }

        if ( exampleElem0 ) {
            exampleElem0.innerHTML = exampleLatex0;
        }

        if ( exampleElem1 ) {
            exampleElem1.innerHTML = exampleLatex1;
        }

        ithExampleRow.addEventListener("click", function() {
            loadPhiPsiExample(exampleLatex0, exampleLatex1);
        });
    }

    var dialogElemScrollContainer = document.getElementById("dialog-element-scroll-container");

    renderAllMathLive(dialogElemScrollContainer);
}

function showDialog(message, noOrOkButtonText, yesButtonText = null, onYes = null, onNo = null, onPreNo = null, cancelable = true) {
    
    var dialogElemScrollContainer = document.createElement("div");

    dialogElemScrollContainer.id = "dialog-element-scroll-container";
    dialogElemScrollContainer.style.width = "100vw";
    dialogElemScrollContainer.style.height = "100vh";
    dialogElemScrollContainer.style.overflow = "hidden";
    dialogElemScrollContainer.style.position = "fixed";
    dialogElemScrollContainer.style.top = "0px";
    dialogElemScrollContainer.style.left = "0px";
    dialogElemScrollContainer.style.zIndex = "100";

    var dialogBackgroundElem = document.createElement("div");
    dialogBackgroundElem.style.position = "absolute";
    dialogBackgroundElem.style.top = "0px";
    dialogBackgroundElem.style.left = "0px";
    dialogBackgroundElem.style.width = "100vw";
    dialogBackgroundElem.style.height = "100vh";
    dialogBackgroundElem.id = "yes-no-dialog-background";
    dialogBackgroundElem.style.backgroundColor = "rgba(0, 0, 0, .7)";

    var dialogElem = document.createElement("table");

    dialogElem.id = "yes-no-dialog-wrapper";

    var row = dialogElem.insertRow(-1);

    var dialogInnerHtml = '<td style="vertical-align:middle; text-align:center;"><div id="yes-no-dialog" style=" max-width: 90vw; position:relative;">';
    dialogInnerHtml += '        <div style="max-height:75vh; overflow:scroll;">';
    dialogInnerHtml += '            ' + message;
    dialogInnerHtml += '        </div>';
    dialogInnerHtml += '        <div style="margin-top:15px;">';
    dialogInnerHtml += '            <button class="large_font small_border_radius bordered general_button" style="margin:5px;" id="dialog-no-button" >'+noOrOkButtonText+'</button>';
    if ( yesButtonText ) {
        dialogInnerHtml += '            <button class="large_font small_border_radius bordered general_button" style="margin:5px;" id="dialog-yes-button">'+yesButtonText+'</button>';
    }
    dialogInnerHtml += '        </div>';
    dialogInnerHtml += '  </div></td>';

    row.innerHTML = dialogInnerHtml;

    dialogElemScrollContainer.appendChild(dialogElem);

    document.body.appendChild(dialogElemScrollContainer);

    document.body.style.overflow = "hidden";

    document.getElementById("yes-no-dialog").parentElement.insertBefore(dialogBackgroundElem, document.getElementById("yes-no-dialog"));

    if ( cancelable ) {
        dialogBackgroundElem.addEventListener("click", function() {
            hideDialog(dialogElemScrollContainer);
        });
    }
    
    document.getElementById("dialog-no-button").addEventListener("click", function() {
        if ( onPreNo ) {
            onPreNo();
        }

        hideDialog(dialogElemScrollContainer);

        if ( onNo ) {
            onNo();
        }
    });
    
    if ( yesButtonText ) {
        document.getElementById("dialog-yes-button").addEventListener("click", function() {
            hideDialog(dialogElemScrollContainer);
    
            if ( onYes ) {
                onYes();
            }
        });
    }
}


function generateRMatrixHtmlSkeleton(baseElemId, title) {
    var isLast = OUTPUT_MATRIX_CONFIGS[OUTPUT_MATRIX_CONFIGS.length-1].base_id == baseElemId;

    var toReturn = "";

    var titleStyle = isLast === true ? "display:none;" : "";

    toReturn += `
        <tr style='${titleStyle}'>
            <td colspan=6 style="white-space:nowrap; width:0%; text-align:left;">
                <span class="result_name" id="${baseElemId}-title">${title}</span>
            </td>
        </tr>`;
        
    toReturn += `<tr style='${titleStyle}'><td colspan=6 class="spacer_row_col_small"></td></tr>`;

    toReturn += `
        <tr>
            <td colspan=6 style="vertical-align:bottom; width:100%;">
                <div class="r_matrix" id="${baseElemId}">
                    
                </div>
            </td>
        </tr>

        <tr>
            <td colspan=6 style="vertical-align:bottom; width:100%;">
                <div class="r_matrix_latex_wrapper" id="${baseElemId}-latex">
                    
                </div>
            </td>
        </tr>

        <tr><td colspan=6 class="spacer_row_col"></td></tr>
    `;

    return toReturn;
}