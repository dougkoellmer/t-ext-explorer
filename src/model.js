


function recalculateOnBackEnd() {
    var phiModule = appState.equation0.value;
    var psiModule = appState.equation1.value;

    setVisibleStatus(VISIBLE_STATUS__CALCULATING);

    var xhr = runCalculationOnServer(phiModule, psiModule, function(responseJson) {

        if ( responseJson.success ) {
            setVisibleStatus(VISIBLE_STATUS__COMPLETE);
        } else {
            if (typeof responseJson === 'string') {
                if ( isResponseCancellation(responseJson) ) {
                    setVisibleStatus(VISIBLE_STATUS__CALCULATION_CANCELLED);
                    showErrorOutput(responseJson);
                } else if ( isResponseNetworkError(responseJson) ) {
                    setVisibleStatus(VISIBLE_STATUS__NETWORK_ERRORS);
                    showErrorOutput(responseJson);
                } else {
                    setVisibleStatus(VISIBLE_STATUS__CALCULATION_ERRORS);
                    showErrorOutput(responseJson);
                }
            } else {
                setVisibleStatus(VISIBLE_STATUS__CALCULATION_ERRORS);
                showErrorOutput(responseJson.error);
            }
            
            return;
        }

        var variant = getVariant();

        if ( variant == VARIANT__PHI_PSI_EMBEDDING ) {

            var intermediateOutputRelatedElems = Array.from(document.getElementsByClassName("intermediate_output_related"));

            if ( getCheckboxOptionValue("show_reduction_steps") === true ) {
                intermediateOutputRelatedElems.forEach(ith => {
                    ith.style.removeProperty("display");
                });

                // setBackEndResultsField(appState.results2, responseJson.aux_matrix_latex);

                var e_a_delimiter = ",\\hspace{5px}";
                setBackEndListResultsField(appState.results5, responseJson.e_matrices_latex, null, null, e_a_delimiter);
                setBackEndListResultsField(appState.resultsAMatrices, responseJson.a_matrices_latex, null, null, e_a_delimiter);
                setBackEndListResultsField(appState.resultsMultiplyByA, responseJson.multiply_by_a_latex, responseJson.max_degrees, responseJson.contains_constant_term_list);

                setBackEndResultsField(appState.unitTestResults, responseJson.unit_test_results);

                // setBackEndListResultsField(appState.resultsMaxDegrees, responseJson.max_degrees);

                setBackEndListResultsField(appState.resultsEMatricesByDegrees, responseJson.Ext_matrix, responseJson.V_matrix, responseJson.W_V_matrix);
            } else {
                intermediateOutputRelatedElems.forEach(ith => {
                    ith.style.display = "none";
                });
            }

            var rightmostBasisTitle = document.getElementById("rightmost-basis-title");

            if ( rightmostBasisTitle ) {
                var rightmostBasisSeconds = shortenRawSeconds(responseJson.seconds_to_generate_rightmost_basis);
                // rightmostBasisTitle.innerHTML = "Rightmost Basis (Execution Time = " + rightmostBasisSeconds + " seconds)"
            }
            
            setBackEndResultsField(appState.results3, responseJson.util_matrix_latex);
            setBackEndResultsField(appState.results4, responseJson.inverse_matrix_latex);

            appState.results3.setAttribute("latex", responseJson.util_matrix_latex);
            appState.results4.setAttribute("latex", responseJson.inverse_matrix_latex);


            var preprocessingTitleElem = document.getElementById("preprocessing-title");

            if ( preprocessingTitleElem ) {
                var secondsForPreprocessing = shortenRawSeconds(responseJson.seconds_to_generate_preprocessing);
                // preprocessingTitleElem.innerHTML = "Preprocessing (Execution Time = " + secondsForPreprocessing + " seconds)"
            }

            var outputInner = document.getElementById("output-inner");
            var metricsRow = document.getElementById("metrics-row");

            for ( var i = 0; i < OUTPUT_MATRIX_CONFIGS.length; i++ ) {
                var ithMatrixConfig = OUTPUT_MATRIX_CONFIGS[i];
                var elemId = ithMatrixConfig.base_id;

                if ( document.getElementById(elemId) ) {
                    // already in dom
                } else {
                    var matrixHtml = generateRMatrixHtmlSkeleton(elemId, ithMatrixConfig.title);

                    // outputInner.insertAdjacentHTML("beforeend", matrixHtml);
                    metricsRow?.insertAdjacentHTML("beforebegin", matrixHtml);
                }

                var isDivTrue = true;
                
                var rMatrix = document.getElementById(elemId)

                setBackEndResultsRMatrix(rMatrix, responseJson[ithMatrixConfig.response_key], isDivTrue);

                resizeRMatrix(rMatrix)
            }

            var metricsHtml = generateMetricsTable(responseJson);
            setBackEndResultsField(appState.metricsElem, metricsHtml, isDivTrue);

            var isCachedResult = responseJson.cached === true;
            var cachedMessage = isCachedResult ? true : false;
            var metricsTitleCachedElem = document.getElementById("metrics-title-cached_status");

            metricsTitleCachedElem.innerHTML = "(cached_response="+cachedMessage+")";

        } else if ( variant == VARIANT__DIAGONAL ) {
            // setBackEndResultsField(appState.reductionSteps, responseJson.final_reduction, isDivTrue);
            setBackEndResultsField(appState.reductionSteps, responseJson.reduction_steps, isDivTrue);
            setBackEndResultsField(appState.finalReduction, responseJson.final_reduction, isDivTrue);
            
        } else if ( variant == VARIANT__RIGHTMOST_BASIS ) {
            // setBackEndResultsField(appState.results2, responseJson.aux_matrix_latex);
            // setBackEndResultsField(appState.results3, responseJson.util_matrix_latex);
            // setBackEndResultsField(appState.results4, responseJson.inverse_matrix_latex);

            // setBackEndResultsField(appState.analysisSteps, responseJson.analysis_steps, isDivTrue);
            setBackEndResultsField(appState.columnLists, responseJson.column_lists, isDivTrue);
            setBackEndResultsField(appState.sortedColumnLists, responseJson.sorted_column_lists, isDivTrue);
            setBackEndResultsField(appState.flagList, responseJson.phi_flag_list, isDivTrue);
            setBackEndResultsField(appState.rightmostBasis, responseJson.rightmost_basis, isDivTrue);

            var seconds = shortenRawSeconds(responseJson.total_seconds_to_generate);
            var isCachedResult = responseJson.cached === true;
            var cachedMessage = isCachedResult ? true : false;
            var metricsHtml = "<table><tr><td style='text-align:right'>Is Cached Response = </td><td>"+cachedMessage+"</td></tr><tr><td style='text-align:right'>Total Execution Time = </td><td>"+seconds+" seconds</td></tr></table>"
            setBackEndResultsField(appState.metricsElem, metricsHtml, isDivTrue);

        }

        if ( FEATURE_FLAG__RENDER_ALL_MATH_OUTPUT_UP_FRONT === true ) {
            renderAllMathLive();
        }
        // setTimeout(resizeRMatrix);


        // setBackEndResultsField(appState.resultsAMatrices, responseJson.psi_module_condensed_latex);
    });

    showLoadingModal(xhr);
}

function deepClone(object) {
    return JSON.parse(JSON.stringify(object));
}

function reversePddMatrix(matrix) {
    matrix = deepClone(matrix);
    var transpose = (mat) => mat[0].map((_, colIndex) => mat.map(row => row[colIndex]));
    matrix = transpose(matrix);

    return matrix;
}

function getWidthWithoutPadding(element) {
  const style = window.getComputedStyle(element);
  const paddingLeft = parseFloat(style.paddingLeft);
  const paddingRight = parseFloat(style.paddingRight);
  
  return element.clientWidth - paddingLeft - paddingRight;
}


// function resizeRMatrix(rMatrixElem) {
//     var allRows = Array.from(rMatrixElem.getElementsByTagName("tr"));

//     var cellSize = "40px"
//     allRows.forEach(ithRow => {
//         var ithRowCells = Array.from(ithRow.getElementsByTagName("td"));

//         ithRowCells.forEach(ithCell => {
//             ithCell.style.width = cellSize;
//             ithCell.style.height = cellSize;

//             if ( ithCell.firstChild ) {
//                 ithCell.firstChild.style.width = cellSize;
//                 ithCell.firstChild.style.height = cellSize;
//             }
//         });
//     });
// }


function resizeRMatrix(rMatrixElem) {
    renderAllMathLive(rMatrixElem)
    var allRows = Array.from(rMatrixElem.getElementsByTagName("tr"));

    if ( allRows.length == 0 ) {
        return;
    }

    var largestCellWidthOrHeight = 0;

    var firstRowCells = Array.from(allRows[0].getElementsByTagName("td"));

    var largestDimensionIsWidth = true;

    firstRowCells.forEach(ithCell => {
        var ithWidth = ithCell.clientWidth;
        var ithHeight = ithCell.clientHeight;

        if ( ithWidth > largestCellWidthOrHeight ) {
            largestCellWidthOrHeight = ithWidth;
            largestDimensionIsWidth = true;
        }

        if ( ithHeight > largestCellWidthOrHeight ) {
            largestCellWidthOrHeight = ithHeight;
            largestDimensionIsWidth = false;
        }
    });

    var cellSize = R_MATRIX_CELL_SIZE_PX;

    allRows.forEach(ithRow => {
        var ithRowCells = Array.from(ithRow.getElementsByTagName("td"));

        ithRowCells.forEach(ithCell => {
            ithCell.style.width = cellSize + "px";
            ithCell.style.height = cellSize + "px";

            if ( ithCell.firstChild ) {
                ithCell.firstChild.style.width = cellSize + "px";
                ithCell.firstChild.style.height = cellSize + "px";
            }
        });
    });

    var rByRTable = rMatrixElem.firstChild;

    // rByRTable.style.tableLayout = "fixed";

    var minimumCellSize = 84;
    var marginBottom = 0; //"40px";

    var tableWidth = rByRTable.clientWidth;
    var tableHeight = rByRTable.clientHeight;
    // var 
    var maxDimension = Math.min(window.innerHeight, getWidthWithoutPadding(document.getElementById("output")));
    // var maxDimension = getWidthWithoutPadding(document.getElementById("output"));
    var scaling = 1;//minimumCellSize / largestCellWidthOrHeight;

    // if ( tableWidth > maxDimension ) {
    //     scaling = maxDimension / tableWidth;
    // } else {
    //     if ( largestCellWidthOrHeight < minimumCellSize ) {
    //          scaling = minimumCellSize / largestCellWidthOrHeight;
    //     } else {
    //         scaling = 1;
    //     }
    // }

    // var redraw = rMatrixElem.clientWidth; // force redraw?

    var additionalMargin = (tableHeight*scaling) - tableHeight;
    // marginBottom *= scaling;
    marginBottom += additionalMargin;

    // rMatrixElem.style.transform = "scale("+scaling+", "+scaling+")";
    // rMatrixElem.style.transformOrigin = "top left";
    rMatrixElem.style.marginBottom = marginBottom + "px";

    // console.log("scaling=" + scaling);

    var zeroBasedRowIndex = 0;
    allRows.forEach(ithRow => {
        var ithRowCells = Array.from(ithRow.getElementsByTagName("td"));

        var zeroBasedColIndex = 0;

        ithRowCells.forEach(ithCell => {
            if ( ithCell.firstChild ) {

                var firstChildWidth = ithCell.firstChild.firstChild.offsetWidth * scaling;
                var inner = null;
                var applyToolTip = true;
                
                if ( ithCell.firstChild.classList.contains("r_by_r_label") ) {
                    applyToolTip = false;
                    inner = ithCell.firstChild.firstChild;
                    if ( inner ) {
                        firstChildWidth = Math.max(inner.clientWidth, inner.clientHeight);
                    }
                } else if ( isMathDivOrSpan(ithCell.firstChild.firstChild) ) {
                    inner = ithCell.firstChild.firstChild
                    var renderElem = getMathFieldRenderElem(inner);
                    if ( renderElem ) {
                        firstChildWidth = Math.max(renderElem.clientWidth, renderElem.clientHeight);
                    }
                } else if ( ithCell.firstChild.classList.contains("r_by_r_zero") ) {
                    inner = ithCell.firstChild.firstChild;
                    if ( inner ) {
                        firstChildWidth = Math.max(inner.clientWidth, inner.clientHeight);
                    }
                }

                var largeAsPossible = false;

                if ( applyToolTip ) {
                    largeAsPossible = true;
                    applyTippy(ithCell, ithCell.firstChild.getAttribute("tooltip-title"));
                }

                // var fontSize = largestCellWidthOrHeight * (firstChildWidth/largestCellWidthOrHeight)
                // console.log(cqi);
                // cqi = cqi > 6 ? 6 : cqi;

                // inner.style.fontSize = fontSize+"px";

                fitText(inner, largeAsPossible);
            }

            zeroBasedColIndex++;
        });

        zeroBasedRowIndex++;
    });

    var defaultBorderSize = 3;

    allRows.forEach(ithRow => {
        var ithRowCells = Array.from(ithRow.getElementsByTagName("td"));

        ithRowCells.forEach(ithCell => {
            ithCell.style.border = (defaultBorderSize / scaling) + "px solid black";
        });
    });

    var latexWrapper = document.getElementById(rMatrixElem.id + "-latex");
    var copyButton = document.getElementById(rMatrixElem.id + "-latex-textarea-copy");
    var rendereredLatex = document.getElementById(rMatrixElem.id + "-latex-textarea");
    renderAllMathLive(rendereredLatex);

    setTimeout(function() {
        var rendereredLatexInnerRender = getMathFieldRenderElem(rendereredLatex);
        latexWrapper.style.width = rendereredLatexInnerRender.clientWidth + "px";
    }, 5);
}

function isErrorMessage(responseMessage) {
    return typeof responseMessage === 'string' && responseMessage.startsWith("Error");
}

function isNoteMessage(responseMessage) {
    return typeof responseMessage === 'string' && responseMessage.startsWith("Note");
}

function setBackEndListResultsField(resultsElem, responseMessage, relatedResponseMessage = null, anotherRelatedResponseMessage = null, delimiter = " ") {
    if ( isErrorMessage(responseMessage) ) {
        setBackEndResultsField(resultsElem, responseMessage);
    } else {
        if ( resultsElem == appState.resultsMultiplyByA ) {
            var SPACER_ROW = '<tr><td class="spacer_row_col"></td> <td style=""></td></tr>';
            // var html = "<table style='border: 2.1px solid black; border-radius:5px; padding:5px;'>";
            var html = "<table style='border-radius:5px; padding:5px;'>";
            var colCount = appState.dForNInput.value;
            for ( var i = 0; i < responseMessage.length; i++ ) {
                var ithMatrix = responseMessage[i];
                var ithDegree = relatedResponseMessage[i];
                var ithContainsConstantTerm = anotherRelatedResponseMessage[i];

                var nRow = Math.floor(i / colCount) + 1;
                var nCol = i % colCount + 1;

                var constantTermMessageStart = "contains constant term = ";
                var constantTermMessage = ithContainsConstantTerm ? constantTermMessageStart + "TRUE" : constantTermMessageStart + "FALSE";
                var nEquals = "<math-span style='border-bottom:1px solid black; font-size:125%'>$N_{"+nRow+","+nCol+"}$ = " + ithDegree + ",\\hspace{5px}\\text{"+constantTermMessage+"}</math-span><div style='height:10px;'></div>";
                
                var ithHtml = "<tr><td><span class='ordinal'>"+(i+1)+".)&nbsp;</span></td><td>"+nEquals+"<math-div>" + ithMatrix + "</math-div></td></tr>";
                
                if ( i != responseMessage.length-1) {
                    ithHtml += SPACER_ROW;
                }
                
                html += ithHtml;
            }
            html += "</table>";

            resultsElem.innerHTML = html;
        } else if ( resultsElem == appState.resultsEMatricesByDegrees ) {
            var html = "<table style='border-collapse:collapse;'>";

            var vrWidthSmall = 2;
            var vrWidthThick = 3.5;
            var generalSpacing = 15;

            var thickVerticalRuleStyle = "border-left: "+vrWidthThick+"px solid black;";

            function spacerRowBreakdownCols(isForW, includeVr, bigSpace) {
                var verticalRuleThick = isForW ? "" : thickVerticalRuleStyle;
                var verticalRuleThin = includeVr ? 'border-left:'+vrWidthSmall+'px solid black;' : "";
                var spacerRowClass = bigSpace ? "spacer_row_col" : "spacer_row_col_small";
                return '<td style="'+verticalRuleThick+'" class="'+spacerRowClass+'"></td><td style=""></td> <td style="'+verticalRuleThin+'"></td> <td></td> <td style="'+verticalRuleThin+'"></td> <td></td>';
            }

            var SPACER_ROW = '<tr>'+spacerRowBreakdownCols(true, true, true)+spacerRowBreakdownCols(false, true, true)+spacerRowBreakdownCols(false, true, true)+'</tr>';
            var SPACER_ROW_NO_VR = '<tr>'+spacerRowBreakdownCols(true, false, false)+spacerRowBreakdownCols(false, false, false)+spacerRowBreakdownCols(false, false, false)+'</tr>';

            var nbsp = "";

            html += "<tr>";
            html += "<td style='border-bottom: 2px solid black; text-align:center' colspan=6><div style='width:unset; margin-bottom:10.932px; display:inline-block; padding-bottom:0px; border-bottom:none;' class='result_name'>For Ext</div></td>";
            html += "<td style='border-bottom: 2px solid black; border-left: "+vrWidthThick+"px solid black; text-align:center' colspan=6><div style='width:unset; margin-bottom:10.932px; display:inline-block; padding-bottom:0px; border-bottom:none;' class='result_name'>For Ext0</div></td>";
            html += "<td style='border-bottom: 2px solid black; border-left: "+vrWidthThick+"px solid black; text-align:center' colspan=6><div style='width:unset; margin-bottom:10.932px; display:inline-block; padding-bottom:0px; border-bottom:none;' class='result_name'>For Embedding</div></td>";
            html += "</tr>";

            html += SPACER_ROW_NO_VR;

            function rMatrixBreakdownHeaders(matrixName) {
                var isForW = matrixName == MATRIX_NAME__W;
                var isForWV = matrixName == MATRIX_NAME__WV;

                var verticalRule = isForW ? "" : thickVerticalRuleStyle;
                var toReturn = "";

                // var firstCol = isForWV ? "$E_{i,j}$\\tau^k (W)" : "$E_{i,j}$\\tau^k";
                // var secondCol = isForWV ? "\\Psi_t($E_{i,j}$c\\tau^k) (V)" : "\\Psi_t($E_{i,j}$c\\tau^k)";

                var firstCol = isForWV ? "$E_{i,j}$\\tau^k (W)" : "$E_{i,j}$\\tau^k";
                var secondCol = isForWV ? "$E_{i,j}$\\tau^k (V)" : "\\Psi_t($E_{i,j}$c\\tau^k)";
                toReturn += "<td style='"+verticalRule+"text-align:center' colspan=2><math-div style='width:unset; display:inline-block;' class='result_name result_name_header'>"+firstCol+"</math-div></td>";
                toReturn += "<td style='text-align:center' colspan=2><math-div style='width:unset; display:inline-block;' class='result_name result_name_header'>"+secondCol+"</math-div></td>";
                toReturn += "<td style='text-align:center' colspan=2><math-div style='width:unset; margin-left:"+generalSpacing+"px; margin-right:"+generalSpacing+"px; display:inline-block;' class='result_name result_name_header'>Reduction</math-div></td>";

                return toReturn;
            }

            var numberOfCols = parseInt(appState.dForMInput.value);

            function rMatrixColumns(rMatrixJson, matrixName) {
                var currentOrdinal = 1;

                var isForW = matrixName == MATRIX_NAME__W;
                var isForV = matrixName == MATRIX_NAME__V;
                var isForWV = matrixName == MATRIX_NAME__WV;

                var toReturn = [];

                var ithListOfMatrices = rMatrixJson.Eij_tau_k_matrices_latex;
                var ithListOfMatricesByPsi = rMatrixJson.Psi_Eij_c_tau_k_matrices_latex;
                var ithListOfReductionMatrices = rMatrixJson.reduction_matrices_latex;
                var ithListOfFinalReductionMatrices = rMatrixJson.final_reduction_matrices_latex;
                var ithListOfReductionCounts = rMatrixJson.reduction_counts;

                for ( var k = 0; k < ithListOfMatrices.length; k++ ) {
                    var kthMatrixAsLatex = ithListOfMatrices[k];
                    var kthMatrixByPsiLatex = ithListOfMatricesByPsi[k];

                    var kthReductionArray = ithListOfReductionMatrices[k];
                    var kthFinalReduction = ithListOfFinalReductionMatrices[k];
                    var kthReductionCount = ithListOfReductionCounts[k];

                    var kthReductionArrayAsLatex = "";
                    var dataReductionLatex = ""
                    if ( kthReductionArray ) {
                        for ( var m = 0; m < kthReductionArray.length; m++ ) {
                            if ( kthReductionArray[m].startsWith("<") ) {
                                kthReductionArrayAsLatex += kthReductionArray[m]; // plain html
                            } else {
                                var leftMargin = "10px";
                                if ( kthReductionArray[m].startsWith("=") ) {
                                    leftMargin = "20px";
                                }
                                kthReductionArrayAsLatex += "<math-div style=\"text-align:left; margin-left:"+leftMargin+";\">" + kthReductionArray[m] + "</math-div>";
                            }
                        }
                    }

                    var fullMathSpanForFinalReduction = "";

                    var clickableReductionLatex = "";

                    if ( kthReductionCount > 0  ) {
                        var reductionCssClass = "clickable_math_elem";

                        var showReductionInline = false;//isForWV;

                        if ( showReductionInline === true ) {
                            fullMathSpanForFinalReduction = "<math-span class='final_reduction "+reductionCssClass+"' data-reduction-latex='"+kthReductionArrayAsLatex+"'>" + kthReductionArrayAsLatex + "</math-span>";
                        } else {
                            clickableReductionLatex = kthFinalReduction;
                            dataReductionLatex = kthReductionArrayAsLatex ? kthReductionArrayAsLatex : "Reduction Count = " + kthReductionCount;
                            fullMathSpanForFinalReduction = "<math-span class='final_reduction "+reductionCssClass+"' data-reduction-latex='"+dataReductionLatex+"'>" + kthFinalReduction + "</math-span>";
                        }
                    } else {
                        fullMathSpanForFinalReduction = "<span>Not Required</span";
                    }

                    var row_Eij = kthMatrixAsLatex.effective_zero_based_row+1;
                    var col_Eij = kthMatrixAsLatex.effective_zero_based_col+1;
                    var degreeK_Eij = kthMatrixAsLatex.degree_of_tau;

                    var row_reduction = kthMatrixByPsiLatex ? kthMatrixByPsiLatex.effective_zero_based_row+1 : -1;
                    var col_reduction = kthMatrixByPsiLatex ? kthMatrixByPsiLatex.effective_zero_based_col+1 : -1;
                    var degreeK_reduction = kthMatrixByPsiLatex ? kthMatrixByPsiLatex.degree_of_tau : -1;

                    var ithLabelBase = "$E_{"+row_Eij+","+col_Eij+"}$\\tau^";

                    var ithLabelByPsiBase = isForWV ? "$E_{"+row_reduction+","+col_reduction+"}$\\tau^" : "\\Psi_t ( $E_{"+row_reduction+","+col_reduction+"}$c\\tau^";

                    var oneBasedOrdinal = currentOrdinal;
                    currentOrdinal++;
                    var ithLabel = ithLabelBase + degreeK_Eij + " = "
                    var ithLabelByPsi = isForWV ? ithLabelByPsiBase + degreeK_reduction + "=" : "";//ithLabelByPsiBase + degreeK_reduction + ") = "
                    var ithLabelForLinearCombination = "";
                    var ithHtml = "";
                    var verticalRule = isForW ? "" : thickVerticalRuleStyle;
                    var kthMatrixByPsiLatex_leftSide = kthMatrixByPsiLatex ? "<math-span>"+ithLabelByPsi+"</math-span>" : "";
                    var kthMatrixByPsiLatex_raw = kthMatrixByPsiLatex ? "<math-span>" + kthMatrixByPsiLatex.raw + "</math-span>" : "";
                    var listOrdinal = isForW ? "<span class='ordinal'>"+oneBasedOrdinal+".)&nbsp;</span>" : nbsp;
                    ithHtml += "<td style='padding-left: 10.12px;"+verticalRule+"'>"+listOrdinal+"<math-span>"+ithLabel+"</math-span></td><td><math-span>" + kthMatrixAsLatex.raw + "</math-span>"+nbsp+"</td>";
                    ithHtml += "<td style='padding-left: 10.12px;border-left:"+vrWidthSmall+"px solid black;'>"+nbsp+kthMatrixByPsiLatex_leftSide+"</td><td>"+kthMatrixByPsiLatex_raw+nbsp+"</td>";

                    var clickJs = clickableReductionLatex ? "onclick='javascript:window.onMathElemClicked(this.dataset.reductionLatex);'" : "";
                    var reductionClass = clickableReductionLatex ? "reduction_td" : "";
                    ithHtml += "<td style='padding-left: 10.12px;border-left:"+vrWidthSmall+"px solid black;'>" + nbsp +ithLabelForLinearCombination+"</td><td "+clickJs+" class='"+reductionClass+"' data-reduction-count='"+kthReductionCount+"' data-reduction-latex='"+dataReductionLatex+"'>"+fullMathSpanForFinalReduction+"</td>";
                    
                    toReturn.push(ithHtml);
                }

                // var currentIndexInRDimension = 0;

                // for ( var i = 0; i < rMatrixJson.degree_k_start_values.length; i++ ) {

                //     var ithDegreeKStartValue = rMatrixJson.degree_k_start_values[i];
                //     var ithDegreeKFinishValue = rMatrixJson.degree_k_finish_values[i];
                //     var ithSliceLength = ithDegreeKFinishValue - ithDegreeKStartValue;

                //     var sliceLimit = currentIndexInRDimension + ithSliceLength

                //     var ithListOfMatrices = rMatrixJson.Eij_tau_k_matrices_latex.slice(currentIndexInRDimension, sliceLimit);
                //     var ithListOfMatricesByPsi = rMatrixJson.Psi_Eij_c_tau_k_matrices_latex.slice(currentIndexInRDimension, sliceLimit);
                //     var ithListOfReductionMatrices = rMatrixJson.reduction_matrices_latex.slice(currentIndexInRDimension, sliceLimit);
                //     var ithListOfFinalReductionMatrices = rMatrixJson.final_reduction_matrices_latex.slice(currentIndexInRDimension, sliceLimit);

                //     currentIndexInRDimension += ithSliceLength;
                    
                //     var col = 1 + (i % numberOfCols);
                //     var row = 1 + (Math.floor(i / (numberOfCols)));

                //     var ithLabelBase = "$E_{"+row+","+col+"}$\\tau^";
                //     var ithLabelByPsiBase = "\\Psi_t ( $E_{"+row+","+col+"}$c\\tau^";

                    
                // }

                return toReturn;
            }

            html += "<tr>";
            html += rMatrixBreakdownHeaders(MATRIX_NAME__W);
            html += rMatrixBreakdownHeaders(MATRIX_NAME__V);
            html += rMatrixBreakdownHeaders(MATRIX_NAME__WV);
            html += "</tr>";

            html+= SPACER_ROW_NO_VR;

            var wMatrixCols = rMatrixColumns(responseMessage, MATRIX_NAME__W);
            var vMatrixCols = rMatrixColumns(relatedResponseMessage, MATRIX_NAME__V);
            var wvMatrixCols = rMatrixColumns(anotherRelatedResponseMessage, MATRIX_NAME__WV);


            for ( var i = 0; i < Math.max(wMatrixCols.length, vMatrixCols.length); i++ ) {
                var ithWMatrixCols = i < wMatrixCols.length ? wMatrixCols[i] : null;
                var ithVMatrixCols = i < vMatrixCols.length ? vMatrixCols[i] : null;
                var ithWVMatrixCols = i < wvMatrixCols.length ? wvMatrixCols[i] : null;

                var ithHtml = "<tr>";
                if ( ithWMatrixCols ) {
                    ithHtml += ithWMatrixCols;
                } else {
                    ithHtml += spacerRowBreakdownCols(true, true, true);
                }

                if ( ithVMatrixCols ) {
                    ithHtml += ithVMatrixCols;
                } else {
                    ithHtml += spacerRowBreakdownCols(false, true, true);
                }

                if ( ithWVMatrixCols ) {
                    ithHtml += ithWVMatrixCols;
                } else {
                    ithHtml += spacerRowBreakdownCols(false, true, true);
                }

                ithHtml += "</tr>";

                ithHtml += SPACER_ROW;
                html += ithHtml;
            }
            
            html += "</table>";

            resultsElem.innerHTML = html;

            var reductionCells = Array.from(resultsElem.getElementsByClassName("reduction_td"));
            reductionCells.forEach(ith => {
                if ( ith.dataset.reductionCount > 0 ) {
                    applyTippy(ith, "Reduction Count = " + ith.dataset.reductionCount);
                }
            });

            var finalReductions = Array.from(document.getElementsByClassName("final_reduction"));
            finalReductions.forEach(ith => {
                ith.addEventListener("click", function() {
                    var reductionLatex = ith.dataset.reductionLatex;

                    if ( reductionLatex ) {
                        // console.log(reductionLatex);
                    }
                });
            });
        } else {
            var as_latex = responseMessage.join(delimiter);
            setBackEndResultsField(resultsElem, as_latex);
        }
    }
}

function setBackEndResultsRMatrix(resultsElem, responseMessage, isDiv = false) {
    setBackEndResultsField(resultsElem, responseMessage, isDiv);
}

function setBackEndResultsField(resultsElem, responseMessage, isDiv = false) {
    if ( resultsElem ) {
        // all good
    } else {
        return;
    }

    if ( isErrorMessage(responseMessage) ) {
        resultsElem.innerHTML = errorHtml(responseMessage);
    } else {

        var innerHTML = "";

        var rMatrixConfig = null;
        var rMatrixConfigIndex = -1;

        rMatrixConfigIndex = OUTPUT_MATRIX_CONFIGS.findIndex(ith => {
            return ith.base_id == resultsElem.id
        });

        if ( rMatrixConfigIndex >= 0 ) {
            rMatrixConfig = OUTPUT_MATRIX_CONFIGS[rMatrixConfigIndex];
        }

        if ( resultsElem == appState.metricsElem ) {
            innerHTML = responseMessage;
        } else if ( resultsElem == appState.rightmostBasis ) {
            var readableOrdinal = 1;

            var ithAnalysisStep = responseMessage;
            
            var ithHtml = "";
            ithHtml += "<tr>";

            ithHtml += "<td style='font-size:150%; font-weight:800;'>"+readableOrdinal+":</td>";

            // ithHtml += "<td><table><tr>";
                // ithHtml += "<td>";

                var rowDimensionFirst = true;
                var matrixTermsLatex = "";
                for ( var k = 0; k < ithAnalysisStep.terms.length; k++ ) {
                    var kthMatrixTerm = ithAnalysisStep.terms[k];
                    var kthMatrixTermLatex = matrixToLatex(kthMatrixTerm, null, null, rowDimensionFirst);
                    matrixTermsLatex += kthMatrixTermLatex
                }

                var kthFlagMatricesU = ithAnalysisStep.flag_matrices_U;
                var kthFlagMatricesPhi = ithAnalysisStep.flag_matrices_Phi;
                var kthFlagMatricesSpanEqualities = ithAnalysisStep.flag_matrices_spans_are_equal;

                // var kthFlagLatex = "<table class=\"flag_breakdown_table\">";
                // kthFlagLatex += "<tr><td class=\"flag_breakdown_title_td flag_breakdown_td_right_border flag_breakdown_td_bottom_border\"><span class=\"flag_breakdown_title\">From U Matrix</span></td><td class=\"flag_breakdown_title_td flag_breakdown_td_right_border flag_breakdown_td_bottom_border\"><span class=\"flag_breakdown_title\">From Phi Module</span></td><td class=\"flag_breakdown_title_td flag_breakdown_td_bottom_border\"><span class=\"flag_breakdown_title\">Span Equality</span></td></tr>";
                // // var kthFlagLatex = "";
                // if ( kthFlagMatricesU ) {
                //     for ( var m = 0; m < kthFlagMatricesU.length; m++ ) {
                //         var areEqualSpans = kthFlagMatricesSpanEqualities[m];
                        
                //         var bottomBorderCssClass = m < kthFlagMatricesU.length -1 ? "flag_breakdown_td_bottom_border" : "";

                //         var equalityStatus = areEqualSpans ? "<span style=\"color:green\">EQUAL</span>" : "<span style=\"color:red\">NOT EQUAL</span>";
                //         var spanEqualityHtml = "<td class=\""+bottomBorderCssClass+"\">"+equalityStatus+"</td>";

                //         // var equalSign = areEqualSpans ? "=" : "\\neq";
                //         // var equalSignTd = "<td><math-div>"+equalSign+"</math-div></td>";

                //         // kthFlagLatex += "<math-div>" + kthFlagMatricesU[m] + "</math-div><math-div>"+equalSign+"</math-div><math-div>" +kthFlagMatricesPhi[m]+ "</math-div>";
                //         kthFlagLatex += "<tr><td class=\"flag_breakdown_td_right_border "+bottomBorderCssClass+"\" style=\"text-align:right;\"><math-div style=\"display:inline-block;\">" + kthFlagMatricesU[m] + "</math-div></td><td class=\"flag_breakdown_td_right_border "+bottomBorderCssClass+"\"><math-div>" +kthFlagMatricesPhi[m]+ "</math-div></td>"+spanEqualityHtml+"</tr>";
                //     }
                // }

                // kthFlagLatex += "</table>";

                // var ordinalHtml = "<span>"+ithAnalysisStep.ordinal+"</span>";
                // var cssClassForU = ithAnalysisStep.is_invertible ? "clickable_math_elem" : "";
                // var clickJs = ithAnalysisStep.is_invertible ? "onclick='javascript:window.onMathElemClicked(this.dataset.flagLatex);'" : "";
                var flagHtml = "";
                
                // if ( ithAnalysisStep.is_invertible ) {
                //     flagHtml += "<td data-flag-latex='"+kthFlagLatex+"' style='cursor:pointer;' "+clickJs+">";
                //     flagHtml += "<span class='"+cssClassForU+"'>" + (ithAnalysisStep.passes_flag_condition ? "Satisfies Flag" : "Fails Flag") + "</span>";
                //     flagHtml += "</td>";
                // } else {
                //     flagHtml += "<td></td>";
                // }
                
                // ithHtml += "<td><math-span id=''>"+matrixTermsLatex+"\\Rightarrow<math-span></td><td><math-span>" + ithAnalysisStep.U_matrix_latex+"\\Rightarrow</math-span></td><td><math-span>" + ithAnalysisStep.U_matrix_inverse_latex + "</math-span></td>" + flagHtml;
                ithHtml += "<td><math-span id=''>"+matrixTermsLatex+"\\Rightarrow<math-span></td><td><math-span>" + ithAnalysisStep.final+"\\Rightarrow</math-span></td><td><math-span>" + ithAnalysisStep.final_inverse + "</math-span></td>";

            // ithHtml += "</tr></table></td>";

            ithHtml += "</tr>";

            ithHtml += "<tr><td class='spacer_row_col_small'></td></tr>"

            innerHTML += ithHtml;
        } else if ( resultsElem == appState.flagList ) {
            innerHTML += "<table>";
            var rowDimensionFirst = true;
            for ( var i = 0; i < responseMessage.length; i++ ) {

                var readableOrdinal = i + 1;

                var ithFlag = responseMessage[i];
                
                var ithHtml = "";
                ithHtml += "<tr>";

                ithHtml += "<td style='font-size:150%; font-weight:800;'>"+readableOrdinal+":</td>";

                ithHtml += "<td><table><tr>";
                    ithHtml += "<td>";

                    // var matrixTermsLatex = "";
                    // for ( var k = 0; k < ithOptimizedStep.matrix_terms.length; k++ ) {
                    //     var kthMatrixTerm = ithOptimizedStep.matrix_terms[k];
                    //     var kthMatrixTermLatex = matrixToLatex(kthMatrixTerm, null, null, rowDimensionFirst);
                    //     matrixTermsLatex += kthMatrixTermLatex
                    // }
                    // var ordinalHtml = "<span>"+ithOptimizedStep.ordinal+"</span>";

                    // var flagLatex = combineMatricesWithVerticalLines(ithFlag.matrix_list);
                    var flagLatex = ithFlag.matrix_list;

                    ithHtml += "<math-span id=''>"+flagLatex+" \\hspace{10px} \\text{Rank} \\hspace{5px} = \\hspace{5px} " + ithFlag.rank + "</math-span>";
                    ithHtml += "</td>";
                    
                ithHtml += "</tr></table></td>";

                ithHtml += "</tr>";

                ithHtml += "<tr><td class='spacer_row_col_small'></td></tr>"

                innerHTML += ithHtml;

            }

            innerHTML += "</table>";
        } else if ( resultsElem == appState.sortedColumnLists ) {
            resultsElem.innerHtml = "<math-span>" + responseMessage + "</math-span>";
        } else if ( resultsElem == appState.columnLists ) {
            innerHTML += "<table>";
            var rowDimensionFirst = true;
            for ( var i = 0; i < responseMessage.length; i++ ) {

                var readableOrdinal = i + 1;

                var ithOptimizedStep = responseMessage[i];
                
                var ithHtml = "";
                ithHtml += "<tr>";

                ithHtml += "<td style='font-size:150%; font-weight:800;'>"+readableOrdinal+":</td>";

                ithHtml += "<td><table><tr>";
                    ithHtml += "<td>";

                    var matrixTermsLatex = "";
                    for ( var k = 0; k < ithOptimizedStep.matrix_terms.length; k++ ) {
                        var kthMatrixTerm = ithOptimizedStep.matrix_terms[k];
                        var kthMatrixTermLatex = matrixToLatex(kthMatrixTerm, null, null, rowDimensionFirst);
                        matrixTermsLatex += kthMatrixTermLatex
                    }
                    var ordinalHtml = "<span>"+ithOptimizedStep.ordinal+"</span>";

                    ithHtml += "<math-span id=''>"+matrixTermsLatex+"\\Rightarrow" + ithOptimizedStep.concatenated_matrix + "</math-span>";
                    ithHtml += "</td>";
                    
                ithHtml += "</tr></table></td>";

                ithHtml += "</tr>";

                ithHtml += "<tr><td class='spacer_row_col_small'></td></tr>"

                innerHTML += ithHtml;

            }

            innerHTML += "</table>";
        } else if ( resultsElem == appState.analysisSteps ) {
            innerHTML += "<table>";
            var rowDimensionFirst = true;
            for ( var i = 0; i < responseMessage.length; i++ ) {

                var readableOrdinal = i + 1;

                var ithAnalysisStep = responseMessage[i];
                
                var ithHtml = "";
                ithHtml += "<tr>";

                ithHtml += "<td style='font-size:150%; font-weight:800;'>"+readableOrdinal+":</td>";

                // ithHtml += "<td><table><tr>";
                    // ithHtml += "<td>";

                    var matrixTermsLatex = "";
                    for ( var k = 0; k < ithAnalysisStep.matrix_terms.length; k++ ) {
                        var kthMatrixTerm = ithAnalysisStep.matrix_terms[k];
                        var kthMatrixTermLatex = matrixToLatex(kthMatrixTerm, null, null, rowDimensionFirst);
                        matrixTermsLatex += kthMatrixTermLatex
                    }

                    var kthFlagMatricesU = ithAnalysisStep.flag_matrices_U;
                    var kthFlagMatricesPhi = ithAnalysisStep.flag_matrices_Phi;
                    var kthFlagMatricesSpanEqualities = ithAnalysisStep.flag_matrices_spans_are_equal;

                    var kthFlagLatex = "<table class=\"flag_breakdown_table\">";
                    kthFlagLatex += "<tr><td class=\"flag_breakdown_title_td flag_breakdown_td_right_border flag_breakdown_td_bottom_border\"><span class=\"flag_breakdown_title\">From U Matrix</span></td><td class=\"flag_breakdown_title_td flag_breakdown_td_right_border flag_breakdown_td_bottom_border\"><span class=\"flag_breakdown_title\">From Phi Module</span></td><td class=\"flag_breakdown_title_td flag_breakdown_td_bottom_border\"><span class=\"flag_breakdown_title\">Span Equality</span></td></tr>";
                    // var kthFlagLatex = "";
                    if ( kthFlagMatricesU ) {
                        for ( var m = 0; m < kthFlagMatricesU.length; m++ ) {
                            var areEqualSpans = kthFlagMatricesSpanEqualities[m];
                            
                            var bottomBorderCssClass = m < kthFlagMatricesU.length -1 ? "flag_breakdown_td_bottom_border" : "";

                            var equalityStatus = areEqualSpans ? "<span style=\"color:green\">EQUAL</span>" : "<span style=\"color:red\">NOT EQUAL</span>";
                            var spanEqualityHtml = "<td class=\""+bottomBorderCssClass+"\">"+equalityStatus+"</td>";

                            // var equalSign = areEqualSpans ? "=" : "\\neq";
                            // var equalSignTd = "<td><math-div>"+equalSign+"</math-div></td>";

                            // kthFlagLatex += "<math-div>" + kthFlagMatricesU[m] + "</math-div><math-div>"+equalSign+"</math-div><math-div>" +kthFlagMatricesPhi[m]+ "</math-div>";
                            kthFlagLatex += "<tr><td class=\"flag_breakdown_td_right_border "+bottomBorderCssClass+"\" style=\"text-align:right;\"><math-div style=\"display:inline-block;\">" + kthFlagMatricesU[m] + "</math-div></td><td class=\"flag_breakdown_td_right_border "+bottomBorderCssClass+"\"><math-div>" +kthFlagMatricesPhi[m]+ "</math-div></td>"+spanEqualityHtml+"</tr>";
                        }
                    }

                    kthFlagLatex += "</table>";

                    // var ordinalHtml = "<span>"+ithAnalysisStep.ordinal+"</span>";
                    var cssClassForU = ithAnalysisStep.is_invertible ? "clickable_math_elem" : "";
                    var clickJs = ithAnalysisStep.is_invertible ? "onclick='javascript:window.onMathElemClicked(this.dataset.flagLatex);'" : "";
                    var flagHtml = "";
                    
                    if ( ithAnalysisStep.is_invertible ) {
                        flagHtml += "<td data-flag-latex='"+kthFlagLatex+"' style='cursor:pointer;' "+clickJs+">";
                        flagHtml += "<span class='"+cssClassForU+"'>" + (ithAnalysisStep.passes_flag_condition ? "Satisfies Flag" : "Fails Flag") + "</span>";
                        flagHtml += "</td>";
                    } else {
                        flagHtml += "<td></td>";
                    }
                    
                    ithHtml += "<td><math-span id=''>"+matrixTermsLatex+"\\Rightarrow<math-span></td><td><math-span>" + ithAnalysisStep.U_matrix_latex+"\\Rightarrow</math-span></td><td><math-span>" + ithAnalysisStep.U_matrix_inverse_latex + "</math-span></td>" + flagHtml;

                // ithHtml += "</tr></table></td>";

                ithHtml += "</tr>";

                ithHtml += "<tr><td class='spacer_row_col_small'></td></tr>"

                innerHTML += ithHtml;

            }

            innerHTML += "</table>";

        } else if ( resultsElem == appState.reductionSteps ) {
            innerHTML += "<table>";
            var rowDimensionFirst = true;
            for ( var i = 0; i < responseMessage.length; i++ ) {
                var elemId_matrixBefore = "matrix-before-" + i;

                var readableOrdinal = i + 1;
                var reductionOrdinal = readableOrdinal; // may change in future.

                var ithReductionStep = responseMessage[i];
                // console.log(ithReductionStep);

                var matrixSize = ithReductionStep.matrix_before.length;
                var UpsilonOrdinalHigher = matrixSize - ithReductionStep.zero_based_col;
                var UpsilonOrdinalLower = matrixSize - ithReductionStep.zero_based_row;
                var kExponent = ithReductionStep.k_exponent_value;
                

                var ithHtml = "";
                ithHtml += "<tr>";

                ithHtml += "<td style='font-size:150%; font-weight:800;'>"+readableOrdinal+":</td>";

                ithHtml += "<td><table><tr>";
                    ithHtml += "<td>";
                    var ithMatrixBeforeLatex = matrixToLatex(ithReductionStep.matrix_before, null, null, rowDimensionFirst);
                    ithHtml += "<math-span id='"+elemId_matrixBefore+"'>"+ithMatrixBeforeLatex+"</math-span>";
                    ithHtml += "</td>";

                    ithHtml += "<td>";
                    ithHtml += "<table>";
                    // ithHtml += "<tr><td>Needs reduction = </td><td>" + ithReductionStep.needed_reduction + "</td></tr>";
                    ithHtml += "<tr><td style='text-align:right'><math-span>\\mu\\tau^k\\cdot\\Upsilon_"+UpsilonOrdinalHigher+" - \\Upsilon_"+UpsilonOrdinalLower+"\\cdot\\mu\\tau^k = </math-span></td><td><math-span>"+ithReductionStep.delta+"</math-span></td></tr>";
                    ithHtml += "<tr><td style='text-align:right'><math-span>\\mu_{"+reductionOrdinal+"}\\tau^"+kExponent+"\\cdot\\Upsilon_"+UpsilonOrdinalHigher+" - \\Upsilon_"+UpsilonOrdinalLower+"\\cdot\\mu_{"+reductionOrdinal+"}\\tau^"+kExponent+" = </math-span></td><td><math-span>"+ithReductionStep.delta_subbed+"</math-span></td></tr>";
                    ithHtml += "<tr><td style='text-align:right'><math-span>f_k(\\mu_{"+reductionOrdinal+"}) = </math-span></td><td><math-span>"+ithReductionStep.mu_function+" = 0</math-span></td></tr>";
                    ithHtml += "<tr><td style='text-align:right'><math-span>"+ithReductionStep.element_before+" \\rightarrow</math-span></td><td><math-span>"+ithReductionStep.element_after+"</math-span></td></tr>";
                    ithHtml += "</table>";
                    ithHtml += "</td>";
                ithHtml += "</tr></table></td>";

                ithHtml += "</tr>";

                ithHtml += "<tr><td class='spacer_row_col_small'></td></tr>"

                innerHTML += ithHtml;

            }

            innerHTML += "</table>";

        } else if ( resultsElem == appState.results2 ) {
            // appState.results2.innerHTML = "<math-span>"+responseMessage+"</math-span>";
            
            // setTimeout(function() {
            //     var mathSpan = appState.results2.getElementsByTagName("math-span")[0];
            //     var verticalSeperators = Array.from(mathSpan.shadowRoot.querySelectorAll(".ML__vertical-separator"));
            //     verticalSeperators.forEach(ith => {
            //         ith.style.borderRight = "2px solid currentColor";
            //     })
            // }, 1000);
        } else if ( resultsElem == appState.unitTestResults && Array.isArray(responseMessage) ) {
            responseMessage = responseMessage.join(",\\hspace{5px}");
        } else if (  rMatrixConfig ) {

            var totalReductionCount = responseMessage.total_reduction_count;
            var totalSecondsToGenerate = responseMessage.total_seconds_to_generate;
            totalSecondsToGenerate = shortenRawSeconds(totalSecondsToGenerate);
            responseMessage = responseMessage.r_matrix

            

            resultsElem.style.cssText = "";

            var resultsElemTitle = document.getElementById(resultsElem.id+"-title");

            var baseTitle = rMatrixConfig.title;
            var baseTooltipTitle = rMatrixConfig.title;
            var titleOrdinal = rMatrixConfigIndex + 1;
            
            var copyButton = '<button style="margin-left:10px;" id="'+(resultsElem.id+'-latex-textarea-copy')+'" type="button">Copy As LaTeX</button>';

            if ( resultsElemTitle ) {
                // resultsElemTitle.innerHTML = baseTitle + " (Reduction Count = "+ totalReductionCount + ", Execution Time = "+totalSecondsToGenerate+" seconds)";
                var matrixSingularOrPlural = rMatrixConfig == OUTPUT_MATRIX_CONFIGS[OUTPUT_MATRIX_CONFIGS.length-2] ? "Matrices" : "Matrix";
                resultsElemTitle.innerHTML = "Output " + titleOrdinal + ": <math-span>" + baseTitle + "</math-span> " + matrixSingularOrPlural;

                resultsElemTitle.parentElement.innerHTML = resultsElemTitle.outerHTML + copyButton;

                var latexWrapper = document.getElementById(resultsElem.id + "-latex");
                if ( latexWrapper ) {
                    addRMatrixLatexToDom(latexWrapper, responseMessage);
                }
            }
            
            innerHTML += "<table class='auto-size-table'>";
            for ( var row = -1; row < responseMessage.length; row++ ) {

                innerHTML += "<tr>";
                
                if ( row === -1 ) {
                    var colCount = responseMessage[0].length;
                    for ( var col = -1; col < colCount; col++ ) {
                        if ( col == -1 ) {
                            innerHTML += "<td style='text-align:center'>";
                            innerHTML += "";
                        } else {
                            innerHTML += "<td style='font-weight:800; text-align:center'>";
                            innerHTML += "<span class='r_by_r_label'><span>"+(col+1) + "</span></span>";
                        }

                        innerHTML += "</td>";
                    }
                } else {
                    var ithRow = responseMessage[row];
                
                    for ( var col = -1; col < ithRow.length; col++ ) {
                        if ( col == -1 ) {
                            innerHTML += "<td style='font-weight:800; text-align:center;'>";
                            innerHTML += "<span class='r_by_r_label'><span>"+(row+1)+"</span></span>";
                        } else {
                            var ithRowColValue = ithRow[col];
                            var oneBasedRow = row + 1;
                            var oneBasedCol = col + 1;
                            var cellTitle = "<math-div>$"+baseTooltipTitle + "_{"+oneBasedRow+"\\text{x}"+oneBasedCol+"}$ = " + ithRowColValue + "</math-div>";

                            
                        
                            if ( ithRowColValue == MATRIX_ELEM__ZERO ) {
                                innerHTML += "<td class='r_by_r_zero_cell' style='text-align:center'>";
                                innerHTML += "<span tooltip-title='"+cellTitle+"' class='r_by_r_entry r_by_r_zero'><span>0</span></span>";
                            } else {
                                innerHTML += "<td class='r_by_r_equation_cell' style='text-align:center'>";
                                innerHTML += '<div tooltip-title="'+cellTitle+'" class="r_by_r_entry_wrapper" >';
                                innerHTML += "<math-span class='r_by_r_entry'>" + ithRowColValue + "</math-span>";
                                innerHTML += "</div>";
                            }
                        }
                        
                        innerHTML += "</td>";
                    }
                }
                
                innerHTML += "</tr>";
            }
            innerHTML += "</table>";
        }

        if ( isMathDivOrSpan(resultsElem) ) {
            resultsElem.textContent = responseMessage;
            resultsElem.render();
        } else {
            if ( isDiv === true ) {
                resultsElem.innerHTML = innerHTML ? innerHTML : "<math-field style='opacity:1; border:none; background:none;' max-matrix-cols=100 disabled>"+responseMessage+"</math-field>"
            } else {
                resultsElem.innerHTML = innerHTML ? innerHTML : "<math-span max-matrix-cols=100>"+responseMessage+"</math-span>"
            }
        }
        
        // responseMessage = responseMessage.replaceAll("begin{matrix}", "begin{pmatrix}");
        // responseMessage = responseMessage.replaceAll("end{matrix}", "end{pmatrix}");

    }
}

function enforceMaxMatrixDimension(equationOrdinal, dimensionInputElem) {
    var equationElem = equationOrdinal == 0 ? appState.equation0 : appState.equation1;
    var latexInputElem = equationOrdinal == 0 ? appState.latex0 : appState.latex1;

    var d = getLargestDimensionOfLatexMatrix(equationElem.value);

    var hadToEnforce = false;

    if ( d > MAX_INPUT_MATRIX_DIMENSION ) {
        resizeMatricesInInputElem(equationOrdinal, MAX_INPUT_MATRIX_DIMENSION);

        hadToEnforce = true;

        showMatrixExceededMaxDimensionsError(equationElem);

        d = MAX_INPUT_MATRIX_DIMENSION;
    }

    setInputValueProgrammaticallyForIntInput(dimensionInputElem, d >= 1 ? d : 1);

    return hadToEnforce;
}

function onInputChanged(forceServerCalculation = false, firstStartUp = false) {

    var variant = getVariant();

    var actuallyChanged = true;

    var staleMessage = forceServerCalculation ? noteHtml("Calculating...") : noteHtml("Stale, click Run in upper right.");
    var inputErrorMessage = errorHtml("Fix input errors then click Run in upper right.");

    if ( variant == VARIANT__PHI_PSI_EMBEDDING ) {
        var currentEquation0 = getMathFieldValueFromLocalStorage(0);
        var currentEquation1 = getMathFieldValueFromLocalStorage(1);
        var firstChanged = true;
        var secondChanged = true;

        enforceMaxMatrixDimension(0, appState.dForNInput);
        enforceMaxMatrixDimension(1, appState.dForMInput);


        if ( appState.equation0.value == currentEquation0 ) {
            firstChanged = false;
        } else {
            // appState.resultsPhiModuleCondensed.innerHTML = staleMessage;
        }

        if ( appState.equation1.value == currentEquation1 ) {
            secondChanged = false;
        } else {
            // appState.resultsPsiModuleCondensed.innerHTML = staleMessage;
        }

        // UPDATE: Now always show stale message, to refresh UI properly from a server error.
        // appState.resultsPhiModuleCondensed.innerHTML = staleMessage;
        // appState.resultsPsiModuleCondensed.innerHTML = staleMessage;

        actuallyChanged = firstChanged || secondChanged;
    } else {
        var currentEquation0 = getMathFieldValueFromLocalStorage(0);
        
        if ( appState.equation0.value == currentEquation0 ) {
            actuallyChanged = false; 
        } else {
            // appState.resultsPhiModuleCondensed.innerHTML = staleMessage;
        }
    }

    // if ( actuallyChanged === false && forceServerCalculation == false && firstStartUp === false ) {
    //     console.log("Runulation not needed, values are the same.");
    //     return;
    // }

    clearInputErrorContainers();
    saveMathFieldValuesToLocalStorage();

    var inputErrorElem = appState.results0;
    var message = "Checking...";
    var checkingHtml = "<div class='output_note'>"+message+"</div>";
    addToOutput(inputErrorElem, checkingHtml);

    if ( forceServerCalculation === true ) {
        recalculateOnBackEnd();
    } else {
        validatePhiPsiInput(appState.equation0.value, appState.equation1.value);
    }

    doResizeRoutine();
}

function setInputValueProgrammaticallyForIntInput(inputElem, value) {
    inputElem.value = value;
    inputElem.setAttribute("old_value", value);
}

function resetEquations() {
    clearMathFieldValuesFromLocalStorage();

    appState.tmodule0 = [];
    appState.tmodule1 = [];

    // if ( appState.nInput ) {
    //     appState.nInput.value = getDefaultForN(0);
    // }

    // if ( appState.mInput ) {
    //     appState.mInput.value = getDefaultForN(1);
    // }

    if ( appState.dForNInput ) {
        setInputValueProgrammaticallyForIntInput(appState.dForNInput, DEFAULT_VALUE__D);
    }
    
    if ( appState.dForMInput ) {
        setInputValueProgrammaticallyForIntInput(appState.dForMInput, DEFAULT_VALUE__D);
    }
    
    loadStartingExample();
}

function loadStartingExample() {
    var variant = getVariant();

    if ( variant == VARIANT__PHI_PSI_EMBEDDING ) {
        loadPhiPsiExample(EXAMPLES__PHI_PSI_EMBEDDDING[0][0], EXAMPLES__PHI_PSI_EMBEDDDING[0][1]);
    }
}

function stripWhitespaceFromAst(astNodeArray) {
    var toReturn = [];

    astNodeArray.forEach(function(ith) {
        if ( ith.type != "whitespace" ) {
            toReturn.push(ith);
        }
    });

    return toReturn;
}

function splitAstArrayByContent(astNodeArray, splitter) {
    var toReturn = [];
    var subArray = [];

    for ( var i = 0; i < astNodeArray.length; i++ ) {
        var ithNode = astNodeArray[i];

        if ( ithNode.content && ithNode.content == splitter ) {
            toReturn.push(subArray);
            subArray = [];
        } else {
            subArray.push(ithNode);
        }
    }

    if ( subArray.length > 0 ) {
        toReturn.push(subArray);
    }

    return toReturn;
}

function astMatrixToPddMatrix(astMatrixNode) {
    var matrixRowNodes = splitAstArrayByContent(astMatrixNode.content, "\\");

    var rowCount = matrixRowNodes.length;

    var pddMatrix = nilPddMatrix();

    for ( var i = 0; i < rowCount; i++ ) {
        var ithRow = matrixRowNodes[i];
        var matrixColNodes = splitAstArrayByContent(ithRow, "&");
        var colCount = matrixColNodes.length;

        pddMatrix = pddMatrix.length > 0 ? pddMatrix : newPddMatrix(colCount, rowCount);

        for ( var k = 0; k < colCount; k++ ) {
            var kthCol = matrixColNodes[k];
            kthCol = stripWhitespaceFromAst(kthCol);
            var colLatex = "";

            if ( kthCol.length > 0 ) {
                colLatex = pdd_node_bundle.toString(kthCol);
            } else {
                if ( FEATURE_FLAG__USE_LATEX_PLACEHOLDERS === true ) {
                    colLatex = "\\placeholder{}";
                } else {
                    colLatex = "0";
                }
            }

            pddMatrix[k][i] = colLatex;
        }
    }

    return pddMatrix;
}

function isAstMatrix(astNode) {
    if ( astNode.env && astNode.env.endsWith("matrix") ) {
        return true;
    } else {
        return false;
    }
}

function isContentPowerOfInt(contentArray, exponent, exponentNotationIndex) {
    
    if ( contentArray[exponentNotationIndex].content.toLowerCase() == "^"+exponent ) {
        return true;
    } else {
        if ( contentArray[exponentNotationIndex].content.toLowerCase() == "^" ) {
            if ( contentArray[exponentNotationIndex+1] && contentArray[exponentNotationIndex+1].type == "group" ) {
                if ( contentArray[exponentNotationIndex+1].content[0] && contentArray[exponentNotationIndex+1].content[0].type == "string" && contentArray[exponentNotationIndex+1].content[0].content == exponent ) {
                    return true;
                } else {
                    return false;
                }
            }
        } else {
            return false;
        }
    }
}


function getPddMatrixColCount(pddMatrix) {
    if ( Array.isArray(pddMatrix) && pddMatrix.length > 0 ) {
        return pddMatrix.length;
    } else {
        return 0;
    }
}

function getPddMatrixRowCount(pddMatrix) {
    if ( Array.isArray(pddMatrix) && pddMatrix.length > 0 ) {
        return pddMatrix[0].length;
    } else {
        return 0;
    }
}

function forEachElemInPddMatrix(pddMatrix, callback) {
    var colCount = getPddMatrixColCount(pddMatrix);
    var rowCount = getPddMatrixRowCount(pddMatrix);

    for( var i = 0; i < rowCount; i++ ) {
        for( var k = 0; k < colCount; k++ ) {
            var returnValue = callback(pddMatrix[k][i], k, i);

            if ( returnValue === false ) {
                return;
            } else {
                if ( returnValue && returnValue !== true ) {
                    pddMatrix[k][i] = returnValue;
                }
            }
        }
    }
}

function newPddMatrix(colCount, rowCount) {
    var toReturn = nilPddMatrix();

    for( var i = 0; i < colCount; i++ ) {
        var col = [];

        for( var k = 0; k < rowCount; k++ ) {

            // if ( ordinal === 0 && i === k ) {
            //     col.push("\\theta");
            // } else {
            //     col.push("");
            // }

            col.push("");
        }

        toReturn.push(col);
    }

    return toReturn;
}

function isElemInUpperRightDiagonal(d, rowIndex, colIndex) {
    if (colIndex > 0 && colIndex < d) {
        var rowLimitExclusive = colIndex;
        if ( rowIndex >= rowLimitExclusive ) {
            return false;
        } else {
            return true;
        }
    } else {
        return false;
    }
}

function isEmptyMatrixElem(latexForMatrixElem) {
    if ( latexForMatrixElem ) {
        if ( latexForMatrixElem.startsWith("\\placeholder{") ) {
            return true;
        } else {
            return false;
        }
    } else {
        return true;
    }
}

function matrixToLatex(matrix, ordinal = null, splitString = null, rowDimensionFirst = false) {
    var latexString = "\\begin{"+DEFAULT_LATEX_MATRIX_ENVIRONMENT+"}";

    if ( splitString ) {
        latexString += splitString;
    }

    var colCount = matrix.length;
    var rowCount = matrix[0].length;

    var defaultEmptyValue = FEATURE_FLAG__USE_LATEX_PLACEHOLDERS === true ? "\\placeholder{}" : "0";

    for (var i = 0; i < rowCount; i++) {
        for (var k = 0; k < colCount; k++) {
            var col = k;
            var row = i;

            var contentsAtPosition = rowDimensionFirst ? matrix[row][col] : matrix[col][row];
            contentsAtPosition = contentsAtPosition ? contentsAtPosition : defaultEmptyValue;
            latexString += contentsAtPosition;

            if ( k < colCount - 1) {
                latexString += " & ";
            } else {
                latexString += " \\\\";
            }
        }
    }

    latexString += "\\end{pmatrix}";

    var variant = getVariant();

    if ( variant == VARIANT__PHI_PSI_EMBEDDING || variant == VARIANT__RIGHTMOST_BASIS ) {
        if ( ordinal != null ) {
            latexString += "\\cdot \\tau^{" + ordinal + "}";
        }
    }
    
    return latexString;
}