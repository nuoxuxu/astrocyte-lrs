switchPlotFromTables <- function(
    ### Required tables
    isoformFeatures,
    exons,
    conditions,

    ### Gene/isoform selection (exactly one must be provided)
    gene = NULL,
    isoform_id = NULL,

    ### Condition selection
    condition1,
    condition2,

    ### Optional annotation tables
    orfAnalysis = NULL,
    domainAnalysis = NULL,
    idrAnalysis = NULL,
    signalPeptideAnalysis = NULL,
    topologyAnalysis = NULL,

    ### Plot arguments (passed through to switchPlot)
    IFcutoff = 0.05,
    dIFcutoff = 0.1,
    alphas = c(0.05, 0.001),
    rescaleTranscripts = TRUE,
    plotTopology = TRUE,
    reverseMinus = TRUE,
    addErrorbars = TRUE,
    logYaxis = FALSE,
    localTheme = theme_bw(base_size = 14),
    additionalArguments = list()
) {
    if (!is.data.frame(isoformFeatures)) {
        stop("'isoformFeatures' must be a data.frame")
    }
    if (class(exons) != 'GRanges') {
        stop("'exons' must be a GRanges object")
    }
    if (!is.data.frame(conditions) ||
        !all(c('condition', 'nrReplicates') %in% colnames(conditions))) {
        stop("'conditions' must be a data.frame with columns 'condition' and 'nrReplicates'")
    }

    switchObj <- new(
        "switchAnalyzeRlist",
        list(
            isoformFeatures = isoformFeatures,
            exons           = exons,
            conditions      = conditions,
            sourceId        = 'fromTables'
        )
    )

    if (!is.null(orfAnalysis))           switchObj$orfAnalysis           <- orfAnalysis
    if (!is.null(domainAnalysis))        switchObj$domainAnalysis        <- domainAnalysis
    if (!is.null(idrAnalysis))           switchObj$idrAnalysis           <- idrAnalysis
    if (!is.null(signalPeptideAnalysis)) switchObj$signalPeptideAnalysis <- signalPeptideAnalysis
    if (!is.null(topologyAnalysis))      switchObj$topologyAnalysis      <- topologyAnalysis

    switchPlot(
        switchAnalyzeRlist  = switchObj,
        gene                = gene,
        isoform_id          = isoform_id,
        condition1          = condition1,
        condition2          = condition2,
        IFcutoff            = IFcutoff,
        dIFcutoff           = dIFcutoff,
        alphas              = alphas,
        rescaleTranscripts  = rescaleTranscripts,
        plotTopology        = plotTopology,
        reverseMinus        = reverseMinus,
        addErrorbars        = addErrorbars,
        logYaxis            = logYaxis,
        localTheme          = localTheme,
        additionalArguments = additionalArguments
    )
}