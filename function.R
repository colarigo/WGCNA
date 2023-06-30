multiData.subset = function(multiData, colIndex = NULL, rowIndex = NULL)
{
  size = checkSets(multiData);
  if (is.null(colIndex)) colIndex = c(1:size$nGenes);
  if (is.null(rowIndex)) rowIndex = lapply(size$nSamples, function(n) {c(1:n)})
  if (length(rowIndex)!=size$nSets) 
    stop("If given, 'rowIndex' must be a list of the same length as 'multiData'.");
  out = list();
  for (set in 1:size$nSets)
    out[[set]] = list(data = multiData[[set]]$data[rowIndex[[set]], colIndex, drop = FALSE]);
  names(out) = names(multiData);
  out;
}

multiData2list = function(multiData)
{
  lapply(multiData, `[[`, 'data');
}

list2multiData = function(data)
{
  out = list();
  for (set in 1:length(data))
    out[[set]] = list(data = data[[set]]);
  names(out) = names(data);
  out;
}

multiData.mapply = function(FUN, ..., MoreArgs = NULL, mdmaSimplify = FALSE, mdma.doCollectGarbage = FALSE,
                            mdma.argIsMultiset = NULL)
{
  dots = list(...);
  if (length(dots)==0) 
    stop("No arguments were specified. Please type ?multiData.mapply to see the help page.");
  dotLengths = sapply(dots, length);
  if (any(dotLengths!=dotLengths[1]))
    stop(spaste("All arguments to vectorize over must have the same length.\n", 
                "Scalar arguments should be put into the 'MoreArgs' argument.\n",
                "Note: lengths of '...' arguments are: ", paste(dotLengths, collapse = ", ")));
  nArgs = length(dots);
  res = list();
  if (is.null(mdma.argIsMultiset)) mdma.argIsMultiset = sapply(dots, isMultiData);

  FUN = match.fun(FUN);
  nSets = dotLengths[1];
  for (set in 1:nSets)
  {
    localArgs = list();
    for (arg in 1:nArgs)
      localArgs[[arg]] = if (mdma.argIsMultiset[arg]) dots[[arg]] [[set]] $ data else dots[[arg]] [[set]];
    names(localArgs) = names(dots);
    res[[set]] = list(data = do.call(FUN, c(localArgs, MoreArgs)));
    if (mdma.doCollectGarbage) collectGarbage();
  }

  names(res) = names(dots[[1]]);

  if (mdmaSimplify)
    return(multiData.simplify(res));

  return(res);
}

multiData.simplify = function(multiData)
{
  len = length(multiData[[1]]$data);
  dim = dim(multiData[[1]]$data);
  simplifiable = TRUE;
  nSets = length(multiData);
  for (set in 1:nSets)
  {
    if (len!=length(multiData[[set]]$data)) simplifiable = FALSE;
    if (!isTRUE(all.equal( dim, dim(multiData[[set]]$data)))) simplifiable = FALSE;
  }
  if (simplifiable)
  {
    if (is.null(dim)) {
       innerDim = len;
       innerNames = names(multiData[[1]]$data);
       if (is.null(innerNames)) innerNames = spaste("X", c(1:len));
    } else {
       innerDim = dim;
       innerNames = dimnames(multiData[[1]]$data);
       if (is.null(innerNames)) 
         innerNames = lapply(innerDim, function(x) {spaste("X", 1:x)})
       nullIN = sapply(innerNames, is.null);
       if (any(nullIN))
         innerNames[nullIN] = lapply(innerDim[nullIN], function(x) {spaste("X", 1:x)})
    }
    setNames = names(multiData);
    if (is.null(setNames)) setNames = spaste("Set_", 1:nSets);
    multiData.s = matrix(NA, prod(innerDim), nSets);
    for (set in 1:nSets)
      multiData.s[, set] = as.vector(multiData[[set]]$data);

    dim(multiData.s) = c(innerDim, nSets);
    if (!is.null(innerNames))
      dimnames(multiData.s) = c (if (is.list(innerNames)) innerNames else list(innerNames), list(setNames));
    return(multiData.s);
  }
  return(multiData);
}

labeledBarplot2 = function(x, names, g = NULL, horiz = FALSE, 
                           srt = if (horiz) 0 else 45, 
                           adj = if (horiz) c(1, 0.5) else c(1, 0.5), 
                           namesColor = 1, cex.names = 1, 
                           addGuide = TRUE, guideColor = "grey30", guideLty = 2, ...)
{
  
  if (!is.null(g)) {
    mp = verboseBarplot(x, g, names.arg = NA, horiz = horiz, ...);
    heights = attr(mp, "height");
  } else {
    mp = barplot(x, names.arg = NA, horiz= horiz, ...);
    heights = x;
    attr(mp, "height") = x;
    attr(mp, "stdErr") = rep(0, length(x));
  }
  
  box = par("usr");
  
  if (horiz)
  {
    yText = mp;
    xText = rep(box[1] - 0.01 * (box[2] - box[1]), length(mp));
  } else {
    xText = mp;
    yText = rep(box[3] - 0.02 * (box[4] - box[3]), length(mp));
  }
  
  text(xText, yText, names, srt = srt, adj = adj, col = namesColor, cex = cex.names, xpd = TRUE);
  if (addGuide) for (i in 1:length(mp))
    if (horiz)
    {
      lines(c(min(0, heights[i]), box[1] - 0.02 * (box[2] - box[1])),  rep(mp[i], 2), 
            col = guideColor, lty = guideLty);
    } else {
      lines(rep(mp[i], 2), c(min(0, heights[i]), box[3] - 0.02 * (box[4] - box[3])), 
            col = guideColor, lty = guideLty);
    }
  
  invisible(mp);
  
}