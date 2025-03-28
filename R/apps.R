
## ================== Custom instantiators and converters for gMatrix ================== ##

#' @name homeology
#' Given a (named) character vector or XStringSet computes a matrix of string distances between
#' all tiles of width + pad.  Returns a gMatrix whose seqnames are the names(seq) and coordinates
#' refer to the provided sequences.  Every GRanges pair ij of the gMatrix represents the edit distance
#' (default levenshtein) between the sequences in the provided sequence tiles (+/- some padding).
#'
#' If flag rc=TRUE, then the distance ij will be the edit distance between the sequence in GRanges i and the reverse
#' complement of sequence in GRanges j. 
#'
#' @param seqs a character vector or XStringSet (e.g. DNAStringSet)
#' @param gr GRanges
#' @param stride specifies both the width and the stride of the sequence bins / tiles around which to measure homeology and which are returned in the output gMatrix
#' @param pad the padding around each tile with which to measure sequence homeology
#' @param rc logical flag specifying whether to provide reverse complement distance
#' @author Marcin Imielinski
#' @export
#' @return new gPair containing the difference between x and y
homeology = function(seqs, gr = NULL, stride = 1, pad = 10, rc = FALSE, verbose = FALSE, method = "levenshtein", ignoreCase = TRUE, ...)
{
  ## short hand for supplying GRanges
  if (is(seqs, 'GRanges'))
  {
    if (is.null(seqs$seq))
      stop('If GRanges provided to homeology, they must have character or DNAStringSet field $seq that matches the width of the corresponding GRanges')

    gr = seqs
    seqs = seqs$seq
    names(seqs) = seqnames(gr)
  }
  
  if (is.null(names(seqs)))
    names(seqs) = 1:length(seqs)

  if (is.character(seqs))
  {
    seqs = DNAStringSet(seqs)
  }

  if (is.null(gr))
    gr = GRanges(names(seqs), IRanges(1, width(seqs)))
  else if (!identical(width(gr), width(seqs)) || !identical(as.character(seqnames(gr)), names(seqs)))
    stop('If gr is provided, gr must be same length as seqs and each gr[i] must have the same width and seqnames as the nchar and name of the corresponding seqs[i]')

  gr = gr.fix(gr)
  if (any(ix <- strand(gr)=='*'))
    strand(gr)[ix] = '+'

  width = stride
  starts = start(gr)
  names(starts) = seqnames(gr)

  grs = gr.fix(GRanges(seqnames(gr), IRanges(1, width(gr))))
  grs$sn.new = dedup(as.character(seqnames(grs)))
  tiles.og  = suppressWarnings(gr.tile(grs, width))
  strand(tiles.og) = strand(gr)[tiles.og$query.id]
  tiles = suppressWarnings(tiles.og + pad)
  end(tiles) =  suppressWarnings(pmin(seqlengths(tiles)[as.character(seqnames(tiles))], end(tiles)))
  start(tiles) = suppressWarnings(pmax(1, start(tiles)))

  ## need to "dedup" seqnames in case the input granges are on the same chromosome
  tiles = GRanges(grs$sn.new[tiles$query.id], ranges(tiles), strand = strand(tiles))
  names(seqs) = dedup(names(seqs))

  if (rc)
  {
    tiles = c(tiles, gr.flipstrand(tiles))
  }

  if (verbose)
    {
      message('Populating ', length(tiles), ' bins with width ', width(tiles)[1], ' and stride ', width)
    }

  tiles$seq = BSgenome::getSeq(seqs, tiles)

  if (verbose)
  {
    message('Computing ', length(tiles), ' by ', length(tiles), ' distance matrix')
  }

  D = as.matrix(Biostrings::stringDist(tiles$seq, method = method, ignoreCase = ignoreCase, upper = TRUE, ...))

  if (rc) ## only keep +- ij pairs
  {
    D = D[which(as.logical(strand(tiles)=='+')), which(as.logical(strand(tiles)=='-'))]
  }

  ## case gMatrix to original (nonoverlapping) tiles (ie without the padding)
  out = GRanges(as.character(seqnames(tiles.og)),
                IRanges(start(tiles.og) + starts[tiles.og$query.id]-1,
                        end(tiles.og) + starts[tiles.og$query.id]-1),
                strand = strand(tiles.og))
  out$seq = tiles$seq
  return(gM(out, D))
}


#' @name straw
#' @title straw
#' @description
#' Instantiates gMatrix from .hic object using straw API https://github.com/theaidenlab/straw/tree/master/R
#' used to extract all of the data in length n query gr (note will do n choose 2 queries since
#' straw only supports pairwise queries)
#'
#' @keywords straw
#' @param hic path to .hic file
#' @param gr granges to query
#' @param norm string specifying normalization to apply ("KR" (default), "NONE", VC")
#' @param res resolution to query (default 10kb)
#' @param type "BP" vs "FRAG"
#' @param mc.cores parallelization to apply for query (useful when length(gr)>2)
#' @return gMatrix 
#' @export
#' @author Marcin Imielinski
straw = function(hic, gr = NULL, norm = "NONE", type = 'BP', res = 1e4, mc.cores = 1, colormap = c('white', 'red', 'black'), ...)
{
  hic = normalizePath(hic)
  if (is.null(gr))
    gr = hg_seqlengths()

  if (is.integer(gr) | is.numeric(gr))
    {
      if (is.null(names(gr))) ## assume this is numeric form of chromosomes
        {
          gr = as.character(gr)
        }
      else ## assume it's a seqlengths
        {
          gr = si2gr(gr)
        }
    }
  
  if (is.character(gr))
    gr = parse.gr(gr, seqlengths = hg_seqlengths())

  if (!inherits(gr, 'GRanges'))
    gr = si2gr(gr)  

  gr = reduce(gr.stripstrand(gr[, c()]))
  grs = as.data.table(gr)[, paste(seqnames, start, end, sep = ":")]
  n = length(gr)
  grs.pairs = chr2.pairs = chr1.pairs = NULL
  grs.singles = paste(grs, grs)
  ## new new new
  r1 = grs
  r2 = grs
  chr1.singles = chr2.singles = as.character(seqnames(gr))
  if (n>1)
    {
      pairs = t(combn(n, 2)) ##pairs
      grs.pairs = paste(grs[pairs[,1]], grs[pairs[,2]])   ## combinations
      #' Tuesday, Nov 07, 2017 04:43:52 PM - Julie: adding as.character()
      r1 = c(r1, grs[pairs[,1]])
      r2 = c(r2, grs[pairs[,2]])
      chr1.pairs = as.character(seqnames(gr)[pairs[,1]])
      chr2.pairs = as.character(seqnames(gr)[pairs[,2]])
    }
  str = paste(norm, hic, c(grs.singles, grs.pairs), type, as.integer(res))
  chr1 = as.character(c(chr1.singles, chr1.pairs))
  chr2 = as.character(c(chr2.singles, chr2.pairs))
  ## use strawR
  out = rbindlist(mcmapply(r1 = r1,
                           r2 = r2,
                           FUN = function(r1, r2)
                           {
                             dt = as.data.table(
                                       strawr::straw(norm = norm,
                                                     unit = type,
                                                     binsize = res,
                                                     fname = hic,
                                                     chr1loc = r1,
                                                     chr2loc = r2))
                             setnames(dt, c("start1", "start2", "counts"))
                             dt[, chr1 := gsub("^([0-9XY]+):.*", "\\1", r1)]
                             dt[, chr2 := gsub("^([0-9XY]+):.*", "\\1", r2)]
                             dt[, end1 := start1 + res - 1]
                             dt[, end2 := start2 + res - 1]
                             return(dt)
                           }, SIMPLIFY = FALSE,  mc.cores = mc.cores), fill = TRUE)  
  ## out = rbindlist(mcmapply(str = str, chr1 = chr1, chr2 = chr2,
  ##                          FUN = function(str, chr1, chr2)
  ##                          {
  ##                            dt = as.data.table(
  ##                              tryCatch(GxG:::straw_R(str),
  ##                                       error = function (e)
  ##                                         data.table()))
  ##                            return(dt)
  ##                          }, SIMPLIFY = FALSE,  mc.cores = mc.cores), fill = TRUE)
  out = out[!is.na(counts), ]
  if (!nrow(out))
    stop('Query resulted in no output, please check .hic file or input coordinates')

  out[, end1 := start1+res-1]
  out[, end2 := start2+res-1]
  out[, str1 := paste0(chr1, ':', start1, '-', end1)]
  out[, str2 := paste0(chr2, ':', start2, '-', end2)]
  gr.out = unique(dt2gr(rbind(out[, .(seqnames = chr1, start = start1, end = end1)],
                              out[, .(seqnames = chr2, start = start2, end = end2)])))
  gr.map = data.table(grs = gr.string(gr.out), ix = 1:length(gr.out), key = 'grs')
  out$i1 = gr.map[.(out$str1), ix]
  out$j1 = gr.map[.(out$str2), ix]  
  out[, i := pmin(i1, j1)]
  out[, j := pmax(i1, j1)]
  out=out[i > 0 & i <= length(gr.out) & j > 0 & j <= length(gr.out)]

  gm = gM(gr.out, out[, .(i, j, value = counts)])
  return(gm) 
}


#' @name hicpro
#' @title hicpro
#' @description
#' Instantiates gMatrix from .matrix and .bed Hi-C pro file 
#'
#' @keywords hicpro
#' @param mat path to .matrix file
#' @param bed path to .bed file
#' @return gMatrix
#' @export
#' @author Marcin Imielinski
hicpro = function(mat, bed = NULL)
{
  if (is.null(bed))
    bed = paste0(gsub('.matrix$', '', mat), '_ord.bed')

  dat = fread(mat)
  setnames(dat, c('i', 'j', 'value'))
  
  if (!file.exists(bed))
    stop('bed file must be provided or one must exist with the same prefix as the .matrix file suffixed with _ord.bed.')
 
  gr = gr.fix(rtracklayer::import(bed))

  if (max(dat$i, dat$j)>length(gr))
    stop('provided .matrix has out of bounds indices relative to provided .bed file, please check to make sure these are compatible')
  return(gM(gr, dat))
}


#' @name cooler
#' @title cooler
#' @description
#' Instantiates gMatrix from .cool and/or .mcool file, +/- at specific locations. 
#'
#' To find viable resolutions (e.g. for mcool files) use info = TRUE and you will get a list
#' of info, including seqlengths. 
#'
#' @keywords cooler
#' @param file path to .cool or m.cool file
#' @param gr optional GRanges of intervals to query
#' @param res optional resolution to pull
#' @param info logical flag (FALSE) specifying whether to return a vector of file info instead of pulling down data, this will include available resolutions (for .mcool format) and fields
#' @return gMatrix or list
#' @export
#' @author Marcin Imielinski
cooler = function(file, gr = NULL, res = NULL, info = FALSE, field = NULL, nochr = FALSE, mc.cores = 1)
{
  contents = subfile = NULL

  ## if mcool need to pick a resolution
  if (grepl('mcool$', file))

  {
    contents = as.data.table(rhdf5::h5ls(file))[grepl('\\/resolutions\\/\\d+/bins$', group), ][!(name %in% c('chrom', 'end', 'start')), ]
    contents[, subfile := gsub('\\/bins$', '', group)]
    contents[, resolution := gsub('(\\/resolutions\\/)|(\\/bins)', '', group) %>% as.integer]
    contents = contents[order(resolution), .(field = name, dim = dim, resolution)]

    if (is.null(res))
      res = contents$resolution[1]

    if (!(res %in% contents$resolution))
      stop(sprintf('resolution not present in file %s: choose one of these available resolutions: %s or use default (%s)', file, paste(contents$resolution %>% unique, collapse = ', '), contents$resolution[1]))

    contents = contents[resolution %in% res, ]

    subfile = paste('resolutions', contents[, resolution[1]], sep = '/')

    ## not sure how to get any field other than count .. eg VC, VC_SQRT
    ## (my guess is these are supposed to get computed on the fly from the marginal
    ## correction factors rather than stored in the file)
    ##
    ## if (!(field %in% contents$field))
    ##   stop(sprintf('resolution not present in file %s: choose one of these available resolutions: %s or use default (%s)', file, paste(contents$resolution %>% unique, collapse = ', '), contents$resolution[1]))

    ## subfile = paste('resolutions', contents[field == field, resolution[1]], sep = '/')
  }

  ## install cooler if not installed on reticulate python virtualenv
  reticulate::py_run_string("try:
    import cooler
    cooler_installed = True
except ImportError:
    cooler_installed = False
")
  if (!reticulate::py$cooler_installed)
    reticulate::py_install("cooler")

  ## construct subfile (ie if we have mcool)
  if (!is.null(subfile))
  {
    uri = paste(file, subfile, sep = '::')
  } else
  {
    uri = file
  }

  ## load data via cooler
  reticulate::py_run_string(sprintf('data = cooler.Cooler("%s")', uri))

  ## get seqlengths from object
  tmp = reticulate::py_run_string('chromsizes = data.chromsizes')$chromsizes
  sl = structure(as.vector(tmp), names = names(tmp))

  if (info)
  {
    info = reticulate::py_run_string(sprintf('info = data.info'))$info
    if (!is.null(contents))
      {
        info$resolutions = contents$resolution %>% unique
        info$default.resolution = contents$resolution[1]
      }
    info$fields = colnames(reticulate::py_run_string(sprintf('pix = data.pixels()[:1]'))$pix)[-c(1:2)]
    info$seqlengths = sl
    return(info)
  }

  ## create our bin GRanges
  ##  bins = reticulate::py_run_string('bins = data.bins()[:]; print(bins)')$bins

  ## takes care of weird reticulate conversion error
  bins = reticulate::py_run_string('bins = data.bins()[:]; bins.chrom= bins.chrom.astype("str")')$bins %>% as.data.table
  bins[, start := start+1]
  setnames(bins, 'chrom', 'seqnames')
  bins = bins %>% dt2gr(seqlengths = sl)

  if (nochr)
    bins = gr.sub(bins, 'chr', '')

  if (!is.null(gr))
  {
    gr.old = gr
    gr = gr.fix(gr, sl, drop = TRUE) %>% gr.stripstrand
    ## more fixing
    end(gr) = pmax(end(gr), 1)
    start(gr) = pmax(start(gr), 1)
    if (length(gr)<length(gr.old))
      warning("Some of the provided ranges had to be dropped / clipped to make compatible with the seqlengths of the .cool / .mcool file")

    if (!length(gr))
      return(gM(GRanges(seqlengths = sl)))

    pix = expand.grid(i = 1:length(gr), j = 1:length(gr))

    res = mclapply(1:nrow(pix), function(k)
      reticulate::py_run_string(sprintf('res = data.matrix(balance = False, join = False, as_pixels = True).fetch("%s", "%s")',
                                        gr[pix$i[k]] %>% gr.string,
                                        gr[pix$j[k]] %>% gr.string))$res %>% as.data.table, mc.cores = mc.cores) %>% rbindlist
                                        #  res = unique(res, by = c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'))
  }
  else
  {
    res = reticulate::py_run_string(sprintf('res = data.matrix(balance = False, join = False, as_pixels = True)[:]'))$res %>% as.data.table
  }

  if (is.null(field))
    field = names(res)[-c(1:2)][1]

  dat = data.table(i = res[[1]]+1, j = res[[2]]+1, value = res[[field]])[, .(i = pmin(i, j), j = pmax(i, j), value = value)]

  gm = gM(bins, dat = dat)
  return(gm)
}


#' @name homeology.wrapper
#' @title homeology.wrapper
#' @param junctions Path to .vcf or bedpe or rds file of junctions or gGraph from which alt edges will be taken
#' @param width width around junction breakpoint around which to search for homeology
#' @param pad number of bases of padding around each sequence position (bin) to use when computing homeology, i.e. we then will be comparing 1 + 2*pad -mer sequences for edit distance
#' @param thresh string distance threshold for calling homeology in a bin
#' @param stride distance in bases between consecutive bins in which we will be measuring homeology
#' @param genome Path to .2bit or ffTrack .rds containing genome sequence
#' @param cores How many cores to use
#' @param flip if flip = FALSE, homeology search for -/- and +/+ junctions is done between a sequence and its reverse complement
#' @param bidirectional adding padding on both sides of each breakpoint (TRUE) or only in the direction of the fused side (FALSE)
#' @param annotate annotate edges in gGraph object and save it in working directory
#' @param savegMatrix save gMatrix object of edit distances
#' @param outdir output directory 
#'
#' @export homeology.wrapper
homeology.wrapper <- function(junctions,
                      width = 50,
                      pad = 0,
                      thresh = 0,
                      stride = 0,
                      genome = "~/DB/GATK/human_g1k_v37.fasta.2bit",
                      cores,
                      flip = FALSE,
                      bidirectional = TRUE,
                      annotate = TRUE,
                      savegMatrix = TRUE,
                      outdir = "./") {

  setDTthreads(1)

  system(paste('mkdir -p', outdir))

  gg = tryCatch(gG(junctions = junctions), error = function(e) NULL)

  if (is.null(gg))
    gg = tryCatch(readRDS(junctions), error = function(e) NULL)

  if (!inherits(gg, "gGraph")) gg = gG(junctions = gg)

  junctions = tryCatch(gg$edges[type == "ALT"]$junctions, error = function(e) NULL)

  if (is.null(junctions) || !inherits(junctions, "Junction"))
    stop('Error reading junctions - input should be either vcf, bedpe, .rds to Junction object, GRangesList, or gGraph')


  standardchr = which(as.logical(as.character(seqnames(junctions$left)) %in% GenomeInfoDb::standardChromosomes(junctions$left) &
                                   as.character(seqnames(junctions$right)) %in% GenomeInfoDb::standardChromosomes(junctions$right)))


  message("junctions with either breakends mapped outside of standard chromosomes will be thrown out")
  message("TRUE = junctions with breakends mapped to standard chromosomes")
  print(NROW(standardchr))

  junctions = junctions[standardchr]

  si = seqinfo(TwoBitFile(genome))

  ll = gr.nochr(junctions$left)
  rr = gr.nochr(junctions$right)

  outofbounds = union(
    GenomicRanges:::get_out_of_bound_index(conform_si(ll, si) + (ceiling(width) / 2)),
    GenomicRanges:::get_out_of_bound_index(conform_si(rr, si) + (ceiling(width) / 2))
  )

  if (NROW(outofbounds) > 0) {
    message(NROW(outofbounds), " out of bound junctions found")
    junctions = junctions[-outofbounds]
  }

  events = data.table(bp1 = gr.string(gr.nochr(junctions$left)),
                      bp2 = gr.string(gr.nochr(junctions$right)))

  print(events)

  #browser()

  cmd=sprintf("homeology.event(events, pad = ceiling(%s/2), thresh = %s, stride = %s, pad2 = %s, genome = '%s', mc.cores = %s, bidirected_search = %s, flip = %s, save_gm = %s)", width, thresh, stride, pad, genome, cores, bidirectional, flip, savegMatrix)

  if (!file.exists(paste0(outdir, "/res.rds"))) {
    message("running ", cmd)
    res = et(cmd)
    message("finished querying for homeology, saving as intermediate")
    saveRDS(res, paste0(outdir, "/res.rds"), compress = FALSE)
  } else {
    message("results already exist...")
    message("reading in to perform post-processing")
    res = readRDS(paste0(outdir, "/res.rds"))
  }


  message("Post-processing")
  stat = res[[3]]

  print(stat)
  stat = as.data.table(cbind(stat, junctions$dt)) %>% dedup.cols

  keep = which(sapply(stat, class) %in% c('character', 'integer', 'numeric', 'factor', 'logical'))
  stat = stat[, keep, with = FALSE]

  if (annotate) {
    message("annotating gGraph edges with homeology features")
    added_cols = c("numfeat", "numfeat2", "numfeat5",
                   "maxfeat", "numfeat10", "numlines5",
                   "maxlines5", "numlines", "maxlines",
                   "maxcor", "numglines", "numglines_pm5",
                   "numglines_p10", "numglines_p20", "hlen")
    if (NROW(stat)) {
      ed.annotate = dplyr::select(khtools::dedup.cols(stat, remove = T), edge.id, matches("numglines"))
      ed.annotate = dplyr::select(khtools::dedup.cols(stat, remove = T), edge.id, one_of(added_cols))
      cols = grep("numglines", colnames(ed.annotate), value = T)
      cols = colnames(ed.annotate)[colnames(ed.annotate) %in% added_cols]
      for (col in cols)
        et(sprintf("gg$edges[as.character(ed.annotate$edge.id)]$mark(%s = ed.annotate$%s)", col, col))
    }
    message("saving annotated gGraph")
    saveRDS(gg, paste0(outdir, "/marked_gGraph.rds"))
  }

  message("Saving junction-level stats")
  fwrite(stat, paste(outdir, 'junctions.txt', sep = '/'), sep = '\t')
  saveRDS(stat, paste(outdir, 'junctions.rds', sep = '/'))

  rawstat = res[[2]]

  if (!inherits(rawstat, "data.table")) setDT(rawstat)

  if (NROW(rawstat)) {
    rawstat = rawstat[
      , which(sapply(rawstat, class) %in%
                c('character', 'integer', 'numeric', 'factor', 'logical')),with = FALSE]
    if (is.null(rawstat$k)) rawstat$k = NA_integer_
    rawstat[, featuid := rleseq(seq, k, na.ignore = T, clump = T)$idx]
    if (is.null(rawstat$featuid)) rawstat$featuid = NA_integer_
    rawstat = merge(rawstat, stat[, .(seq, edge.id)], by = "seq")
  }

  message("Saving junction feature-level stats")
  fwrite(rawstat, paste(outdir, 'rawstats.txt', sep = "/"), sep = '\t')

  saveRDS(rawstat, paste(outdir, 'rawstats.rds', sep = '/'))

  if (savegMatrix) {
    message("Saving gMatrix")
    saveRDS(res[[1]], paste(outdir, 'gMatrix.rds', sep = '/'))
  }

  return(list(gMList = res[[1]], rawstat = rawstat, stat = stat))
}

#' @name homeology.event
#' @title homeology.event
#'
#' @param event data.table of junctions 
#' @param pad width around junction breakpoint around which to search for homeology
#' @param thresh string distance threshold for calling homeology in a bin
#' @param stride distance in bases between consecutive bins in which we will be measuring homeology
#' @param pad2 number of bases of padding around each sequence position (bin) to use when computing homeology
#' @param flip if flip = FALSE, homeology search for -/- and +/+ junctions is done between a sequence and its reverse complement
#' @param mc.cores number of cores to use
#' @param anchor if TRUE, anchor the junctions to the genome 
#' @param deanchor_gm if TRUE, deanchor the gMatrix
#' @param mat if TRUE, return the gMatrix as a matrix
#' @param genome path to .2bit or ffTrack .rds containing genome sequence
#' @param bidirected_search if TRUE, add padding on both sides of each breakpoint, if FALSE, add padding only in the direction of the fused side
#' @param save_gm if TRUE, save the gMatrix
#'
#' @return
#' @export homeology.event
#'
#' @examples
homeology.event = function (event, 
                pad = 100, 
                thresh = 2, 
                stride = 1, 
                pad2 = 5,
                flip = FALSE, 
                mc.cores = 1, 
                anchor = TRUE, 
                deanchor_gm = TRUE,
                mat = FALSE,
                genome = "~/DB/GATK/human_g1k_v37.fasta.2bit",
                bidirected_search = TRUE,
                save_gm = TRUE)
{
  if (!NROW(event)) {
    return(list(gm = list(), rawres = data.table(), res = data.table()))
  }
  event = copy2(event)
  ## event = data.table(bp1 = gr.string(event$left), bp2 = gr.string(event$right))
  
  gmfeat = function(gm, thresh = 3, op = "<=") {
    if (is(gm, "gMatrix")) {
      mat = gm$mat
      mat = (mat + t(mat))/2
    }
    else mat = gm
    ## dehom = deanchor(jJ(GRangesList(c(dg(evbp1)[dg(i)], dg(evbp2)[dg(i)]))),  dg(hom))
    ## demat = dehom$mat[dg(ix1), dg(ix2)]
    im = imager::as.cimg(as.matrix(mat))
    cmd = sprintf("im %s %s", op, thresh)
    px = eval(parse(text = cmd))
    if (sum(px) == 0)
      return(data.table())
    sp = imager::split_connected(px)
    if (NROW(sp) == 0)
      return(data.table())
    relmat = (1 - (mat / max(mat, na.rm = T)))

    res = rbindlist(lapply(
      1:length(sp),
      function(k) as.data.table(which(as.matrix(sp[[k]]),
                                      arr.ind = TRUE))[
                                        ,.(k = k, i = pmin(row, col),
                                           j = pmax(row, col), levd = mat[cbind(row, col)],
                                           relsim = relmat[cbind(row, col)],
                                           iw = rownames(mat)[row],
                                           jw = colnames(mat)[col])]), fill = TRUE)
    if (NROW(res) > 0) {
      griw = gr2dt(with(parse.gr2(res$iw, meta = res), {
        SHIFT = abs(min(pmin(start, end))) + 1;
        GenomicRanges::shift(dynGet("data"), SHIFT) %>% gr.spreduce(k)
      }))[, .(iwid = sum(width)), by = .(k = as.integer(k))]
      grjw = gr2dt(with(parse.gr2(res$jw, meta = res), {
        SHIFT = abs(min(pmin(start, end))) + 1;
        GenomicRanges::shift(dynGet("data"), SHIFT) %>% gr.spreduce(k)
      }))[, .(jwid = sum(width)), by = .(k = as.integer(k))]
      res = merge.data.table(merge.data.table(res, griw, by = 'k', all.x = TRUE),
                             grjw, by = "k", all.x = TRUE)[, minpx := pmin(iwid, jwid)]
      res[, `:=`(N = .N, r = cor(i, j)), by = k]
    }
    return(res)
  }
  
  #stats function for gmatrix
  gmstats = function(res) {
    if (NROW(res) > 0) {
      res[!duplicated(k)][
        ,.(numfeat = sum(N > 0), numfeat2 = sum(N > 2), numfeat5 = sum(N > 5),
           numfeat10 = sum(N > 10), maxfeat = max(N),
           numlines5 = sum(r > 0.5, na.rm = TRUE),
           maxlines5 = max(c(0, N[r > 0.5]), na.rm = TRUE),
           numlines = sum(r > 0.9, na.rm = TRUE),
           maxlines = max(c(0, N[r > 0.9]), na.rm = TRUE),
           maxcor = max(r, na.rm = TRUE),
           numglines = sum(na2false(r > 0.9) &
                             na2false(N >= 24) &
                             na2false(floor(N / minpx) <= 4)),
           numglines_pm5 = sum(na2false(r > 0.9) & na2false(minpx >= 8)),
           numglines_p10 = sum(na2false(r > 0.9) & na2false(minpx >= 15)),
           numglines_p20 = sum(na2false(r > 0.9) & na2false(minpx >= 25)),
           hlen = max(ifelse(na2false(r > 0.9), minpx, 0L)))]
    }
    else data.table(numfeat = 0, maxfeat = 0)
  }
  
  evbp1 = gr.end(gr.flipstrand(parse.gr(event$bp1)))
  if (flip)
    evbp2 = gr.flipstrand(gr.end(parse.gr(event$bp2)))
  else evbp2 = gr.end(parse.gr(event$bp2))


  if (isTRUE(bidirected_search)) {
    win = c(GRanges("Left:0") + pad, GRanges("Right:0") + pad)
    query.bp1 = evbp1 + pad
    query.bp2 = evbp2 + pad
  } else {
    win = c(gr.resize(GRanges("Left:0"), fix = "end", pad * 2, pad = F),
            gr.resize(GRanges("Right:0"), fix = "start", pad * 2, pad = F))
    query.bp1 = gr.resize(evbp1, pad * 2, pad = FALSE, fix = "end") # fix at end because was flipped
    bp2dir = if (flip) "end" else "start"
    query.bp2 = gr.resize(evbp2, pad * 2, pad = FALSE, fix = bp2dir)
  } # querying by both sides of junction, or just looking in the direction of the fused side of junction
  
  event$query.bp1 = gr.string(query.bp1)
  event$query.bp2 = gr.string(query.bp2)
  
  seq1 = ffTrack::get_seq(genome,
                          query.bp1)
  seq2 = ffTrack::get_seq(genome,
                          query.bp2)


  ifun = function(i) {
    tryCatch({
      message(i)
      this.env = environment()
      if (!anchor)
        win = c(evbp1[i], evbp2[i]) + pad
      win$seq = c(seq1[i], seq2[i])
      hom = homeology(win, stride = stride, pad = pad2)
      ix1 = which(hom$gr %^% win[1])
      ix2 = which(hom$gr %^% win[2])
      gm = gmfeat(hom$mat[ix1, ix2], thresh = thresh)
      gms = gmstats(gm)
      gm[, seq := this.env$i]
      gms[, seq := this.env$i]
      if (anchor && deanchor_gm)
        hom = deanchor(jJ(GRangesList(grbind(evbp1[i], evbp2[i]))), hom)
      return(list(hom = hom, gm = gm, gms = gms))
    }, error = function(e) printerr(i))
  }

  lst = mclapply(1:length(seq1), ifun, mc.cores = mc.cores)
  lst = purrr::transpose(lst)

  rawres = merge.data.table(event[, seq := seq_len(.N)], as.data.table(rbindlist(lst[[2]], fill = T)), by = "seq", all.x = TRUE)
  res = cbind(event, rbindlist(lst[[3]], fill = TRUE))

  if (save_gm)
    return(list(gm = lst[[1]], rawres = rawres, res = res))
  else
    return(list(gm = NULL, rawres = rawres, res = res))
}
