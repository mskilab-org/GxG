#' @name  GxG
#' @title GxG
#' @description
#' Matrices and pairwise features in a GenomicRanges / data.table framework
#'
#' Copyright (C) 2018 Marcin Imielinski
#'
#'    This program is free software: you can redistribute it and/or modify
#'    it under the terms of the GNU General Public License as published by
#'    the Free Software Foundation, either version 3 of the License, or
#'    (at your option) any later version.
#'
#'    This program is distributed in the hope that it will be useful,
#'    but WITHOUT ANY WARRANTY; without even the implied warranty of
#'    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#'    GNU General Public License for more details.
#'
#'    You should have received a copy of the GNU General Public License
#'    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#'
#'    Github: https://github.com/mskilab/GxG
#'
#' @import methods
#' @import R6
#' @import data.table
#' @import Matrix
#' @import GenomicRanges
#' @import igraph
#' @importFrom reshape2 melt
#' @import gUtils
#' @import gTrack
"_PACKAGE"


## ================== gMatrix class definition ================== ##
#' @export
gMatrix = setClass("gMatrix")
gMatrix = R6::R6Class("gMatrix",
                    public = list(

                      #' @name GMATRIX class constructor
                      #' @description 
                      #' Builds gPair object from length n granges and dat nxn matrix or data.table with columns
                      #' $i $j and $value, where $i and $j are integers between 1 and n indexing gr
                      #' @param gr GRanges 
                      #' @param dat length(gr) matrix or data.table / data.frame with field $i, $j, and $value indexing gr
                      #' @param fill fill value (default 0) for missing entries
                      #' @param full whether to explicitly set missing entries to fill (warning setting full = FALSE may distort certain calculations, though safe for summing etc)
                      #' @author Marcin Imielinski                         
                      initialize = function(gr = NULL, dat = NULL, fill = 0, full = NULL, na.rm = TRUE, agg.fun = sum)
                      {
                        private$pdat = data.table()
                        private$pgr = NULL
                        private$pfill = fill

                        if (!is.null(gr))
                          private$pgr = gr.stripstrand(gr[, c()])

                        if (is.null(full))
                          private$pfull = is.na(fill) || fill != 0 ## full by default if fill not NA and != 0
                        else
                          private$pfull = full

                        private$pna.rm = na.rm

                        if (!is.null(dat))
                        {
                          if (inherits(dat, 'data.frame'))
                          {
                            if (!all(c('i', 'j', 'value') %in% colnames(dat)))
                            {
                              stop('missing column values in data.table dat, either provide matrix or data.table with fields $i, $j, and $value')
                            }

                            if (any(dat$i>length(gr))| any(dat$j>length(gr)))
                            {
                              stop("out of bounds indices in dat data.table / data.frame, indices $i and $j are integers that index gr, and thus can't exceed the length")
                            }
                          } else
                          {
                            if (!is.matrix(dat) | !is.array(dat))
                              stop('dat must be either matrix, data.table, or data.frame')

                            if (!identical(length(gr), unique(dim(dat))))
                              stop('if dat is matrix it must be square with dimensions that are equal to length(gr)')

                            dat = melt(dat)
                            setnames(dat, c("i", "j", "value"))
                            private$pdat = pdat
                          }

                          dat = dat[, .(i = pmin(i, j), j = pmax(i, j), value)]

                          ## coerce fill (eg if NA) to proper classe
                          private$pfill = as(private$pfill, class(dat$value))

                          setkeyv(dat, c('i', 'j'))
                        }

                        if (!is.null(gr))
                          {
                            if (!is(gr, 'GRanges'))
                              stop('gr must be null or GRanges')

                            gr = gr.stripstrand(private$pgr)
                            dgr = disjoin(gr)
                            dgr = dgr[order(gr.match(dgr, gr))] ## resort in order of gr
                            is.diff = !identical(gr, dgr)
                            is.dup = any(duplicated(dat, by = c("i", "j")))
                            if (is.diff | is.dup) ## then we need to aggregate or collapse
                            {
                              if (is.diff)
                              {
                                warning('input intervals had to be disjoined, aggregating via agg.fun')
                                private$pgr = dgr[, c()]
                              }

                              if (!is.null(dat))
                              {
                                if (is.diff)
                                  gr$index = gr.match(gr, private$pgr)
                                else
                                  gr$index = 1:length(gr)

                                if (is.dup)
                                {
                                  warning('input data had some duplicates, aggregating via agg.fun')
                                }


                                dat = dat[, .(i = gr$index[i], j = gr$index[j], value)][, .(value = private$pagg.fun(value)), keyby = .(i, j)][, .(i = pmin(i,j), j = pmax(i,j), value = value)]
                                setkeyv(dat, c('i', 'j'))
                              }
                            }

                            ## now we "fill in" dat with fill value (ie we use full representation of matrix
                            ## but in data.table format) marking missing entries with "missing" column

                            if (private$pfull) ## fill out missing entries only only if fill != 0
                              {
                                private$pdat = as.data.table(expand.grid(i = 1:length(gr), j = 1:length(gr)))[i<=j, ][, value := private$pfill]
                                setkeyv(private$pdat, c('i', 'j'))
                                private$pdat[.(dat$i, dat$j), value := dat$value]
#                                private$pdat[.(dat$j, dat$i), value := dat$value] ## symmetrify
                                private$pdat[, id := 1:.N]
                              }
                            else if (!is.null(dat))
                              {
                                ##                                private$pdat = rbind(dat, dat[i!=j, .(j = i, i = j, value)])
                                private$pdat = dat[, id := 1:.N]
                                ##                                setkeyv(private$pdat, c('i', 'j'))
                              }
                          }

                        return(self)
                      },

                      #' @name disjoin
                      #' @description
                      #' disjoins current gMatrix definition with a set of input
                      #' ranges, distributing values appropriately 
                      #' @param gmats gMatix or list of gMatrices
                      #' @author Marcin Imielinski                         
                      disjoin = function(gr = NULL)
                      {
                        if (is.null(gr))
                          return(invisible(self))
                        dgr = disjoin(grbind(private$pgr, gr))
                        gr.new = dgr %*% private$pgr
                        gr.new$index = 1:length(gr.new)
                        grdt = gr2dt(gr.new)
                        dt.new = merge(
                          merge(private$pdat, grdt[, .(i = subject.id, inew = index)], by = 'i',
                                allow.cartesian = TRUE),
                          grdt[, .(j = subject.id, jnew = index)], by = 'j',
                          allow.cartesian = TRUE)[, .(
                                        i = inew,
                                        j = jnew,
                                        value = value)][, id := 1:.N]

                        tmp = gMatrix$new(gr.new, dt.new, fill = self$fill, full = self$full,
                                          agg.fun = self$agg.fun, na.rm = self$na.rm)
                        private$pgr = tmp$gr[, c()]
                        private$pdat = tmp$dat
                        return(invisible(self))
                      },

                      #' @name sweep
                      #' @description
                      #' sweeps the vector with the provided operation from all values
                      #' (no margin to specify since square matrix, so applied to both margins)
                      #' TODO: maybe need to do in ? both directions
                      #' @param vec vector of length 1 or length(self$gr)
                      #' @param MARGIN to apply to, 1 (i) vs 2 (j) (will symmetrify anyway)
                      #' @param FUN 
                      #' @author Marcin Imielinski                         
                      sweep = function(vec, FUN = "-", MARGIN = c(1))
                      {
                        if (length(vec)!=1 & length(vec)!=length(private$pgr))
                          stop('length of vec must be 1 or length(mat)')

                        if (length(vec) == 1)
                          vec = rep(vec, length(gr))

                        if (1 %in% MARGIN)
                          private$pdat[, value := eval(parse(text = sprintf('"%s"(value, vec[i])', FUN)))]

                        if (2 %in% MARGIN)
                          private$pdat[, value := eval(parse(text = sprintf('"%s"(value, vec[j])', FUN)))]

                        return(invisible(self))
                      },

                      #' @name set
                      #' @description
                      #' sets matrix "patch" defined by gr pair to given (scalar) value
                      #' and if either gr1 or gr2 not specified then fills in the entire
                      #' "column", disjoining gMatrix if need be
                      #' @param gmats gMatrix or list of gMatrices
                      #' @author Marcin Imielinski                         
                      set = function(val, gr1 = NULL, gr2 = NULL)
                      {                        
                        if (is.null(gr1))
                          gr1 = private$pgr

                        if (is.null(gr2))
                          gr2 = private$pgr

                        ## we disjoin matrix to accomodate these new values
                        self$disjoin(grbind(gr1, gr2))
                        
                        ij = expand.grid(i = (gr1 %*% private$pgr)$subject.id,
                                         j = (gr2 %*% private$pgr)$subject.id)
                        
                        newval = private$pdat[ij, ]

                        ## these values already exist in the dat
                        if (any(!is.na(newval$id)))
                        {
                          oldval = newval[!is.na(id), .(i, j)]
                          private$pdat[oldval, value := val]
                        }

                        ## these values we need to add to pdat
                        if (any(is.na(newval$id)))
                        {
                          newval = newval[is.na(id), ]
                          newval[, value := val]
                          newval[, id := nrow(private$pdat)+1:nrow(newval)]
                          private$pdat = rbind(private$pdat, newval)
                        }

                        return(invisible(self))
                      },

                      #' @name agg
                      #' @description
                      #' aggregates, ie collapses or refactors gr to the
                      #' provided intervals, aggregating via agg.fun
                      #' and dividing by the total "area" (ie doing a weighted
                      #' average)
                      #' @param gmats gMatix or list of gMatrices
                      #' @author Marcin Imielinski                         
                      agg = function(gr = NULL, FUN = private$pagg.fun)
                      {
                        if (is.null(gr))
                          return(invisible(self))
                        gr.new = gr
                        ov = gr2dt(gr.new %*% private$pgr)
                        wid = as.numeric(width(private$pgr))
                        dt.tmp = merge(
                          merge(private$pdat,
                                ov[, .(i = subject.id, inew = query.id)], by = 'i'),
                          ov[, .(j = subject.id, jnew = query.id)], by = 'j')
                        dt.tmp[, area := wid[i]*wid[j]]
                        dt.new = dt.tmp[, .(value = FUN(value*area, na.rm = private$pna.rm) /
                                              FUN(area, na.rm = private$pna.rm)),
                                        by = .(i = inew, j = jnew)]

                        tmp = gMatrix$new(gr.new, dt.new, fill = self$fill, full = self$full,
                                          agg.fun = self$agg.fun)

                        private$pgr = tmp$gr
                        private$pdat = tmp$dat
                        return(invisible(self))
                      },

                      #' @name merge
                      #' @description
                      #' merge with (list of) other gmats
                      #' and applies vector input, scalar valued function (default sum)
                      #' to values of interval pairs that overlap
                      #' @param gmats gMatix or list of gMatrices
                      #' @author Marcin Imielinski                         
                      merge = function(gmats = NULL, FUN = self$agg.fun)
                      {
                        if (is.null(gmats))
                          return(self)
                        
                        if (!is.list(gmats))
                          gmats = c(list(self), list(gmats))

                        if (!all(sapply(gmats, function(x) class(x)[1])=='gMatrix'))
                          stop('gmats must be gMatrix or list of gMatrices')
                        
                        dts = rbindlist(lapply(1:length(gmats), function(x) gmats[[x]]$dat[, source := x]))
                        grs = do.call(grbind, lapply(1:length(gmats), function(i) {gr = gmats[[i]]$gr; gr$source = i; gr$old = 1:length(gr); gr}))
                        gr.new = disjoin(grs)
                        grs$index = match(grs, gr.new)
                        tmp.dt = data.table(index = grs$index, source = grs$source, old = grs$old)
                        setkeyv(tmp.dt, c('source', 'old'))
                        dts$i = tmp.dt[.(dts$source, dts$i), index]
                        dts$j = tmp.dt[.(dts$source, dts$j), index]
                        dt.new = dts[, .(value = FUN(value)), by = .(i, j)]
                        dt.new[, id := 1:.N]
                       
                        tmp = gMatrix$new(grs, dt.new, fill = self$fill,
                                           full = self$full, agg.fun = self$agg.fun,
                                           na.rm = self$na.rm
                                          )

                        private$pgr = tmp$gr
                        private$pdat = tmp$dat
                        return(invisible(self))
                      },

                      #' @name subset
                      #' @description
                      #' subsets in place
                      #' Allows subseting of the gPair object using bracket notation
                      #' @param i integer or self$length logical vector specifying subset
                      #' @author Marcin Imielinski                         
                      subset = function(gr1 = NULL, gr2 = NULL)
                      {
                        grb = grbind(gr1, gr2)

                        gr.new = disjoin(grb) %*% private$pgr
                        gr.new = gr.new[order(gr.match(gr.new, grb))] ## sort gr.new with respect to grb
                        gr.new$index = 1:length(gr.new)
                        grdt = gr2dt(gr.new)
                        dt.new = merge(
                          merge(private$pdat, grdt[, .(i = subject.id, inew = index)], by = 'i',
                                allow.cartesian = TRUE),
                          grdt[, .(j = subject.id, jnew = index)], by = 'j',
                          allow.cartesian = TRUE)[, .(
                                        i = inew,
                                        j = jnew,
                                        value = value)][, id := 1:.N]

                      
                        ikeep = gr.new$index[gr.new %^% gr1]
                        jkeep = gr.new$index[gr.new %^% gr2]

                        ## we want dt.new that either have ikeep jkeep
                        ## or jkeep ikeep combo (but not jkeep jkeep
                        ## combo or ikeep ikeep combo)
                        dt.new = dt.new[(i %in% ikeep & j %in% jkeep) |
                               (i %in% jkeep & j %in% ikeep), ]
                        
                        return(gMatrix$new(gr.new[, c()], dt.new, full = self$full, na.rm = self$na.rm, agg.fun = self$agg.fun, fill = self$fill))
                      },
                      
                      #' @name gtrack
                      #' @description
                      #' Outputs a gTrack of this matrix
                      #' @param colormap colors to use with this colormap, default = c('white', 'red', 'black')
                      #' @param clim length 2 vector color limits to which to assign the min and max color
                      #' @author Marcin Imielinski                         
                      gtrack = function(colormap = c('white', 'red', 'black'),
                                        clim = NA, quantile = 0.01, ...)
                      {
                        dat = private$pdat
                        cmap.min = clim[1]
                        cmap.max = clim[2]

                        if (is.na(cmap.min))
                          cmap.min = quantile(dat$value, pmin(quantile, 1-quantile))

                        if (is.na(cmap.max))
                          cmap.max = quantile(dat$value, pmax(1-quantile, quantile))

                        mdata = Matrix::sparseMatrix(dat$i, dat$j, x = dat$value,
                                                     dims = c(length(self$gr), length(self$gr)),
                                                     symmetric = TRUE)
                        return(gTrack(private$pgr, mdata = mdata, colormap = colormap,
                                      cmap.min = cmap.min, cmap.max = cmap.max, ...))
                      },
                      
                      #' @name peaks
                      #' @description
                      #'
                      #' Returns gMatrix with TRUE value where there is a peak
                      #' (ie each pixel is >= than it's neighbors,
                      #' by default full = FALSE (ie only TRUE entries are returned in the dat
                      #' structure
                      #'
                      #' @param reverse will look at valleys instead of peaks
                      #' @return na.rm variable associated with this gMatrix object
                      #' @author Marcin Imielinski                         
                      peaks = function(reverse = FALSE, full = FALSE)
                      {
                        pdat = self$dat
                        ## test all neighbors
                        .nb = function(pdat, ii= 0, jj = 0)
                        {
                          nb = pdat[.(pdat$i+ii, pdat$j + jj), ]
                          if (any(ix <- is.na(nb$id)))
                            nb[ix, value := self$fill]

                          ## check seqnames matching between pdat and nb, NA if not matching
                          sn = data.table(id = 1:length(self$gr), sn = as.character(seqnames(self$gr)), key = 'id')
                          to.na = sn[.(pdat$i), sn] != sn[.(nb$i), sn] | sn[.(pdat$j), sn] != sn[.(nb$j), sn]
                          if (any(to.na, na.rm = TRUE))
                            nb[to.na, value := NA]
                          nb$value
                        }

                        ## only 5 neighbors after we take into account diagonal symmetry
                        if (reverse)
                          pdat$is.peak = pdat$value <= .nb(pdat, 1,0) &
                            pdat$value <= .nb(pdat, 1,1) &
                            pdat$value <= .nb(pdat, 1,-1) &
                            pdat$value <= .nb(pdat, -1,-1) &
                            pdat$value <= .nb(pdat, -1,0) 
                        else
                          pdat$is.peak = pdat$value >= .nb(pdat, 1,0) &
                            pdat$value >= .nb(pdat, 1,1) &
                            pdat$value >= .nb(pdat, 1,-1) &
                            pdat$value >= .nb(pdat, -1,-1) &
                            pdat$value >= .nb(pdat, -1,0) 
                        
                        dat.new = pdat[is.peak == TRUE, .(i, j, value = is.peak)]
                        
                        tmp = gMatrix$new(self$gr, dat.new, fill = self$fill,
                                          full = full, agg.fun = self$agg.fun)

                        private$pgr = tmp$gr
                        private$pdat = tmp$dat
                        return(invisible(self))
                      },

                      #' @name pack
                      #' @description
                      #'
                      #' "packs" gMatrix onto a single chromosome by concatenating
                      #' GRanges (even if they have different seqnames)
                      #' onto a single seqname, specified by the seqname arg
                      #' @param seqname seqname to recast the ranges into 
                      #' @author Marcin Imielinski                         
                      pack = function(seqname = '1')
                      {
                        private$pgr = dt2gr(gr2dt(private$pgr)[, .(seqnames = seqname, start = cumsum(width)-width+1,
                                                             end = cumsum(width))])
                        return(invisible(self))
                      },

                      #' @name mat
                      #' @description
                      #'
                      #' Returns matrix form of dat, i.e. of granges pairs -> values in sparseMatrix format (self$fill == 0)
                      #' full matrix if (self$fill != 0) or
                      #  or if (if baselevel == TRUE) of bp pairs (warning this matrix can be huge)
                      #' 
                      #' @param reverse will look at valleys instead of peaks
                      #' @author Marcin Imielinski                         
                      mat = function(baselevel = FALSE, upperdiag = TRUE)
                      {
                        if (length(private$pgr)==0)
                          return(sparseMatrix())

                        if (!baselevel)
                        {
                          if (self$fill ==0)
                            out = sparseMatrix(1, 1, x = 0, dims = rep(length(private$pgr),2))
                          else
                            out = matrix(self$fill, nrow = length(private$pgr), ncol = length(private$pgr))

                          if (!is.null(private$pdat))
                          {
                            tmp = private$pdat[value!=0, ]
                            if (nrow(tmp)>0)
                              {
                                out[cbind(tmp$i, tmp$j)] = tmp$value
                                if (!upperdiag)
                                  out[cbind(tmp$j, tmp$i)] = tmp$value
                              }
                            
                            tmp = private$pdat[is.na(value), ]
                            if (nrow(tmp)>0)
                              {
                                out[cbind(tmp$i, tmp$j)] =NA

                                if (!upperdiag)
                                  out[cbind(tmp$j, tmp$i)] = tmp$value
                              }
                          }
                          colnames(out) = rownames(out) = gr.string(private$pgr)
                        }
                        else
                        {
                          stop('not yet supported!')
                        }
                        
                        return(out)
                      },

                      #' @name print
                      #' @description
                      #' 
                      #' Prints out gMatrix Object. Prints the length and the GRangesList of the junctions.
                      #' @author Marcin Imielinski                         
                      print = function()
                      {
                        ngr = length(private$pgr)
                        ndt = length(private$pdat)                           
                        message(sprintf('gMatrix with %s intervals and %s entries',
                                        ngr, ndt), appendLF = FALSE)
                        
                        if (length(private$pgr) > 0){
                          message(' comprising:\n')
                          print(self$gr)

                          if (!is.null(private$pdat))
                            print(self$dat)
                        }
                        else
                          message('\n')
                      }
                    ),

                    private = list(
                      pgr = NULL, ## non-overlapping ranges
                      pdat = NULL, ## data.table
                      pfill = 0,
                      pfull = FALSE, ## whether matrix is full ie has all values represented
                      pagg.fun = sum,
                      pna.rm = TRUE
                    ),

                    active = list(
                      #' @name fill
                      #' @description
                      #'
                      #' Returns fill
                      #'
                      #' @return fill variable associated with this gMatrix object
                      #' @author Marcin Imielinski                         
                      fill = function() private$pfill,

                      #' @name footprint
                      #' @description
                      #'
                      #' Returns GRanges footprint of gMatrix
                      #'
                      #' @return fill variable associated with this gMatrix object
                      #' @author Marcin Imielinski                         
                      footprint = function() reduce(private$pgr),

                      #' @name gt
                      #' @description
                      #'
                      #' Returns default gTrack
                      #'
                      #' @return default gTrack
                      #' @author Marcin Imielinski                         
                      gt = function() self$gtrack(),

                      #' @name full
                      #' @description
                      #'
                      #' Returns full
                      #'
                      #' @return full variable associated with this gMatrix object
                      #' @author Marcin Imielinski                         
                      full = function() private$pfull,

                      #' @name agg.fun
                      #' @description
                      #'
                      #' Returns agg.fun
                      #'
                      #' @return agg.fun variable associated with this gMatrix object
                      #' @author Marcin Imielinski                         
                      agg.fun = function() private$pagg.fun,

                      #' @name na.rm
                      #' @description
                      #'
                      #' Returns na.rm
                      #'
                      #' @return na.rm variable associated with this gMatrix object
                      #' @author Marcin Imielinski                         
                      na.rm = function() private$pna.rm,
                      
                      #' @name gr
                      #' @description
                      #'
                      #' Returns or sets the GRanges associated with this gMatrix object,
                      #' useful for "refactoring" coordinates e.g. around a common
                      #' coordinate space or "squishing" matrices into a contiguous set of data.
                      #' even distorting the width of the coordinates
                      #'
                      #' @param gr must be the same length (but not the same width)
                      #' @return gets or sets the GRanges associated with this gMatrix object
                      #' @author Marcin Imielinski                         
                      gr = function(gr)
                      {
                        if (!missing(gr))
                        {
                          if (length(gr) != length(private$pgr))
                            stop('Length mismatch with inputted gr and gr of this objects: lengths (But not widths) must match')
                          tmp = gMatrix$new(gr, private$pdat, fill = self$fill,
                                            full = self$full, agg.fun = self$agg.fun, na.rm = self$na.rm)
                          private$pgr = tmp$gr
                          private$pdat = tmp$dat

                          return(invisible(self))
                        }
                        else
                        {
                          return(private$pgr)
                        }
                      },
                      
                      #' @name length
                      #' @description
                      #' 
                      #' Returns the number of junctions in this gMatrix Object.
                      #' 
                      #' @return Number of junctions in this gMatrix Object
                      #' @author Marcin Imielinski                         
                      length = function()
                      {                               
                        return(length(private$pgr))
                      },

                      copy = function() self$clone(),

                      #' @name dat
                      #' @description
                      #'
                      #' Returns the data.table
                      #'
                      #' @return data.table associated with this object.
                      #' @author Marcin Imielinski                         
                      dat = function()
                      {
                        return(copy(private$pdat))
                      },

                      value = function(val)
                      {
                        if(!missing(val)){
                          private$pdat[, value := val]
                        } else{
                          return(private$pdat$value)
                        }
                      }

                    )
                    )


## ================== Non-Member Functions for gMatrix ================== ##

#' @name gM
#' @title gM
#' @description
#'
#' gMatrix instantiator 
#' 
#' @param gr GRanges around which to build a gMatrix
#' @param dat data.table of $i and $j indexing gr and field $value OR a vector valued expression involving terms i and j that will be applied to expand.grid
#' @param full logical flag whether to explicitly "fill" missing entries in provided data.table (default = FALSE)
#' @param fill value with which to fill missing values (default 0)
#' @param agg.fun default function with which to aggregate values, this function should take a vector and na.rm = TRUE and return a scalar
#' @param na.rm default na.rm argument to agg.fun
#' @return A new gMatrix object
#' @author Marcin Imielinski                         
#' @export
gM = function(gr = NULL, dat = NULL, full = FALSE, fill = 0, agg.fun = sum, na.rm = TRUE)
{
  if (deparse(substitute(dat)) != "NULL")
  {
    is.dat = tryCatch(is.data.table(dat) | inherits(dat, 'matrix') | inherits(dat, 'Matrix') | inherits(dat, 'array'), error = function(e) FALSE)
    if (!is.dat)
    {
        dat.new = as.data.table(expand.grid(i = 1:length(gr), j = 1:length(gr)))[i<=j, ]
        val = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(dat)))), parent.frame()), dat.new, parent.frame(2)), error = function(e) NULL)
        dat.new$value = val
        dat = dat.new
    }
    else
    {
      if (inherits(dat, 'matrix') | inherits(dat, 'Matrix') | inherits(dat, 'array'))
      {
        ij = which(dat!=fill, arr.ind = TRUE)
        dat = cbind(as.data.table(ij)[, .(i = row, j = col)], data.table(value = dat[ij]))
      }
    }

    if (is.null(dat))
    {
      stop('dat argument is malformed, should either be data.table with fields $i, $j, and $value OR expression involving $i and $j that will be evaluated on the expanded grid of ij pairs')
    }
  }

  return(gMatrix$new(gr, dat, full = full, fill = fill, agg.fun = agg.fun, na.rm = na.rm))
}


#' @name [.gMatrix
#' @title gMatrix
#' @description
#'
#' Overloads subset operator for gMatrix. Allows subsetting of gMatrix via index.
#' Also, allows for data.table style queries on gMatrix GRangesList and
#' corresponding metadata.
#'
#' @param obj gMatrix object This is the gMatrix object to be subset
#' @param i integer, logical, or expression in gMatrix metadata used to subset gMatrix
#' @return A new gMatrix object that contains only the given id's
#' @author Marcin Imielinski                         
#' @export
'[.gMatrix' = function(obj, i = NULL, j = NULL, full = FALSE){

  if (is.null(j) & is.null(i))
    return(clone(obj))

  if (inherits(i, 'GRanges') | inherits(j, 'GRanges'))
   {
    if (!is.null(j) && !inherits(j, 'GRanges'))
      stop('if i is GRanges then j must be NULL or GRanges')

    if (!is.null(i) && !inherits(i, 'GRanges'))
      stop('if j is GRanges then i must be NULL or GRanges') 

    if (is.null(j))
      j = obj$gr

    if (is.null(i))
      i = obj$gr

    ret = obj$subset(i, j)
  } else
  {
    if (is.logical(i))
      i = which(i)
    
    if (is.logical(j))
      j = which(j)

    if (is.array(i))
      i = as.matrix(i)

    if (is.matrix(i))
      i = as.data.table(i)[, .(i = V1, j = V2)]

    if (is.data.frame(i))
      i = as.data.table(i)

    if (!is.data.table(i))
    {
      if (is.null(i))
        i = seq_along(obj$gr)
      
      if (is.null(j))
        j = seq_along(obj$gr)
      
      if (max(c(i, j))>length(obj$gr))
        stop('indices out of bounds')

      uind = unique(c(i, j))
      dt.new = merge(
        merge(obj$dat, data.table(i = uind, inew = 1:length(uind)),
              by = 'i', allow.cartesian = TRUE),              
        data.table(j = uind, jnew = 1:length(uind)), by = 'j', allow.cartesian = TRUE)[, .(
                      i = inew, j = jnew, value = value)][, id := 1:.N]

      ikeep = which(uind %in% i)
      jkeep = which(uind %in% j)

      dt.new = dt.new[(i %in% ikeep & j %in% jkeep) |
                               (i %in% jkeep & j %in% ikeep), ]

      gr.new = obj$gr[uind]
      ret = gMatrix$new(gr.new, dt.new, fill = obj$fill, full = obj$full)    
    }
    else
    {
      if (!is.null(j))
        warning('j value ignored if matrix or data.table input provided')

      if (max(c(i$i, i$j))>length(obj$gr))
        stop('indices out of bounds')
      
      uind = unique(i$i, i$j)
    
      if(is.null(i$i) || is.null(i$j))
        stop('indices must be integer or data.table / data.frame with fields $i and $j or two column matrix / array')

      ij = i
      
      ij = ij[, .(i = pmin(i, j), j = pmax(i, j))]
      
      ## map old indices to new "unique indices"
      map = data.table(old = uind)[, new := 1:.N]
      setkey(map, old)
      
      ## create new indices and dt.new
      gr.new = obj$gr[uind]
      dt.new = obj$dat[.(ij$i, ij$j), .(i = map[.(i), new], j = map[.(j), new], value, id)]
      
      if (!obj$full)
        dt.new = dt.new[!is.na(id), ]
      
      ret = gMatrix$new(gr.new, dt.new, fill = obj$fill, full = obj$full)    
    }
  }
  
  if (full)
    ret = ret$full

  return(ret)
}

#' @name grunif
#' @title grunif
#' @description
#'
#' gMatrix instantiator creates random matrix with uniform distribution 
#' 
#' @param gr GRanges around which to build a gMatrix
#' @param full logical flag whether to explicitly "fill" missing entries in provided data.table (default = FALSE)
#' @param fill value with which to fill missing values (default 0)
#' @param agg.fun default function with which to aggregate values, this function should take a vector and na.rm = TRUE and return a scalar
#' @param na.rm default na.rm argument to agg.fun
#' @export
grunif = function(gr, full = FALSE, fill = 0, agg.fun = sum, na.rm = TRUE)
{
  return(gM(gr = gr, dat = runif(length(i)), fill = fill, full = full, agg.fun = agg.fun, na.rm = na.rm))
}



#' @name grnorm
#' @title grnorm
#' @description
#'
#' gMatrix instantiator creates random matrix with normal distribution 
#' 
#' @param gr GRanges around which to build a gMatrix
#' @param mean mean parameter funtion of i and j
#' @param sd scalar sd parameter to rnorm
#' @param full logical flag whether to explicitly "fill" missing entries in provided data.table (default = FALSE)
#' @param fill value with which to fill missing values (default 0)
#' @param agg.fun default function with which to aggregate values, this function should take a vector and na.rm = TRUE and return a scalar
#' @param na.rm default na.rm argument to agg.fun
#' @export
grnorm = function(gr, mean = 0, sd = 1, full = FALSE, fill = 0, agg.fun = sum, na.rm = TRUE)
{
  init = gM(gr, dat = mean, fill = fill, full = full, agg.fun = agg.fun, na.rm = na.rm)
  new.dat = init$dat[, value := rnorm(value, sd = sd)]
  return(gM(gr, dat = new.dat, fill = fill, full = full, agg.fun = agg.fun, na.rm = na.rm))
}

#' @name grpois
#' @title grpois
#' @description
#'
#' gMatrix instantiator creates random matrix with poisson
#' 
#' @param gr GRanges around which to build a gMatrix
#' @param full logical flag whether to explicitly "fill" missing entries in provided data.table (default = FALSE)
#' @param fill value with which to fill missing values (default 0)
#' @param agg.fun default function with which to aggregate values, this function should take a vector and na.rm = TRUE and return a scalar
#' @param na.rm default na.rm argument to agg.fun
#' @export
grpois = function(gr, lambda = 1, full = FALSE, fill = 0, agg.fun = sum, na.rm = TRUE)
{
  init = gM(gr, dat = lambda, fill = fill, full = full, agg.fun = agg.fun, na.rm = na.rm)
  new.dat = init$dat[, value := rpois(n = length(value), lambda = value)]
  return(gM(gr, dat = new.dat, fill = fill, full = full, agg.fun = agg.fun, na.rm = na.rm))
}

#' @name cocount
#' @title cocount
#' @description
#'
#' gMatrix instantiator creates contact "counts" of "events", which are GRanges grouped into
#' pairs / triples / groups (e.g. paired end reads, linked readsa, SPRITE contacts) via "by"
#' column.
#' 
#' @param events GRanges of events
#' @param by column of events GRanges that allows co-counting
#' @param bins bins to count to, default is disjoin(events)
#' @param frac whether to count fractional overlap of event with bin as 1 (frac == FALSE) or based on width
#' @export
cocount = function(events, bins = disjoin(events), by = names(values(events))[1], frac = FALSE, full = FALSE, fill = 0, na.rm = TRUE)
{
  if (length(events)==0)
    return(gM(gr = bins, full = full, fill = fill, agg.fun = agg.fun, na.rm = na.rm))
  
  if (is.na(by))
    stop('by must be specified and a metadata column of events GRanges')
           
  if (!(by %in% names(values(events))))
    stop('by must be a metadata column of events GRanges')

  events$group = values(events)[[by]]
  tmp = gr2dt(events[, c("group")] %*% bins)
  
  if (frac)
    tmp[, weight := width/sum(width), by = group]
  else
    tmp[, weight := 1, by = group]

  tmp = tmp[, .(bid = subject.id, group = as.integer(factor(group)), weight)]

  ## sum weights inside bin pairs that share a group
  dat = merge(tmp, tmp, by = c('group'), allow.cartesian = TRUE)[, .(value = sum(weight.x * weight.y)), by = .(i = bid.x, j = bid.y)]

  return(gM(bins, dat, full = full, fill = fill, na.rm = na.rm, agg.fun = sum))
}




#' @name gmatalign
#' @description
#'
#' Utility function "aligns" two gMatrices around their disjoin
#' 
#' @param mat1 gMatrix object
#' @param mat2 gMatrix object
#' @return a new gMatrix object whose values are the sum of the inputs of the
#' @author Marcin Imielinski                         
gmatalign = function(mat1, mat2)
{

  if (identical(mat1$gr, mat2$gr))
    gr.new = mat1$gr
  else
    {
      gr.new = disjoin(grbind(mat1$gr, mat2$gr))
      
      ## these will have identical indices
      mat1  = mat1$copy$disjoin(gr.new)
      mat2  = mat2$copy$disjoin(gr.new)     
    }

  dat1 = mat1$dat[, .(i, j, id1 = id, val1 = value)]
  dat2 = mat2$dat[, .(i, j, id2 = id, val2 = value)]
 
  dat.new = merge(dat1, dat2, by = c('i', 'j'), all = TRUE)

  if (any(ix <- dat.new[, is.na(id1)]))
    dat.new[ix, val1 := mat1$fill]

  if (any(ix <- dat.new[, is.na(id2)]))
    dat.new[ix, val2 := mat2$fill]

  return(list(gr = gr.new, dat = dat.new))
}


#' @name +.gMatrix
#' @description
#'
#' Allows addition of gMatrices
#' 
#' @param mat1 gMatrix object
#' @param mat2 gMatrix object
#' @return a new gMatrix object whose values are the sum of the inputs of the
#' @author Marcin Imielinski                         
#' @export
'+.gMatrix' = function(mat1, mat2)
{
  if (is(mat1, 'gMatrix') & is(mat2, 'gMatrix'))
    {
      tmp = gmatalign(mat1, mat2)
      gr.new = tmp$gr
      dat.new = tmp$dat
    }
  else if (is(mat2, 'gMatrix'))
  {
    val1 = mat1
    mat1 = mat2
    gr.new = mat2$gr
    dat.new = mat2$dat
    dat.new[, val2 := value]
  }
  else 
  {
    val2 = mat2
    mat2 = mat1
    gr.new = mat1$gr
    dat.new = mat1$dat
    dat.new[, val1 := value]
  }

  dat.new[, value := val1 + val2]
  
  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg.fun,
                 na.rm = mat1$na.rm))
}

#' @name -.gMatrix
#' @description
#'
#' Allows subtraction of gMatrices
#' 
#' @param mat1 gMatrix object
#' @param mat2 gMatrix object
#' @return a new gMatrix object whose values are the sum of the inputs of the
#' @author Marcin Imielinski                         
#' @export
'-.gMatrix' = function(mat1, mat2)
{
  if (is(mat1, 'gMatrix') & is(mat2, 'gMatrix'))
  {
    tmp = gmatalign(mat1, mat2)
    gr.new = tmp$gr
    dat.new = tmp$dat
  }
  else if (is(mat2, 'gMatrix'))
  {
    val1 = mat1
    mat1 = mat2
    gr.new = mat2$gr
    dat.new = mat2$dat
    dat.new[, val2 := value]
  }
  else 
  {
    val2 = mat2
    mat2 = mat1
    gr.new = mat1$gr
    dat.new = mat1$dat
    dat.new[, val1 := value]
  }

  dat.new[, value := val1 - val2]
  
  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg.fun,
                 na.rm = mat1$na.rm))
}


#' @name *.gMatrix
#' @description
#'
#' Allows subtraction of gMatrices
#' 
#' @param mat1 gMatrix object
#' @param mat2 gMatrix object
#' @return a new gMatrix object whose values are the sum of the inputs of the
#' @author Marcin Imielinski                         
#' @export
'*.gMatrix' = function(mat1, mat2)
{
  if (is(mat1, 'gMatrix') & is(mat2, 'gMatrix'))
  {
    tmp = gmatalign(mat1, mat2)
    gr.new = tmp$gr
    dat.new = tmp$dat
  }
  else if (is(mat2, 'gMatrix'))
  {
    val1 = mat1
    mat1 = mat2
    gr.new = mat2$gr
    dat.new = mat2$dat
    dat.new[, val2 := value]
  }
  else 
  {
    val2 = mat2
    mat2 = mat1
    gr.new = mat1$gr
    dat.new = mat1$dat
    dat.new[, val1 := value]
  }

  dat.new[, value := val1 * val2]
  
  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg.fun,
                     na.rm = mat1$na.rm))
  gr.new = disjoin(grbind(mat1$gr, mat2$gr))
}

#' @name /.gMatrix
#' @description
#'
#' Allows division of gMatrices
#' 
#' @param mat1 gMatrix object
#' @param mat2 gMatrix object
#' @return a new gMatrix object whose values are the sum of the inputs of the
#' @author Marcin Imielinski                         
#' @export
'/.gMatrix' = function(mat1, mat2)
{
  if (is(mat1, 'gMatrix') & is(mat2, 'gMatrix'))
  {
    tmp = gmatalign(mat1, mat2)
    gr.new = tmp$gr
    dat.new = tmp$dat
  }
  else if (is(mat2, 'gMatrix'))
  {
    val1 = mat1
    mat1 = mat2
    gr.new = mat2$gr
    dat.new = mat2$dat
    dat.new[, val2 := value]
  }
  else 
  {
    val2 = mat2
    mat2 = mat1
    gr.new = mat1$gr
    dat.new = mat1$dat
    dat.new[, val1 := value]
  }

  dat.new[, value := val1 / val2]
  
  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg.fun,
                     na.rm = mat1$na.rm))
}

#' @name ==.gMatrix
#' @description
#'
#' Allows comparison of matrices
#' 
#' @param mat1 gMatrix object
#' @param mat2 gMatrix object
#' @return a new gMatrix object whose values are the sum of the inputs of the
#' @author Marcin Imielinski                         
#' @export
'==.gMatrix' = function(mat1, mat2)
{
  if (is(mat1, 'gMatrix') & is(mat2, 'gMatrix'))
  {
    tmp = gmatalign(mat1, mat2)
    gr.new = tmp$gr
    dat.new = tmp$dat
  }
  else if (is(mat2, 'gMatrix'))
  {
    val1 = mat1
    mat1 = mat2
    gr.new = mat2$gr
    dat.new = mat2$dat
    dat.new[, val2 := value]
  }
  else 
  {
    val2 = mat2
    mat2 = mat1
    gr.new = mat1$gr
    dat.new = mat1$dat
    dat.new[, val1 := value]
  }

  dat.new[, value := val1 == val2]
  
  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg.fun,
                     na.rm = mat1$na.rm))
}

#' @name <=.gMatrix
#' @description
#'
#' Allows comparison of matrices
#' 
#' @param mat1 gMatrix object
#' @param mat2 gMatrix object
#' @return a new gMatrix object whose values are the sum of the inputs of the
#' @author Marcin Imielinski                         
#' @export
'<=.gMatrix' = function(mat1, mat2)
{
  if (is(mat1, 'gMatrix') & is(mat2, 'gMatrix'))
  {
    tmp = gmatalign(mat1, mat2)
    gr.new = tmp$gr
    dat.new = tmp$dat
  }
  else if (is(mat2, 'gMatrix'))
  {
    val1 = mat1
    mat1 = mat2
    gr.new = mat2$gr
    dat.new = mat2$dat
    dat.new[, val2 := value]
  }
  else 
  {
    val2 = mat2
    mat2 = mat1
    gr.new = mat1$gr
    dat.new = mat1$dat
    dat.new[, val1 := value]
  }

  dat.new[, value := val1 <= val2]
  
  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg.fun,
                     na.rm = mat1$na.rm))
}

#' @name >=.gMatrix
#' @description
#'
#' Allows comparison of matrices
#' 
#' @param mat1 gMatrix object
#' @param mat2 gMatrix object
#' @return a new gMatrix object whose values are the sum of the inputs of the
#' @author Marcin Imielinski                         
#' @export
'>=.gMatrix' = function(mat1, mat2)
{
  if (is(mat1, 'gMatrix') & is(mat2, 'gMatrix'))
  {
    tmp = gmatalign(mat1, mat2)
    gr.new = tmp$gr
    dat.new = tmp$dat
  }
  else if (is(mat2, 'gMatrix'))
  {
    val1 = mat1
    mat1 = mat2
    gr.new = mat2$gr
    dat.new = mat2$dat
    dat.new[, val2 := value]
  }
  else 
  {
    val2 = mat2
    mat2 = mat1
    gr.new = mat1$gr
    dat.new = mat1$dat
    dat.new[, val1 := value]
  }

  dat.new[, value := val1 >= val2]
  
  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg.fun,
                     na.rm = mat1$na.rm))
}

#' @name >.gMatrix
#' @description
#'
#' Allows comparison of matrices
#' 
#' @param mat1 gMatrix object
#' @param mat2 gMatrix object
#' @return a new gMatrix object whose values are the sum of the inputs of the
#' @author Marcin Imielinski                         
#' @export
'>.gMatrix' = function(mat1, mat2)
{
  if (is(mat1, 'gMatrix') & is(mat2, 'gMatrix'))
  {
    tmp = gmatalign(mat1, mat2)
    gr.new = tmp$gr
    dat.new = tmp$dat
  }
  else if (is(mat2, 'gMatrix'))
  {
    val1 = mat1
    mat1 = mat2
    gr.new = mat2$gr
    dat.new = mat2$dat
    dat.new[, val2 := value]
  }
  else 
  {
    val2 = mat2
    mat2 = mat1
    gr.new = mat1$gr
    dat.new = mat1$dat
    dat.new[, val1 := value]
  }

  dat.new[, value := val1 > val2]
  
  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg.fun,
                     na.rm = mat1$na.rm))
}



#' @name <.gMatrix
#' @description
#'
#' Allows comparison of matrices
#' 
#' @param mat1 gMatrix object
#' @param mat2 gMatrix object
#' @return a new gMatrix object whose values are the sum of the inputs of the
#' @author Marcin Imielinski                         
#' @export
'<.gMatrix' = function(mat1, mat2)
{
  if (is(mat1, 'gMatrix') & is(mat2, 'gMatrix'))
  {
    tmp = gmatalign(mat1, mat2)
    gr.new = tmp$gr
    dat.new = tmp$dat
  }
  else if (is(mat2, 'gMatrix'))
  {
    val1 = mat1
    mat1 = mat2
    gr.new = mat2$gr
    dat.new = mat2$dat
    dat.new[, val2 := value]
  }
  else 
  {
    val2 = mat2
    mat2 = mat1
    gr.new = mat1$gr
    dat.new = mat1$dat
    dat.new[, val1 := value]
  }

  dat.new[, value := val1 < val2]
  
  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg.fun,
                     na.rm = mat1$na.rm))
}

#' @name !=.gMatrix
#' @description
#'
#' Allows comparison of matrices
#' 
#' @param mat1 gMatrix object
#' @param mat2 gMatrix object
#' @return a new gMatrix object whose values are the sum of the inputs of the
#' @author Marcin Imielinski                         
#' @export
'!=.gMatrix' = function(mat1, mat2)
{
  if (is(mat1, 'gMatrix') & is(mat2, 'gMatrix'))
  {
    tmp = gmatalign(mat1, mat2)
    gr.new = tmp$gr
    dat.new = tmp$dat
  }
  else if (is(mat2, 'gMatrix'))
  {
    val1 = mat1
    mat1 = mat2
    gr.new = mat2$gr
    dat.new = mat2$dat
    dat.new[, val2 := value]
  }
  else 
  {
    val2 = mat2
    mat2 = mat1
    gr.new = mat1$gr
    dat.new = mat1$dat
    dat.new[, val1 := value]
  }

  dat.new[, value := val1 != val2]
  
  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg.fun,
                     na.rm = mat1$na.rm))
}


#' @name &.gMatrix
#' @description
#'
#' Allows comparison of matrices
#' 
#' @param mat1 gMatrix object
#' @param mat2 gMatrix object
#' @return a new gMatrix object whose values are the sum of the inputs of the
#' @author Marcin Imielinski                         
#' @export
'&.gMatrix' = function(mat1, mat2)
{
  if (is(mat1, 'gMatrix') & is(mat2, 'gMatrix'))
  {
    tmp = gmatalign(mat1, mat2)
    gr.new = tmp$gr
    dat.new = tmp$dat
  }
  else if (is(mat2, 'gMatrix'))
  {
    val1 = mat1
    mat1 = mat2
    gr.new = mat2$gr
    dat.new = mat2$dat
    dat.new[, val2 := value]
  }
  else 
  {
    val2 = mat2
    mat2 = mat1
    gr.new = mat1$gr
    dat.new = mat1$dat
    dat.new[, val1 := value]
  }

  dat.new[, value := val1 & val2]
  
  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg.fun,
                     na.rm = mat1$na.rm))
}

#' @name |.gMatrix
#' @description
#'
#' Allows comparison of matrices
#' 
#' @param mat1 gMatrix object
#' @param mat2 gMatrix object
#' @return a new gMatrix object whose values are the sum of the inputs of the
#' @author Marcin Imielinski                         
#' @export
'|.gMatrix' = function(mat1, mat2)
{
  if (is(mat1, 'gMatrix') & is(mat2, 'gMatrix'))
  {
    tmp = gmatalign(mat1, mat2)
    gr.new = tmp$gr
    dat.new = tmp$dat
  }
  else if (is(mat2, 'gMatrix'))
  {
    val1 = mat1
    mat1 = mat2
    gr.new = mat2$gr
    dat.new = mat2$dat
    dat.new[, val2 := value]
  }
  else 
  {
    val2 = mat2
    mat2 = mat1
    gr.new = mat1$gr
    dat.new = mat1$dat
    dat.new[, val1 := value]
  }

  dat.new[, value := val1 | val2]
  
  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg.fun,
                     na.rm = mat1$na.rm))
}

#' @name !.gMatrix
#' @description
#'
#' Allows negation of matrix
#' 
#' @param mat gMatrix object
#' @return a new gMatrix object whose values are the sum of the inputs
#' @author Marcin Imielinski                         
#' @export
'!.gMatrix' = function(mat)
{
  dat.new = mat$dat[, value := !value]
  gr.new = mat$gr
  
  return(gMatrix$new(gr.new, dat.new, fill = mat$fill, full = mat$full, agg.fun = mat$agg.fun,
                     na.rm = mat$na.rm))
}

#' @name sign.gMatrix
#' @description
#'
#' Allows sign of matrix
#' 
#' @param mat gMatrix object
#' @return a new gMatrix object whose values are the sign of inputs
#' @author Marcin Imielinski                         
#' @export
'sign.gMatrix' = function(mat)
{
  dat.new = mat$dat[, value := sign(value)]
  gr.new = mat$gr
  
  return(gMatrix$new(gr.new, dat.new, fill = mat$fill, full = mat$full, agg.fun = mat$agg.fun,
                     na.rm = mat$na.rm))
}


#' @name log.gMatrix
#' @description
#'
#' Allows log of matrix values
#' 
#' @param mat gMatrix object
#' @return a new gMatrix object whose values are the log of the inputs
#' @author Marcin Imielinski                         
#' @export
'log.gMatrix' = function(mat)
{
  dat.new = mat$dat[, value := log(value)]
  gr.new = mat$gr
  
  return(gMatrix$new(gr.new, dat.new, fill = mat$fill, full = mat$full, agg.fun = mat$agg.fun,
                     na.rm = mat$na.rm))
}

#' @name all.gMatrix
#' @description
#'
#' Returns TRUE if all entries are TRUE
#' 
#' @param mat gMatrix object
#' @return a new gMatrix object whose values are the log of the inputs
#' @author Marcin Imielinski                         
#' @export
'all.gMatrix' = function(mat, na.rm = mat$na.rm)
{
  return(all(mat$value, na.rm = na.rm))
}

#' @name any.gMatrix
#' @description
#'
#' Returns TRUE if any entries are TRUE
#' 
#' @param mat gMatrix object
#' @return a new gMatrix object whose values are the log of the inputs
#' @author Marcin Imielinski                         
#' @export
'any.gMatrix' = function(mat, na.rm = mat$na.rm)
{
  return(any(mat$value, na.rm = na.rm))
}

#' @name dim.gMatrix
#' @description
#'
#' Dimensiox
#' 
#' @param mat gMatrix object
#' @return dimensions of matrix
#' @author Marcin Imielinski                         
#' @export
'dim.gMatrix' = function(mat)
{
  if (length(mat)==0)
    return(c(0, 0))

  rep(sum(as.numeric(width(mat$gr))), 2)
}


#' @name length.gMatrix
#' @description
#'
#' Dimensiox
#' 
#' @param mat gMatrix object
#' @return dimensions of matrix
#' @author Marcin Imielinski                         
#' @export
'length.gMatrix' = function(mat)
{
  length(mat$gr)
}


#' @name ^.gMatrix
#' @description
#'
#' Allows powers of values in matrices
#' 
#' @param mat1 gMatrix object
#' @param power numeric
#' @return a new gMatrix object whose values are the sum of the inputs of the
#' @author Marcin Imielinski                         
#' @export
'^.gMatrix' = function(mat1, power)
{
  val2 = power
  gr.new = mat1$gr
  dat.new = mat1$dat
  dat.new[, val1 := value]
  dat.new[, value := val1^power]
  
  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg.fun,
                     na.rm = mat1$na.rm))
}


#' @name refresh
#' @description
#' 
#' Updates gMatrix object to reflect changes in source code
#' 
#' @param gMatrix object
#' @return gMatrix object
#' @author Marcin Imielinski                         
#' @exportMethod refresh
#' @export
setGeneric("refresh", function(x) standardGeneric("refresh"))
setMethod("refresh", "gMatrix",
          function(x) {
            return(gMatrix$new(gr = x$gr, dat = x$dat, na.rm = x$na.rm, fill = x$fill, full = x$full, agg.fun = x$agg.fun))
          })



## ================== gPair class definition ================== ##
#' @export
gPair = setClass("gPair")
gPair = R6::R6Class("gPair",
                    public = list(

                      #' @name gPair class constructor
                      #' @description 
                      #' Builds gPair class, grl must be GRangesList with each GRanges of length 2 JUNCTION
                      #' Empty junctions are removed
                      #' If grl is empty GRangesList, gPairs class is empty
                      #' @param gr1 GRanges gr1 
                      #' @param gr2 GRanges gr2
                      #' @param gr1 optional metadata
                      #' @author Marcin Imielinski                         
                      initialize = function(gr1 = NULL, gr2 = NULL,
                                            meta = NULL, ...)
                      {

                        if (length(gr1) != length(gr1))
                          stop('Inputs gr1 and gr2 must be the same length')
                        
                        
                        if (length(gr1)>0)
                        {
                          if (nrow(meta)!=length(gr1))
                            stop('meta must have the same number of rows as gr1 and gr2')
                        }
                        
                        private$pgr1 = gr1
                        private$pgr2 = gr2
                        private$pmeta = meta                           
                        return(self)
                      },
                      
                      #' @name removeDups
                      #' @description
                      #'
                      #' Removes duplicate junctions in this gPair Object
                      #'
                      #' this gPair Object has no duplicate junctions
                      #' @author Marcin Imielinski
                      removeDups = function()

                      {
                        uix = !duplicated(cbind(as.data.table(gr1), as.data.table(gr2), meta))
                        private$pgr1 = private$pgr1[uix]
                        private$pgr2 = private$pgr2[uix]
                        private$pmeta = private$pmeta[uix,]
                        return(self)
                      },

                      #' @name subset
                      #' @description 
                      #' Allows subseting of the gPair object using bracket notation
                      #' @param i integer or self$length logical vector specifying subset 
                      subset = function(i)
                      {
                        if (is.null(i)){
                          i = integer()
                        }

                        if (is.logical(i)){
                          i = which(i)
                        }
                        if (length(i)>0 && (is.numeric(i) | is.integer(i))) {
                          if (max(i, na.rm = TRUE)>self$length) {
                            stop('index out of bounds')
                          }
                          
                          if (any(i<0))
                          {
                            if (!all(i<0)){
                              stop('cannot mix positive with negative subscripts for gPair object')
                            }
                            
                            i = setdiff(1:self$length, abs(i))
                          }
                        }
                        
                        private$pairs = private$pairs[i]
                        
                        return(self)
                      },
                      
                      #' @name print
                      #' @description
                      #' 
                      #' Prints out the gPair Object. Prints the length and the GRangesList of the junctions.
                      print = function()
                      {
                        message("gPair Object with ", self$length, " junctions\n")
                        if (self$length>0)
                        {
                          HEAD = 1:pmin(4, self$length)
                          pairdt = data.table(pair = paste(gr.string(private$pgr1[HEAD])), '<->',
                                              gr.string(private$pgr1[HEAD]))

                          if (ncol(self$dt)>0)
                          {
                            pairdt = cbind(bpdt[, "pair", with = FALSE], self$dt[HEAD, ])
                          }
                          print(pairdt)
                          more = self$length-HEAD[length(HEAD)]
                          if (more>0)
                            message('... (', more,' additional pairs)')
                        }
                      }
                    ),

                    private = list(
                      gr1 = NULL, ## second interval
                      gr2 = NULL, ## first interval
                      meta = NULL ## data.table
                    ),
                    
                    active = list(
                      #' @name grl
                      #' @description
                      #'
                      #' Returns the GRangesList of the junctions in this gPair Object
                      #'
                      #' @return GRangesList of the junctions in this gPair Object
                      #' @author Marcin Imielinski                         
                      grl = function()
                      {
                        grl = grl.pivot(list(private$pgr1, private$pgr2))
                        values(grl) = private$pgr1
                        return(grl)
                      },
                      
                      #' @name length
                      #' @description
                      #' 
                      #' Returns the number of junctions in this gPair Object.
                      #' 
                      #' @return Number of junctions in this gPair Object
                      length = function()
                      {                               
                        return(length(private$pairs))
                      },

                      copy = function() self$clone(),

                      #' @name dt
                      #' @description
                      #'
                      #' Returns the GRangesList of the junctions in the gPair Object as a data.table.
                      #'
                      #' @return data.table GRangesList of the junctions coverted to a data.table
                      dt = function()
                      {
                        return(private$pdt)
                      },

                      #' @name span
                      #' @description
                      #' 
                      #' Returns the distance between breakpoint pairs on the genome
                      #' 
                      #' @return vector of spans of all junctions in junction object
                      span = function()
                      {                               
                        ifelse(seqnames(private$pgr1) == seqnames(private$pgr2),
                               pmin(
                                 abs(start(private$pgr1)-end(private$pgr2)),
                                 abs(start(private$pgr1)-start(private$pgr2)),
                                 abs(end(private$pgr1)-end(private$pgr2)),
                                 abs(end(private$pgr1)-start(private$pgr2))),
                               Inf)                                 
                      },

                      #' @name sign
                      #' @description
                      #' 
                      #' Sign is the product of strands, i.e. ++ and -- is 1
                      #' and "+-" and "-+" is -1
                      #' 
                      #' @return vector of orientations of each junction in junction object
                      sign = function()
                      {                               
                        return(grl.eval(private$pairs, sign((sign(strand[1]=='+')-0.5)*(sign(strand[2]=='+')-0.5))))
                      }
                    )
                    )

## ================== Non-Member Functions for gPair ================== ##

#' @name length
#' @title length.gPair
#' @description
#' 
#' The number of junctions in this gPair Object
#'
#' @param gPair a gPair Object
#' @return the number of junctions in the gPair Object
#' @author Marcin Imielinski                         
#' @export
`length.gPair` = function(gPair)
{
  return(gPair$length)
}


#' @name c
#' @title c.gPair
#' @description
#'
#' Concatenates gPair Objects
#'
#' @param gPair object
#'
#' @return a new concatenated gPair Object
#' @author Marcin Imielinski                         
#' @author Marcin Imielinski
#' @export
`c.gPair` = function(...)
{                            
  pairs.list=list(...)

  isg = sapply(pairs.list, inherits, 'gPair')
  
  if(any(!isg)){
    stop('Error: All inputs must be of class gPair.')
  }

  ## Get all the pairs to create new gPair Object
  grll = lapply(pairs.list, function(x) x$grl)
  grlm = lapply(grll, function(x) as.data.table(values(x)))
  newgrl = dodo.call(grl.bind, lapply(grll, function(x) {values(x) = NULL; x}))
  values(newgrl) = rbindlist(grlm, fill = TRUE)
  return (gPair$new(newgrl))
}


#' @name union
#' @title union.gPair
#' @description
                                        #o' 
#' Returns a new gPair Object which is the union of x and y.
#' 
#' @param x a gPair Object
#' @param y a gPair Object
#' @author Marcin Imielinski
#' @return new gPair Object containing the union of x and y
#' @author Marcin Imielinski                         

#' @export
setMethod("union", c('gPair', "gPair"), function(x, y, pad = 0, ignore.strand = FALSE, ...)
{
  newJunc=c(x, y)
  return(unique(newJunc, pad, ignore.strand))
})

#' @name unique
#' @title unique.gPair
#' @description
#' 
#' Returns the subset of gPair object that it is unique
#' 
#' @param x a gPair Object
#' @author Marcin Imielinski
#' @return new gPair Object containing the union of x and y
#' @author Marcin Imielinski                         
#' @export
"unique.gPair" = function(x, pad = 0, ignore.strand = FALSE)
{
  if (pad==0)
    return(x$copy$removeDups())
  else
    return(x[!ra.duplicated(x$grl)])
}
setMethod("unique", c('gPair'), unique.gPair)


#' @name setdiff
#' @title setdiff.gPair
#' @description
#' 
#' Returns a new gPair Object which is the difference between x and y.
#'
#' @param x a gPair Object
#' @param y a gPair Object
#' @author Marcin Imielinski
#' @exportMethod setdiff
#' @return new gPair Object containing the difference between x and y
#' @author Marcin Imielinski                         
#' @export
setMethod("setdiff", c('gPair', "gPair"), function(x, y, pad = 0, ...)
{  
  ov = ra.overlaps(x$grl, y$grl, pad = pad)
  ix = setdiff(1:length(x$grl), ov[,1])
  return(x[ix])
})


#' @name intersect
#' @title intersect.gPair
#' @description
#' 
#' Returns a new gPair Object which is the intersection of x and y.
#' 
#' @param x a gPair Object
#' @param y a gPair Object
#' @author Marcin Imielinski  
#' @return new gPair Object containing the intersection of x and y
#' @author Marcin Imielinski                         
#' @export
setMethod("intersect", c('gPair', 'gPair'), function(x, y, pad = 0, ...) {

  ov = ra.overlaps(x$grl, y$grl, pad = pad)
  return(unique(x[ov[, 'ra1.ix']], pad = pad))
})


#' @name [.gPair
#' @title gPair
#' @description
#'
#' Overloads subset operator for gPair. Allows subsetting of gPair via index.
#' Also, allows for data.table style queries on gPair GRangesList and
#' corresponding metadata.
#'
#' @param obj gPair object This is the gPair object to be subset
#' @param i integer, logical, or expression in gPair metadata used to subset gPair
#' @return A new gPair object that contains only the given id's
#' @author Marcin Imielinski                         
#' 
#' @export
'[.gPair' = function(obj, i, j, with = TRUE){
  pairs = obj$clone()
  if (!missing(i))
  {
    if (with)
    {
      inew = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(i)))), parent.frame()),pairs$dt, parent.frame(2)), error = function(e) NULL)
      if (is.null(inew))
        inew = i ## just give up
    }
    else
    {
      inew = i
    }
    pairs$subset(inew)
  }

  if (!missing("j"))
  {
    grl = pairs$grl
    values(grl) = values(grl)[, j, drop = FALSE]
    pairs = jJ(grl)
  }
  return(pairs)
}


#' @name +.gPair
#' @description
#'
#' Allows padding of junctions. The rpad will be added to the left of a "+" junction
#' and to the right of "-" junction.
#'
#' @param jj gPair Object
#' @param pad Positive number representing amount to pad the gPair Object.
#' @return a new gPair class with the padding applied
#' @export
'+.gPair' = function(jj, pad)
{
  new.grl = resize(jj$grl, pad+1, fix="end")
  values(new.grl) = values(jj$grl)
  return(gPair$new(new.grl))
}

#' @name union
#' Returns a new gPair object which is the union of x and y.
#' 
#' @param x a gPair Object
#' @param y a gPair Object
#' @author Rick Mortenson
#' @exportMethod union
#' @return new gPair Object containing the union of x and y
setMethod("union", c("gPair", "gPair"),
          function(x, y) {
            newJunc=c(x, y)
            newJunc$removeDups()
            return(newJunc)
          })

#' @name setdiff
#' Returns a new gPair object which is the difference between x and y (id's).
#'
#' @param x a gPair Object
#' @param y a gPair Object
#' @author Rick Mortenson
#' @exportMethod setdiff
#' @return new gPair containing the difference between x and y
setMethod("setdiff", c("gPair", "gPair"),
          function(x, y) {
            ## Make sure that both come from the same graph                                         
            overlaps=ra.overlaps(x$grl, y$grl)
            overlaps=overlaps[, "ra1.ix"]
            all.ix=c(1:length(x$grl))
            dif.ix=setdiff(all.ix, overlaps)             
            return(gPair$new(x$grl[dif.ix]))
            
          })

#' @name refresh
#' @description
#' 
#' Updates gPair object to reflect changes in source code
#' 
#' @param gPair object
#' @return gPair object
#' @exportMethod refresh
#' @export
setMethod("refresh", "gPair",
          function(x) {
            return(gPair$new(x$pairs))
          })

