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
                      #' @param gr GRanges or gPair
                      #' @param dat length(gr) matrix or data.table / data.frame with field $i, $j, and $value indexing gr
                      #' @param fill fill value (default 0) for missing entries
                      #' @param field if gPair has a field then this will be used as a weight when adding / superimposing gPairs
                      #' @param full whether to explicitly set missing entries to fill (warning setting full = FALSE may distort certain calculations, though safe for summing etc)
                      #' @author Marcin Imielinski                         
                      initialize = function(gr = NULL, dat = NULL, field = NULL, fill = 0, full = NULL, na.rm = TRUE, agg.fun = sum, lower.tri = FALSE)
                      {
                        private$pdat = data.table()
                        private$pgr = NULL
                        private$pfill = fill

                        if (is.null(full))
                          private$pfull = is.na(fill) || fill != 0 ## full by default if fill not NA and != 0
                        else
                          private$pfull = full

                        private$pna.rm = na.rm

                        if (is(gr, 'gPair'))
                        {
                          gp = gr
                          if (is.null(field))
                            weight = rep(1, length(gp))
                          else
                            weight = gp$meta[[field]]

                          ## these will be merged via aggregation below
                          gr = grbind(gp$gr1, gp$gr2)
                          dat = data.table(i = 1:length(gp),
                                           j = length(gp) + 1:length(gp),
                                           value = weight)
                        }

                        if (!is.null(gr))
                          private$pgr = gr.stripstrand(gr)
                        
                        if (!is.null(dat) && nrow(dat)>0)
                        {
                          if (inherits(dat, 'data.frame'))
                          {
                            if (!all(c('i', 'j', 'value') %in% colnames(dat)))
                            {
                              stop('missing column $value in data.table dat, either provide matrix or data.table with fields $i, $j, and $value')
                            }

                            if (any(dat$i>length(gr))| any(dat$j>length(gr)))
                            {
                              stop("out of bounds indices in dat data.table / data.frame, indices $i and $j are integers that index gr, and thus can't exceed the length")
                            }
                          } else
                          {
                            if (!is(dat, 'Matrix') & !is.matrix(dat) & !is.array(dat))
                              stop('dat must be either Matrix, matrix, data.table, or data.frame')

                            if (!identical(length(gr), unique(dim(dat))))
                              stop('if dat is matrix it must be square with dimensions that are equal to length(gr)')

                            if (is(dat, 'Matrix'))
                              dat = as.data.table(which(dat != private$pfill, arr.ind = TRUE))[, .(i = row, j = col)][, value := dat[cbind(i, j)]]
                            else
                              dat = melt(dat)

                            setnames(dat, c("i", "j", "value"))
                            nr = nrow(dat)
                            if (!lower.tri)
                              dat = dat[i<=j, ] ## only take upper triangle
                            if (nrow(dat)<nr)
                              warning('Lower triangular values from input matrix were ignored: when using matrix input please put all data values in diagonal and upper triangle or use lower.tri = TRUE')                            
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

                            if ('index' %in% names(values(gr)))
                            {
                              warning('Metadata column $index is reserved for gMatrix operations, removing / ignoring from provided GRanges')
                              gr$index = NULL
                            }

                            gr = gr.stripstrand(private$pgr)
                            dgr = gr.disjoin(gr)
                            dgr = dgr[order(gr.match(dgr, gr))] ## resort in order of gr
                            is.diff = !identical(gr[, c()], dgr[, c()])
                            if (is.diff) ## then we need to aggregate or collapse
                            {
                              warning('input intervals had to be disjoined, aggregating via agg.fun')
                              private$pgr = dgr
                              
                              map = gr2dt(gr %*% private$pgr[, c()])[, .(old = query.id, new = subject.id)]
                              dat = merge(merge(dat, map, by.x = 'i', by.y = 'old', allow.cartesian = TRUE),
                                          map, by.x = 'j', by.y = 'old', allow.cartesian = TRUE)[, .(i = new.x, j = new.y, value = value)]
                              
                            }

                            if (!is.null(dat))
                            {
                              is.dup = any(duplicated(dat, by = c("i", "j")))
                              if (is.dup)
                              {
                                warning('input data had some duplicates, aggregating via agg.fun')

                                dat = dat[, .(i = pmin(i,j), j = pmax(i,j), value = value)][, .(value = private$pagg.fun(value)), keyby = .(i, j)]
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
                            else if (!is.null(dat) && nrow(dat)>0)
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
                        dgr = sort(gr.fix(gr.fix(disjoin(grbind(private$pgr, gr)), private$pgr), gr))
                        private$pgr$query.id = NULL
                        grdt = gr2dt(gr.merge(dgr, gr.fix(private$pgr, dgr), all = TRUE)) 
                        setkey(grdt, query.id)
                        grdt$index = 1:nrow(grdt)

                        if (nrow(private$pdat)>0)
                        {
                            dti = merge(private$pdat, grdt[, .(i = subject.id, inew = index)],
                                        by = 'i',
                                        allow.cartesian = TRUE)
                            setkey(dti, j)
                            dtj = grdt[, .(j = subject.id, jnew = index)]
                            setkey(dtj, j)
                            dt.new = merge(
                              dti,
                              dtj, by = 'j',
                              allow.cartesian = TRUE)[, .(
                                            i = pmin(inew, jnew),
                                            j = pmax(inew, jnew),
                                            value = value)][, id := 1:.N]
                          }
                        else
                          dt.new = data.table()

                        dt.new = unique(dt.new, by = c('i', 'j'))
                        gr.new = dt2gr(grdt, seqlengths = seqlengths(gr))
                        gr.new$index = NULL ## strip index

                        return(gMatrix$new(gr.new, dt.new, fill = self$fill, full = self$full,
                                          agg.fun = self$agg.fun, na.rm = self$na.rm))
                        ## private$pgr = tmp$gr[, c()]
                        ## private$pdat = tmp$dat
                        ## return(invisible(self))
                      },

                      #' @name sweep
                      #' @description
                      #' sweeps the vector with the provided operation from all values
                      #' (no margin to specify since square matrix, so applied to both margins)
                      #' TODO: maybe need to do in ? both directions
                      #' NOTE: modifies gMatrix in place
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

                      #' @name metaset
                      #' @description
                      #'
                      #' sets metadata of gMatrix GRanges, either adding columns or changing
                      #' existing values
                      #' @param ... name value pairs where value is a vector / scalar that broadcasts to $gr
                      #' @author Marcin Imielinski
                       metaset = function(...)
                       {
                         args = list(...)
                         
                         for (arg in names(args)){
                           values(private$pgr)[[arg]] = args[[arg]]
                         }
                         return(invisible(self))
                       },

                      #' @name metaset
                      #' @description
                      #'
                      #' sets metadata of gMatrix GRanges, either adding columns or changing
                      #' existing values
                      #' @param ... name value pairs where value is a vector / scalar that broadcasts to $gr
                      #' @author Marcin Imielinski
                      annotate = function(y)
                      {
                        gr.new = gr.val(self$gr, y, val = names(values(y)))
                        do.call(self$metaset, as.list(values(gr.new)))
                        return(invisible(self))
                      },

                      #' @name clusters
                      #' @description
                      #'
                      #' finds clusters or communities of bins in data based on
                      #' by building a graph that uses all nonzero pixels as edge weights
                      #' on a graph
                      #' 
                      #' @param mode clustering mode either 'weak', 'louvain', 'rwlaktrap')finding weakly connected components 
                      #' @author Marcin Imielinski
                      clusters = function(mode = 'louvain', gpair = TRUE)
                      {
                        G = graph.adjacency(self$mat!=0)

                        if (mode %in% c("weak")){
                          membership = igraph::clusters(G, mode)$membership
                        }
                        else if (mode== 'walktrap'){
                          membership = igraph::cluster_walktrap(G)$membership
                        }
                        else if (mode== 'louvain'){
                          membership = igraph::cluster_walktrap(G)$membership
                        }
                        else
                        {
                          stop(paste('clustering mode', mode, 'not supported'))
                        }

                        ## resort so most popular clusters are named first
                        rename = data.table(names(rev(sort(table(c(membership))))))[, structure(1:.N, names = V1)]
                        membership = rename[as.character(membership)]

                        if (gpair)
                        {
                          tmp.gr = private$pgr
                          tmp.gr$cluster = membership
                          tmp.dat = private$pdat[value!=0,]
                          tmp.dat$cl = membership[tmp.dat$i]
                          tmp.dat$cl.2 = membership[tmp.dat$j]
                          setkey(tmp.dat, cl) ## sort by cluster                          
                          gr1 = tmp.gr[tmp.dat$i]
                          gr2 = tmp.gr[tmp.dat$j]
                          gp = gP(gr1, gr2, meta = data.table(cluster = tmp.dat$cl, cluster2 = tmp.dat$cl.2, value = tmp.dat$value))
                          return(gp)
                        }

                        return(membership)
                      },

                      #' @name transform
                      #' @description
                      #' transform gMatrix values by some function
                      #' @param FUN
                      #' @param ... other args to FUN can be scalar or gMatrix
                      #' @author Marcin Imielinski                         
                      transform = function(FUN, ...)
                      {
                        dt.new = self$dat                        
                        args = list(...)

                        grd = self$gr
                        if (any(ix <- sapply(lapply(args, class), '[', 1) == 'gMatrix'))
                          {
                            grd = sort(disjoin(do.call(grbind, c(list(grd), lapply(args[ix], function(x) x$gr)))))
                          }

                        selfd = self$disjoin(grd)+0
                        dt.new = selfd$dat
                        gr.new = selfd$gr
                        funcall = 'FUN(value)'

                        if (length(args)>0)
                        {
                          argnm = paste0(names(args), '')
                          unnamed = nchar(argnm)==0
                          vname = paste0('value', 1:length(args))
                          argval = lapply(1:length(args), function(x) NA)

                          for (i in 1:length(args))
                          {
                            if (is(args[[i]], 'gMatrix'))
                            {
                              tmp = merge(selfd, args[[i]]$disjoin(grd), gmatrix = TRUE)[[2]]
                              if (!identical(selfd$gr[, c()], tmp$gr[,c()]))
                                stop('error merging self with provided gMatrix arguments')
                              dt.new[[vname[i]]] = tmp$dat[.(dt.new$i, dt.new$j), ]$value
                            }
                            else
                            {                             
                              argval[[i]] = args[[i]] ## replace argval with what was provided
                              vname[i] = sprintf('argval[[%s]]', i) ## a bit round about hack to get the FUN call to access argval for non gMatrix
                            }
                          }

                          ## build function call using named an unnamed arguments
                          funcall = paste0('FUN(value,', paste(ifelse(unnamed, vname, paste(argnm, vname, sep = '=')), collapse  = ','), ')')
                        }

                        dt.new = eval(parse(text = paste('dt.new[, value := ', funcall, ']')))

                        return(gMatrix$new(gr.new, dt.new, fill = self$fill, full = self$full,
                                          agg.fun = self$agg.fun))
                      },

                      #' @name set
                      #' @description
                      #' sets matrix "patch" defined by gr pair to given (scalar) value
                      #' and if either gr1 or gr2 not specified then fills in the entire
                      #' "column", disjoining gMatrix if need be
                      #' @param val scalar value
                      #' @param gr2 GRanges to set
                      #' @param gr1 GRanges to set
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
                      #' provided intervals, aggregating via agg.fun.
                      #' if weighted = TRUE, computes an area weighted area of each function
                      #' i.e. sum_ij (w_i*w_j*FUN(value))
                      #' average), returns new gMatrix
                      #' @param gr GRanges
                      #' @param weighted logical determining whether to compute a width weighted average
                      #' @author Marcin Imielinski                         
                      agg = function(gr = NULL, FUN = private$pagg.fun, weighted = FALSE)
                      {
                        if (is.null(gr))
                          return(invisible(self))
                        gr.new = gr

                        pgr = private$pgr
                        pdat = private$pdat

                        if (weighted)
                        {
                          tmpself = self[gr, gr]+0
                          pgr = tmpself$gr
                          pdat = tmpself$dat
                        }

                        if (!inherits(gr, 'GRanges'))
                          stop('input to $agg should be a GRanges')
                        ov = gr2dt(gr.new %*% pgr)
                        wid = as.numeric(width(pgr))
                        if (nrow(private$pdat)>0 && nrow(ov)>0)
                          {
                            dt.tmp = merge(
                              merge(pdat,
                                    ov[, .(i = subject.id, inew = query.id)], by = 'i', allow.cartesian = TRUE),
                              ov[, .(j = subject.id, jnew = query.id)], by = 'j', allow.cartesian = TRUE)[inew<=jnew, ]
                            if (weighted)
                            {
                              dt.tmp[, area := wid[i]*wid[j]]
                              dt.new = dt.tmp[, .(value = FUN(value*area, na.rm = private$pna.rm) /
                                                    FUN(area, na.rm = private$pna.rm)),
                                              by = .(i = inew, j = jnew)]
                            }
                            else
                              dt.new = dt.tmp[, .(value = FUN(value, na.rm = private$pna.rm)),
                                              by = .(i = inew, j = jnew)]
                          }
                        else
                          dt.new = pdat

                        gm = gMatrix$new(gr.new, dt.new, fill = self$fill, full = self$full,
                                         agg.fun = self$agg.fun)
                        return(gm)

                      },

                      #' @name merge
                      #' @description
                      #' merge with (list of) other gmats
                      #' and applies vector input, scalar valued function (default sum)
                      #' to values of interval pairs that overlap
                      #' returns new gMatrix 
                      #' @param gmats gMatix or list of gMatrices
                      #' @author Marcin Imielinski                         
                      merge = function(gmats = NULL, FUN = self$agg.fun, simplify = FALSE, mc.cores = 1, verbose = FALSE)
                      {
                        if (is.null(gmats))
                          return(self)
                        
                        if (!is.list(gmats))
                          gmats = c(list(self), list(gmats))

                        if (!all(sapply(gmats, function(x) class(x)[1])=='gMatrix'))
                          stop('gmats must be gMatrix or list of gMatrices')

                        if (simplify)
                          gmats = mclapply(1:length(gmats), function(x)
                          {
                            if (verbose)
                              message('Simplifying ', x)
                            refresh(gmats[[x]])$simplify
                          }, mc.cores = mc.cores)
                       
                        dts = rbindlist(mclapply(1:length(gmats), function(x)
                        {
                          dat = gmats[[x]]$dat
                          if (nrow(dat)>0)
                            dat[, source := x]
                          dat
                        }, mc.cores = mc.cores), fill = TRUE)

                        ## concatenate grs from input gMatrices
                        grs = do.call(grbind, lapply(1:length(gmats), function(i) {gr = gmats[[i]]$gr; gr$source = i; gr$old = 1:length(gr); gr}))
                        ## disjoin and match 
                        gr.new = disjoin(grs)
                        dt.final = NULL
                        if (nrow(dts))
                        {
                          grs = grs %*% gr.new[, c()]
                          ## map old x source --> new
                          map.dt = data.table(index = grs$subject.id, source = grs$source, old = grs$old)
                          
                          
                          ## also want to keep track of all "missing" ie source x index combos that don't have
                          ## a corresponding value in old - i.e. left join of cartesian product of grs and source
                          ## with map.dt
                          other.dt = merge(as.data.table(
                            expand.grid(index = 1:length(gr.new), source = 1:length(gmats))),
                            map.dt, by = c('source', 'index'), allow.cartesian = TRUE, all.x = TRUE)[is.na(old), ]
                          map.dt = rbind(map.dt, other.dt)
                          setkeyv(map.dt, c('source', 'old'))
                          
                          ## count.grs 
                          ngrs = sapply(gmats, length)
                          fills = sapply(gmats, function(x) x$fill)
                          
                        ## we keep it "simple" if all fills are 0 and function is sum
                        keep.it.simple = all(fills == 0) & identical(FUN, sum)

                        ## first merge
                        dt1 = merge(dts[, .(i,j, value, source)], map.dt,
                                    by.x = c('i', 'source'), by.y = c('old', 'source'),
                                    allow.cartesian = TRUE, all.y = !keep.it.simple)

                        ## any NA in j represent right join rows that are in map.dt but missing in dts
                        ## we fill these in with all possible j values
                        if (!keep.it.simple)
                        {
                          dt1 = rbind(dt1[!is.na(j), ],
                                      dt1[is.na(j), .(j = 1:ngrs[source], fill = fills[source]),
                                          by = .(source, i, value, index)], fill = TRUE)
                          if (any(!is.na(dt1$fill)))
                            dt1[!is.na(fill), value := fill]
                        }

                        if (identical(FUN, sum)) ## can get rid of 0 val rows
                          dt1 = dt1[value!=0, ]

                          dt2 = merge(
                            dt1,
                            map.dt, all.y = !keep.it.simple, 
                            by.x = c('j', 'source'), by.y = c('old', 'source'),
                            allow.cartesian = TRUE)

                          ## we do similar things for missing i in dt2
                          if (!keep.it.simple)
                          {
                            if (any(is.na(dt2$i)))
                              dt2 = rbind(dt2[!is.na(i), ],
                                          dt2[is.na(i), .(i = 1:ngrs[source], fill = fills[source]),
                                              by = .(source, j, value, index.x, index.y)], fill = TRUE)
                            
                            ## apply fill
                            if (any(!is.na(dt2$fill)))
                              dt2[!is.na(fill), value := fill]                                                       
                          }

                          ## if i == j merge will create duplicate pixels which can cause over-counting
                          ## remove these
                          dt2 = dt2[index.x <= index.y, ]                        
                          dt2[, i := index.x]
                          dt2[, j := index.y]

                          dt.final = dt2[, .(value = FUN(value, na.rm = self$na.rm)), by = .(i, j)]
                          dt.final[, id := 1:.N]
                        }
                        ret = gMatrix$new(gr.new, dt.final, fill = self$fill,
                                           full = self$full, agg.fun = self$agg.fun,
                                           na.rm = self$na.rm
                                          )
                        return(ret)
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
                        gr.new$index = NULL

                        ## we want dt.new that either have ikeep jkeep
                        ## or jkeep ikeep combo (but not jkeep jkeep
                        ## combo or ikeep ikeep combo)
                        dt.new = dt.new[(i %in% ikeep & j %in% jkeep) |
                               (i %in% jkeep & j %in% ikeep), ]
                        
                        return(gMatrix$new(gr.new, dt.new, full = self$full, na.rm = self$na.rm, agg.fun = self$agg.fun, fill = self$fill))
                      },
                      
                      #' @name gtrack
                      #' @description
                      #' Outputs a gTrack of this matrix
                      #' @param colormap colors to use with this colormap, default = c('white', 'red', 'black')
                      #' @param clim length 2 vector color limits to which to assign the min and max color
                      #' @author Marcin Imielinski                         
                      gtrack = function(name = '', colormap = c('white', 'red', 'black'), 
                                        clim = NA, cmap.min = NA, cmap.max = NA, quantile = 0.01, max.ranges = 3e3,
                                        ...)
                      {
                        dat = private$pdat
                        if (is.na(cmap.min))
                          cmap.min = clim[1]

                        if (is.na(cmap.max))
                          cmap.max = clim[2]
                        
                        if (is.na(cmap.min))
                          cmap.min = pmin(self$fill, quantile(dat$value, pmin(quantile, 1-quantile), na.rm = TRUE))

                        if (is.na(cmap.max))
                          cmap.max = pmin(max(setdiff(dat$value, Inf), na.rm = TRUE), pmax(self$fill, quantile(dat[value>cmap.min, value], pmax(1-quantile, quantile), na.rm = TRUE)))

                        mdata = Matrix::sparseMatrix(dat$i, dat$j, x = dat$value,
                                                     dims = c(length(self$gr), length(self$gr)),
                                                     symmetric = TRUE)
                        return(gTrack(private$pgr, name = name, mdata = mdata, colormap = colormap,
                                      max.ranges = max.ranges,
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

                      #' @name hicpro
                      #' @description
                      #'
                      #' Writes gMatrix to disk in 'hicpro' format
                      #' ie prefix_ord.bed file and i j value "sparse" text file prefix.matrix
                      #' where prefix is file path input with .matrix stripped
                      #' @param prefix
                      #' @author Marcin Imielinski                         
                      hicpro = function(file, verbose = TRUE)
                      {
                        prefix = gsub('\\.matrix$', '', file)
                        out.matrix = paste0(prefix, '.matrix')
                        out.bed = paste0(prefix, '_ord.bed')

                        if (!file.exists(dirname(file)))
                        {
                          if (verbose)
                            ggmessage('Making dir' , dirname(file))
                          system(paste('mkdir -p', dirname(file)))
                        }                       
                        
                        if (verbose)
                          ggmessage('Dumping hicpro style .matrix file to ' , out.matrix)
                        fwrite(self$dat[, .(i, j, value)], out.matrix, sep = '\t', col.names = FALSE)
                        if (verbose)
                          ggmessage('Dumping hicpro style .bed file to ', out.bed)
                        rtracklayer::export(self$gr, out.bed)
                        
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

                      #' @name print
                      #' @description
                      #' 
                      #' Prints out gMatrix Object. Prints the length and the GRangesList of the junctions.
                      #' @author Marcin Imielinski                         
                      print = function()
                      {
                        ngr = length(private$pgr)
                        ndt = length(private$pdat)                           
                        ggmessage(sprintf('gMatrix with %s intervals and %s entries',
                                        ngr, ndt), appendLF = FALSE)
                        
                        if (length(private$pgr) > 0){
                          ggmessage(' containing:\n')
                          print(self$gr)

                          if (!is.null(private$pdat))
                            print(self$dat)
                        }
                        else
                          ggmessage('\n')
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
                      fill = function(value) {
                        if (missing(value))
                        {
                          private$pfill
                        }
                        else
                        {
                          if (length(value) == 1 && (is.numeric(value) | is.integer(value)))
                          {                            
                            private$pfill = value
                          }
                          else
                          {
                            stop('value must be a scalar numeric')                            
                          }
                          return(invisible(self))
                        }
                      },
                        
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
                      full = function(value)
                        {
                          if (missing(value))
                          {
                            private$pfull
                          }
                          else
                          {
                            if (length(value) == 1 && (is.logical(value)))
                            {
                              private$pfull = value
                              if (private$pfull) ## fill by adding 0
                                private$pdat = (self + 0)$dat
                              else ## remove all non pfill rows from pdat
                                private$pdat = private$pdat[value != private$pfill, ]
                            }
                            else
                            {
                              stop('value must be a scalar logical')                            
                            }
                            return(invisible(self))
                          }
                        },

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

                      #' @name mat
                      #' @description
                      #'
                      #' Returns matrix form of dat, i.e. of granges pairs -> values in sparseMatrix format (self$fill == 0)
                      #' full matrix if (self$fill != 0) or
                                        #  or if (if baselevel == TRUE) of bp pairs (warning this matrix can be huge)
                      #'
                      #' @param reverse will look at valleys instead of peaks
                      #' @author Marcin Imielinski
                      mat = function()
                      {
                        upperdiag = TRUE
                        if (length(private$pgr)==0)
                          return(sparseMatrix())

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
                        
                        return(out)
                      },

                      #' @name simplify
                      #' @description
                      #' simplifies gMatrix by merging all reference adjacent intervals
                      #' with identical value profiles 
                      #' @author Marcin Imielinski                         
                      simplify = function()
                      {
                        mat = self$mat
                        ## symmetrify
                        mat[which(lower.tri(mat))] = t(mat)[which(lower.tri(mat))]
                        grdt = gr2dt(self$gr)[, grsig := as.integer(factor(apply(mat, 1, paste, collapse = ' ')))]
                        grdt[, same := c(diff(grsig)==0, NA), by = seqnames]
                        grdt[, run := label.runs(same), by = seqnames]
                        grdt[, run := ifelse(is.na(run) & !is.na(shift(run)), shift(run), run)]
                        grdt[, nruns := max(c(0L,run), na.rm = TRUE), by = seqnames]
                        grdt[is.na(run), run := 1:.N + nruns, by = seqnames]
                        grdt[, row := 1:.N]
                        newgr = dt2gr(grdt[  , .(start = min(start), end = max(end), row = row[1]), by = .(seqnames, run)], seqlengths = seqlengths(self$gr))                        
                        newmat = self$mat[newgr$row, newgr$row, drop = FALSE]
                        newgm = gMatrix$new(newgr[, names(values(self$gr))], dat = newmat)
                        private$pgr = newgm$gr
                        private$pdat = newgm$dat
                        return(invisible(self))
                      },
                      
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
#' @param gr GRanges around which to build a gMatrix or gPair object to sum
#' @param dat data.table of $i and $j indexing gr and field $value OR a vector valued expression involving terms i and j that will be applied to expand.grid
#' @param outer character scalar field of gr that will be used to generate an "outer product" of gr with itself, which given vector x for numerics will make each gmatrix entry equal to x_i*x_j and for character vector or factor fields will make each vector equal to x_i == x_j
#' @param field if gr is a gPair then field can be a character that specifies a gPair metadata numeric or integer metadata field that will be used as weights to sum the gPair objects
#' @param full logical flag whether to explicitly "fill" missing entries in provided data.table (default = FALSE)
#' @param fill value with which to fill missing values (default 0)
#' @param agg.fun default function with which to aggregate values, this function should take a vector and na.rm = TRUE and return a scalar
#' @param na.rm default na.rm argument to agg.fun
#' @return A new gMatrix object
#' @author Marcin Imielinski                         
#' @export
gM = function(gr = NULL, dat = NULL, outer = NULL, field = NULL, full = FALSE, fill = 0, agg.fun = sum, na.rm = TRUE, lower.tri = FALSE)
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
        ijdt = as.data.table(ij)[, .(i = row, j = col)]
        if (!lower.tri)
          ijdt = ijdt[i<=j, ] ## only take upper triangle
#        if (nrow(ijdt)<nrow(ij))
 #         warning('Lower triangular values from input matrix were ignored, to include (and aggregate) lower triangle values please use lower.tri = TRUE')
        dat = cbind(ijdt, data.table(value = dat[ijdt[, cbind(i, j)]]))
      }
    }

    if (is.null(dat))
    {
      stop('dat argument is malformed, should either be data.table with fields $i, $j, and $value OR expression involving $i and $j that will be evaluated on the expanded grid of ij pairs')
    }
  } else if (!is.null(outer)) ## generate dat corresponding outer product
  {
    if (!(outer %in% names(values(gr))))
      stop('outer should be field of gr from which outer product will be computed')

    dat = as.data.table(expand.grid(i = 1:length(gr), j = 1:length(gr)))[i <= j, ]

    val = values(gr)[[outer]]
    if (is.factor(val) | is.character(val))
    {
      dat[, value := sign(val[i] == val[j])]
    }
    else
    {
      dat[, value := val[i]*val[j]]
    }
    dat = dat[value!=fill, ]
  }

  return(gMatrix$new(gr, dat, field = field, full = full, fill = fill, agg.fun = agg.fun, na.rm = na.rm, lower.tri = lower.tri))
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
'[.gMatrix' = function(obj, i = NULL, j = NULL, with = TRUE, drop = FALSE, full = FALSE){

  if (with)
  {
    dt = gr2dt(obj$gr)
    inew = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(i)))), parent.frame()), dt, parent.frame(2)), error = function(e) NULL)
    if (is.null(inew))
    {
      inew = i ## just give up
    }
    i = inew

    jnew = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(j)))), parent.frame()), dt, parent.frame(2)), error = function(e) NULL)
    if (is.null(jnew))
    {
      jnew = j ## just give up
    }
    j = jnew
  }

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

      if (drop == TRUE)  ## drop output to GRanges, only works if single i bin is specified
      {
        if (is.vector(i) && length(unique(i))==1)
        {  
          ui = i[1]
          ret = obj$gr
          ret$value = obj$fill
          tmp = rbind(obj$dat[i == ui, .(ind = j, value)],
                      obj$dat[j == ui, .(ind = i, value)])
          ret$value[tmp$ind] = tmp$value
          return(ret)
        }
        else
          warning('drop = TRUE only drops to GRanges when a single index or single GRanges bin is provided as an argument')
      }
        
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
grnorm = function(gr, mean = NULL, sd = NULL, full = FALSE, fill = 0, agg.fun = sum, na.rm = TRUE)
{
  if (is(gr, 'gMatrix'))
  {
    if (is.null(mean))
      mean = gr + 0 ## adding 0 will make it "full"
    gr = gr$gr 
  }

  if (is.null(mean))
  {
    mean = 0
  }

  if (is.null(sd))
  {
    sd = 1
  }
  
  if (!is(mean, 'gMatrix'))
  {
    mean = gM(gr, mean, fill = fill, full = full, agg.fun = agg.fun, na.rm = na.rm)
  }

  if (!is(sd, 'gMatrix'))
  {
    sd = gM(gr, sd, fill = fill, full = full, agg.fun = agg.fun, na.rm = na.rm)
  }

  if (!(identical(mean$gr[, c()], gr)))
    mean = mean$disjoin(gr)

  if (!(identical(sd$gr[, c()], gr)))
    sd = sd$disjoin(mean$gr)

  ## now mean and sd should be both aligned
  new.dat = mean$dat[, value := rnorm(.N, value, sd = sd$dat$value)]
  return(gM(gr, dat = new.dat, fill = fill, full = full, agg.fun = agg.fun, na.rm = na.rm))
}

#' @name grpois
#' @title grpois
#' @description
#'
#' gMatrix instantiator creates random matrix with poisson
#' 
#' @param gr GRanges around which to build a gMatrix or gMatrix from which to draw "lambda" parameter
#' @param lambda scalar lambda parameter of poisson or gMatrix from which to draw lambda parameter)
#' @param full logical flag whether to explicitly "fill" missing entries in provided data.table (default = FALSE)
#' @param fill value with which to fill missing values (default 0)
#' @param agg.fun default function with which to aggregate values, this function should take a vector and na.rm = TRUE and return a scalar
#' @param na.rm default na.rm argument to agg.fun
#' @export
grpois = function(gr, lambda = NULL, full = FALSE, fill = 0, agg.fun = sum, na.rm = TRUE)
{
  if (is(gr, 'gMatrix'))
  {
    if (is.null(lambda))
      lambda = gr
    gr = gr$gr 
  }

  if (is.null(lambda))
  {
    lambda = 1
  }
  else if (is(lambda, 'gMatrix'))
  {
    init = lambda
  }
  else  
    init = gM(gr, lambda, fill = fill, full = full, agg.fun = agg.fun, na.rm = na.rm)

  if (!(identical(init$gr[, c()], gr[, c()])))
    init = init$disjoin(gr)

  new.dat = init$dat[, value := rpois(n = length(value), lambda = value)]
  return(gM(init$gr, dat = new.dat, fill = fill, full = full, agg.fun = agg.fun, na.rm = na.rm))
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

  if (nrow(tmp)>0)
    tmp = tmp[!is.na(group), ]

  if (!nrow(tmp))
    {
      if (length(bins)>0)
        return(gM(bins))
      else
        return(gM())
    }
 
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
  if (identical(mat1$gr[, c()], mat2$gr[, c()]))
    gr.new = mat1$gr
  else
    {
      gr.new = sort(disjoin(grbind(mat1$gr, mat2$gr)))
      
      ## these will have identical indices
      mat1  = mat1$copy$disjoin(gr.new)
      mat2  = mat2$copy$disjoin(gr.new)
      gr.new = mat1$gr ## contents may have shifted during disjoin
    }

  dat1 = data.table(i = as.integer(NA), j = as.integer(NA),
                    id1 = as.integer(NA), val1 = as.numeric(NA))
  dat2 = data.table(i = as.integer(NA), j = as.integer(NA),
                    id1 = as.integer(NA), val1 = as.numeric(NA))

  if (nrow(mat1$dat)>0)
    dat1 = mat1$dat[, .(i, j, id1 = id, val1 = value)]

  if (nrow(mat2$dat)>0)
    dat2 = mat2$dat[, .(i, j, id2 = id, val2 = value)]

  dat.new = merge(dat1, dat2, by = c('i', 'j'), all = TRUE)[!is.na(i), ]

  if (any(ix <- dat.new[, is.na(id1)]))
    dat.new[ix, val1 := mat1$fill]

  if (any(ix <- dat.new[, is.na(id2)]))
    dat.new[ix, val2 := mat2$fill]

  return(list(gr = gr.new, dat = dat.new))
}


#' @name merge.gMatrix
#' @description
#'
#' merges two gMatrix objects and returns a data.table of intervals and
#' and values for each of the inputs or (if gmatrix = TRUE) a length 2 list of
#' gMatrix objects.
#'
#' Note: only disjoins bins, does not aggregate any bins
#' 
#' @param mat1 gMatrix object
#' @param mat2 gMatrix object
#' @param gmatrix logical flag (default FALSE) whether to return a data.table (if FALSE) or a length 2 list of gMatrix objects (if TRUE)
#' @param meta logical flag (default FALSE) whether to return the metadata of the ranges as data.table columns
#' @return a data.table or length 2 list of gMatrix objects
#' @author Marcin Imielinski                         
#' @export
'merge.gMatrix' = function(x, y, gmatrix = FALSE, meta = TRUE)
{
  if (!is(x, 'gMatrix'))
    x = gM(y$gr, x)

  if (!is(y, 'gMatrix'))
    y = gM(x$gr, y)
  
  tmp = gmatalign(x, y)
  gr1 = gr2 = tmp$gr[, c()]
  dat.new = tmp$dat

  if (meta && ncol(values(x$gr)))
    values(gr1) = values(x$gr)[gr.match(gr1, x$gr), , drop = FALSE]
  
  if (meta && ncol(values(y$gr)))
    values(gr2) = values(y$gr)[gr.match(gr2, y$gr), , drop = FALSE]

  miss1 = !(tmp$gr %^% x$gr)
  miss2 = !(tmp$gr %^% y$gr)

  dat.new[, missing1 := miss1[i]]
  dat.new[, missing2 := miss2[j]]

  if (gmatrix)
  {
    return(list(
      gMatrix$new(gr1, dat.new[, .(i, j, value = val1)], fill = x$fill, full = x$full, agg.fun = x$agg, na.rm = x$na.rm),
      gMatrix$new(gr2, dat.new[, .(i, j, value = val2)], fill = y$fill, full = y$full, agg.fun = y$agg, na.rm = y$na.rm)
    ))
  }
  else
  {
    grdt1 = data.table(gr1 = gr.string(gr1))
    if (ncol(values(gr1)))
      grdt1 = cbind(grdt1, as.data.table(values(gr1)))

    grdt2 = data.table(gr2 = gr.string(gr2))
    if (ncol(values(gr2)))
      grdt2 = cbind(grdt2, as.data.table(values(gr2)))

    return(cbind(grdt1[dat.new$i, ], grdt2[dat.new$j, ], dat.new[, .(i, j, missing1, missing2, val1, val2)]))
  }
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
  if (!is(mat1, 'gMatrix'))
    mat1 = gM(mat2$gr, mat1)

  if (!is(mat2, 'gMatrix'))
    mat2 = gM(mat1$gr, mat2)
  
  tmp = gmatalign(mat1, mat2)
  gr.new = tmp$gr
  dat.new = tmp$dat

  dat.new[, value := val1 + val2]

  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg, na.rm = mat1$na.rm))
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
  if (missing(mat2))
    return(mat1*-1)

  if (!is(mat1, 'gMatrix'))
    mat1 = gM(mat2$gr, mat1)

  if (!is(mat2, 'gMatrix'))
    mat2 = gM(mat1$gr, mat2)
  
  tmp = gmatalign(mat1, mat2)
  gr.new = tmp$gr
  dat.new = tmp$dat

  dat.new[, value := val1 - val2]

  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg, na.rm = mat1$na.rm))
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
  if (!is(mat1, 'gMatrix'))
    mat1 = gM(mat2$gr, mat1)

  if (!is(mat2, 'gMatrix'))
    mat2 = gM(mat1$gr, mat2)
  
  tmp = gmatalign(mat1, mat2)
  gr.new = tmp$gr
  dat.new = tmp$dat

  dat.new[, value := val1 * val2]

  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg, na.rm = mat1$na.rm))
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
  if (!is(mat1, 'gMatrix'))
    mat1 = gM(mat2$gr, mat1)

  if (!is(mat2, 'gMatrix'))
    mat2 = gM(mat1$gr, mat2)
  
  tmp = gmatalign(mat1, mat2+0)
  gr.new = tmp$gr
  dat.new = tmp$dat

  dat.new[, value := val1 / val2]

  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg, na.rm = mat1$na.rm))
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
  if (!is(mat1, 'gMatrix'))
    mat1 = gM(mat2$gr, mat1)

  if (!is(mat2, 'gMatrix'))
    mat2 = gM(mat1$gr, mat2)
  
  tmp = gmatalign(mat1, mat2)
  gr.new = tmp$gr
  dat.new = tmp$dat

  dat.new[, value := val1 == val2]

  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg, na.rm = mat1$na.rm))
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
  if (!is(mat1, 'gMatrix'))
    mat1 = gM(mat2$gr, mat1)

  if (!is(mat2, 'gMatrix'))
    mat2 = gM(mat1$gr, mat2)
  
  tmp = gmatalign(mat1, mat2)
  gr.new = tmp$gr
  dat.new = tmp$dat

  dat.new[, value := val1 <= val2]

  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg, na.rm = mat1$na.rm))
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
  if (!is(mat1, 'gMatrix'))
    mat1 = gM(mat2$gr, mat1)

  if (!is(mat2, 'gMatrix'))
    mat2 = gM(mat1$gr, mat2)
  
  tmp = gmatalign(mat1, mat2)
  gr.new = tmp$gr
  dat.new = tmp$dat

  dat.new[, value := val1 >= val2]

  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg, na.rm = mat1$na.rm))
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
  if (!is(mat1, 'gMatrix'))
    mat1 = gM(mat2$gr, mat1)

  if (!is(mat2, 'gMatrix'))
    mat2 = gM(mat1$gr, mat2)
  
  tmp = gmatalign(mat1, mat2)
  gr.new = tmp$gr
  dat.new = tmp$dat

  dat.new[, value := val1 > val2]

  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg, na.rm = mat1$na.rm))
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
  if (!is(mat1, 'gMatrix'))
    mat1 = gM(mat2$gr, mat1)

  if (!is(mat2, 'gMatrix'))
    mat2 = gM(mat1$gr, mat2)
  
  tmp = gmatalign(mat1, mat2)
  gr.new = tmp$gr
  dat.new = tmp$dat

  dat.new[, value := val1 < val2]

  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg, na.rm = mat1$na.rm))
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
  if (!is(mat1, 'gMatrix'))
    mat1 = gM(mat2$gr, mat1)

  if (!is(mat2, 'gMatrix'))
    mat2 = gM(mat1$gr, mat2)
  
  tmp = gmatalign(mat1, mat2)
  gr.new = tmp$gr
  dat.new = tmp$dat

  dat.new[, value := val1 != val2]

  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg, na.rm = mat1$na.rm))
}


#' @name %$%
#' @title annotate metadata GRanges $gr of a gMatrix with those of another
#' @description
#'
#' gm %$% gr
#'
#' will populate the GRanges metadata of gMatrix gm (i.e. gm$gr)
#' with the metadata of gr
#'
#' @return gMatrix with GRanges metadata  populated
#' @rdname gr.val
#' @docType methods
#' @param x \code{gMatrix} object
#' @aliases %$%,gMatrix-method
#' @exportMethod %$%
#' @author Marcin Imielinski
setMethod("%$%", signature(x = "gMatrix"), function(x, y) {
  gr.new = gr.val(x$gr, y, val = names(values(y)))
  do.call(x$metaset, as.list(values(gr.new)))
})





#' @name gm2dat
#' @description
#'
#' Runs a glm on a gMatrix treating node (GRanges) metadata as (1D, marginal) covariates and
#' allows optional named gMatrix arguments which are treated as "bivariates" (ie 2D covariates) in a
#' glm that has the following form:
#'
#' value_ij ~ sum_k(alpha_k*cov_ik + beta_k*cov_jk) + sum_k(gamma_k*cov_jk*cov_jk) + sum_k(delta_k*biv_ijk)
#'
#' where alpha_k, beta_k, gamma_k, and delta_k are the inferred coefficients on the 1D covariates (alpha, beta),
#' the bivariate interactions of the marginal covariates (gamma_k) and any provided bivariates (delta_k)
#'
#' 
#' @param gm gMatrix object whose value is to be modeled
#' @param covariates 1D covariate metadata fields (default names(values(gm$gr)))
#' @param interactions logical flag (default TRUE) whether to model interaction
#' @param offset gMatrix to use as an offset variable (default NULL)
#' @param ... named arguments specifying gMatrix objects to use as "bivariates"
#' @return data.table
#' @author Marcin Imielinski                         
gm2dat = function(gm, covariates = NULL, interactions = TRUE, offset = NULL, family = gaussian, verbose = FALSE, ...)
{
  if (is(gm, 'GRanges'))
    gm = gM(gm)

  if (!is.null(offset) && !is(offset, 'gMatrix'))
    stop('offset argument must be a gMatrix')

  if (is.null(covariates))
    covariates = c()
#    covariates = names(values(gm$gr))
  bivariates = list(...)

  if (length(bivariates)>0 | !is.null(offset)) ## harmonize gm with all bivariates and offsets via agg
  {
    if (verbose)
    {
      ggmessage('harmonizing primary and offset / bivariate gMatrix objects')
    }

    if (!all(sapply(lapply(bivariates, 'class'), '[', 1)=='gMatrix'))
      stop('please check your arguments - all bivariates provided as named arguments must be gMatrix objects')

    bivoff = NULL
    if (length(bivariates)>0)
      bivoff = bivariates

    if  (!is.null(offset))
      bivoff = c(bivoff, list(offset))

    ## disjoin and remerge with gm$gr to get gr metadata (also limiting to footprint of gm)
    dgr = disjoin(do.call(grbind, c(list(gm$gr),lapply(bivoff, function(x) x$gr)))) %*% gm$gr

    if (verbose)
    {
      ggmessage('disjoined intervals')
    }

    ## refactor this gm and all bivariates around dgr, and make each full
    if (length(bivariates)>0)
    {
      bivariates = lapply(bivariates, function(x) x$agg(dgr)+0)
      if (verbose)
      {
        ggmessage('aggregated bivariate gMatrix objects')
      }
    }

    if (!is.null(offset))
    {
      offset = offset$agg(dgr)+0
      if (verbose)
      {
        ggmessage('aggregated offset gMatrix')
      }
    }

    gm = gm$agg(dgr)
    if (verbose)
    {
      ggmessage('aggregated primary gMatrix')
    }
  }

  gm = gm+0 ## make gm full
  dat = gm$dat
  dat$id = NULL

  if (length(covariates)>0)
  {
    covnmsi = paste0('covi_', covariates)
    covnmsj = paste0('covj_', covariates)
    
    for (k in 1:length(covariates))
    {
      dat[[covnmsi[k]]] = values(gm$gr)[[covariates[k]]][dat$i]
      dat[[covnmsj[k]]] = values(gm$gr)[[covariates[k]]][dat$j]
    }
  }

  if (length(bivariates)>0) ## bivariates should be synced with gm at this point, ie dat has identical dimensions
  {
    bivnms= paste0('biv_', names(bivariates))
    for (k in 1:length(bivariates)) ## populate 
    {
      dat[[bivnms[k]]]= bivariates[[k]]$dat$value
    }
  }

  if (!is.null(offset)) ## offset is either a field of gr or a gMatrix
  {
    dat[['off']] = offset$dat$value
  }

  if (verbose)
  {
    ggmessage('populated glm data matrix')
  }

  return(list(dat = dat, gr = gm$gr))
}

#' @name gglm
#' @description
#'
#' Runs a glm on a gMatrix treating node (GRanges) metadata as (1D, marginal) covariates and
#' allows optional named gMatrix arguments which are treated as "bivariates" (ie 2D covariates) in a
#' glm that has the following form:
#'
#' value_ij ~ sum_k(alpha_k*cov_ik + beta_k*cov_jk) + sum_k(gamma_k*cov_jk*cov_jk) + sum_k(delta_k*biv_ijk)
#'
#' where alpha_k, beta_k, gamma_k, and delta_k are the inferred coefficients on the 1D covariates (alpha, beta),
#' the bivariate interactions of the marginal covariates (gamma_k) and any provided bivariates (delta_k)
#'
#' 
#' @param data gMatrix object whose value is to be modeled
#' @param covariates 1D covariate metadata fields (default names(values(data$gr)))
#' @param interactions logical flag (default TRUE) whether to model interaction
#' @param family family to use (default gaussian)
#' @param nb logical flag whether to use glm.nb (FALSE)
#' @param zinb logical flag whether to use zeroinfl negative binomial model (FALSE)
#' @param zinp logical flag whether to use zeroinfl poisson model (FALSE)
#' @param offset gMatrix to use as an offset variable (default NULL)
#' @param subsample integer number of data points to subsample (default NULL)
#' @param ... named arguments specifying gMatrix objects to use as "bivariates"
#' @return model
#' @author Marcin Imielinski                         
#' @export
gglm = function(data, covariates = NULL, interactions = TRUE, offset = NULL,
                nb = FALSE, zinb = FALSE, zinp = FALSE, 
                family = gaussian, na.action = na.pass, subsample = NULL,
                verbose = FALSE,
                ...)
{
  ## make data data.table from provided gMatrix object data
  dat = gm2dat(data, covariates = covariates, interactions = interactions, offset = offset,
               verbose = verbose,
               ...)$dat

  .clean = function(x) gsub('\\W', "_", x)

  vars = setdiff(names(dat), c('value', 'id','i','j'))

  ## convert vars to formula (a bit of a hack using the names)
  covi = .clean(grep('^covi', vars, value = TRUE))
  covj = .clean(gsub('^covi', 'covj',covi))

  covij = c()
  if (interactions & length(covi))
    covij = paste0(covi, '*', covj)

  bivj = c()
  bivariates = list(...)
  if (length(bivariates)>0)
    bivj = paste0('biv_', .clean(names(bivariates)))

  off  = c()
  if (!is.null(offset))
    off = 'offset(off)'

  fmstring = paste('value ~', paste(c(covi,covj,covij, bivj, off), collapse = ' + '))
  fm = formula(fmstring)
  if (verbose)
  {
    ggmessage('running glm with formula:\n\t', fmstring)
  }

  ## if subsample integer provided then will run on subset of data
  if (!is.null(subsample) && subsample<nrow(dat))
    dat = dat[sample(.N, subsample), ]

  setnames(dat, .clean(names(dat)))
  if (nb)
    model = MASS::glm.nb(dat, formula = fm, na.action = na.action)
  else if (zinb)
    model = pscl::zeroinfl(dat, formula = fm, dist = 'negbin', na.action = na.action)
  else if (zinp)
    model = pscl::zeroinfl(dat, formula = fm, dist = 'poisson', na.action = na.action)
  else
    model = glm(dat, formula = fm, family = family, na.action = na.action)

  ## keep track of covariate names etc, will be useful for predict phase
  model$covariates = covariates
  model$bivariates = names(bivariates)
  model$has.offset = !is.null(offset)
  model$interactions = interactions
  model$type = 'gglm'
  return(model)
}


#' @name gpredict
#' @description
#'
#' Applies gglm "gmodel" (ie output of gglm) to a new gMatrix +/- covariates treating node (GRanges)
#' metadata as (1D, marginal) covariates and
#' allows optional named gMatrix arguments which are treated as "bivariates" (ie 2D covariates) the 
#' glm that has the following form:
#'
#' value_ij ~ sum_k(alpha_k*cov_ik + beta_k*cov_jk) + sum_k(gamma_k*cov_jk*cov_jk) + sum_k(delta_k*biv_ijk)
#'
#' where alpha_k, beta_k, gamma_k, and delta_k are the inferred coefficients on the 1D covariates (alpha, beta),
#' the bivariate interactions of the marginal covariates (gamma_k) and any provided bivariates (delta_k)
#'
#' The data must be compatible with the provided model (i.e. the same named covariates must be provided 
#' 
#' @param gmodel must have been created by gglm or gglm.nb, has fields $covariates, $bivariates, $offset, $interactions, $type = 'gglm'
#' @param newdata GRanges with metadata
#' @param offset gMatrix offset
#' @param type character specifying linnk, response, or terms
#' @param ... named arguments specifying gMatrix objects to use as "bivariates"
#' @return gMatrix with prediction values which are (based on type var) 
#' @author Marcin Imielinski                         
#' @export
gpredict = function(gmodel, newdata, offset = NULL, type = 'response', terms = NULL, na.action = na.pass, verbose = FALSE, na.rm = TRUE, agg.fun = sum, fill = 0, full = FALSE, ...)
{
  ## make data matrix from provided data
  if (!is.null(gmodel$type) && gmodel$type != 'gglm')
    stop('provided model must be output of gglm, gglm.nb, or other related model fitting funtions in GxG')

  if (gmodel$has.offset && is.null(offset))
    stop('provided model has offset term, please provide offset input')

  bivariates = list(...)

  if (length(gmodel$bivariates) && length(setdiff(gmodel$bivariates, names(bivariates))))
    stop('provided model has bivariate terms for which values were not provided as named arguments to gpredict, see ?gpredict documentation')

 if (length(gmodel$covariates) && length(setdiff(gmodel$covariates, names(values(newdata)))))
    stop('provided model has covariate terms for which columns do not exist in provided gMatrix object, please check inputs and see ?gpredict documentation')

  covariates = gmodel$covariates
  interactions = gmodel$interactions
  if (length(bivariates)>0)
    bivariates = bivariates[gmodel$bivariates]

  ## formulate data matrix using provided gm
  tmp = gm2dat(newdata, covariates = covariates, interactions = interactions, offset = offset,
               verbose = verbose,
               ...)
  dat = tmp$dat
  .clean = function(x) gsub('\\W', "_", x)
  setnames(dat, .clean(names(dat)))
  val = predict(gmodel, dat, type = type, terms = terms, na.action = na.action)
  dat[, value := val]
  return(gMatrix$new(tmp$gr, dat, fill = fill, full = full, agg.fun = agg.fun, na.rm = na.rm))
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
  if (!is(mat1, 'gMatrix'))
    mat1 = gM(mat2$gr, mat1)

  if (!is(mat2, 'gMatrix'))
    mat2 = gM(mat1$gr, mat2)
  
  tmp = gmatalign(mat1, mat2)
  gr.new = tmp$gr
  dat.new = tmp$dat

  dat.new[, value := val1 & val2]

  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg, na.rm = mat1$na.rm))
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
  if (!is(mat1, 'gMatrix'))
    mat1 = gM(mat2$gr, mat1)

  if (!is(mat2, 'gMatrix'))
    mat2 = gM(mat1$gr, mat2)
  
  tmp = gmatalign(mat1, mat2)
  gr.new = tmp$gr
  dat.new = tmp$dat

  dat.new[, value := val1 | val2]

  return(gMatrix$new(gr.new, dat.new, fill = mat1$fill, full = mat1$full, agg.fun = mat1$agg, na.rm = mat1$na.rm))
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
  mat = mat+0 ## make full
  dat.new = mat$dat[, value := log(value)]
  gr.new = mat$gr
  
  return(gMatrix$new(gr.new, dat.new, fill = mat$fill, full = mat$full, agg.fun = mat$agg.fun,
                     na.rm = mat$na.rm))
}


#' @name exp.gMatrix
#' @description
#'
#' Allows exp of matrix values
#' 
#' @param mat gMatrix object
#' @return a new gMatrix object whose values are the exp of the inputs
#' @author Marcin Imielinski                         
#' @export
'exp.gMatrix' = function(mat)
{
  mat = mat+0 ## make full
  dat.new = mat$dat[, value := exp(value)]
  gr.new = mat$gr
  
  return(gMatrix$new(gr.new, dat.new, fill = mat$fill, full = mat$full, agg.fun = mat$agg.fun,
                     na.rm = mat$na.rm))
}


#' @name round.gMatrix
#' @description
#'
#' Round matrix values
#' 
#' @param mat gMatrix object
#' @return a new gMatrix object whose values are rounded
#' @author Marcin Imielinski                         
#' @export
'round.gMatrix' = function(mat, digits = 0)
{
  mat = mat+0 ## make full
  dat.new = mat$dat[, value := round(value, digits = digits)]
  gr.new = mat$gr
  
  return(gMatrix$new(gr.new, dat.new, fill = mat$fill, full = mat$full, agg.fun = mat$agg.fun,
                     na.rm = mat$na.rm))
}


#' @name signif.gMatrix
#' @description
#'
#' Round matrix values to signif digits
#' 
#' @param mat gMatrix object
#' @return a new gMatrix object whose values are rounded
#' @author Marcin Imielinski                         
#' @export
'signif.gMatrix' = function(mat, digits = 6)
{
  mat = mat+0 ## make full
  dat.new = mat$dat[, value := signif(value, digits = digits)]
  gr.new = mat$gr
  
  return(gMatrix$new(gr.new, dat.new, fill = mat$fill, full = mat$full, agg.fun = mat$agg.fun,
                     na.rm = mat$na.rm))
}


#' @name floor.gMatrix
#' @description
#'
#' Round matrix values
#' 
#' @param mat gMatrix object
#' @return a new gMatrix object whose values are rounded
#' @author Marcin Imielinski                         
#' @export
'floor.gMatrix' = function(mat)
{
  mat = mat+0 ## make full
  dat.new = mat$dat[, value := floor(value)]
  gr.new = mat$gr
  
  return(gMatrix$new(gr.new, dat.new, fill = mat$fill, full = mat$full, agg.fun = mat$agg.fun,
                     na.rm = mat$na.rm))
}



#' @name ceiling.gMatrix
#' @description
#'
#' Round matrix values
#' 
#' @param mat gMatrix object
#' @return a new gMatrix object whose values are rounded
#' @author Marcin Imielinski                         
#' @export
'ceiling.gMatrix' = function(mat)
{
  mat = mat+0 ## make full
  dat.new = mat$dat[, value := ceiling(value)]
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

#' @name max.gMatrix
#' @description
#'
#' Returns max value 
#' 
#' @param mat gMatrix object
#' @return maximum value of gMatrix
#' @author Marcin Imielinski                         
#' @export
'max.gMatrix' = function(mat, na.rm = mat$na.rm)
{
  return(max(mat$value, na.rm = na.rm))
}

#' @name min.gMatrix
#' @description
#'
#' Returns min value 
#' 
#' @param mat gMatrix object
#' @return maximum value of gMatrix
#' @author Marcin Imielinski                         
#' @export
'min.gMatrix' = function(mat, na.rm = mat$na.rm)
{
  return(min(mat$value, na.rm = na.rm))
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


#' @name seqinfo
#' @title seqinfo
#' @description
#'
#' @param x a gMatrix
#'
#' @return the seqinfo of this gMatrix
#' @exportMethod seqinfo
#' @export
setMethod("seqinfo", c("gMatrix"),
          function(x) {
            return(seqinfo(x$gr))
          })


#' @name seqlengths
#' @title seqlengths
#' @description
#'
#' @param x a gMatrix
#'
#' @return the seqlengths of this graph
#' @exportMethod seqlengths
#' @export
setMethod("seqlengths", c("gMatrix"),
          function(x) {
            return(seqlengths(x$gr))
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
                        if (is.null(gr2))
                          gr2 = gr1

                        if (is.null(gr1))
                          gr1 = gr2

                        if (is.null(gr1))
                          {
                            private$pgr1 = GRanges()
                            private$pgr2 = GRanges()
                            private$pmeta = data.table()
                            return(self)
                          }

                        if (!inherits(gr1, 'GRanges') |
                            !inherits(gr2, 'GRanges'))
                          stop('gPair can only be initialized from GRanges')

                        grix = data.table(i = 1:length(gr1),
                                          j = 1:length(gr2))
                        gr1 = gr1[grix$i]
                        gr2 = gr2[grix$j]

                        if (is.null(meta))
                          meta = as.data.table(values(gr1))

                        if (length(gr1) != length(gr1))
                          stop('Inputs gr1 and gr2 must be the same length')
                                                
                        if (length(gr1)>0)
                        {
                          if (nrow(meta)>0 & nrow(meta)!=length(gr1))
                            stop('meta must have the same number of rows as gr1 and gr2')
                        }
                        
                        private$pgr1 = gr.stripstrand(gr1)
                        private$pgr2 = gr.stripstrand(gr2)
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

                      #' @name gtrack
                      #' @description
                      #' Outputs a gTrack of this gPair object where each interval is a GRangesList
                      #' @param colormap colors to use with this colormap, default = c('white', 'red', 'black')
                      #' @param ... gTrack arguments
                      #' @author Marcin Imielinski                         
                      gtrack = function(name = '', stack.gap = 1e9, cex.label = 0.3, col = alpha('gray', 0.5), ...)                                        
                      {
                        return(gTrack(name = name, data = self$grl, cex.label = cex.label, stack.gap = stack.gap, col = col, ...))
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

                        private$pgr1 = private$pgr1[i]
                        private$pgr2 = private$pgr2[i]
                        private$pmeta = private$pmeta[i]
                        
                        return(self)
                      },

                      #' @name set
                      #' @description
                      #'
                      #' sets metadata of gPair, either adding columns or changing
                      #' existing values
                      #' @param ... name value pairs where value is a vector / scalar that broadcasts to $gr
                      #' @author Marcin Imielinski
                      set = function(...)
                       {
                         args = list(...)
                         
                         for (arg in names(args)){
                           private$pmeta[[arg]] = args[[arg]]
                         }
                         return(invisible(self))
                       },

                      #' @name print
                      #' @description
                      #' 
                      #' Prints out the gPair Object. Prints the length and the GRangesList of the junctions.
                      print = function()
                      {
                        ggmessage("gPair Object containing  ", self$length, " pairs.")
                        if (self$length>0)
                        {
                          HEAD = 1:pmin(4, self$length)
                          pairdt = self$dt

                          if (ncol(self$dt)>0)
                          {
                            pairdt = self$dt[HEAD, ]
                          }
                          print(pairdt)
                          more = self$length-HEAD[length(HEAD)]
                          if (more>0)
                            ggmessage('... (', more,' additional pairs)')
                        }
                      }
                    ),

                    private = list(
                      pgr1 = NULL, ## second interval
                      pgr2 = NULL, ## first interval
                      pmeta = NULL ## data.table
                    ),
                    
                    active = list(
                      #' @name grl
                      #' @description
                      #'
                      #' Returns the GRangesList of the pairs in this gPair Object
                      #'
                      #' @return GRangesList of the pairs in this gPair Object
                      #' @author Marcin Imielinski                         
                      grl = function()
                      {
                        grl = split(grbind(private$pgr1, private$pgr2),
                                    rep(1:length(private$pgr1), 2))
                        if (nrow(private$pmeta)>0)
                          values(grl) = private$pmeta
                        return(grl)
                      },

                      #' @name gr1
                      #' @description
                      #'
                      #' Returns the GRangesList of the "first" granges in each pair
                      #'
                      #' @return GRanges of "first" granges in each pair
                      #' @author Marcin Imielinski                         
                      gr1 = function()
                      {
                        return(private$pgr1)
                      },

                      #' @name gr2
                      #' @description
                      #'
                      #' Returns the GRangesList of the "first" granges in each pair
                      #'
                      #' @return GRanges of "first" granges in each pair
                      #' @author Marcin Imielinski                         
                      gr2 = function()
                      {
                        return(private$pgr2)
                      },

                      #' @name gt
                      #' @description
                      #' Outputs a gTrack of this gPair object where each interval is a GRangesList
                      #' @param colormap colors to use with this colormap, default = c('white', 'red', 'black')
                      #' @param ... gTrack arguments
                      #' @author Marcin Imielinski                         
                      gt = function()                                     
                      {
                        return(self$gtrack())
                      },
                                            
                      #' @name length
                      #' @description
                      #' 
                      #' Returns the number of pairs in this gPair Object.
                      #' 
                      #' @return Number of pairs in this gPair Object
                      length = function()
                      {                               
                        return(length(private$pgr1))
                      },

                      copy = function() self$clone(),

                      #' @name dt
                      #' @description
                      #'
                      #' Returns the GRangesList of the pairs in the gPair Object as a data.table.
                      #'
                      #' @return data.table GRangesList of the pairs coverted to a data.table
                      dt = function()
                      {
                        out = data.table(gr1 = gr.string(private$pgr1),
                                         gr2 = gr.string(private$pgr2))
                        if (ncol(private$pmeta))
                          out = cbind(out, copy(private$pmeta))
                        return(out)
                      },

                      #' @name meta
                      #' @description
                      #'
                      #' Returns the metadata of this 
                      #'
                      #' @return data.table GRangesList of the pairs coverted to a data.table
                      meta = function()
                      {
                        return(private$pmeta)
                      },

                      #' @name span
                      #' @description
                      #' 
                      #' Returns the distance between breakpoint pairs on the genome
                      #' 
                                        #
                      #' @return vector of spans of all pairs in junction object
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

#' @name gP
#' @title gP
#' @description
#'
#' gPair instantiator, creates a pair with some metadata
#' 
#' @param gr1 GRanges around which to build a gPair, also can be a GRAngesList, in which case all intra-list item pairs are created.
#' @param gr2 optional second granges same length as the first (but possibly different width), by default is set to gr1
#' @param meta data.table of metadata (by default pulled from values(gr1)
#' @return A new gPair object
#' @author Marcin Imielinski                         
#' @export
gP = function(gr1 = NULL, gr2 = NULL,  meta = NULL)
{
  if (is(gr1, 'GRangesList'))
  {
    tmpgr = grl.unlist(grl)
    tmpgr$ix = 1:length(tmpgr)
    tmpdt = gr2dt(tmpgr)[, as.data.table(expand.grid(i = ix, j = ix))[i<=j, ], by = grl.ix]

    gr1 = tmpgr[tmpdt$i]
    gr2 = tmpgr[tmpdt$j]
  }
  return(gPair$new(gr1, gr2, meta))
}


#' @name length
#' @title length.gPair
#' @description
#' 
#' The number of pairs in this gPair Object
#'
#' @param gPair a gPair Object
#' @return the number of pairs in the gPair Object
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
  gr1l = lapply(pairs.list, function(x) x$gr1)
  gr2l = lapply(pairs.list, function(x) x$gr2)
  grlm = lapply(pairs.list, function(x) x$meta)
  gr1.new = dodo.call(grbind, lapply(gr1l, function(x) {values(x) = NULL; x}))
  gr2.new = dodo.call(grbind, lapply(gr2l, function(x) {values(x) = NULL; x}))
  meta.new = rbindlist(grlm, fill = TRUE)
  return (gPair$new(gr1.new, gr2.new, meta.new))
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
      inew = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(i)))), parent.frame()), pairs$meta, parent.frame(2)), error = function(e) NULL)
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
    pairs = gP(obj$gr1, obj$gr2, obj$meta[, j])
  }
  return(pairs)
}


#' @name +.gPair
#' @description
#'
#' Allows padding of pairs. The rpad will be added to the left of a "+" junction
#' and to the right of "-" junction.
#'
#' @param jj gPair Object
#' @param pad Positive number representing amount to pad the gPair Object.
#' @return a new gPair class with the padding applied
#' @export
'+.gPair' = function(gp, pad)
{
  return(gPair$new(gp$gr1 + pad, gp$gr2 + pad, meta = gp$meta))
}

#' @name -.gPair
#' @description
#'
#' Allows padding of pairs. The rpad will be added to the left of a "+" junction
#' and to the right of "-" junction.
#'
#' @param jj gPair Object
#' @param pad Positive number representing amount to pad the gPair Object.
#' @return a new gPair class with the padding applied
#' @export
'-.gPair' = function(gp, pad)
{
  return(gPair$new(gp$gr1 - pad, gp$gr2 - pad, meta = gp$meta))
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

#' @name seqinfo
#' @title seqinfo
#' @description
#'
#' @param x a gPair
#'
#' @return the seqinfo of this gPair
#' @exportMethod seqinfo
#' @export
setMethod("seqinfo", c("gPair"),
          function(x) {
            return(seqinfo(x$gr1))
          })


#' @name seqlengths
#' @title seqlengths
#' @description
#'
#' @param x a gPair
#'
#' @return the seqlengths of this graph
#' @exportMethod seqlengths
#' @export
setMethod("seqlengths", c("gPair"),
          function(x) {
            return(seqlengths(x$gr1))
          })



#' @name ggmessage
#' @title ggmessage
#' @rdname internal
#' @return timestamped GxG message
ggmessage = function(..., pre = 'GxG')
    message(pre, ' ', paste0(as.character(Sys.time()), ': '), ...)

