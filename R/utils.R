#' @name alpha
#' @title alpha
#' @description
#' Give transparency value to colors
#'
#' Takes provided colors and gives them the specified alpha (ie transparency) value
#'
#' @author Marcin Imielinski
#' @param col RGB color
#' @keywords internal
#' @export
alpha = function(col, alpha)
{
  col.rgb = col2rgb(col)
  out = rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
  names(out) = names(col)
  return(out)
}

#' @name dodo.call
#' @title dodo.call
#' @description
#' do.call implemented using eval parse for those pesky (e.g. S4) case when do.call does not work
dodo.call = function(FUN, args)
{
    if (!is.character(FUN))
        FUN = substitute(FUN)
    cmd = paste(FUN, '(', paste('args[[', 1:length(args), ']]', collapse = ','), ')', sep = '')
    return(eval(parse(text = cmd)))
}

