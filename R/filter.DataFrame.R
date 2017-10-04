#this makes filter work with DataFrame
filter.DataFrame <- function(.data,...) {
	    DataFrame(filter(as.data.frame(.data),...))
}
