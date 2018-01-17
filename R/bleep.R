#' bleep
#' 
#' Rings a terminal bell after long running processes. It just calls a system
#' command called "bleep" which should look something like:
#' #!/bin/bash
#' paplay /usr/share/sounds/freedesktop/stereo/complete.oga &
#'
#' @return A noise of your choosing
#' @export
bleep <- function(){
  system('bleep &')
# this puts up a notification in the gnome notification tray
  system('notify-send -t 3000 "tcrSeqR" "Done"')
}
