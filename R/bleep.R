# an R terminal bell (require a shell script or alias which makes a noise and exits)
bleep <- function(){
  system('bleep &')
  system('notify-send -t 3000 "immunoSeqR" "Done"')
}
