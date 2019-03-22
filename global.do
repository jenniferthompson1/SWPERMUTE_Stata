 adopath + "H:\My Documents\PhD\Stata_command\progs"

  global ProgDir "H:\My Documents\PhD\Stata_command\progs"

 
 global DataDir "H:\My Documents\Data\XpertMTBRIF\Mar2017"

window menu append submenu "stUser" "&Cluster RCTs"
window menu append item "Cluster RCTs" "Permute for stepped wedge trials (&swpermute)" "db swpermute"
window menu refresh

