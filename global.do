 adopath + "H:\My Documents\Stepped Wedge Trials\PhD\Stata_swpermute\Code"

 adopath + "C:\Users\lsh304314\Filr\My Files\My Documents\Stepped Wedge Trials\PhD\Stata_swpermute\Code"
 
  global ProgDir "H:\My Documents\Stepped Wedge Trials\PhD\Stata_swpermute\Code"

 
 global DataDir "H:\My Documents\Data\XpertMTBRIF\Mar2017"

window menu append submenu "stUser" "&Cluster RCTs"
window menu append item "Cluster RCTs" "Permute for stepped wedge trials (&swpermute)" "db swpermute"
window menu refresh

