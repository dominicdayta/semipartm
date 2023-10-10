# codes for combining multiple results into one file
partial_file_list = list(
  "3000doc_3500word_0.7spar_m1.RDS" = c(
    "3000doc3500word_partial/3000doc_3500word_0.7spar_m1_1-20.RDS"
  ),
  "3000doc_3500word_0.7spar_m2.RDS" = c(
    "3000doc3500word_partial/3000doc_3500word_0.7spar_m2_1-20.RDS"
  ),
  "3000doc_3500word_0.9spar_m1.RDS" = c(
    "3000doc3500word_partial/3000doc_3500word_0.9spar_m1_1-20.RDS"
  ),
  "3000doc_3500word_0.9spar_m2.RDS" = c(
    "3000doc3500word_partial/3000doc_3500word_0.9spar_m2_1-20.RDS"
  ),
  "3000doc_3500word_0.99spar_m1.RDS" = c(
    "3000doc3500word_partial/3000doc_3500word_0.99spar_m1_1-20.RDS"
  ),
  "3000doc_3500word_0.99spar_m2.RDS" = c(
    "3000doc3500word_partial/3000doc_3500word_0.99spar_m2_1-20.RDS"
  ),
  "3000doc_1500word_0.7spar_m1.RDS" = c(
    "3000doc1500word_partial/3000doc_1500word_0.7spar_m1_1-20.RDS"
  ),
  "3000doc_1500word_0.7spar_m2.RDS" = c(
    "3000doc1500word_partial/3000doc_1500word_0.7spar_m2_1-20.RDS"
  ),
  "3000doc_1500word_0.9spar_m1.RDS" = c(
    "3000doc1500word_partial/3000doc_1500word_0.9spar_m1_1-20.RDS"
  ),
  "3000doc_1500word_0.9spar_m2.RDS" = c(
    "3000doc1500word_partial/3000doc_1500word_0.9spar_m2_1-20.RDS"
  ),
  "3000doc_1500word_0.99spar_m1.RDS" = c(
    "3000doc1500word_partial/3000doc_1500word_0.99spar_m1_1-20.RDS"
  ),
  "3000doc_1500word_0.99spar_m2.RDS" = c(
    "3000doc1500word_partial/3000doc_1500word_0.99spar_m2_1-20.RDS"
  )
)

# file directories
base_dir = function(x){
  return(paste0("/Users/dominic/Documents/semipartm_thesis/Notebooks/Data/",x))
}

for(current_file_name in names(partial_file_list)){
  partial_files = partial_file_list[[current_file_name]]
  file = 0
  
  for(partial_file in partial_files){
    file = file + 1
    file_contents = readRDS(base_dir(partial_file))
    
    # now combine results
    if(file == 1){
      class.stat = file_contents$class.stat
      method.stat = file_contents$method.stat
      twolevel.B = file_contents$twolevel.B
      cv.list = file_contents$cv.list
      nword = file_contents$nword
      ndoc = file_contents$ndoc
      spar = file_contents$spar
      m = file_contents$m
      comp.time = file_contents$comp.time
      
    }else{
      class.stat = rbind(file_contents$class.stat,class.stat)
      method.stat = rbind(file_contents$method.stat, method.stat)
      twolevel.B = rbind(file_contents$twolevel.B, twolevel.B)
      cv.list = c(file_contents$cv.list,cv.list)
      comp.time = file_contents$comp.time + comp.time
      
    }
  }
  
  saveRDS(list(
    "class.stat" = class.stat,
    "method.stat" = method.stat,
    "twolevel.B" = twolevel.B,
    "cv.list" = cv.list,
    "nword" = nword,
    "ndoc" = ndoc,
    "spar" = spar,
    "m" = m,
    "comp.time" = comp.time
  ), base_dir(current_file_name))
  
}
