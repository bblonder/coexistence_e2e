# assumes input has the initial abundances on lines 2-10 (9 species total)
n_sp = 9
offset = 2
input = readLines('GMF with correct weibell dispersal and output.xml')

dir.create("inputs")

combinations = expand.grid(data.frame(t(matrix(nrow=n_sp,ncol=2,data=c(0,1),byrow=TRUE))))

for (i in 1:nrow(combinations))
{
  # copy over all the text
  output_this = input
  
  # modulate the abundances from 25 to 0 for each species
  ids_to_edit = which(combinations[i,]==0)
  if (length(ids_to_edit)>0)
  {
    output_this[offset+ids_to_edit-1] = gsub("25.0","0.0",output_this[offset+ids_to_edit-1],fixed=TRUE)
    message('.')
  }
  
  # rewrite output file
  file_id = paste(combinations[i,],collapse="-")
  output_this[length(output_this)] = gsub("GMFstr1.out",sprintf("outputs/GMF_%s.out",file_id), output_this[length(output_this)],fixed=TRUE)

  output_this_final = paste(output_this,collapse="\n")
  writeLines(text=output_this_final,con = sprintf("inputs/input_GMF_%s.xml",file_id))
  
  # put all input files on the virtual machine, then in SORTIE create a new batch file, running each input 5x
}