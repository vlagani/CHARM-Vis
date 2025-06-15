# Function for converting gene names between mouse and chicken

# input: 
# - names_to_convert: names to convert
# - organism_from: mouse/chicken
# - organism_to: mouse/chicken. Can be the same as organism_from
# - type_fromm: ensembl/name/mixed. The type of gene names in names_to_convert. 
# - type_to: ensembl/name/mixed. The type of gene names to return
# - annot_table: table with the information for the conversion
# output:
# - converted_names: vector with the converted gene names. 
# - converted_ids: Boolean vector indicating which names are converted. 
# notes:
# - names that cannot be converted are left unchanged

# test
# rm(list = ls())
# load('annot_table.RData')
# names_to_convert <- annot_table$Gene.name.chicken[1:10]
# #names_to_convert[c(3, 7)] <- annot_table$Gene.stable.ID.chicken[c(3, 7)]
# names_to_convert[c(4, 8)] <- NA
# organism_from = 'chicken'
# organism_to = 'chicken'
# #type_from = 'mixed'
# type_from = 'name'
# #type_to = 'mixed'
# type_to = 'ensembl'
# annot_table = NULL

# test
# names_to_convert = rownames(test_count_values)
# organism_from = 'chicken'
# organism_to = 'mouse' 
# type_from = 'mixed'
# type_to = 'ensembl' 
# annot_table = annot_table

#test
# names_to_convert = gs_list[[1]][[5]]
# organism_from = 'human'
# organism_to = 'mouse'
# type_from = 'name'
# type_to = 'name' 
# annot_table = annot_table

# function
name_conversion <- function(names_to_convert, organism_from, 
                            organism_to, 
                            type_from = 'name', type_to = 'name',
                            annot_table = NULL){
  
  # check
  if(type_to == 'mixed' && type_from != 'mixed'){
    stop('if "type_to" is mixed, also "type_from" must be mixed.')
  }
  
  # doouble conversion in case of mixed
  if(type_from == 'mixed' && type_to == 'mixed'){
    
    # first conversion
    tmp <- name_conversion(names_to_convert = names_to_convert, 
                           organism_from = organism_from,
                           organism_to = organism_to,
                           type_from = 'name', 
                           type_to = 'name',
                           annot_table = annot_table)
    converted_ids <- tmp$converted_ids
    
    # second conversion 
    tmp <- name_conversion(names_to_convert = tmp$converted_names, 
                           organism_from = organism_from, 
                           organism_to = organism_to,
                           type_from = 'ensembl', 
                           type_to = 'ensembl',
                           annot_table = annot_table)
    converted_ids <- converted_ids | tmp$converted_ids
    
    # return
    return(list(converted_names = tmp$converted_names, 
                converted_ids = converted_ids))
    
  }
  if(type_from == 'mixed' && !(type_to == 'mixed')){
    
    # first conversion
    tmp <- name_conversion(names_to_convert = names_to_convert, 
                           organism_from = organism_from, 
                           organism_to = organism_to,
                           type_from = 'name', 
                           type_to = type_to,
                           annot_table = annot_table)
    converted_ids <- tmp$converted_ids
    
    # second conversion 
    tmp <- name_conversion(names_to_convert = tmp$converted_names, 
                           organism_from = organism_from, 
                           organism_to = organism_to, 
                           type_from = 'ensembl', 
                           type_to = type_to,
                           annot_table = annot_table)
    converted_ids <- converted_ids | tmp$converted_ids
    
    # return
    return(list(converted_names = tmp$converted_names, 
                converted_ids = converted_ids))
    
  }
  
  # loading the annotation table
  if(is.null(annot_table)){
    #load('annot_table.RData')
    print('Annotation table not provided')
    return()
  }
  
  # organism to
  #organism_to <- ifelse(organism_from == 'chicken', 'mouse', 'chicken')
  
  # names from
  if(type_from == 'name'){
    names_from = annot_table[[paste0('Gene.name.', organism_from)]]
  }
  if(type_from == 'ensembl'){
    names_from = annot_table[[paste0('Gene.stable.ID.', organism_from)]]
  }
  
  # names to
  if(type_to == 'name'){
    names_to = annot_table[[paste0('Gene.name.', organism_to)]]
  }
  if(type_to == 'ensembl'){
    names_to = annot_table[[paste0('Gene.stable.ID.', organism_to)]]
  }
  
  # named vector of names
  names(names_to) <- names_from
  
  # which names can be converted?
  converted_ids <- (names_to_convert %in% names_from) &
                      !is.na(names_to_convert) &
                      names_to_convert != ''
  
  # converted names
  converted_names <- names_to[names_to_convert]
  converted_names[!converted_ids] <- names_to_convert[!converted_ids]
  names(converted_names) <- NULL
  
  # return
  return(list(converted_names = converted_names, 
              converted_ids = converted_ids))
  
}

# test
# name_conversion(names_to_convert, organism_from, organism_to,
#                 type_from, type_to, annot_table)
# load('annot_table.RData')
# annot_table[annot_table$Gene.stable.ID.chicken == 'ENSGALG00000036909', ]
# annot_table[annot_table$Gene.stable.ID.mouse == 'ENSMUSG00000064345', ]
