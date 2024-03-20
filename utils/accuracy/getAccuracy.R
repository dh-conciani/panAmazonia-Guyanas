## get accuracy for test version in the mapbiomas collection 8.0
## --- --- --- EXPORT ACCURACY METRICS PER ECORREGION

## get libraries
library(rgee)
library(reticulate)
library(caret)
library(reshape2)

## initialize
ee_Initialize()

## set directories
file_path <- 'projects/mapbiomas-workspace/public/'

##  define files to be computed
file_name <- c( #'collection8/mapbiomas_collection80_integration_v1',
  # 'CERRADO_col8_gapfill_incidence_temporal_frequency_geomorphology_spatial_v13',
  # 'CERRADO_col8_gapfill_incidence_temporal_frequency_geomorphology_spatial_v5',
  # 'CERRADO_col8_gapfill_incidence_temporal_frequency_geomorphology_spatial_v6',
  # 'CERRADO_col8_gapfill_incidence_temporal_frequency_geomorphology_spatial_v7',
  # 'CERRADO_col8_gapfill_incidence_temporal_frequency_geomorphology_spatial_v8',
  # 'CERRADO_col8_gapfill_incidence_temporal_frequency_geomorphology_spatial_v9',
  # 'CERRADO_col8_gapfill_incidence_temporal_frequency_geomorphology_spatial_v10',
  #'CERRADO_col8_gapfill_incidence_temporal_frequency_geomorphology_spatial_v11',
  'collection7_1/mapbiomas_collection71_integration_v1'
)

## set output path (local)
output <- './tables/'

## import validation points
validation_points = ee$FeatureCollection('projects/mapbiomas-workspace/VALIDACAO/mapbiomas_85k_col2_points_w_edge_and_edited_v1')

## set classes to select from validation dataset
selectClasses = c(
  'LAVOURA TEMPORÁRIA',            
  'LAVOURA PERENE',           
  'CANA',      
  'PASTAGEM',       
  'FORMAÇÃO FLORESTAL',
  'RIO, LAGO E OCEANO',
  'FORMAÇÃO CAMPESTRE',
  'FORMAÇÃO SAVÂNICA',
  'OUTRA ÁREA NÃO VEGETADA',
  'CAMPO ALAGADO E ÁREA PANTANOSA'
)

## get classification regions
regions <- ee$FeatureCollection('users/dh-conciani/collection7/classification_regions/vector_v2')
regions_list <- sort(regions$aggregate_array('mapb')$getInfo())

## set years to be processed
years <- seq(from=1985, to=2018)

## set dictionary
classes <- ee$Dictionary(list(
  'LAVOURA TEMPORÁRIA'= 21,            
  'LAVOURA PERENE'= 21,           
  'CANA'= 21,      
  'PASTAGEM'= 21,       
  'FORMAÇÃO FLORESTAL'= 3,
  'RIO, LAGO E OCEANO'= 33,
  'FORMAÇÃO CAMPESTRE'= 12,
  'FORMAÇÃO SAVÂNICA'= 4,
  'OUTRA ÁREA NÃO VEGETADA'= 25,
  'CAMPO ALAGADO E ÁREA PANTANOSA'= 11
))

## for each file
for (i in 1:length(unique(file_name))) {
  print(paste0('processing file --> ', file_name[i]))
  ## read file [i]
  collection <- ee$Image(paste0(file_path, file_name[i]))
  
  ## set recipes
  recipe_metrics <- as.data.frame(NULL)
  recipe_table <- as.data.frame(NULL)
  
  ## for each year
  for(j in 1:length(unique(years))) {
    ## select year [j] and remap
    print(paste0('processing year -->', years[j]))
    
    collection_ij <- collection$select(paste0('classification_', years[j]))$
      remap(c(3, 4, 5, 11, 12, 29, 15, 19, 39, 20, 40, 41, 46, 47, 48, 21, 23, 24, 30, 25, 33, 31),
            c(3, 4, 3, 11, 12, 12, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 25, 25, 25, 25, 33, 33))$
      rename(paste0('classification_', years[j]))
    
    ## get validation points
    validation_ij <- validation_points$
      filterMetadata('POINTEDITE', 'not_equals', 'true')$
      filter(ee$Filter$inList(paste0('CLASS_', years[j]), selectClasses))$
      filterBounds(regions)$
      map(function(feature) {
        return(feature$set('year', years[j])$
                 set('reference', classes$get(feature$get(paste0('CLASS_', years[j])))))
      })
    
    ## for each region
    for(k in 1:length(unique(regions_list))) {
      print(paste0('processing region --> ', regions_list[k]))
      
      ## clip classification for the region
      classification_ijk <- collection_ij$clip(
        regions$filterMetadata('mapb', 'equals', regions_list[k]))
      
      ## clip val
      validation_ijk <- validation_ij$filterBounds(
        regions$filterMetadata('mapb', 'equals', regions_list[k])$geometry())
      
      ## extract classification value for each point and pair it with the reference data
      paired_data <- classification_ijk$sampleRegions(
        collection= validation_ijk, 
        scale= 30,
        geometries= FALSE)
      
      ## get map classes
      mapped <- paired_data$aggregate_array(paste0('classification_', years[j]))$getInfo()
      ref <- paired_data$aggregate_array('reference')$getInfo()
      
      ## transform lists into factors and merge them 
      toCompute <- as.data.frame(cbind(reference= mapped,
                                       predicted= ref))
      
      ## subset by considering classes that have reference points
      toCompute <- subset(toCompute, predicted %in% unique(toCompute$predicted)[
        which(unique(toCompute$predicted) %in% unique(toCompute$reference))
      ]
      )
      
      ## compute confusion matrix
      ## Check for classes
      unique_predicted <- unique(toCompute$predicted)
      unique_reference <- unique(toCompute$reference)
      
      ## Check for missing classes in the reference
      missing_classes <- unique_predicted[!(unique_predicted %in% unique_reference)]
      if (length(missing_classes) > 0) {
        warning("Classe presentes na variável 'predicted', mas não estão presentes na variável 'reference': ", paste(missing_classes, collapse = ", "))
      }
      
      ## Check for missing classes in the 'predicted' variable
      missing_classes <- unique_reference[!(unique_reference %in% unique_predicted)]
      if (length(missing_classes) > 0) {
        warning("Classe presentes na variável 'reference', mas não estão presentes na variável 'predicted': ", paste(missing_classes, collapse = ", "))
      }
      
      ## Filter only classes present in both variables
      toCompute <- subset(toCompute, predicted %in% unique_reference & reference %in% unique_predicted)
      
      confusion <- confusionMatrix(data = as.factor(toCompute$predicted),
                                   reference = as.factor(toCompute$reference))
      
      ## get metrics
      metrics <- rbind(melt(confusion$overall), 
                       melt(confusion$byClass[,11]))
      
      ## build results 
      ## insert variables name
      metrics$variable <- row.names(metrics)
      metrics$region <- unique(regions_list)[k]
      metrics$year <- years[j]
      metrics$file <- file_name[i]
      
      ## get confusion table
      confusionTable <- as.data.frame(confusion$table)
      confusionTable$region <- unique(regions_list)[k]
      confusionTable$year <- years[j]
      confusionTable$file <- file_name[i]
      
      ## bind data 
      recipe_metrics <- rbind(recipe_metrics, metrics)
      recipe_table <- rbind(recipe_table, confusionTable)
    }
  }
  
  ## save file results
  print('exporting results')
  
  # export path
  metrics_file <- file.path(output, paste0('metrics_', strsplit(file_name[i], "/")[[1]][length(strsplit(file_name[i], "/")[[1]])], '.csv'))
  table_file <- file.path(output, paste0('table_', strsplit(file_name[i], "/")[[1]][length(strsplit(file_name[i], "/")[[1]])], '.csv'))
  
  print('exporting results')
  write.csv(recipe_metrics, file = metrics_file)
  write.csv(recipe_table, file = table_file)
}

print('done, enjoy :)')
