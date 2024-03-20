## get accuracy for guyanas regions, in the scope of pan-amazonia mapbiomas project 
## --- --- --- EXPORT ACCURACY METRICS PER COUNTRY
## dhemerson.costa@ipam.org.br

## get libraries
library(rgee)
library(reticulate)
library(caret)
library(reshape2)

## initialize
ee_Initialize()

## set directories
file_path <- 'projects/mapbiomas-raisg/public/'

##  define files to be computed
file_name <- c(
  'collection5/mapbiomas_raisg_panamazonia_collection5_integration_v1',
  'collection4/mapbiomas_raisg_panamazonia_collection4_integration_v1',
  'collection3/mapbiomas_raisg_panamazonia_collection3_integration_v1',
  'collection2/mapbiomas_raisg_panamazonia_collection2_integration_v1'
)

## set output path (local)
output <- './table/'

## import validation points
validation_points = ee$FeatureCollection('users/vieiramesquita/MAPBIOMAS/mapbiomas_amazonia_50K_RAISG_plus_Brasil_v6')

## set classes to be validated 
selectClasses = c(
  'Formación Forestal',            
  'Formação Florestal',           
  'Formación Natural No Forestal Inundable',      
  'Otra Formación Natural No Forestal',       
  'Outra Formação Natural Não Florestal"',
  'Mosaico de Agricultura y/o Pasto',
  'Área sin Vegetación',
  'Río, Lago u Océano'
)

## set dictionary for each class
classes <- ee$Dictionary(list(
  'Formación Forestal'= 3,            
  'Formação Florestal'= 3,           
  'Formación Natural No Forestal Inundable'= 11,      
  'Otra Formación Natural No Forestal'= 12,       
  'Outra Formação Natural Não Florestal'= 12,
  'Mosaico de Agricultura y/o Pasto'= 21,
  'Área sin Vegetación'= 25,
  'Río, Lago u Océano'= 33
))

## get classification regions
regions <- ee$FeatureCollection('projects/mapbiomas-raisg/DATOS_AUXILIARES/VECTORES/paises-5')

## list regions to be computed (use name to filter)
regions_list <- c(
  'Guyane Française',
  'Suriname',
  'Guyana'
)

## set years to be processed (bcz collections have different years, i moved to inner loop)
#years <- seq(from=1985, to=2018)

## for each file
for (i in 1:length(unique(file_name))) {
  print(paste0('processing file --> ', file_name[i]))
  ## read file [i]
  collection <- ee$Image(paste0(file_path, file_name[i]))
  
  ## get years
  years <- as.numeric(sapply(strsplit(collection$bandNames()$getInfo(), "_"), `[`, 2))
  
  ## set recipes
  recipe_metrics <- as.data.frame(NULL)
  recipe_table <- as.data.frame(NULL)
  
  ## for each year
  for(j in 1:length(unique(years))) {
    ## select year [j] and remap
    print(paste0('processing year --> ', years[j]))
    
    ## reamp to match validation points 
    collection_ij <- collection$select(paste0('classification_', years[j]))$
      remap(c(3, 4, 5, 6, 11, 12, 29, 13, 15, 18, 9,  35, 21, 23, 24, 30, 25, 33, 34, 27),
            c(3, 3, 3, 3, 11, 12, 12, 12, 21, 21, 21, 21, 21, 25, 25, 25, 25, 33, 33, 27))$
      rename(paste0('classification_', years[j]))
    
    ## get validation points
    validation_ij <- validation_points$
      filterMetadata('POINTEDITE', 'not_equals', 'TRUE')$
      filter(ee$Filter$inList(paste0('CLASS_', years[j]), selectClasses))$
      map(function(feature) {
        return(feature$set('year', years[j])$
                 set('reference', classes$get(feature$get(paste0('CLASS_', years[j])))))
      })

        ## for each region
    for(k in 1:length(unique(regions_list))) {
      print(paste0('processing region --> ', regions_list[k]))
      
      ## clip classification for the region
      classification_ijk <- collection_ij$clip(
        regions$filterMetadata('name', 'equals', regions_list[k]))
      
      ## clip val
      validation_ijk <- validation_ij$filterBounds(
        regions$filterMetadata('name', 'equals', regions_list[k])$geometry())
      
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
        warning("classes occurs in the variable 'predicted' but not in 'reference': ", paste(missing_classes, collapse = ", "))
      }
      
      ## Check for missing classes in the 'predicted' variable
      missing_classes <- unique_reference[!(unique_reference %in% unique_predicted)]
      if (length(missing_classes) > 0) {
        warning("classes occurs in the variable 'reference', but not in 'predicted': ", paste(missing_classes, collapse = ", "))
      }
      
      ## Filter only classes present in both variables
      toCompute <- subset(toCompute, predicted %in% unique_reference & reference %in% unique_predicted)
      
      ## get confusion matrix 
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
