#' A prepLearnSets Function
#'
#' Genera datos de entrenamiento a partir de una base de datos mediante diferentes técnicas de cross-validation
#' @param Y.tr vector de datos de clase "factor" donde se indican las condiciones experimentales (a comparar) de cada individuo
#' @param learnSetNames método a utilizar para la generación de "learningsSets". Los diferentes métodos (LOOCV, CV, MCCV, bootstrap) se explican en \code{\link{GenerateLearningsets}} 
#' @param compName nombre de la comparación, se usara como identificador y para el nombre de archivos resultantes
#' @param resultsDir carpeta donde se guardaran los resultados. p.e. "results/". Solo necesario cuando saveLearnSet = TRUE
#' @param fold Gives the number of CV-groups. Used only when method="CV"
#' @param niter Number of iterations.
#' @param saveLearnSet valor logico que indica si se guardan los resultados en disco
#' @export prepLearnSets
#' @import CMA
#' @author Miriam Mota \email{mmota.foix@@gmail.com}
#' @seealso Text with \code{\link{GenerateLearningsets}} 
#' @examples
#' y <- factor(c(rep("A",10),rep("B",10)))
#' lSet <- prepLearnSets(Y.tr = y , learnSetNames = "LOOCV", compName = "AvsB", saveLearnSet = FALSE) 
#' @return learningSets datos resultantes, "datos de entrenamiento". Objeto con clase "learningsets".
#' @return learningSetsFileName nombre del archivo guardado con los "datos de entrenamiento"
#' @keywords cma predictor biomarcador clasificador learningsets
#' @references M Slawski, M Daumer and A-L Boulesteix, CMA – a comprehensive Bioconductor package for supervised classification with high dimensional data

prepLearnSets <- function(Y.tr, learnSetNames, compName,resultsDir,  fold = 5 , niter = 100, saveLearnSet = T){
  #if (!grepl("LOOCV|bootstrap|CV|MCCV",learnSetNames))  
  #  cat("Error in match.arg(learnSetNames, c('LOOCV', 'CV', 'MCCV', 'bootstrap')): \n'arg' should be one of 'LOOCV', 'CV', 'MCCV', 'bootstrap'")
  if(saveLearnSet) if (!dir.exists(resultsDir)) dir.create(resultsDir)
  if (learnSetNames == "LOOCV") loo <- GenerateLearningsets(y = Y.tr, method = learnSetNames)
  if (learnSetNames == "bootstrap") loo <- GenerateLearningsets(y = Y.tr, method = learnSetNames,
                                                                niter = niter,ntrain = floor(2/3*length(Y.tr)))
  if (learnSetNames == "CV") loo <- GenerateLearningsets(y = Y.tr, method = learnSetNames,fold = fold, niter = niter)
  if (learnSetNames == "MCCV") loo <- GenerateLearningsets(y = Y.tr, method = learnSetNames,
                                                           niter = niter,ntrain = floor(2/3*length(Y.tr)))
  learningSets <- list(oneFold = loo)
  learningSetsFileName <- paste0(compName, "_learningSets_", if (learnSetNames != "LOOCV") {paste0(niter, "iter")},".Rda")
  if (saveLearnSet) {
    save(learningSets, file = file.path(resultsDir, learningSetsFileName))
  }
  return(list(learningSets = learningSets,learningSetsFileName = learningSetsFileName))
}








