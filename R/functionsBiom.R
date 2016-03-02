#' A prepLearnSets Function
#'
#' Genera datos de entrenamiento a partir de una base de datos mediante diferentes técnicas de cross-validation
#' @param Y.tr vector de datos de clase "factor" donde se indican las condiciones experimentales (a comparar) de cada individuo
#' @param learnSetNames método a utilizar para la generación de "learningsSets". Los diferentes métodos (LOOCV, CV, MCCV, bootstrap) se explican en \code{\link{GenerateLearningsets}} 
#' @param compName nombre de la comparación, se usara como identificador y para el nombre de archivos resultantes
#' @param resultsDir carpeta donde se guardaran los resultados. p.e. "results/". Solo necesario cuando saveLearnSet = TRUE
#' @param fold Gives the number of CV-groups. Used only when method="CV"
#' @param niter Number of iterations (s.details).
#' @param saveLearnSet valor logico que indica si se guardan los resultados en disco
#' @export prepLearnSets
#' @import CMA
#' @author Miriam Mota <mmota.foix@@gmail.com>
#' @seealso Text with \code{\link{GenerateLearningsets}} 
#' @examples
#' y <- factor(c(rep("A",10),rep("B",10)))
#' lSet <- prepLearnSets(Y.tr = y , learnSetNames = "LOOCV", compName = "AvsB", saveLearnSet = F) 
#' @return learningSets datos resultantes, "datos de entrenamiento". Objeto con clase "learningsets".
#' @return learningSetsFileName nombre del archivo guardado con los "datos de entrenamiento"
#' @keywords cma predictor biomarcador clasificador learningsets
#' @references M Slawski, M Daumer and A-L Boulesteix, CMA – a comprehensive Bioconductor package for supervised classification with high dimensional data

prepLearnSets <- function(Y.tr, learnSetNames, compName,resultsDir,  fold = 5 , niter = 100, saveLearnSet = T){
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



#' A createClassif Function
#'
#' Construcción y validación externa de un biomarcador.
#' @param X.tr matriz númerica con los valores de los datos. Donde las columnas son los "features" y las filas los individuos. Estos datos serán usados para la construcción del biomarcador
#' @param Y.tr vector de datos de clase "factor" donde se indican las condiciones experimentales (a comparar) de cada individuo indicado en las filas de X.tr
#' @param learningSets An object of class learningsets. May be missing, then the complete datasets is used as learning set.
#' @param learnSetNames método a utilizar para la generación de "learningsSets". Los diferentes métodos (LOOCV, CV, MCCV, bootstrap) se explican en "See Also"
#' @param selMethodNames vector caracter indicando métodos de seleccion (máximo 3?), ver posibles métodos en \code{\link{GeneSelection}} . (p.e: c("t.test"))
#' @param numGenes2Sel vector númeric donde se indican el número de genes incluidos en cada predictor.
#' @param classifierNames vector caracter indicando el método de classificador usado .(p.ej. c("CMA::dldaCMA", "CMA::rfCMA", "CMA::pnnCMA", "CMA::plrCMA") ), parámetro classifier \code{\link{GeneSelection}}
#' @param resultsDir carpeta donde se guardaran los resultados. p.e. "results/". Solo necesario cuando saveLearnSet = TRUE
#' @param compName nombre de la comparación, se usara como identificador y para el nombre de archivos resultantes
#' @param fold Gives the number of CV-groups. Used only when method="CV"
#' @param niter Number of iterations (s.details).
#' @param ntoplist número de features que se muestran en la lista de candidatos
#' @param X.new  Solo si validation = TRUE. matriz númerica con los valores de los datos. Donde las columnas son los "features" y las filas los individuos. Estos datos serán usados para la validación del biomarcador
#' @param Y.new Solo si validation = TRUE. vector de datos de clase "factor" donde se indican las condiciones experimentales (a comparar) de cada individuo indicado en las filas de X.new
#' @param validation valor lógico que indica si se incluye una base de datos de validación.
#' @param cond3 valor lógico que indica si se la comparación, es decir, los factores son igual a 3. 
#' @param saveResults valor logico que indica si se guardan los resultados en disco
#' @keywords cma predictor biomarcador clasificador learningsets
#' @export createClassif
#' @import CMA xlsx
#' @author Miriam Mota <mmota.foix@@gmail.com>
#' @seealso Text with \code{\link{GeneSelection}} 
#' @examples
#'  x <- cbind(matrix(rnorm(400,500),nrow = 40),matrix(rnorm(400,5),nrow = 40))
#'  colnames(x) <- paste0("a",1:ncol(x))
#'  rownames(x) <- paste0("aa",1:nrow(x))
#'  y <- factor(c(rep("A",20),rep("B",20)))
#'  lSet <- prepLearnSets(Y.tr = y , learnSetNames = "LOOCV", compName = "AvsB", saveLearnSet = F)  
#'  resF <- createClassif(X.tr = x,
#'                        Y.tr = y,
#'                       learningSets  = lSet$learningSets,
#'                       learnSetNames = c("LOOCV"),
#'                       selMethodNames = c("t.test"),
#'                       numGenes2Sel = c(3,5,10),
#'                       classifierNames = c("CMA::dldaCMA", "CMA::rfCMA", "CMA::pnnCMA") , 
#'                       resultsDir = "example/",
#'                       compName = "AvsB",
#'                       validation=FALSE,
#'                       ntoplist = 10)



createClassif <- function(X.tr, Y.tr, learningSets, learnSetNames, selMethodNames, numGenes2Sel = c(3,5,10), classifierNames,
                          cond3=FALSE, isTunable, resultsDir, compName, niter,ntoplist = 25, X.new, Y.new, validation = TRUE,
                          saveResults = TRUE)
{
  classifs <- list();   geneSels <- list() ;   results <- list(); misscls <- list()
  for (i in 1:length(learningSets)) {
    for (j in 1:length(selMethodNames)) {
      print(paste("sel method: ", selMethodNames[j], sep = ""))
      selected  <- CMA::GeneSelection(X.tr, Y.tr, learningsets = learningSets[[i]], method = selMethodNames[j])
      itemName <- paste(learnSetNames[i], selMethodNames[j], sep = ".")
      geneSels[[itemName]] <- selected
      for (numGenes in numGenes2Sel) { 
        print(paste("num genes: ", numGenes, sep = ""))
        for (k in 1:length(classifierNames)) {
          print(paste("classif: ", classifierNames[k]), sep = "")
          myClassifier <- eval(parse(text = classifierNames[k]))
          # if (isTunable[k]) {     
          # tuneVals <- CMA::tune(X = X.tr, y = Y.tr, learningsets = learningSets[[i]],
          # genesel = selected, nbgene = numGenes, 
          # classifier = myClassifier, grids = list( mtry = 2, nodesize = 3))
          # classif <- CMA::classification(X = X.tr, y = Y.tr, learningsets = learningSets[[i]], 
          # genesel = selected, nbgene = numGenes, 
          # classifier = myClassifier, 
          # tuneres = tuneVals)#, scheme = "observationwise")
          
          # }else{
          classif <- CMA::classification(X = X.tr, y = Y.tr, learningsets = learningSets[[i]], 
                                         genesel = selected, nbgene = numGenes, 
                                         classifier = myClassifier)#, scheme = "observationwise")
          itemName <- paste(learnSetNames[i], selMethodNames[j], numGenes, classifierNames[k], sep = ".")
          if (validation) {
            pred <- CMA::prediction(X.tr = X.tr, 
                                    y.tr = Y.tr,
                                    X.new = X.new,
                                    classifier = myClassifier,
                                    genesel = selected,
                                    nbgene = numGenes)
            predvalues <- factor(show(pred),c(levels(Y.tr)))
            Y.new <- factor(Y.new, c(levels(Y.tr)))
            tabpred <- table(predvalues, Y.new)
            misscl <- 1 - (sum(diag(tabpred))/sum(tabpred))
            misscls[[itemName]] <- misscl
          }
          classifs[[itemName]] <- classif
        }
      }
    }
  }
  
  ## save results features seleccionats i classificacio
  selectedGenesFileName <- paste0(compName, "_selectedGenes_",if (learnSetNames != "LOOCV") {paste0(niter, "iter")}, ".Rda")
  save(geneSels, file = file.path(resultsDir, selectedGenesFileName))
  classifsFileName <- paste0(compName, "_classifs_", if (learnSetNames != "LOOCV") {paste0(niter, "iter")}, ".Rda")
  save(classifs, file = file.path(resultsDir, classifsFileName))
  
  # corbes ROC
  resClass <- lapply(classifs, join)
  
  if (cond3) {
    if (!dir.exists(paste0(resultsDir,"/contrastsMatrix"))) dir.create(paste0(resultsDir,"/contrastsMatrix"))
    missclsINT <- NULL
    for (i in 1:length(resClass)) {
      txt_name <- file.path(resultsDir,paste0("contrastsMatrix/",compName,names(resClass)[[i]],"contrastMatrix.txt"))  
      sink(txt_name)
      ftable(resClass[[i]])
      sink()
      ftab.comp <- readLines(txt_name)
      missclsINT[i] <- as.numeric(strsplit(ftab.comp,":",fixed = T)[[2]][2])
      names(missclsINT)[i] <- names(resClass)[[i]]
    }
  }
  
  if (!cond3) { # no funciona be, cal mirar-ho
    pdf(file.path(resultsDir,paste0(compName, "_ROC", if (learnSetNames != "LOOCV") {paste0(niter, "iter")}, ".pdf")),onefile = T)
    par(mfrow = c(2,2))
    for (i in seq(along = resClass)) try(roc(eval(resClass[[i]]), main = names(resClass)[[i]]),TRUE)
    dev.off()  
  }
  
  # features seleccionats
  topLists <- lapply(geneSels, toplist, ntoplist)
  res25O <- as.data.frame(topLists)
  results[["table_all"]] <- x <- res25O
  
  genIdxs <- x[,1]
  geneNames <- colnames(X.tr)
  myGeneNames <- geneNames[as.integer(genIdxs)]
  results[["selectedTable"]] <- data.frame(feature = myGeneNames, IndexImportance = x[,2])#, unlist(misscls))
  
  
  # guardem features candidats
  biomarkersFileName <- paste0(compName, "Results",if (learnSetNames != "LOOCV")  {paste0(niter, "iter")} , ".xls")
  if(saveResults) write.xlsx(results[["selectedTable"]], file = file.path(resultsDir, biomarkersFileName), row.names = FALSE , sheetName = "candidateBiomarkers")
  
  
  ## COMPARE RESULTS 
  compMeasures <-  c("misclassification", "sensitivity", "specificity")
  s1 <- c(rep(learnSetNames[1], length(selMethodNames)*length(numGenes2Sel)*length(classifierNames)))
  s2 <- c(rep(selMethodNames[1], length(numGenes2Sel)*length(classifierNames)))
  s3 <- rep(c(rep(numGenes2Sel[1], length(classifierNames)), 
              rep(numGenes2Sel[2], length(classifierNames)), 
              rep(numGenes2Sel[3], length(classifierNames))), length(selMethodNames))
  s4 <- rep(classifierNames, length(selMethodNames)*length(numGenes2Sel))
  s <- paste(s1, s2, s3, s4, sep = ".")
  
  
  if (!cond3) {
    compClassifs <- CMA::compare(classifs, measure = compMeasures)
    resultsClassif <- data.frame(CrossVal = s1, VarSel = s2, numGenes = s3, Classif = s4, compClassifs)
  }
  else{
    compClassifs <- missclsINT
    resultsClassif <- data.frame(CrossVal = s1, VarSel = s2, numGenes = s3, Classif = s4, misclassification = compClassifs)
  }
  
  ## SAVE COMPARED RESULTS 
  if(saveResults) 
  {
    write.xlsx(resultsClassif, 
               file = file.path(resultsDir,
                                paste0(compName, "Results", if (learnSetNames != "LOOCV")
                                {paste0(niter, "iter")}, ".xls")),sheetName = "classif",append = TRUE)
  }
  if (validation) {
    resValid <- cbind(resultsClassif,misclassifTEST = unlist(misscls))
    if(saveResults) 
    {
    write.xlsx(resValid, file = file.path(resultsDir,
                                          paste0(compName, "Results", if (learnSetNames != "LOOCV")
                                          {paste0(niter, "iter")}, ".xls")),sheetName = "validation", append = TRUE)
    }
  }
  
  if (validation) {
    return(list( candFeat = results, cl = resClass, ResultsClassif = resultsClassif,misscls = misscls ))
  }
  else{
    return(list( candFeat = results, cl = resClass, ResultsClassif = resultsClassif))
  }
}




