% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepLearnSets.R
\name{prepLearnSets}
\alias{prepLearnSets}
\title{A prepLearnSets Function}
\usage{
prepLearnSets(Y.tr, learnSetNames, compName, resultsDir, fold = 5,
  niter = 100, saveLearnSet = T)
}
\arguments{
\item{Y.tr}{vector de datos de clase "factor" donde se indican las condiciones experimentales (a comparar) de cada individuo}

\item{learnSetNames}{método a utilizar para la generación de "learningsSets". Los diferentes métodos (LOOCV, CV, MCCV, bootstrap) se explican en \code{\link{GenerateLearningsets}}}

\item{compName}{nombre de la comparación, se usara como identificador y para el nombre de archivos resultantes}

\item{resultsDir}{carpeta donde se guardaran los resultados. p.e. "results/". Solo necesario cuando saveLearnSet = TRUE}

\item{fold}{Gives the number of CV-groups. Used only when method="CV"}

\item{niter}{Number of iterations.}

\item{saveLearnSet}{valor logico que indica si se guardan los resultados en disco}
}
\value{
learningSets datos resultantes, "datos de entrenamiento". Objeto con clase "learningsets".

learningSetsFileName nombre del archivo guardado con los "datos de entrenamiento"
}
\description{
Genera datos de entrenamiento a partir de una base de datos mediante diferentes técnicas de cross-validation
}
\examples{
y <- factor(c(rep("A",10),rep("B",10)))
lSet <- prepLearnSets(Y.tr = y , learnSetNames = "LOOCV", compName = "AvsB", saveLearnSet = FALSE) 
}
\references{
M Slawski, M Daumer and A-L Boulesteix, CMA – a comprehensive Bioconductor package for supervised classification with high dimensional data
}
\seealso{
Text with \code{\link{GenerateLearningsets}}
}
\author{
Miriam Mota \email{mmota.foix@gmail.com}
}
\keyword{biomarcador}
\keyword{clasificador}
\keyword{cma}
\keyword{learningsets}
\keyword{predictor}
