#' A prepComp Function
#'
#' Prepara los datos de un comparación de m niveles para usar en las funciones prepLearnSets i createClassif.
#' @param BD data.frame donde las filas son las muestras y las columnas los features además de una columna con el grupo o condición experimental 
#' @param cond vector de strings que indique las diferentes condiciones experimentales en el orden que se desee. (p.e. c("bb","a"))
#' @param Y nombre de la columna en formato string donde se indica la condición experimental
#' @param appLog valor lógico que indica si se aplica logaritmo en base 2 a los datos.
#' @export prepComp
#' @import stringr
#' @author Miriam Mota \email{mmota.foix@@gmail.com}
#' @examples
#' dat.prov <- data.frame(matrix(rnorm(1000),ncol=20))
#' colnames(dat.prov) <- paste0("gen",1:dim(dat.prov)[2])
#' rownames(dat.prov) <- paste0("sample",1:dim(dat.prov)[1])
#' dim(dat.prov)
#' dat.prov$condition <- c(rep("at.3",20),rep("b",20),rep("cc.4",10))
#' BvsA <- prepComp(BD = dat.prov, cond = c("b","a"), Y = "condition")
#' @return data.comp: data.frame original, donde las filas son las muestras (de las condiciones experimentales indicadas) y las columnas los "features". IMPORTANTE: no son necesariamente en formato númerico. 
#' @return X: matriz numérica, donde las filas son las muestras (de las condiciones experimentales indicadas) y las columnas los "features" 
#' @return Y: vector "factor", donde se indica la condición experimental correspondiente a cada una de las muestras
#' @keywords prepare comparison



prepComp <- function(BD, cond, Y, appLog = FALSE){
  str.search <- substr(paste0(cond,"|", collapse = ""), 1, nchar(paste0(cond,"|", collapse = ""))-1) 
  data.comp <- BD[grep(str.search,BD[,Y]),]
  condY <- data.comp[,Y]
  data.comp <- data.comp[,-which(names(data.comp) == Y)]
  
  #factor al que correspon cada mostra
  Y.tr <- as.factor(str_extract(condY,str.search))
  Y.tr <- factor(Y.tr, cond)
  
  #matriu dades sense condicio experimental
  X.tr <- as.matrix(data.comp) # comprobar que sigui una matriu numerica
  colX <- rownames(X.tr) 
  X.tr <- apply(X.tr,2, function(x) as.numeric(as.character(x) ) )
  rownames(X.tr)  <- colX
  if (appLog) {X.tr <- log2(X.tr)}
  return(list(data.comp = data.comp, X = X.tr, Y = Y.tr ))
}