\name{hSol}
\alias{hSol}

\title{Calculo anchos de banda}
\description{
Calcular mediante diferentes metodos anchos de banda para usar en el kernel
}
\usage{
hSol(data,kernel="gaus",plot=FALSE,method="sirverman",...)
}

\arguments{
  \item{data}{Observaciones de la variable con las que se estima la densidad kernel}
  \item{kernel}{Tipo de kernel para usar en la estimacion ("gaus","unif","Epkov","Triang","Bwgt","Trgt")}
  \item{plot}{Grafico para metodos "CV.MV" y "CV.LSCV"}
  \item{method}{Metodo para calcular el ancho de banda ("sirverman","Ref.Normal","SobreSuav","CV.LSCV",
  "CV.MV")}
}

\details{
  Detalles de los metodos usados en la funcion corespondiente
}
\value{
  Valores que retorna la funcion
}
\references{
  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.

  Fox, J. and Weisberg, S. (2019)
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.

  Atkinson, A. C. (1985)
  \emph{Plots, Transformations, and Regression.} Oxford.
}
\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[stats]{qqplot}}, \code{\link[stats]{qqnorm}},
  \code{\link[stats]{qqline}}, \code{\link{showLabels}}}

\examples{
  poner ejemplo de uso de funcion con datos reales
}

\keyword{distribution}
\keyword{univar}
\keyword{regression}

