\name{KerEs}
\alias{KerEs}

\title{Calculo estimación kernel}
\description{
Calcula y grafica la estimación por método kernel para una variable tipo continua univariada especificando tipo de kernel, ancho de banda. Incluye métodos kernel-acotado y kernel con ventana variable.
}
\usage{
KerES(data,kernel="gaus",h=NULL,tipo="usual",point=NULL,alphaB=0.05,alphaKV=0.5,...)
summary(KerES)
autoplot(KerES)
}

\arguments{
  \item{data}{Observaciones de la variable con las que se estima la densidad kernel}
  \item{kernel}{Tipo de kernel para usar en la estimacion ("gaus","unif","Epkov","Triang","Bwgt","Trgt")}
  \item{h}{Ancho de banda seleccionado para la estimacion kernel}
  \item{tipo}{Tipo de kernel seleccionado ("usual","variable","acotado")}
  \item{point}{Puntos donde se calcula la estimacion kernel por defecto}
  \item{alphaB}{Nivel de confianza para las bandas estimadas}
  \item{alphaKV}{Alpha correspondiente a el kernel con ancho de banda variable}
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
