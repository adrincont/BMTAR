\name{tsregim}
\alias{tsregim}

\title{Creation of class ``tsregim'' for some data}
\description{
Creation and organization of data
}
\usage{
tsregim(Yt, Zt = NULL, Xt = NULL, r = NULL, ...)
autoplot(tsregim)
print(tsregim)
}

\arguments{
  \item{Yt}{matrix type object, observed process (admit NA values). Not NULL}
  \item{Zt}{matrix type object, threshold process (admit NA values). Default NULL}
  \item{Xt}{matrix type object, covariate process (admit NA values). Default NULL}
  \item{r}{numeric type, threshold value if known. Default NULL}
}

\details{
  Detalles de los metodos usados en la funcion corespondiente
}
\value{
  Return a list type object of class ``tsregim''
  $Yt observed process
  $Zt if enter
  $Xt if enter
  $N number of observations
  $k number of variables

  If r known
  $r threshold value
  $Ind numeric type, number of the regime each observation belongs
  $Summary_r data.frame type, number and proportion of observations in each regime

}
\references{
  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.
}
\author{Valeria Bejarano \email{vbejaranos@unal.edu.co}}
\author{Andrey Rincon \email{adrincont@unal.edu.co}}

\examples{
  poner ejemplo de uso de funcion con datos reales
}

\keyword{Observed process}

