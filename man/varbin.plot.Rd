\name{varbin.plot}
\alias{varbin.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{varbin.plot
%%  ~~function to do ... ~~
}
\description{Generate simple plot to visualize binning results}
\usage{
varbin.plot(ivTable)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ivTable}{Output from either varbin, varbin.factor, varbin.monotonic or varbin.kink.}
%%     ~~Describe \code{x} here~~

}
%%\details{
%%%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
The command varbin.plot generates a simple plot with the Weight of Evidence values on the y-axis and the cutpoints/binnings on the x-axis. Gives a nice overview of the functional form and the relatioship between the explanatory variable and the dependent variable.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%%\references{
%%%% ~put references to the literature/web site here ~
%%}
%%\author{
%%%%  ~~who you are~~
%%}
%%\note{
%%%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%%%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
# Set seed and generate data
set.seed(1337)
target <- as.numeric(runif(10000, 0, 1)<0.2)
age <- round(rnorm(10000, 40, 15), 0)
age[age<20] <- round(rnorm(sum(age<20), 40, 5), 0)
age[age>95] <- round(rnorm(sum(age>95), 40, 5), 0)
inc <- round(rnorm(10000, 100000, 10000), 0)
educ <- sample(c("MSC", "BSC", "SELF", "PHD", "OTHER"), 10000, replace=TRUE)
df <- data.frame(target=target, age=age, inc=inc, educ=educ)

# Perform restricted binning
result <- varbin.kink(df, "inc", "target")
# Plot result
varbin.plot(result)
}
