\name{varbin}
\alias{varbin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{varbin
%%  ~~function to do ... ~~
}
\description{Optimal binning of numerical variable}
\usage{
varbin(df, x, y, p=0.05, custom_vec=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{A data frame}
  \item{x}{String. Name of continuous variable in data frame.}
  \item{y}{String. Name of binary response variable (0,1) in data frame.}
  \item{p}{Percentage of records per bin. Default 5 pct. (0.05). This parameter only accepts
values greater than 0.00 (0 pct.) and lower than 0.50 (50 pct.).}
  \item{custom_vec}{Numerical input vector with custom cutpoints. E.g. custom_vec=c(20, 50, 75) for a variable
  representing age, will result in the cutpoints [<20, <50, <75, >=75]. NA results in the default unrestricted (most optimal) binning.}
%%     ~~Describe \code{x} here~~

}
%%\details{
%%%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
The command varbin generates a data frame with necessary info and utilities for binning.
The user should save the output result so it can be used with e.g. varbin.plot,  or varbin.convert.
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

# Perform unrestricted binning
result <- varbin(df, "age", "target")

# Perform custom binning
result2 <- varbin(df, "age", "target", custom_vec=c(30,40,60,75))
}
