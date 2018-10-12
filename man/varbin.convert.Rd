\name{varbin.convert}
\alias{varbin.convert}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{varbin.convert
%%  ~~function to do ... ~~
}
\description{Generate new variable based on constructed binnings}
\usage{
varbin.convert(df, ivTable, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{A data frame}
  \item{ivTable}{Output from either varbin, varbin.factor, varbin.monotonic or varbin.kink.}
  \item{x}{String. Name of variable in data frame for which binninngs should be applied.}
%%     ~~Describe \code{x} here~~

}
%%\details{
%%%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
The command varbin.convert appends a new variable named "WoE_[x]" to the data frame. The new variable consist of the Weight of Evidence values from the resulting binning.
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

# Split train/test
df_train <- df[1:5000, ]
df_test <- df[5001:nrow(df), ]

# Perform unrestricted binnings
result <- varbin.factor(df_train, "educ", "target")
result2 <- varbin(df_train, "age", "target")

# Convert test data
df_new <- varbin.convert(rbind(df_train, df_test), result,"educ")
df_new <- varbin.convert(df_new, result2, "age")
}
