\dontrun{

##############################################################################################
## Boundary value problem: Ozon formation with fixed ozon/oxygen ratio at 
## final time point
##############################################################################################

library(bvpSolve)

# O2 + O <-> O3
# diff = O2 - O3
# build_O3 = const.
f <- c(
  O3 = " build_O3 * O2 * O - decay_O3 * O3",
  O2 = "-build_O3 * O2 * O + decay_O3 * O3",
  O  = "-build_O3 * O2 * O + decay_O3 * O3",
  diff = "-2 * build_O3 * O2 * O + 2 * decay_O3 * O3", 
  build_O3 = "0"
)

bound <- data.frame(
    name = names(f),
    yini = c(0, 3, 2, 3, NA),
    yend = c(NA, NA, NA, 0, NA)
)

# Generate ODE function
func <- funC(f, jacobian="full", boundary = bound)

# Initialize times, states, parameters and forcings
times <- seq(0, 15, by = .1)
pars <- c(decay_O3 = .1)
xguess <- times
yguess <- matrix(c(1, 1, 1, 1, 1), ncol=length(times), nrow = length(f))

# Solve BVP
out <- bvptwpC(x = times, func = func, parms = pars, xguess = xguess, yguess = yguess)

# Solve BVP for different ini values, end values and parameters
yini <- c(O3 = 2)
yend <- c(diff = 0.2)
pars <- c(decay_O3 = .01)
out <- bvptwpC(yini = yini, yend = yend, x = times, func = func, 
	       parms = pars, xguess = xguess, yguess = yguess)

# Plot solution
par(mfcol=c(1,2))
t <- out[,1]
M1 <- out[,2:5]
M2 <- cbind(out[,6], pars)

matplot(t, M1, type="l", lty=1, col=1:4, xlab="time", ylab="value", main="states")
legend("topright", legend = c("O3", "O2", "O", "O2 - O3"), lty=1, col=1:4)
matplot(t, M2, type="l", lty=1, col=1:2, xlab="time", ylab="value", main="parameters")
legend("right", legend = c("build_O3", "decay_O3"), lty=1, col=1:2)

}
