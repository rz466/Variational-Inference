logit <- function(x) {
	log(x) - log(1 - x)
}

logistic <- function(x) {
	1 / (1 + exp(-x))
}

log_sum_exp <- function(x) {
	x.max <- max(x);
	log(sum(exp(x - x.max))) + x.max
}
