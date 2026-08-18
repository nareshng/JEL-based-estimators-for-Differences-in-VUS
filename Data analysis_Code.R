

## =========================
## Real Data Analysis Code
## =========================

library(readxl)
library(emplik)

start.time <- Sys.time()

## ---- Load data
setwd("/Users/")
dat <- read_excel("Data.xlsx", sheet="Sheet1")



s <- split(dat, dat$DX)
s1 <- s$CN
s2 <- s$LMCI
s3 <- s$AD

## ---- Build samples (2D)
x <- cbind(ABETA1 = s1$ABETA1, PTAU = s1$PTAU)
y <- cbind(ABETA1 = s2$ABETA1, PTAU = s2$PTAU)
z <- cbind(ABETA1 = s3$ABETA1, PTAU = s3$PTAU)

n1 <- nrow(x); n2 <- nrow(y); n3 <- nrow(z)
n  <- n1 + n2 + n3

stopifnot(n1 > 1, n2 > 1, n3 > 1)  # jackknife needs this


## ==========================================================
## 1) Fast computation of S = sum_j L(y_j)*G(y_j)
##    where L(y)=#(x<y), G(y)=#(z>y)   (strict inequalities)
## ==========================================================
order_stat_aux <- function(xv, yv, zv) {
  xs <- sort(xv)
  zs <- sort(zv)
  n3 <- length(zv)
  
  # L(y) = #(x < y)  (strict)
  L <- findInterval(yv, xs, left.open = TRUE)
  
  # G(y) = #(z > y)  (strict) = n3 - #(z <= y)
  G <- n3 - findInterval(yv, zs, rightmost.closed = TRUE)
  
  term <- L * G
  S <- sum(term)
  
  # For jackknife updates we want y sorted, with L and G aligned.
  o <- order(yv)
  y_sorted <- yv[o]
  Ls <- L[o]
  Gs <- G[o]
  
  # prefix sums of L (for fast sum_{y<z} L(y))
  prefL <- c(0L, cumsum(Ls))
  
  # suffix sums of G (for fast sum_{y>x} G(y))
  suffG <- c(rev(cumsum(rev(Gs))), 0L)
  
  list(
    S = S,
    term = term,
    y_sorted = y_sorted,
    L_sorted = Ls,
    G_sorted = Gs,
    prefL = prefL,
    suffG = suffG
  )
}

## Fast helper: sum_{y > t} G(y) using y_sorted + suffG
sumG_gt <- function(t, y_sorted, suffG) {
  # idx = #(y <= t)
  idx <- findInterval(t, y_sorted, rightmost.closed = TRUE)
  # y > t are positions (idx+1 .. n2)
  suffG[idx + 1L]
}

## Fast helper: sum_{y < t} L(y) using y_sorted + prefL
sumL_lt <- function(t, y_sorted, prefL) {
  # idx = #(y < t)
  idx <- findInterval(t, y_sorted, left.open = TRUE)
  prefL[idx + 1L]
}


## ---- Compute theta_hat = P(x1<y1<z1) - P(x2<y2<z2)
aux1 <- order_stat_aux(x[,1], y[,1], z[,1])
aux2 <- order_stat_aux(x[,2], y[,2], z[,2])

theta_hat <- (aux1$S - aux2$S) / (as.numeric(n1) * as.numeric(n2)* as.numeric(n3))
theta_hat


## ==========================================================
## 2) jackknife theta_{-l} for all observations
## ==========================================================
theta_minus <- numeric(n)

## ---- Leave-one-out from group X (size n1)
den_x <- (n1 - 1) * n2 * n3
for (i in seq_len(n1)) {
  S1_i <- aux1$S - sumG_gt(x[i,1], aux1$y_sorted, aux1$suffG)
  S2_i <- aux2$S - sumG_gt(x[i,2], aux2$y_sorted, aux2$suffG)
  theta_minus[i] <- (S1_i - S2_i) / den_x
}

## ---- Leave-one-out from group Y (size n2)
den_y <- n1 * (n2 - 1) * n3
for (j in seq_len(n2)) {
  S1_j <- aux1$S - aux1$term[j]
  S2_j <- aux2$S - aux2$term[j]
  theta_minus[n1 + j] <- (S1_j - S2_j) / den_y
}

## ---- Leave-one-out from group Z (size n3)
den_z <- n1 * n2 * (n3 - 1)
for (k in seq_len(n3)) {
  S1_k <- aux1$S - sumL_lt(z[k,1], aux1$y_sorted, aux1$prefL)
  S2_k <- aux2$S - sumL_lt(z[k,2], aux2$y_sorted, aux2$prefL)
  theta_minus[n1 + n2 + k] <- (S1_k - S2_k) / den_z
}

## ---- Standard jackknife pseudo-values for JEL
v <- n * theta_hat - (n - 1) * theta_minus
mean(v)  # should be close to theta_hat (jackknife identity)


## ==========================================================
## 3) JEL confidence intervals (emplik)
## ==========================================================
myfun6 <- function(theta, v) el.test(v, mu = theta)

el.test(v, mu = theta_hat)

CI90 <- findUL(step = 0.01, initStep = 0, fun = myfun6, v = v, MLE = theta_hat, level = 2.705543)
CI95 <- findUL(step = 0.01, initStep = 0, fun = myfun6, v = v, MLE = theta_hat, level = 3.84146)
CI99 <- findUL(step = 0.01, initStep = 0, fun = myfun6, v = v, MLE = theta_hat, level = 5.023886)

CI90[c("Low","Up")]
CI95[c("Low","Up")]
CI99[c("Low","Up")]


## ==========================================================
## 4) Normal CI variance estimate 
## ==========================================================
v_x <- n1 * theta_hat - (n1 - 1) * theta_minus[1:n1]
v_y <- n2 * theta_hat - (n2 - 1) * theta_minus[(n1 + 1):(n1 + n2)]
v_z <- n3 * theta_hat - (n3 - 1) * theta_minus[(n1 + n2 + 1):n]

Var_hat <- var(v_x)/n1 + var(v_y)/n2 + var(v_z)/n3
Var_hat

c(
  "90% L" = theta_hat - 1.645 * sqrt(Var_hat),
  "90% U" = theta_hat + 1.645 * sqrt(Var_hat),
  "95% L" = theta_hat - 1.96  * sqrt(Var_hat),
  "95% U" = theta_hat + 1.96  * sqrt(Var_hat),
  "99% L" = theta_hat - 2.576 * sqrt(Var_hat),
  "99% U" = theta_hat + 2.576 * sqrt(Var_hat)
)


## ==========================================================
## 5) Bootstrap kernel-based method
## ==========================================================
silverman_h <- function(v) {
  0.9 * min(sd(v), IQR(v)/1.34) * length(v)^(-0.2)
}

kernel_theta_fast <- function(x, y, z) {
  n1 <- nrow(x); n2 <- nrow(y); n3 <- nrow(z)
  
  h1 <- silverman_h(x[,1]); h2 <- silverman_h(y[,1]); h3 <- silverman_h(z[,1])
  h4 <- silverman_h(x[,2]); h5 <- silverman_h(y[,2]); h6 <- silverman_h(z[,2])
  
  s12_1 <- sqrt(h1^2 + h2^2)
  s23_1 <- sqrt(h2^2 + h3^2)
  s12_2 <- sqrt(h4^2 + h5^2)
  s23_2 <- sqrt(h5^2 + h6^2)
  
  # For each y_j:
  # A_j = mean_i Phi((y_j - x_i)/s12)
  # B_j = mean_k Phi((z_k - y_j)/s23)
  A1 <- rowMeans(pnorm( outer(y[,1], x[,1], "-") / s12_1 ))
  B1 <- colMeans(pnorm( outer(z[,1], y[,1], "-") / s23_1 ))
  contrib1 <- mean(A1 * B1)
  
  A2 <- rowMeans(pnorm( outer(y[,2], x[,2], "-") / s12_2 ))
  B2 <- colMeans(pnorm( outer(z[,2], y[,2], "-") / s23_2 ))
  contrib2 <- mean(A2 * B2)
  
  contrib1 - contrib2
}

set.seed(9)
B <- 1000
kstat <- numeric(B)

for (b in seq_len(B)) {
  xb <- x[sample.int(n1, n1, replace = TRUE), , drop = FALSE]
  yb <- y[sample.int(n2, n2, replace = TRUE), , drop = FALSE]
  zb <- z[sample.int(n3, n3, replace = TRUE), , drop = FALSE]
  kstat[b] <- kernel_theta_fast(xb, yb, zb)
}

quantile(kstat, c(0.05, 0.95))    # 90%
quantile(kstat, c(0.025, 0.975))  # 95%
quantile(kstat, c(0.005, 0.995))  # 99%


end.time <- Sys.time()
end.time - start.time
