# R6-6 (ADDENDUM_R6_6.md md5 c21147e3, POST HOC): ancestral constraint (cross-clade, eggNOG-OG-bridged background
# omega) x deep-specific relaxation. Run from ~/deepsea_conv/work. Read-only on R6_3 / R6_2b outputs.
set.seed(20261105)
CL <- c("Actinopteri_m","Annelida","Arthropoda","Chondrichthyes","Cnidaria","Echinodermata","Mollusca","Porifera")
OUT <- "addon/R6/R6_6"; dir.create(OUT, showWarnings = FALSE)
rt <- function(p) read.delim(p, colClasses = "character")
F <- rt("addon/R6/R6_3/r6_3_fits.tsv"); F <- F[F$bph_ok == "1", ]
for (v in c("relax_p","relax_K","wpur_T","wpur_B")) F[[v]] <- as.numeric(F[[v]])
A <- do.call(rbind, lapply(CL, function(c) { a <- rt(sprintf("func_v2/%s/hog_annot.tsv", c)); a <- a[a$og != "", ]
  data.frame(clade = c, hog = a$hog, og = a$og) }))
F <- merge(F, A, by = c("clade","hog"), all.x = TRUE)
lv <- function(w) log10(w + .001)
# donor table: real fits, mean per (og, clade)
R <- F[F$label == "real" & !is.na(F$og), ]
D <- aggregate(list(v = lv(R$wpur_B)), list(og = R$og, clade = R$clade), mean)
Sraw <- function(clade, og) mapply(function(c, o) { x <- D$v[D$og == o & D$clade != c]; if (length(x)) -mean(x) else NA },
                                   clade, og)
z <- function(x) (x - mean(x)) / sd(x)
## ---------------- primary ----------------
P <- F[!is.na(F$og), ]
pair <- intersect(paste(P$clade, P$hog)[P$label == "real"], paste(P$clade, P$hog)[P$label == "null1"])
P <- P[paste(P$clade, P$hog) %in% pair, ]
U <- unique(P[, c("clade","hog","og")]); U$S0 <- Sraw(U$clade, U$og); U <- U[!is.na(U$S0), ]; U$S <- z(U$S0)
P <- merge(P, U[, c("clade","hog","S0","S")], by = c("clade","hog"))
P$relaxed <- as.integer(P$relax_p < .05 & P$relax_K < 1); P$real <- as.integer(P$label == "real")
nNA <- sum(is.na(P$relaxed)); P <- P[!is.na(P$relaxed), ]
P <- P[order(P$clade, P$hog, P$label), ]
cat("R6-6 primary: HOGs", length(unique(paste(P$clade, P$hog))), " rows", nrow(P), " dropped NA", nNA, "\n")
X0 <- model.matrix(~ S * real + factor(clade), P); ji <- which(colnames(X0) == "S:real"); js <- which(colnames(X0) == "S")
bint <- function(S, X = X0, y = P$relaxed) { X[, js] <- S; X[, ji] <- S * P$real
  coef(glm.fit(X, y, family = binomial()))[ji] }
key <- paste(P$clade, P$hog); uk <- unique(key); uc <- sub(" .*", "", uk); Su <- P$S[match(uk, key)]
b0 <- bint(P$S)
permS <- function() { s <- Su; for (c in CL) { i <- which(uc == c); s[i] <- s[i][sample.int(length(i))] }; s[match(key, uk)] }
bp <- replicate(10000, bint(permS()))
p_prim <- (1 + sum(abs(bp) >= abs(b0))) / 10001
idxl <- split(seq_along(key), key)
bb <- replicate(2000, { k <- unlist(lapply(CL, function(c) { i <- uk[uc == c]; sample(i, length(i), TRUE) }))
  idx <- unlist(idxl[k], use.names = FALSE); Q <- P[idx, ]
  X <- model.matrix(~ S * real + factor(clade), Q); coef(glm.fit(X, Q$relaxed, family = binomial()))[ji] })
ci <- quantile(bb, c(.025, .975), na.rm = TRUE)
cat(sprintf("beta_int (S:real, logOR/SD) = %.4f  perm P(2s) = %.4f  boot95 %.4f..%.4f\n", b0, p_prim, ci[1], ci[2]))
pc <- sapply(CL, function(c) { Q <- P[P$clade == c, ]; if (nrow(Q) < 20) return(NA)
  coef(glm(relaxed ~ S * real, binomial, Q))["S:real"] })
cat("per-clade beta_int:\n"); print(round(pc, 4))
cat("clades with pooled sign:", sum(sign(pc) == sign(b0), na.rm = TRUE), "of", sum(!is.na(pc)), "\n")
## ---------------- secondary ----------------
H <- reshape(P[, c("clade","hog","S","label","wpur_T","wpur_B")], idvar = c("clade","hog","S"), timevar = "label",
             direction = "wide")
H$d <- (lv(H$wpur_T.real) - lv(H$wpur_B.real)) - (lv(H$wpur_T.null1) - lv(H$wpur_B.null1))
H <- H[order(H$clade, H$hog), ]; H <- H[!is.na(H$d), ]   # fix 1: HOGs with one side dropped (NA relax) have NA d
Xs <- model.matrix(~ S + factor(clade), H); bs <- function(S) { Xs[, "S"] <- S; qr.coef(qr(Xs), H$d)["S"] }
s1 <- bs(H$S); hk <- paste(H$clade, H$hog)
s1p <- replicate(10000, { s <- H$S; for (c in CL) { i <- which(H$clade == c); s[i] <- s[i][sample.int(length(i))] }; bs(s) })
cat(sprintf("s1 Dpur(real)-Dpur(null1) ~ S: slope %.4f  perm P(2s) = %.4f (secondary, raw)\n", s1,
            (1 + sum(abs(s1p) >= abs(s1))) / 10001))
P$Sq <- cut(P$S, quantile(P$S, 0:4 / 4), include.lowest = TRUE, labels = paste0("Q", 1:4))
cat("s2 relaxed rate by S quartile (Q4 = strongest ancestral constraint):\n")
print(round(tapply(P$relaxed, list(P$Sq, P$label), mean), 4))
print(table(P$Sq, P$label))
own <- -lv(F$wpur_B[match(paste(P$clade, P$hog, "real"), paste(F$clade, F$hog, F$label))])  # fix 2: look up in F, not P
cat(sprintf("s3 (calibration, coupled same-fit baseline) beta_int = %.4f\n", bint(z(own))))
## ---------------- simulation calibration ----------------
W <- merge(D, D, by = "og"); W <- W[W$clade.x < W$clade.y, ]
rp <- sapply(split(W, paste(W$clade.x, W$clade.y)), function(w) if (nrow(w) >= 10) cor(w$v.x, w$v.y) else NA)
r <- mean(rp, na.rm = TRUE)
cat(sprintf("reliability r = %.4f (mean of %d clade-pair correlations, n>=10 shared OGs; range %.3f..%.3f)\n",
            r, sum(!is.na(rp)), min(rp, na.rm = TRUE), max(rp, na.rm = TRUE)))
n <- nrow(H); r <- max(min(r, 1), 0)
sim <- replicate(2000, { t <- rnorm(n); e <- function() sqrt(1 - r) * rnorm(n)
  own <- sqrt(r) * t + e(); crs <- sqrt(r) * t + e(); fg <- sqrt(r) * t + e(); d <- fg - own
  c(own = coef(lm(d ~ own))[2], cross = coef(lm(d ~ crs))[2]) })
cat(sprintf("sim null slope of Dpur on baseline: own-fit mean %.4f (sd %.4f); cross-clade mean %.4f (sd %.4f)\n",
            mean(sim[1, ]), sd(sim[1, ]), mean(sim[2, ]), sd(sim[2, ])))
## ---------------- blind direction check (R6-2b) ----------------
B <- rt("addon/R6/R6_2b/r6_2_fits.tsv"); B <- B[B$fit == "ok", ]
B$logK <- log10(pmin(pmax(as.numeric(B$K), .001), 50))
B <- B[!(paste(B$clade, B$hog) %in% paste(F$clade, F$hog)), ]
M <- aggregate(list(m = B$logK), list(clade = B$clade, hog = B$hog, set = B$set), mean)
M <- merge(M[M$set == "real", -3], M[M$set == "null1", -3], by = c("clade","hog"), suffixes = c(".r",".n"))
M$y <- M$m.r - M$m.n; M <- merge(M, A, by = c("clade","hog")); M$S0 <- Sraw(M$clade, M$og); M <- M[!is.na(M$S0), ]
M$S <- z(M$S0); M <- M[order(M$clade, M$hog), ]
Xb <- model.matrix(~ S + factor(M$clade), M); bb2 <- function(S) { Xb[, "S"] <- S; qr.coef(qr(Xb), M$y)["S"] }
br <- bb2(M$S)
brp <- replicate(10000, { s <- M$S; for (c in unique(M$clade)) { i <- which(M$clade == c); s[i] <- s[i][sample.int(length(i))] }; bb2(s) })
agree <- sign(br) == -sign(b0)
p1 <- if (b0 > 0) (1 + sum(brp <= br)) / 10001 else (1 + sum(brp >= br)) / 10001
cat(sprintf("blind R6-2b: HOGs %d (clades %d); slope of dlogK on S = %.4f; agrees in direction: %s; one-sided P = %.4f\n",
            nrow(M), length(unique(M$clade)), br, agree, p1))
cat("blind per-clade slope:\n"); print(round(sapply(split(M, M$clade), function(m) coef(lm(y ~ S, m))[2]), 4))
## ---------------- decision ----------------
sup <- p_prim < .05 && (ci[1] > 0 || ci[2] < 0) && sum(sign(pc) == sign(b0), na.rm = TRUE) >= 6 && agree && p1 < .05
cat("DECISION (ADDENDUM_R6_6):", if (sup) "SUPPORTED" else "NOT SUPPORTED", "\n")
write.table(U, file.path(OUT, "r6_6_units.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(M[, c("clade","hog","og","S0","S","y")], file.path(OUT, "r6_6_blind_units.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
