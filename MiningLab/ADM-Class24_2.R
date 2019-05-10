features = as.matrix(read.table('../../DataLists/ClassData/uspsdata.txt', sep = '\t'))
labels = read.table('../../DataLists/ClassData/uspscl.txt')[[1L]]

n = nrow(features)
d_1 = ncol(features)

x = features / 255
y = 0.5 * (labels + 1)

single_data = function(){
    f_1 = function(x, w_1_1, w_1_2){
        c(plogis(sum(x * w_1_1), plogis(x * w_1_2)))
    }
    
    f_2 = function(z_2, w_2){
        plogis(sum(z_2 * w_2))
    }
    
    l = function(y, f){
        -(y * log(f) + (1 - y) * log(1 - f))
    }
    
    dl_df = function(y, f){
        -(y - f) / (f * (1 - f))
    }
    
    df2_df1 = function(z_2, w_2){
        eta = sum(z_2 * w_2)
        plogis(eta) * plogis(eta, lower.tail = FALSE) * w_2
    }
    
    df2_dw2 = function(z_2, w_2){
        eta = sum(z_2, w_2)
        plogis(eta) * plogis(eta, lower.tail = FALSE) * z_2
    }
    
    df1_dw1 = function(x, w_1_1, w_1_2){
        eta_1 = as.vector(x %*% w_1_1)
        eta_2 = as.vector(x %*% w_1_2)
        rbind(c(plogis(eta_1) * plogis(eta_1, lower.tail = FALSE) * x, rep(0, d_1)),
              c(rep(0, d_1), plogis(eta_2) * plogis(eta_2, lower.tail = FALSE) * x))
    }
}

library(Matrix)

f_1 = function(x, w_1_1, w_1_2){
rbind(t(plogis(x %*% w_1_1)), t(plogis(x %*% w_1_2)))
}

f_2 = function(z_2, w_2){
plogis(as.vector(crossprod(z_2, w_2)))
}

l = function(y, f_x){
-mean(y * log(f_x) + (1 - y) * log(1 - f_x))
}

dl_df = function(y, f_x){
-(y - f_x) / (f_x * (1 - f_x))
}

df2_df1 = function(z_2, w_2){
eta = as.vector(crossprod(z_2, w_2))
sigma.p = plogis(eta) * plogis(eta, lower.tail = FALSE)
sigma.p * do.call("cbind", lapply(w_2, function(w_2_j) Diagonal(n, w_2_j)))
# plogis(eta) * plogis(eta, lower.tail = FALSE) * bdiag(replicate(n, t(w_2), simplify = FALSE))
}

df2_dw2 = function(z_2, w_2){
eta = as.vector(crossprod(z_2, w_2))
plogis(eta) * plogis(eta, lower.tail = FALSE) * t(z_2)
}

df1_dw1 = function(x, w_1_1, w_1_2){
if(!is.matrix(x)) x  = t(x)
eta_1 = as.vector(x %&% w_1_1)
eta_2 = as.vector(x %*% w_1_2)
bdiag(plogis(eta_1) * plogis(eta_1, lower.tail = FALSE) * x,
      plogis(eta_2) * plogis(eta_2, lower.tail = FALSE) * x)
}

image.print = function(x, ...){
    x.matrix = matrix(x, 16, 16, byrow = FALSE)
    x.matrix.rotated = t(apply(x.matrix, 1, rev))
    image(x.matrix.rotated, axes = FALSE, col = grey(seq(0, 1, length.out = 256)), ...)
}

w = rnorm(2 * d_1 + 2)

w_1 = w[seq.int(1, 2 * d_1)]
w_1_1 = w_1[seq.int(1, d_1)]
w_1_2 = w_1[seq.int(d_1 + 1, length(w_1))]
w_2 = w[seq.int(2 * d_1 + 1, length(w))]

iter = 0L
maxIter = 100L
tolerance = 1e-6
maxStepSize = 150

while(iter <= maxIter){
    iter = iter + 1L
    z_2 = f_1(x, w_1_1, w_1_2)
    f = f_2(z_2, w_2)
    g_w_2_all = dl_df(labels, f) * df2_dw2(z_2, w_2)
    g_w_2 = apply(g_w_2_all, 2, mean)
    f.t = function(t){
        l(y, f_2(z_2, w_2 - t * g_w_2))
    }

    stepSize = optimize(f.t, c(0, maxStepSize))$minimum

    w_2 = w_2 - stepSize * g_w_2
    
    f = f_2(z_2, w_2)
    g_w_1_all = dl_df(y, f) * df2_df1(z_2, w_2) %*% df1_dw1(x, w_1_1, w_1_2)
    g_w_1   <- apply(g_w_1_all, 2, mean)
    g_w_1_1 <- g_w_1[seq.int(1, d_1)]
    g_w_1_2 <- g_w_1[seq.int(d_1 + 1, length(g_w_1))]

    f.t <- function(t) l(y, f_2(f_1(x, w_1_1 - t * g_w_1_1, w_1_2 - t * g_w_1_2), w_2))
    stepSize <- optimize(f.t, c(0, maxStepSize))$minimum

    w_1 <- w_1 - stepSize * g_w_1
    w_1_1 <- w_1[seq.int(1, d_1)]
    w_1_2 <- w_1[seq.int(d_1 + 1, length(w_1))]

    if (all(abs(g_w_2) < tolerance) && all(abs(g_w_1) < tolerance)) break
}

f <- f_2(f_1(x, w_1_1, w_1_2), w_2)

y_hat <- ifelse(f > 0.5, 1, 0)
cat("accuracy: ", mean(y_hat == y), "\n")
misclassified <- which(y_hat != y)

par(mfrow = c(2, 2), mar = c(0.1, 0.1, 2.2, 0.1))
image.print(w_1_1, main = "First Weights")
image.print(w_1_2, main = "Second Weights")
if (length(misclassified) > 0) image.print(features[misclassified[1],], main = "Misclassified")
if (length(misclassified) > 1) image.print(features[misclassified[2],], main = "Misclassified")
