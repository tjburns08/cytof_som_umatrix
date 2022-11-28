umatrixForEsom_adapted <- function (Weights, 
                                    Lines = NULL, 
                                    Columns = NULL, 
                                    Toroid = T, 
                                    dist_measure = c("euclidean", "chebyshev"),
                                    neighbor_agg = c("mean", "max")) 
{
    EsomNeurons = Umatrix:::ListAsEsomNeurons(Weights, Lines, Columns)
    nachbarn <- function(k, EsomNeurons, Toroid = FALSE) {
        M <- dim(EsomNeurons)[1]
        N <- dim(EsomNeurons)[2]
        if (Toroid) {
            pos1 = c(k[1] - 1, k[1] - 1, k[1] - 1, k[1], k[1], 
                     k[1] + 1, k[1] + 1, k[1] + 1)%%M
            pos1[which(pos1 == 0)] = M
            pos2 = c(k[2] - 1, k[2], k[2] + 1, k[2] - 1, k[2] + 
                         1, k[2] - 1, k[2], k[2] + 1)%%N
            pos2[which(pos2 == 0)] = N
            nb = cbind(pos1, pos2)
        }
        else {
            if (k[1] == 1) {
                if (k[2] == 1) {
                    nb = rbind(c(1, 2), c(2, 2), c(2, 1))
                }
                else {
                    if (k[2] == N) {
                        nb = rbind(c(1, (N - 1)), c(2, (N - 1)), 
                                   c(2, N))
                    }
                    else {
                        nb = rbind(c(1, (k[2] - 1)), c(1, (k[2] + 
                                                               1)), c(2, (k[2] - 1)), c(2, k[2]), c(2, 
                                                                                                    (k[2] + 1)))
                    }
                }
            }
            if (k[1] == M) {
                if (k[2] == 1) {
                    nb = rbind(c((M - 1), 1), c((M - 1), 2), c(M, 
                                                               2))
                }
                else {
                    if (k[2] == N) {
                        nb = rbind(c((M - 1), (N - 1)), c((M - 1), 
                                                          N), c(M, (N - 1)))
                    }
                    else {
                        nb = rbind(c((M - 1), (k[2] - 1)), c((M - 
                                                                  1), k[2]), c((M - 1), (k[2] + 1)), c(M, 
                                                                                                       (k[2] - 1)), c(M, (k[2] + 1)))
                    }
                }
            }
            if (k[1] != 1 && k[1] != M) {
                if (k[2] == 1) {
                    nb = rbind(c((k[1] - 1), 1), c((k[1] - 1), 
                                                   2), c(k[1], 2), c((k[1] + 1), 1), c((k[1] + 
                                                                                            1), 2))
                }
                else {
                    if (k[2] == N) {
                        nb = rbind(c((k[1] - 1), (N - 1)), c((k[1] - 
                                                                  1), N), c(k[1], (N - 1)), c((k[1] + 1), 
                                                                                              (N - 1)), c((k[1] + 1), N))
                    }
                    else {
                        nb = rbind(c((k[1] - 1), (k[2] - 1)), c((k[1] - 
                                                                     1), k[2]), c((k[1] - 1), (k[2] + 1)), 
                                   c(k[1], (k[2] - 1)), c(k[1], (k[2] + 1)), 
                                   c((k[1] + 1), (k[2] - 1)), c((k[1] + 1), 
                                                                k[2]), c((k[1] + 1), (k[2] + 1)))
                    }
                }
            }
        }
        return(nb)
    }
    if (is.list(EsomNeurons)) 
        EsomNeurons = EsomNeurons$EsomNeurons
    k = dim(EsomNeurons)[1]
    m = dim(EsomNeurons)[2]
    Umatrix = matrix(0, k, m)
    d = dim(EsomNeurons)[3]
    if (is.null(d)) {
        stop("EsomNeurons wts has to be an array[1:Lines,1:Columns,1:Weights], use ListAsEsomNeurons")
    }
    for (i in 1:k) {
        for (j in 1:m) {
            nbs = nachbarn(c(i, j), EsomNeurons, Toroid)
            wij = EsomNeurons[i, j, ]
            n.nbs = dim(nbs)[1]
            for (l in 1:n.nbs) {
                nij = EsomNeurons[nbs[l, 1], nbs[l, 2], ]
                
                # Adding a new series of if-else statements for distance
                #   and neighbor aggregate measures
                if(dist_measure == "euclidean") {
                    # L2 norm distance for neighbors
                    if(neighbor_agg == "max") {
                        Umatrix[i, j] = max(Umatrix[i, j], sqrt(sum((wij - nij)^2)))
                    } else {
                        Umatrix[i, j] = Umatrix[i, j] + sqrt(sum((wij - nij)^2))
                    }
                } else if(dist_measure == "chebyshev") {
                    if(neighbor_agg == "max") {
                        # Max of neighbor cheb distance
                        Umatrix[i, j] = max(Umatrix[i, j],  max(abs(wij - nij))) # Max of neighbor cheb distances
                    } else {
                        # Cheb distance for neighbors
                        Umatrix[i, j] = Umatrix[i, j] + max(abs(wij - nij)) # Cheb distance
                    }
                }
            }
            # Average value of whatever we've got
            if(neighbor_agg == "mean") {
                Umatrix[i, j] = Umatrix[i, j]/n.nbs
            }
        }
    }
    return(Umatrix)
}
