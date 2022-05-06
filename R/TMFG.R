#' Triangulated Maximally Filtered Graph
#'
#' @description Applies the Triangulated Maximally Filtered Graph (TMFG) filtering method
#' (\strong{please see and cite Massara et al., 2016}). The TMFG method uses a structural
#' constraint that limits the number of zero-order correlations included in the network
#' (3n - 6; where \emph{n} is the number of variables). The TMFG algorithm begins by
#' identifying four variables which have the largest sum of correlations to all other
#' variables. Then, it iteratively adds each variable with the largest sum of three
#' correlations to nodes already in the network until all variables have been added to
#' the network. This structure can be associated with the inverse correlation matrix
#' (i.e., precision matrix) to be turned into a GGM (i.e., partial correlation network)
#' by using Local-Global Inversion Method (see Barfuss et al., 2016 for more details).
#' See Details for more information on this network estimation method.
#'
#' @param cormat A correlation matrix
#'
#' @return Returns a list containing:
#'
#' \item{A}{The filtered adjacency matrix}
#'
#' \item{separators}{The separators (3-cliques) in the network}
#'
#' \item{cliques}{The cliques (4-cliques) in the network}
#'
#' @details The TMFG method applies a structural constraint on the network,
#' which restrains the network to retain a certain number of edges (3\emph{n}-6, where \emph{n}
#' is the number of nodes; Massara et al., 2016). The network is also composed of 3- and 4-node
#' cliques (i.e., sets of connected nodes; a triangle and tetrahedron, respectively). The
#' TMFG method constructs a network using zero-order correlations and the resulting network
#' can be associated with the inverse covariance matrix
#' (yielding a GGM; Barfuss, Massara, Di Matteo, & Aste, 2016).
#' Notably, the TMFG can use any association measure and thus does not assume the data is multivariate normal.
#'
#' Construction begins by forming a tetrahedron of the four nodes that have
#' the highest sum of correlations that are greater than the average correlation in the
#' correlation matrix. Next, the algorithm iteratively identifies the node that maximizes
#' its sum of correlations to a connected set of three nodes (triangles) already included
#' in the network and then adds that node to the network. The process is completed once
#' every node is connected in the network. In this process, the network automatically
#' generates what's called a planar network. A planar network is a network that could be
#' drawn on a sphere with no edges crossing (often, however, the networks are depicted
#' with edges crossing; Tumminello, Aste, Di Matteo, & Mantegna, 2005).
#'
#' @examples
#' # Pearson's correlation only for CRAN checks
#' A <- TMFG(cor(wmt2[,7:24]))$A
#'
#' @references
#' Barfuss, W., Massara, G. P., Di Matteo, T., & Aste, T. (2016).
#' Parsimonious modeling with information filtering networks.
#' \emph{Physical Review E}, \emph{94}, 062306.
#'
#' Christensen, A. P., Kenett, Y. N., Aste, T., Silvia, P. J., & Kwapil, T. R. (2018).
#' Network structure of the Wisconsin Schizotypy Scales-Short Forms: Examining psychometric network filtering approaches.
#' \emph{Behavior Research Methods}, \emph{50}, 2531-2550.
#'
#' Massara, G. P., Di Matteo, T., & Aste, T. (2016).
#' Network filtering for big data: Triangulated maximally filtered graph.
#' \emph{Journal of Complex Networks}, \emph{5}, 161-178.
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
# TMFG Filtering Method----
# Updated 04.05.2022
TMFG <-function (cormat)
{
    # Number of nodes
    n <- ncol(cormat)

    # Signed correlations
    tcormat <- cormat
    cormat <- abs(cormat)

    # Let user know matrix is too small for TMFG estimation
    # It is still okay to proceed
    if(n < 9)
    {print("Matrix is too small")}

    # Initialize sparse edge matrix
    nodeTO <- sort(c(rep(1:n,n)))
    nodeFROM <- c(rep(1:n,n))
    nodeWEIGHT <- as.vector(cormat)

    # Initialize matrices
    M <- cbind(nodeTO, nodeFROM, nodeWEIGHT) # sparse node-weight matrix
    in_v <- matrix(nrow=nrow(cormat), ncol=1) # inserted vertices matrix
    ou_v <- matrix(nrow=nrow(cormat), ncol=1) # not yet inserted vertices matrix
    tri <- matrix(nrow=((2*n)-4), ncol=3) # triangles matrix
    separators <- matrix(nrow=n-4, ncol=3) # matrix of 3-cliques (non-face triangles)

    # Find 3 vertices with largest strength
    s <- rowSums(cormat*(cormat > mean(matrix(unlist(cormat), nrow=1)))*1)

    # Order vertices with largest strength
    # and grab the top 4
    in_v[1:4] <- order(s,decreasing=TRUE)[1:4]

    # Set aside nodes that are not in the top 4
    ou_v <- setdiff(1:nrow(in_v),in_v)

    # Build tetrahedron with the largest strength
    ## Insert triangles
    tri[1,]<-in_v[1:3,]
    tri[2,]<-in_v[2:4,]
    tri[3,]<-in_v[c(1,2,4),]
    tri[4,]<-in_v[c(1,3,4),]

    # Initialize sparse TMFG matrix
    S <- matrix(nrow=(3*nrow(cormat)-6),ncol=3)

    # Algorithm for traditional network
    S[1,] <- c(in_v[1],in_v[2],1)
    S[2,] <- c(in_v[1],in_v[3],1)
    S[3,] <- c(in_v[1],in_v[4],1)
    S[4,] <- c(in_v[2],in_v[3],1)
    S[5,] <- c(in_v[2],in_v[4],1)
    S[6,] <- c(in_v[3],in_v[4],1)

    #build initial gain table
    gain <- matrix(-Inf,nrow=n,ncol=(2*(n-2)))
    gain[ou_v,1] <- rowSums(cormat[ou_v,(tri[1,])])
    gain[ou_v,2] <- rowSums(cormat[ou_v,(tri[2,])])
    gain[ou_v,3] <- rowSums(cormat[ou_v,(tri[3,])])
    gain[ou_v,4] <- rowSums(cormat[ou_v,(tri[4,])])

    ntri <- 4 #number of triangles
    gij <- matrix(nrow=1,ncol=ncol(gain))
    v <- matrix(nrow=1,ncol=ncol(gain))
    ve <- array()
    tr <- 0
    for(e in 5:n)
    {
        if(length(ou_v)==1){
            ve<-ou_v
            v<-1
            w<-1
            tr<-which.max(gain[ou_v,])
        }else{
            for(q in 1:ncol(gain))
            {
                gij[,q] <- max(gain[ou_v,q])
                v[,q] <- which.max(gain[ou_v,q])
                tr <- which.max(gij)
            }

            ve <- ou_v[v[tr]]
            w <- v[tr]
        }

        #update vertex lists
        ou_v<-ou_v[-w]
        in_v[e]<-ve

        #update adjacency matrix
        for(u in 1:length(tri[tr,]))
        {
            cou<-6+((3*(e-5))+u)
            S[cou,]<-cbind(ve,tri[tr,u],1)
        }

        #update 3-clique list
        separators[e-4,]<-tri[tr,]
        #update triangle list replacing 1 and adding 2 triangles
        tri[ntri+1,]<-cbind(rbind(tri[tr,c(1,3)]),ve)
        tri[ntri+2,]<-cbind(rbind(tri[tr,c(2,3)]),ve)
        tri[tr,]<-cbind(rbind(tri[tr,c(1,2)]),ve)
        #update gain table
        gain[ve,]<-0
        gain[ou_v,tr]<-rowSums(cormat[ou_v,tri[tr,],drop=FALSE])
        gain[ou_v,ntri+1]<-rowSums(cormat[ou_v,tri[ntri+1,],drop=FALSE])
        gain[ou_v,ntri+2]<-rowSums(cormat[ou_v,tri[ntri+2,],drop=FALSE])

        #update triangles
        ntri<-ntri+2
    }
    cliques<-rbind(in_v[1:4],(cbind(separators,in_v[5:ncol(cormat)])))

    L<-S
    L[,1]<-S[,2]
    L[,2]<-S[,1]
    K<-rbind(S,L)

    x <- matrix(0, nrow = ncol(cormat), ncol = ncol(cormat))

    for(i in 1:nrow(K))
    {
        x[K[i,1], K[i,2]] <- 1
        x[K[i,2], K[i,1]] <- 1
    }

    diag(x)<-1

    for(r in 1:nrow(x))
        for(z in 1:ncol(x))
        {if(x[r,z]==1){x[r,z]<-tcormat[r,z]}}

    colnames(x)<-colnames(cormat)
    x <- as.data.frame(x)
    row.names(x)<-colnames(x)
    x <- as.matrix(x)

    return(list(A=x, separators=separators, cliques=cliques))
}
