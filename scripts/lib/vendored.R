## Vendored helper functions.
##
## Verbatim copies of three functions from David Duchene's shared ~/scripts
## collection, brought into the repository so that results here do not depend on
## an unversioned directory outside it. ~/scripts remains the upstream and is
## untouched; this is a copy, not a move.
##
## Provenance (md5 of the upstream file at the time of copying, 2026-08-06):
##   pathnode   ~/scripts/pathnode.R   c7bb3cf5bf429d52bcb2327983023698
##   rtt.cov    ~/scripts/rtt.cov.R    45026112cfe8b28f633933f0115c7045
##   stemmy     ~/scripts/stemmy.R     774e683142eafdc94b74b8e1521730aa
##
## Function bodies are byte-identical to upstream so behaviour cannot drift on
## copying. The only additions are this header and the library() call below --
## pathnode() calls require(phangorn) internally for Ancestors(), and rtt.cov()
## needs phangorn::midpoint, so the dependency is declared once at file level
## rather than re-checked on every one of 112,199 calls.
##
## If you improve one of these upstream, re-copy it here and update the md5.

suppressPackageStartupMessages(library(phangorn))

# ---------------------------------------------------------- ~/scripts/pathnode.R
pathnode <- function (phylo, tipsonly = T)
{
    require(phangorn)
    di.tr <- dist.nodes(phylo)
    root.tr <- phylo$edge[, 1][!(phylo$edge[, 1] %in% phylo$edge[,
        2])][1]
    tr.depth <- max(di.tr[as.numeric(colnames(di.tr)) == root.tr,
        ])
    if (tipsonly == TRUE) {
        roottotippath <- di.tr[as.numeric(rownames(di.tr)) ==
            root.tr, 1:length(phylo$tip.label)]
        nodesinpath <- sapply(1:length(phylo$tip.label), function(x) length(Ancestors(phylo,
            x)))
    }
    else {
        roottotippath <- di.tr[as.numeric(rownames(di.tr)) ==
            root.tr, ]
        nodesinpath <- sapply(1:(length(phylo$tip.label) + phylo$Nnode),
            function(x) length(Ancestors(phylo, x)))
    }
    #plot(roottotippath, nodesinpath, xlab = "Root-to-tip path length", ylab = "Number of parent nodes", pch = 20)
    return(list(roottotippath = roottotippath, nodesinpath = nodesinpath))
}

# ----------------------------------------------------------- ~/scripts/rtt.cov.R
rtt.cov <- function(tre){
        if(!is.rooted(tre)) tre <- midpoint(tre)
        rtts <- pathnode(tre)[[1]]
        rttcov <- sd(rtts, na.rm = T) / mean(rtts, na.rm = T)
        return(rttcov)
}

# ------------------------------------------------------------ ~/scripts/stemmy.R
stemmy <- function(tre) sum(tre$edge.length[which(tre$edge[,2] > Ntip(tre))], na.rm = T) / sum(tre$edge.length, na.rm = T)
