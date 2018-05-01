lp = function(trainingset, testset, column, upperBound = 1e16){

    matrixA = subset(trainingset, select = column)
    effort = trainingset$Effort
    nV = length(column)

    nConstrains = 2*length(column) + 4*dim(trainingset)[1]
    nVar = length(column) + 2*dim(trainingset)[1]
    mat1 = matrix(NA, ncol = nVar, nrow = dim(trainingset)[1])
    mat2 = matrix(NA, ncol = nVar, nrow = dim(trainingset)[1])
    mat3 = matrix(NA, ncol = nVar, nrow = dim(trainingset)[1])
    mat4 = matrix(0, ncol = nVar, nrow = nV)
    mat5 = matrix(0, ncol = nVar, nrow = nV)
    mat6 = matrix(0, ncol = nVar, nrow = dim(trainingset)[1])
    
    tzip = (1:dim(mat1)[1]*2+nV-1)
    cVec = rep(0, nVar)
    cVec[tzip] = 1
    for(i in 1:dim(matrixA)[1]){
      for(j in 1:dim(matrixA)[2]){
        mat1[i,j] = matrixA[i,j]
      }
    }
    mat1[is.na(mat1)] = 0
    b1 = rep(0, dim(mat1)[1])
    mat2 = mat1
   
    tzip= cbind((1:dim(mat2)[1]), (1:dim(mat2)[1]*2+nV-1))
    mat2[tzip] = -1
    effzip= cbind((1:dim(mat2)[1]), (1:dim(mat2)[1]*2+nV))
    mat2[effzip] = -1
    b2 = rep(0, dim(mat2)[1])
    
    mat3 = mat2
    tzip= cbind((1:dim(mat3)[1]), (1:dim(mat3)[1]*2+nV-1))
    mat3[tzip] = 1
    b3 = rep(0, dim(mat3)[1])
    
    diag(mat4) = 1
    b4 = rep(0, dim(mat4)[1])
    diag(mat5) = 1
    
    b5 = rep(upperBound, dim(mat5)[1])
    effzip= cbind((1:dim(mat6)[1]), (1:dim(mat6)[1]*2+nV))
    mat6[effzip] = 1
    b6 = effort
    A = rbind(mat1,mat2,mat3,mat4,mat5,mat6)
    constr = c(rep(">", dim(mat1)[1]), rep("<=", dim(mat2)[1]), rep(">=", dim(mat3)[1]), rep(">=", dim(mat4)[1]), rep("<=", dim(mat5)[1]), rep("=", dim(mat6)[1]))
    B = c(b1,b2,b3,b4,b5,b6)
    solSimplex = solveLP(cVec, B, A, const.dir = constr,lpSolve = TRUE, zero = 1e-16, tol = 1E-3)
    variablesValue = solSimplex$solution[1:nV]
    matrixATest = subset(testset, select = column)
    measured = testset$Effort
    
    prediction = rowSums(sweep(matrixATest, MARGIN = 2, variablesValue, `*`))

    return(list(data.frame(prediction,measured), variablesValue))
}