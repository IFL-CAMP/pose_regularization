function D = vec_matmult(A,transA,B,transB)
    if transA
        A = permute(A,[2,3,1]);
    else
        A = permute(A,[1,3,2]);
    end
    if transB
        B = permute(B,[2,1,3]);
    end
    D = reshape(sum(bsxfun(@times,reshape(A,[size(A,1),size(A,2),1,size(A,3)]),reshape(B,[1,size(B)])),2),size(A,1),size(B,2),[]);
end