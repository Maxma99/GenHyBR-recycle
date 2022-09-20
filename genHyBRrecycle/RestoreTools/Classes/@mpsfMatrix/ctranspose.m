function A = ctranspose(A);
%
%  CTRANSPOSE The transpose of the mspfMatrix is
%             needed for matrix-vector multiply in iterative
%             restoration methods.  
%

%  J. Nagy  2/16/12

if A.transpose == 0
  A.transpose = 1;
else
  A.transpose = 0;
end
