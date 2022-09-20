function A = ctranspose(A);
%
%  CTRANSPOSE The transpose of the mspfMatrix0 is
%             needed for matrix-vector multiply in iterative
%             restoration methods.  
%

%  J. Nagy  2/23/12

if A.transpose == 0
  A.transpose = 1;
else
  A.transpose = 0;
end
