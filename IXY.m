function I = IXY(py_x, px)
% px is  a row vector
% py_x is the channel matrix, rows are input, colmns are outputs 
% (each row, one input, ==> vector pf p(y|.) as a row)
[nx,ny]=size(py_x);

[a,b]=size(px);

if b==1
    px=px';
end

py=px*py_x;
Hy=-sum(py.*log(py+eps));

Hy_x=-px*sum(py_x.*log(py_x+eps),2); % H(Y|X)=E_x (H(Y|x))

I=Hy-Hy_x;





