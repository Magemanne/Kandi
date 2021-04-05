datamatrix = zeros(6,1270)
datamatrix(2,1:length(u)) = u;
datamatrix(3,1:length(u)) = u;
b = u( 1:find( u > 100000, 1) );
