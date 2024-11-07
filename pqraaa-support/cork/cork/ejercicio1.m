A = 2*eye(10,10);
A = A + diag(-ones(9,1), 1) ;
A = A + diag(-ones(9,1), -1);

v = [0:9]' ;
b = sin( pi*v/9 ) ;

x = A \ b ;

plot( x, 'r*' )