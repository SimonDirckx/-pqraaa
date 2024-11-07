pow = 17:31;

Gflops = [ 7.417955 0.491959 0.483068 0.358251 0.249366
10.076897 0.862598 0.827435 0.521475 0.418573
0.000000 1.473994 1.424338 0.961421 0.685409
18.266415 2.229162 2.129302 1.683350 0.962846
25.060345 2.780442 2.891769 2.482300 1.169514
31.929608 3.330811 3.269664 3.122634 1.496362
21.334047 3.535270 3.441769 3.683177 1.672910
30.204368 4.009108 4.196560 4.491795 1.883472
0.000000 3.824490 3.847874 4.952118 2.042917
29.209406 1.395131 2.450325 3.765521 1.977939
0.000000 3.341329 3.854319 4.531418 2.317750
0.000000 3.535953 4.107432 5.082386 2.507971
0.000000 3.721223 4.164557 5.170724 2.528914
28.578635 4.191889 4.616543 5.694363 2.594143
29.257216 4.463114 5.137354 6.092524 2.852739 ];

Vars = [ 0.793798 22.328185 22.823241 30.052837 44.469741
0.928256 26.188612 27.746124 44.681241 56.118818
1.250224 31.459064 33.600854 51.088721 73.020287
2.200090 44.045533 47.123332 60.224594 107.004087
4.355269 76.094461 73.898164 85.335953 183.775587
8.807876 134.778699 137.858811 146.404226 306.133079
37.552932 268.002757 274.126703 257.598752 573.960957
53.109794 494.932512 476.521808 450.796118 1068.727941
220.254169 1092.045406 1094.584071 854.911451 2067.953972
284.512154 6318.145840 3597.373228 2341.677832 4447.206404
663.567074 5473.153597 4750.624144 4033.943640 7873.952926
1669.920485 10748.327124 9249.718005 7462.104954 15017.422910
3539.529820 21131.117536 18887.224697 15217.121666 30927.331112
5672.441850 38801.342315 35195.029563 28575.956076 62662.151220
11466.069369 75189.620734 65302.380079 55161.371918 117436.198071 ];

figure(1)
relGflops = zeros( size( Gflops,1 ), size( Gflops, 2 )-2 );
relVars   = relGflops;
for i=1:size(Gflops,1)
	for j=1:size(relGflops,2)
		relGflops(i,j)=Gflops(i,j+2) ./ Gflops(i,2);
		relVars(i,j)=Vars(i,j+2) ./ Vars(i,2);
	end
end
bar( pow, relGflops );
xlabel( 'problem size (log_2n)' );
ylabel( 'performance compared to 1x1x64' );
legend( {'1x2x32','1x8x8','2x4x8'}, 'location', 'northwest' );
axis( [16 32 0 2.8] );
title( 'Relative performance (higher is better)' );
print( "performance.eps", "-deps", "-color", "-F:17" );

figure(2)
bar( pow, relVars );
%bar( pow, log10(Vars(:,2:end)) );
xlabel( 'problem size (log_2n)' );
%ylabel( 'log_{10}(stddev) in ms.' );
ylabel( 'stddev compared to 1x1x64' );
%legend( {'1x1x64','1x2x32','1x8x8','2x4x8'}, 'location', 'northwest' );
legend( {'1x2x32','1x8x8','2x4x8'}, 'location', 'northeast' );
axis( [16 32 0 2.5] );
title( 'Relative performance stability (lower is better)' );
print( "stability.eps", "-deps", "-color", "-F:17" );

