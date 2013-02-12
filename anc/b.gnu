set term post enhanced color "Helvetica" 24 eps
set size square

set output "multiBCs.eps"
set title "CMSSM10.1.1, SOFTSUSY3.3.6"
set yrange [0.3:2]
set logscale x
set xrange [10:3e16]
set xtics ("M_Z" 91.187, "M_{S}" 2.e3, "" 0.1e16, "" 0.5e16, "" 1e16, "" 1.5e16, "M_{GUT}" 1.841953e+16)
plot "twoCond" i 0 u ($1):($5/500) tit "M_i/(500 GeV)" w l ls 1 lw 4 lc 2, \
"twoCond" i 0 u ($1):($6/500) notit w l ls 1 lw 4 lc 2, \
"twoCond" i 0 u ($1):($7/500) notit w l ls 1 lw 4 lc 2, \
"twoCond" i 0 u ($1):2 tit "g_{1,2}" w l ls 1 lw 4, \
"twoCond" i 0 u ($1):3 notit w l ls 1 lc 1 lw 4, \
"twoCond" i 0 u ($1):4 tit "h_t" w l lw 4 lc 3, \
"twoCond" i 1 u ($1):2 tit "boundary condition" w p  pt 7, \
"twoCond" i 1 u ($1):3 notit w p pt 7 lc 0, \
"twoCond" i 1 u ($1):4 notit w p pt 7 lc 0, \
"twoCond" i 2 u ($1):($5/500) notit w p pt 7 lc 0
