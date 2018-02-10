set terminal postscript eps enhanced color "Helvetica" 30
set xlabel "M_{QEDxQCD}/TeV"
set ylabel "m/GeV"
#set title "95% CL limits, narrow width"
#set logscale x
set mxtics 5
set termoption dash
set for [i=1:5] linetype i dt i
set style line 1 lt 2 lc rgb "red" lw 3
set style line 2 lt 2 lc rgb "orange" lw 2
set style line 3 lt 2 lc rgb "yellow" lw 3
set style line 4 lt 2 lc rgb "green" lw 2
set label 1 "Allanach, Bhatia, Iyer" at 0,1416 font "Helvetica,12"
set key font "Helvetica,20"
set title "m_0=M_{1/2}=5 TeV, A_0=-14 TeV, tan{/Symbol b}=20"
set output "mh.eps"
plot "outputTest" u ($1/1000):2 tit "h^0/GeV" w l lw 4 

set output "mH.eps"
plot "outputTest" u ($1/1000):3 tit "A^0/GeV" w l lw 4, \
"outputTest" u ($1/1000):4 tit "H^0/GeV" w l lw 4, \
"outputTest" u ($1/1000):5 tit "H^+/GeV" w l lw 4
