set term post enhanced color "Helvetica" 24

set ylabel "M_Z^2(pred)/M_Z^2(exp)"
#set yrange [-40:40]
set title "m_0=350+(k*500) GeV, M_{1/2}=300 GeV, tan{/Symbol b}=10, A_0=0"

set output "largeSlope.eps"
set xlabel "[{/Symbol m}(M_{SUSY})/TeV]^2"
plot "cmssmScans" i 0 u (($1)*($1)/1000000):2 tit "default" w l, 1 notit w l ls 1 lc 2

set output "multi.eps"
set xlabel "{/Symbol m}(M_{SUSY})/GeV"
plot "cmssmScans" i 0 u 1:2 tit "default" w l, 1 notit w l ls 1 lc 2

set output "zoom1.eps"
set title "m_0=2360 GeV, M_{1/2}=300 GeV, tan{/Symbol b}=10, A_0=0"
plot "cmssmZoom" i 0 u 1:2  tit "default" w l, 1 notit w l ls 1 lc 2

set output "zoom2.eps"
set yrange [8.5:9.1]
set xrange [50:55]
plot "cmssmZoom" i 0 u 1:2  tit "default" w l, 1 notit w l ls 1 lc 2
set autoscale
set title "m_0=350+(k*500) GeV, M_{1/2}=300 GeV, tan{/Symbol b}=10, A_0=0"

set output "multi1.eps"
plot "cmssmScans" i 0 tit "default" w l, 1 notit w l ls 1 lc 2, \
"cmssmScans" i 1 tit "m_t=173.5-1" w l ls 2 lc 3

set output "multi2.eps"
plot "cmssmScans" i 0 tit "default" w l, 1 notit w l ls 1 lc 2, \
"cmssmScans" i 2 tit "m_t=173.5+1" w l ls 2 lc 3

set output "multi3.eps"
plot "cmssmScans" i 0 tit "default" w l, 1 notit w l ls 1 lc 2, \
"cmssmScans" i 3 tit "{/Symbol a}_s(M_Z)=0.1187-0.001" w l ls 2 lc 3

set output "multi4.eps"
plot "cmssmScans" i 0 tit "default" w l, 1 notit w l ls 1 lc 2, \
"cmssmScans" i 4 tit "{/Symbol a}_s(M_Z)=0.1187+0.001" w l ls 2 lc 3

set output "multi7.eps"
plot "cmssmScans" i 0 tit "default" w l, 1 notit w l ls 1 lc 2, \
"cmssmScans" i 7 tit "A_0=500 GeV" w l ls 2 lc 3

set output "multi8.eps"
plot "cmssmScans" i 0 tit "default" w l, 1 notit w l ls 1 lc 2, \
"cmssmScans" i 8 tit "A_0=-500 GeV" w l ls 2 lc 3

set output "multi9.eps"
plot "cmssmScans" i 0 tit "default" w l, 1 notit w l ls 1 lc 2, \
"cmssmScans" i 9 tit "tan {/Symbol b}=40" w l ls 2 lc 3

set output "multi10.eps"
plot "cmssmScans" i 0 tit "default" w l, 1 notit w l ls 1 lc 2, \
"cmssmScans" i 10 tit "M_{1/2}=200 GeV" w l ls 2 lc 3

set output "multi11.eps"
plot "cmssmScans" i 0 tit "default" w l, 1 notit w l ls 1 lc 2, \
"cmssmScans" i 11 tit "M_{1/2}=400 GeV" w l ls 2 lc 3

set title "m_0=2350 GeV, M_{1/2}=300 GeV, tan{/Symbol b}=10, A_0=0"
unset ylabel
set xtics 0.1
set output "reason.eps"
plot "zoomEwsb" i 0 u 1:(4*($2-($7-$8))) tit \
"M_Z-(m_{{/Symbol c}@_1^0}+m_{{/Symbol c}@_2^0})" w l, 0 notit w l, \
"zoomEwsb" i 0 u 1:(10*($3-8.98)-0.5) tit "f(M_Z^2(pred)/M_Z^2(exp))" w l



set output "reason2.eps"
plot "zoomEwsb" i 0 u 1:11 tit "{/Symbol P}_{WW}^{T}(M_W)" w l
set xtics 1
#plot "zoomEwsb" i 1 u 1:($7+$6), "zoomEwsb" i 1 u 1:5
#plot "zoomEwsb" i 1 u 1:(($5-($7+$6))) tit \
#"M_Z-2m_{{/Symbol c}@_2^0}" w l, 0 notit w l, \
#"zoomEwsb" i 1 u 1:(10*($3-8.4)) tit "f(M_Z^2(pred)/M_Z^2(exp))" w l
unset xtics


!gs -q -sPAPERSIZE=a4 -dNOPAUSE -dBATCH -sDEVICE=pswrite -sOutputFile=concat.ps *.eps
!mpage -2 concat.ps > b.ps
!ps2pdf b.ps 
#!rm *.ps


