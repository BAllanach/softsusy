set term post enhanced color "Helvetica" 24

set ylabel "M_Z^2(pred)/M_Z^2(exp)"
set xlabel "{/Symbol m}(M_{SUSY})/GeV"

set output "muPoint.eps"
set title "{/Symbol m}(BC)=(8.7*k+0.1) GeV,  m_0=3.1 TeV, M_{1/2}=300 GeV, tan{/Symbol b}=10" 
plot "muPoint" i 0 u 1:2 notit w l lw 4, 1 notit w l lc 2 lw 4


set output "multi.eps"
set title "M_{1/2}=300 GeV, tan{/Symbol b}=10, A_0=0"
set xlabel "{/Symbol m}(M_{SUSY})/GeV"
plot "cmssmScans" i 0 u 1:3 tit "m_0=0.6 TeV" w l lw 4, 1 notit w l ls 0, \
"cmssmScans" i 1 u 1:3 tit "m_0=2.6 TeV" w l lw 4, \
"cmssmScans" i 2 u 1:3 tit "m_0=3.1 TeV" w l lw 4, \
"cmssmScans" i 3 u 1:3 tit "m_0=3.6 TeV" w l lw 4, \
"cmssmScans" i 4 u 1:2 notit w l ls 0, \
"cmssmScans" i 5 u 1:2 notit w l ls 0, \
"cmssmScans" i 6 u 1:2 notit w l ls 0

set output "zoom.eps"
set title "m_0=3.1 TeV, M_{1/2}=300 GeV, tan{/Symbol b}=10, A_0=0"
set yrange [-8:15]
unset ylabel
set key spacing 1.5
plot "cmssmScans" i 2 u 1:3 tit "M@_Z^2(pred)/M@_Z^2(exp)" w l lw 4, \
"cmssmScans" i 2 u 1:(($2-abs($6)-abs($7))/10) tit "(M_Z(exp)-|M_{{/Symbol c}@_1^0}|-|M_{{/Symbol c}@_2^0|})/10 GeV" w l lw 4, \
"cmssmScans" i 2 u 1:(($2-abs($5)-abs($5))/10) tit "(M_Z(exp)-2|M_{{/Symbol c}@_1^+}|)/10 GeV" w l lw 4, \
0 notit w l lt 0, \
"cmssmScans" i 7 u 1:2 notit w l ls 0, \
"cmssmScans" i 8 u 1:2 notit w l ls 0, \
"cmssmScans" i 9 u 1:2 notit w l ls 0, \
"cmssmScans" i 2 u 1:((($5))/10) tit "M_{{/Symbol c}@_1^+}/10 GeV" w l lw 4 lt 4, \
"cmssmScans" i 2 u 1:(5000*$8/($2*$2)) tit "5x10^3 {/Symbol P}@_{Z}^{T}/M@_Z^2(exp)" w l lt 5 lw 4
set output 

set autoscale
set output "inv.eps"
plot "oneScan" u 1:($6/1.6e6) tit "m@_{H_2}^2(M_S)/1.6e+06 GeV^2" w l lw 4
 
set output "inv2.eps"
plot "oneScan" u 1:($9) tit "h_t(M_Z)" w l lw 4
 


!gs -q -sPAPERSIZE=a4 -dNOPAUSE -dBATCH -sDEVICE=pswrite -sOutputFile=concat.ps *.eps
!mpage -2 concat.ps > b.ps
!ps2pdf b.ps 
!rm *.ps



