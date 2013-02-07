set term post enhanced color "Helvetica" 24
set xlabel "{/Symbol m}(M_{SUSY})/GeV"
set ylabel "M_Z^2(pred)/M_Z^2(exp)"
set yrange [-40:40]
set title "m_0=350+(k*500) GeV, M_{1/2}=300 GeV, tan{/Symbol b}=10, A_0=0"

set output "multi.eps"
plot "cmssmScans" i 0 tit "default" w l, 1 notit w l ls 1 lc 2

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

