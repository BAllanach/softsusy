# WARNING: get the grid size EXACTY right or there'll be trouble!
set term post enhanced color "Helvetica" 24
set style line 2 lt 4 lw 4
set size square
min(a,b) = (a < b) ? a : b

set pm3d corners2color c1
set view map

set dgrid3d 21,21,8

# should be max number of solutions+1 (for zero solutions)

min(a,b) = (a<b) ? a : b

set label 1 "Allanach, Bednyakov, de Ruiz Austri, 2013" at 0,1.03 font \
"Helvetica,12" 
set xlabel "m_0/TeV" 
set ylabel "M_{1/2}/TeV" 

# mu<0 tb=10
set title "{/Symbol m}>0, tan{/Symbol b}=30, A_0=-2m_0"

set label 2 "{/Symbol D}_{m_h}/m_h" at 6,1.06
set output "atlasScanMh.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):6 notit w pm3d

min(a, b) = a<b ? a : b

#set output "omega50.eps"
#splot "dat.out" u ($1/1000):($2/1000):(min($6,100)) notit w pm3d

set label 2 "{/Symbol D}_{m_A}/m_A" at 6,1.06
set output "atlasScanMA.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):8 notit w pm3d

set label 2 "{/Symbol D}_{m_H}/m_H" at 6,1.06
set output "atlasScanMH.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):10 notit w pm3d

set label 2 "{/Symbol D}_{m_{H^+}}/m_{H^+}" at 6,1.06
set output "atlasScanMHp.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):12 notit w pm3d

set label 2 "{/Symbol D}_{m_g}/m_g" at 6,1.08
set output "atlasScanMg.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):14 notit w pm3d

set label 2 "{/Symbol D}_{m_q}/m_q" at 6,1.08
set output "atlasScanMq.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):16 notit w pm3d

set label 2 "{/Symbol D}_{m_{eL}}/m_{eL}" at 6,1.06
set output "atlasScanMeL.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):18 notit w pm3d

set label 2 "{/Symbol D}_{m_{eR}}/m_{eR}" at 6,1.06
set output "atlasScanMeR.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):20 notit w pm3d

set label 2 "{/Symbol D}_{m_{{/Symbol c}@_1^0}}/m_{{/Symbol c}@_1^0}" at 6,1.08
set output "atlasScanMneut1.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):22 notit w pm3d

set label 2 "{/Symbol D}_{m_{{/Symbol c}@_2^0}}/m_{{/Symbol c}@_2^0}" at 6,1.08
set output "atlasScanMneut2.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):24 notit w pm3d

set label 2 "{/Symbol D}_{m_{{/Symbol c}@_3^0}}/m_{{/Symbol c}@_3^0}" at 6,1.08
set output "atlasScanMneut3.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):26 notit w pm3d

set label 2 "{/Symbol D}_{m_{{/Symbol c}@_4^0}}/m_{{/Symbol c}@_4^0}" at 6,1.08
set output "atlasScanMneut4.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):28 notit w pm3d

set label 2 "{/Symbol D}_{m_{tL}}/m_{tL}" at 6,1.08
set output "atlasScanMtL.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):30 notit w pm3d

set label 2 "{/Symbol D}_{m_{tR}}/m_{tR}" at 6,1.08
set output "atlasScanMtR.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):32 notit w pm3d

set label 2 "{/Symbol D}_{m_{bL}}/m_{bL}" at 6,1.08
set output "atlasScanMbL.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):34 notit w pm3d

set label 2 "{/Symbol D}_{m_{bR}}/m_{bR}" at 6,1.08
set output "atlasScanMbR.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):36 notit w pm3d

set label 2 "{/Symbol D}_{m_{tauL}}/m_{tauL}" at 6,1.08
set output "atlasScanMtauL.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):38 notit w pm3d

set label 2 "{/Symbol D}_{m_{tauR}}/m_{tauR}" at 6,1.08
set output "atlasScanMtauR.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):40 notit w pm3d

set label 2 "{/Symbol D}_{/Symbol m}/{/Symbol m}" at 6,1.08
set output "atlasScanMu.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):48 notit w pm3d

set label 2 "{/Symbol D}_{m_{{/Symbol c}@_1^+}}/m_{{/Symbol c}@_1^+}" at 6,1.08
set output "atlasScanMch1.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):50 notit w pm3d

set label 2 "{/Symbol D}_{m_{{/Symbol c}@_2^+}}/m_{{/Symbol c}@_2^+}" at 6,1.08
set output "atlasScanMch2.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):52 notit w pm3d

set label 2 "{/Symbol D}{/Symbol W}_{CDM}h^2" at 6,1.08
set title "{/Symbol m}>0, tan{/Symbol b}=50, A_0=0"
set output "hiTbScanOm.eps"
splot "hiTb.dat" u ($1/1000):($2/1000):($53-$54) notit w pm3d

unset title
unset xlabel
unset ylabel
unset xtics
unset ytics
set xrange [0:6]
set yrange [0.1:1]

set label 1 "123" at 1,0.5 tc rgb "white" font "Helvetica, 30" rotate by 280
set label 2 "125" at 1.9,0.5 tc rgb "white" font "Helvetica, 30" rotate by 280
set label 3 "126" at 3,0.5 tc rgb "white" font "Helvetica, 30" rotate by 280
set label 4 "128" at 4.8,0.5 tc rgb "white" font "Helvetica, 30" rotate by 280
set output "../atlasScanMh2.eps"
plot "mh.data" notit w l lw 4 lc rgb "white"

set label 1 "1000" at 0.95,0.43tc rgb "white" font "Helvetica, 30" rotate by 280
set label 2 "2000" at 2.05,0.43 tc rgb "white" font "Helvetica, 30" rotate by 280
set label 3 "3000" at 3.1,0.43 tc rgb "white" font "Helvetica, 30" rotate by 280
set label 4 "4000" at 4.2,0.43 tc rgb "white" font "Helvetica, 30" rotate by 280
set output "../atlasScanMA2.eps"
plot "ma.data" notit w l lw 4 lc rgb "white"

set label 1 "1000" at 3,0.4 tc rgb "white" font "Helvetica, 30" rotate by -5
set label 2 "2000" at 3,0.87 tc rgb "white" font "Helvetica, 30" rotate by -5
unset label 3
unset label 4
set output "../atlasScanMg2.eps"
plot "mg.data" notit w l lw 4 lc rgb "white"

set label 1 "100" at 3,0.2 tc rgb "white" font "Helvetica, 30" rotate by -5
set label 2 "200" at 3,0.43 tc rgb "white" font "Helvetica, 30" rotate by -5
set label 3 "300" at 3,0.67 tc rgb "white" font "Helvetica, 30" rotate by -5
set label 4 "400" at 3,0.93 tc rgb "white" font "Helvetica, 30" rotate by -5
set output "../atlasScanMneut12.eps"
plot "mneut1.data" notit w l lw 4 lc rgb "white"

set label 1 "600" at 1.9,0.47tc rgb "white" font "Helvetica, 30" rotate by 320
set label 2 "1000" at 3.1,0.63 tc rgb "white" font "Helvetica, 30" rotate by 310
set label 3 "1400" at 5.2,0.8 tc rgb "white" font "Helvetica, 30" rotate by 300
unset label 4 
set output "../atlasScanMtR2.eps"
plot "mtR.data" notit w l lw 4 lc rgb "white"

set label 1 "1000" at 0.95,0.43tc rgb "white" font "Helvetica, 30" rotate by 280
set label 2 "2000" at 2.05,0.43 tc rgb "white" font "Helvetica, 30" rotate by 280
set label 3 "3000" at 3.1,0.43 tc rgb "white" font "Helvetica, 30" rotate by 280
set label 4 "4000" at 4.2,0.43 tc rgb "white" font "Helvetica, 30" rotate by 280
set output "../atlasScanMq2.eps"
plot "mq.data" notit w l lw 4 lc rgb "white"

set label 1 "0.09" at 0.14,0.43tc rgb "green" font "Helvetica, 30" rotate by 280
set label 2 "0.11" at 1.4,0.3 tc rgb "green" font "Helvetica, 30" rotate by 70
set label 3 ""
set label 4 ""
set output "../hiTbScanOm2.eps"
plot "omtb.data" notit w l lw 4 lc 2, \
"ewsb3.data" notit w filledcurve y1=0.12 lc 0, \
"ewsb2.data" notit w filledcurve y1=0.12 lc 9

unset label 1
unset label 2
set output "../atlasExcl.eps"
plot "ATLAS_SUSY_MSUGRA_201308.txt" u ($1/1000):($2/1000) notit w l lw 16 lc 2 \
lt 3
