# WARNING: get the grid size EXACTY right or there'll be trouble!
set term post enhanced color "Helvetica" 24
set style line 2 lt 4 lw 4
set size square
min(a,b) = (a < b) ? a : b

set pm3d corners2color c1
set view map

set dgrid3d 20,19,8

# should be max number of solutions+1 (for zero solutions)

min(a,b) = (a<b) ? a : b

set label 1 "Allanach, Bednyakov, Ruiz de Austri, 2014" at 0,1.03 font \
"Helvetica,12" 
set xlabel "m_0/TeV" 
set ylabel "M_{1/2}/TeV" 

# mu<0 tb=10
set title "{/Symbol m}>0, tan{/Symbol b}=30, A_0=-2m_0"


set label 2 "{/Symbol D}_{m_h}/GeV" at 6,1.06
set output "atlasScanMh.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):($9-$5) notit w pm3d

min(a, b) = a<b ? a : b

set label 2 "{/Symbol D}_{m_A}/m_A" at 6,1.08
set output "atlasScanMA.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):(1-$10/$14) notit w pm3d

set label 2 "{/Symbol D}_{m_g}/m_g" at 6,1.08
set output "atlasScanMg.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):(1-$25/$29) notit w pm3d

set label 2 "{/Symbol D}_{m_q}/m_q" at 6,1.08
set output "atlasScanMq.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):(1-$30/$34) notit w pm3d

set label 2 "{/Symbol D}_{m_{{/Symbol c}@_1^0}}/m_{{/Symbol c}@_1^0}" at 6,1.08
set output "atlasScanMneut1.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):(1-$45/$49) notit w pm3d

set label 2 "{/Symbol D}_{m_{{/Symbol c}@_3^0}}/m_{{/Symbol c}@_3^0}" at 6,1.08
set output "atlasScanMneut3.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):26 notit w pm3d

set label 2 "{/Symbol D}_{m_{tR}}/m_{tR}" at 6,1.08
set output "atlasScanMtR.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):(1-$70/$74) notit w pm3d

unset cbrange 
set label 2 "{/Symbol Ds}/{/Symbol s}" at 6,1.08
set output "atlasScanDs.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):(($139-$135)/$139) notit w pm3d


unset cbrange
set label 2 "{/Symbol D}{/Symbol W}_{CDM}h^2" at 6,1.08
set title "{/Symbol m}>0, tan{/Symbol b}=50, A_0=0"
set output "hiTbScanOm.eps"
splot "hiTb.dat" u ($1/1000):($2/1000):($53-$54) notit w pm3d
unset label 2

set xrange [-4:4]
set yrange [2:59]
set pm3d corners2color c1
set dgrid3d 21,20,8
set xlabel "A_0/TeV"
set ylabel "tan{/Symbol b}"
set label 1 "Allanach, Bednyakov, Ruiz de Austri, 2014" at -4,62 font "Helvetica,12" 
set title "{/Symbol m}>0, m_0=2 TeV, M_{1/2}=600 GeV"


set dgrid3d 21,21,8

set label 2 "{/Symbol D}Y@_{bt}^{(All)}-{/Symbol D}Y@_{bt}^{(None)}" at 4.7,65
set output "atlasScanYpNone.eps"
splot "tbScan.dat" u ($3/1000):($4):($155-$151) notit w pm3d

set label 2 "{/Symbol D}Y@_{bt}^{(All)}" at 4.7,65
set output "atlasScanYpAll.eps"
splot "tbScan.dat" u ($3/1000):4:155 notit w pm3d



#set cbrange [-0.2:0.6]
set label 2 "{/Symbol D}Y@_{b{/Symbol t}}^{(All)}-{/Symbol D}Y@_{b{/Symbol t}}^{(None)}" at 4.7,65
set output "atlasScanYbNone.eps"
splot "tbScan.dat" u ($3/1000):4:($150-$146) notit w pm3d

set label 2 "{/Symbol D}Y@_{b{/Symbol t}}^{(All)}" at 4.7,65
set output "atlasScanYbAll.eps"
splot "tbScan.dat" u ($3/1000):4:150 notit w pm3d



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
set output "atlasScanMh2.eps"
plot "mh.data" notit w l lw 4 lc rgb "white"

set label 1 "1000" at 0.95,0.43tc rgb "white" font "Helvetica, 30" rotate by 280
set label 2 "2000" at 2.05,0.43 tc rgb "white" font "Helvetica, 30" rotate by 280
set label 3 "3000" at 3.1,0.43 tc rgb "white" font "Helvetica, 30" rotate by 280
set label 4 "4000" at 4.2,0.43 tc rgb "white" font "Helvetica, 30" rotate by 280
set output "atlasScanMA2.eps"
plot "ma.data" notit w l lw 4 lc rgb "white"

set label 1 "1000" at 3,0.4 tc rgb "white" font "Helvetica, 30" rotate by -5
set label 2 "2000" at 3,0.87 tc rgb "white" font "Helvetica, 30" rotate by -5
unset label 3
unset label 4
set output "atlasScanMg2.eps"
plot "mg.data" notit w l lw 4 lc rgb "white"

set label 1 "100" at 3,0.2 tc rgb "white" font "Helvetica, 30" rotate by -5
set label 2 "200" at 3,0.43 tc rgb "white" font "Helvetica, 30" rotate by -5
set label 3 "300" at 3,0.67 tc rgb "white" font "Helvetica, 30" rotate by -5
set label 4 "400" at 3,0.93 tc rgb "white" font "Helvetica, 30" rotate by -5
set output "atlasScanMneut12.eps"
plot "mneut1.data" notit w l lw 4 lc rgb "white"

set label 1 "600" at 1.9,0.47tc rgb "white" font "Helvetica, 30" rotate by 320
set label 2 "1000" at 3.1,0.63 tc rgb "white" font "Helvetica, 30" rotate by 310
set label 3 "1400" at 5.2,0.8 tc rgb "white" font "Helvetica, 30" rotate by 300
unset label 4 
set output "atlasScanMtR2.eps"
plot "mtR.data" notit w l lw 4 lc rgb "white"

set label 1 "1000" at 0.95,0.43tc rgb "white" font "Helvetica, 30" rotate by 280
set label 2 "2000" at 2.05,0.43 tc rgb "white" font "Helvetica, 30" rotate by 280
set label 3 "3000" at 3.1,0.43 tc rgb "white" font "Helvetica, 30" rotate by 280
set label 4 "4000" at 4.2,0.43 tc rgb "white" font "Helvetica, 30" rotate by 280
set output "atlasScanMq2.eps"
plot "mq.data" notit w l lw 4 lc rgb "white"

set label 1 "0.10" at 0.14,0.43tc rgb "green" font "Helvetica, 30" rotate by 280
set label 2 "0.15" at 1.4,0.3 tc rgb "green" font "Helvetica, 30" rotate by 70
set label 3 ""
set label 4 ""
set output "hiTbScanOm2.eps"
plot "omtb.data" notit w l lw 4 lc 2, \
"ewsb3.data" notit w filledcurve y1=0.12 lc 0, \
"ewsb2.data" notit w filledcurve y1=0.12 lc 9

unset label 1
unset label 2
set output "atlasExcl.eps"
plot "ATLAS_SUSY_MSUGRA_201308.txt" u ($1/1000):($2/1000) notit w l lw 16 lc 2 \
lt 3
set autoscale

set label 1 "Allanach, Bednyakov, Ruiz de Austri, 2014" at 0,0.01 font \
"Helvetica,12" 
set output "tbScanDy.eps"
set title "{/Symbol m}>0, m_0=3 TeV, M_{1/2}=2 TeV, A_0=-6 TeV"
set xlabel "tan{/Symbol b}"
set ylabel "Y_b(M_{GUT})-Y_{/Symbol t}(M_{GUT})"
set xtics 
set ytics
plot "tbScan1d.dat" u 4:146 tit "None" w l lc 3 lw 4, \
"tbScan1d.dat" u 4:150 tit "{/Symbol D}All" w l lw 4 lc 1

set label 1 "Allanach, Bednyakov, Ruiz de Austri, 2014" at 0,0.82 font \
"Helvetica,12" 
set output "tbScanDyt.eps"
set ylabel "Y_t(M_{GUT})-Y_b(M_{GUT})"
plot "tbScan1d.dat" u 4:151 tit "None" w l lc 3 lw 4, \
"tbScan1d.dat" u 4:155 tit "{/Symbol D}All" w l lw 4 lc 1

set autoscale
set label 1 "Allanach, Bednyakov, Ruiz de Austri, 2014" at 0,0.00024 font \
"Helvetica,12" 
set key bottom
set title "{/Symbol m}>0, tan{/Symbol b}=30, A_0=-2m_0"
set output "mScanAs.eps"
set ylabel "{/Symbol a}_3(M_{GUT})-{/Symbol a}_1(M_{GUT})"
set xlabel "m_0/TeV=M_{1/2}/TeV"
plot "massScan.dat" u ($1/1000):141 tit "None" w l lc 3 lw 4, \
"massScan.dat" u ($1/1000):145 tit "{/Symbol D}All" w l lw 4 lc 1, \
0 notit w l lc 0

set cbrange [-0.0024:0]
set label 2 "{/Symbol Da}(None)" at 6,1.06
set output "atlasScanDaNone.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):141 notit w pm3d

set label 2 "{/Symbol Da}(All)" at 6,1.06
set output "atlasScanDaAll.eps"
splot "atlas_scan.dat" u ($1/1000):($2/1000):145 notit w pm3d
