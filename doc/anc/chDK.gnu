# gnuplot 5.0 patchlevel 3
set term post color enhanced "Helvetica" 24
set size square

set key spacing 1.3
set output "chiDK.eps"
set logscale x
set for [i=1:9] linetype i dt i
set xlabel "{/Symbol D}m/GeV"
set ylabel "{/Symbol c}@^+_1 branching ratio"
#set ytics (0,0.2,0.4,0.6,0.8,1)
set label 1 at 0.15,0.07 "{/Symbol n}_e{e}^+{/Symbol c}@_1^0" textcolor rgb "red"
set label 2 at 1.8,0.2 "{/Symbol n}_{/Symbol m}{/Symbol m}^+{/Symbol c}@_1^0" textcolor rgb "blue"
set label 3 at 1.8,0.07 "{/Symbol n}_{/Symbol t}{/Symbol t}^+{/Symbol c}@_1^0" textcolor rgb "magenta"
set label 4 at 0.15,0.92 "{/Symbol p}^+{/Symbol c}@_1^0" textcolor rgb "black"
set label 5 at 0.4,0.25 "{/Symbol p}^+{/Symbol p}^0{/Symbol c}@_1^0" textcolor rgb "green"
set label 6 at 1.8,0.65 "{/Symbol c}@_1^0jj" textcolor rgb "orange"
plot [0.1:10] [0:1] "decay.out" u 3:5 notit w l lt 1 \
lc rgb "black" lw 4, \
"decay.out" u 3:6 notit w l lt 2 \
lc rgb "green" lw 4,  \
"decay.out" u 3:7 notit w l lt 3 \
lc rgb "red" lw 4, \
"decay.out" u 3:8 notit w\
 l lt 4 lc rgb "blue" lw 4, \
"decay.out" u 3:9 notit \
w l lt 5 lc rgb "magenta" lw 4, \
"decay.out" u 3:4 notit w l lt 6 lc rgb "orange" lw 4
unset ytics
unset label 1
unset label 2
unset label 3
unset label 4
unset label 5

set ytics
set output "chiLifetime.eps"
set ylabel "log_{10}({/Symobl t}/s)"
plot "decay.out" u 3:(log10(0.658e-24/$10)) notit w l lt 1 lc rgb "black" lw 4

!epstopdf --autorotate=All chiDK.eps; epstopdf --autorotate=All chiLifetime.eps; rm chiDK.eps chiLifetime.eps
