set dgrid3d 20,19,8
set contour 
set nosurface
set cntrparam bspline 
set noclip

set xrange [0:6]
set yrange [0:1]

set table "mh.data"
set cntrparam levels discrete 123,125,126,128
splot "atlas_scan.dat" u ($1/1000):($2/1000):5 

#set xrange [0:3.5]
#set yrange [0.1:1]
set table "omtb.data"
set cntrparam levels discrete 0.10,0.15
splot "hiTb.dat" u ($1/1000):($2/1000):129

set table "ewsb2.data"
set cntrparam levels discrete 0
splot "hiTb.dat" u ($1/1000):($2/1000):110

set table "ewsb3.data"
set cntrparam levels discrete 0
splot "hiTb.dat" u ($1/1000):($2/1000):114

set xrange [0:6]

set table "ma.data"
set cntrparam levels discrete 1000,2000,3000,4000
splot "atlas_scan.dat" u ($1/1000):($2/1000):14

set table "mg.data"
set cntrparam levels discrete 1000,2000,3000,4000
splot "atlas_scan.dat" u ($1/1000):($2/1000):29

set table "mneut1.data"
set cntrparam levels discrete 100,200,300,400
splot "atlas_scan.dat" u ($1/1000):($2/1000):49

set table "mtR.data"
set cntrparam levels discrete 600,1000,1400
splot "atlas_scan.dat" u ($1/1000):($2/1000):74

set table "mq.data"
set cntrparam levels discrete 1000,2000,3000,4000
splot "atlas_scan.dat" u ($1/1000):($2/1000):34
