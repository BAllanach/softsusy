set dgrid3d 21,21,8
set contour 
set nosurface
set cntrparam bspline 
set noclip

set table "mh.data"
set cntrparam levels discrete 123,125,126,128
splot "atlas_scan.dat" u ($1/1000):($2/1000):5 

set table "omtb.data"
set cntrparam levels discrete 0.09,0.13
splot "hiTb.dat" u ($1/1000):($2/1000):53

set table "ewsb2.data"
set cntrparam levels discrete 0
splot "hiTb.dat" u ($1/1000):($2/1000):48

set table "ewsb3.data"
set cntrparam levels discrete 0
splot "hiTb.dat" u ($1/1000):($2/1000):47

set table "ma.data"
set cntrparam levels discrete 1000,2000,3000,4000
splot "atlas_scan.dat" u ($1/1000):($2/1000):7

set table "mg.data"
set cntrparam levels discrete 1000,2000,3000,4000
splot "atlas_scan.dat" u ($1/1000):($2/1000):13

set table "mneut1.data"
set cntrparam levels discrete 100,200,300,400
splot "atlas_scan.dat" u ($1/1000):($2/1000):21

set table "mtR.data"
set cntrparam levels discrete 600,1000,1400
splot "atlas_scan.dat" u ($1/1000):($2/1000):31

set table "mq.data"
set cntrparam levels discrete 1000,2000,3000,4000
splot "atlas_scan.dat" u ($1/1000):($2/1000):15
