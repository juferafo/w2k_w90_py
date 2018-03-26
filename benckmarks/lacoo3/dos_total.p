# Gnuplot script file for plotting data in file "case.dosievup"
# This file is called   dos.p
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "DOS of LaCoO3"
set xlabel "ENERGY [eV]"
set ylabel "DOS"
set zeroaxis
plot    "lacoo3.dos1evup" using 1:2 title 'total DOS' with line, \
        "lacoo3.dos1evup" using 1:3 title 'total d Co' with line, \
        "lacoo3.dos1evup" using 1:4 title 'total La' with line
set term png
set output "dos.png"
replot
