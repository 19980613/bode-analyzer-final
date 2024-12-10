set terminal png size 800,500
set output 'bode_plot.png'
set logscale x
set grid
set xlabel 'Frequency (Hz)'
set ylabel 'Amplitude (dB)'
set y2label 'Phase (Degrees)'
set ytics nomirror
set y2tics
plot 'bode_plot_data.csv' using 1:2 with lines title 'Amplitude (dB)', \
     'bode_plot_data.csv' using 1:3 axes x1y2 with lines title 'Phase (Degrees)'
