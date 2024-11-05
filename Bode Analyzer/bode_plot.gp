set terminal png size 900,600
set output 'bode_plot.png'
set logscale x
set grid
set xlabel 'Frequency (Hz)'
set ylabel 'Amplitude (dB)'
set y2label 'Phase (Degrees)'
set ytics nomirror
set y2tics
set lmargin at screen 0.1
set rmargin at screen 0.9
set tmargin at screen 0.8
set bmargin at screen 0.1
set label 1 'G(s) = 5 / 5s^2 + 3s^1 + 3'at screen 0.25, 0.9 center
set label 2 'Phasenreserve (Grad) = 43.6058'at screen 0.6, 0.93 left
set label 3 'Amplitudenreserve (dB) = 0'at screen 0.6, 0.89 left
plot 'bode_plot_data.csv' using 1:2 with lines title 'Amplitude (dB)', \
     'bode_plot_data.csv' using 1:3 axes x1y2 with lines title 'Phase (Degrees)'
