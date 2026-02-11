# Live 2D colormap of FT2.dat

# Use an interactive terminal (wxt or qt)
set terminal wxt size 800,600 enhanced font 'Arial,14'
# OR: set terminal qt

# Title and labels
set title "2D Colormap of ew 1.5.dat"
set xlabel "z"
set ylabel "x"

set size ratio 1

# Color map settings
set pm3d map                # top-down view
set palette rgb 33,13,10    # smooth palette
set colorbox                # show color scale

# Optional: set axis and color scale ranges
# set xrange [0.75:]
# set yrange [40:60]
#set cbrange [-0.5:1.8]

# Plot
splot '../data/newrho2DNeps1.000000eps_i1.403153ew1.000000dx0.025000rhob_solvent_l0.603936rhob_solute_l0.019845lambdaB5.000000.dat' using ($1):($2):6 with pm3d notitle
pause -1
