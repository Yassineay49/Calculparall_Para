# Paramètres de sortie
set terminal pngcairo size 1024,768 enhanced
set output "solution_with_proc_colors.png"

# Titres et grilles
set title "Solution numérique pour cas test 3"
set xlabel "x"
set ylabel "y"
set zlabel "u"
set grid

# Désactivation de la palette
unset colorbox

# Tracer les fichiers avec des couleurs fixes
splot \
    "Sol.3.1.0.txt" using 1:2:3 with linespoints lc rgb "blue" title "Proc 0", \
    "Sol.3.1.1.txt" using 1:2:3 with linespoints lc rgb "green" title "Proc 1", \
    "Sol.3.1.2.txt" using 1:2:3 with linespoints lc rgb "red" title "Proc 2" 

    