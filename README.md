# Calculparall_Para

Pour génerer le Makefile :

cmake CMakeLists.txt 
------------------------------
Pour compiler le code et l'exécuter : 

make
mpirun -n  6 ./CHP  ( ou /usr/lib64/openmpi/bin/mpirun -n 4 ./CHP si le chemin de mpirun est différent)

Modifiez 6 par le nombre de processeurs que vous souhaitez activer
--------------------------------------

Pour regrouper les résultats :

g++ -o run settings.cpp Concatrenation.cpp
./run

"Puis entrez le nombre de procs activé après la demande (par exemple 6)"

--------------------------------------
Pour afficher les résultats dans gnuplot :

load "plot.txt"

--------------------------------------
Si vous souhaitez tracer les résultats de chaque processeur avec une couleur différente, utilisez :

load "plot_solutions.gp"

Attention : Le fichier plot_solutions.gp doit être modifié manuellement en fonction du nombre de processeurs utilisés. 

--------------------------------------

Les paramètres du problème, tels que : Nx, Ny, nr, tf, xmax, ymax  d'autres configurations, sont définisdans le fichier data.txt.
Vous pouvez modifier ce fichier pour ajuster les paramètres selon vos besoins avant d'exécuter le programme.

