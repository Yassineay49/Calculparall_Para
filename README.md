# Rapport du Projet Calcul Parallèle

## Étude et mise en œuvre des méthodes de partitionnement et de décomposition de domaine pour le calcul parallèle

### Introduction du problème

La résolution de l'équation de conduction instationnaire est un problème fondamental en physique et en mathématiques appliquées, représentant la diffusion de la chaleur dans un domaine donné au fil du temps. Dans cette étude, nous considérons l'équation suivante définie sur le domaine [0, Lₓ] × [0, Lᵧ] en ℝ² :

∂ₜu(x, y, t) - D∆u(x, y, t) = f(x, y, t),
u|ᵧ₀ = g(x, y, t),
u|ᵧ₁ = h(x, y, t),

où u(x, y, t) représente la température en un point (x, y) du domaine à l'instant t, D est le coefficient de diffusion thermique, et f(x, y, t) est une source thermique.

Pour valider la résolution numérique de cette équation, plusieurs cas tests seront utilisés en fixant Lₓ = Lᵧ = 1 et D = 1. Ces cas incluent des solutions stationnaires et instationnaires, afin de vérifier la précision et la stabilité des solutions calculées :

#### Solutions stationnaires
* **Cas test 1:** f = 2 × (x - x² + y - y²), g = 0, et h = 0, la solution reste constante aux frontières.
* **Cas test 2:** f = sin(x) + cos(y), g = sin(x) + cos(y), et h = sin(x) + cos(y), la solution respecte une dépendance sinusoïdale en bordure du domaine.

#### Solution instationnaire périodique
* **Cas test 3:** f = e^(-(x - Lₓ/2)²) e^(-(y - Lᵧ/2)²) cos(πt/2), g = 0, et h = 1.

## Calculparall_Para

### Compilation et exécution

#### Pour générer le Makefile :
```bash
cmake CMakeLists.txt
```

Pour compiler le code et l'exécuter :

```bash
make
mpirun -n 6 ./CHP
```
Ou si le chemin de mpirun est différent :

```bash
/usr/lib64/openmpi/bin/mpirun -n 6 ./CHP
```
Modifiez 6 par le nombre de processeurs que vous souhaitez activer
### Pour regrouper les résultats :

```bash
g++ -o run settings.cpp Concatrenation.cpp
./run
```

### Pour afficher les résultats dans gnuplot :
```bash 
load "plot.txt"
```
### Pour tracer les résultats de chaque processeur avec une couleur différente :
```bash
load "plot_solutions.gp"
```

Attention : Le fichier plot_solutions.gp doit être modifié manuellement en fonction du nombre de processeurs utilisés.

### Configuration
Les paramètres du problème, tels que : Nx, Ny, nr, tf, xmax, ymax et d'autres configurations, sont définis dans le fichier data.txt. Vous pouvez modifier ce fichier pour ajuster les paramètres selon vos besoins avant d'exécuter le programme.
