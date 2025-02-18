# Tests du Solveur Parallèle pour Équations de Chaleur

Ce fichier décrit les tests unitaires utilisés pour valider le solveur parallèle d'équations de chaleur implémenté avec MPI.

## Structure des Tests

Les tests utilisent le framework Google Test (gtest) avec une configuration personnalisée pour l'affichage des résultats.
Une classe `TestInfo` a été créée pour fournir des informations détaillées sur l'exécution de chaque test.

## Tests Implémentés

### Test1_SizeProc
**Description**: Test de la fonction size_proc pour le calcul de la taille des domaines
- Vérifie la fonctionnalité de base de la fonction `size_proc`
- Teste différentes configurations de domaines et de processus
- S'assure que les tailles calculées sont cohérentes et positives

### Test2_VectorOperations
**Description**: Test des opérations vectorielles (somme et produit scalaire)
- Teste la fonction `somme_vecteur` qui additionne deux vecteurs
- Teste la fonction `produit_scalaire` qui calcule le produit scalaire de deux vecteurs
- Vérifie la précision numérique des résultats

### Test3_ErrorCalculation
**Description**: Test de la fonction de calcul d'erreur
- Vérifie la fonction `erreur` qui calcule la différence relative entre deux vecteurs
- Assure la précision du calcul d'erreur pour des valeurs connues

### Test4_MatrixVectorProduct
**Description**: Test du produit matrice-vecteur séquentiel
- Teste la fonction `produitmatvectSEQ` qui calcule le produit d'une matrice tridiagonale avec un vecteur
- Vérifie que le résultat a la dimension attendue

### Test5_BICGstab
**Description**: Test du solveur BICGstab
- Teste l'algorithme BiCGStab en version séquentielle
- Vérifie que la solution a la dimension attendue

### Test6_ParallelOperations
**Description**: Test des opérations parallèles avec MPI
- Initialise l'environnement MPI pour tester les fonctions parallèles
- Teste la fonction `produitmatvect1` en parallèle
- Vérifie que chaque processus traite une partie non vide du domaine

### Test7_BoundaryConditions
**Description**: Test des conditions aux limites de Robin
- Teste l'implémentation des conditions aux limites de Robin dans le produit matrice-vecteur
- Vérifie que les conditions aux limites modifient correctement les valeurs aux bords

### Test9_Functions
**Description**: Test des fonctions source et conditions aux limites
- Vérifie les fonctions `f`, `g` et `h` pour différents cas
- Teste le cas stationnaire (cas 1)
- Teste le cas sinusoïdal (cas 2)
- Teste le cas de la gaussienne dépendante du temps (cas 3)
- Vérifie la précision numérique des calculs

### Test10_DomainDecomposition
**Description**: Test de la décomposition de domaine et équilibrage de charge
- Teste la fonction `charge_a` qui calcule la distribution du domaine entre processus
- Vérifie que les sous-domaines couvrent exactement le domaine global
- Teste l'équilibrage de charge pour différents nombres de processus
- Vérifie l'absence de chevauchement ou de trous entre sous-domaines adjacents

### Test11_BICGstab
**Description**: Test de la convergence du solveur BICGstab
- Évalue le taux de convergence du solveur BiCGStab
- Vérifie que le résidu diminue avec l'augmentation du nombre d'itérations
- Teste pour différentes valeurs du paramètre kmax (10, 100, 1000)
- Prend en compte la précision machine pour éviter les faux négatifs

### Test12_MPICommunication
**Description**: Test de la communication MPI
- Teste l'envoi et la réception de messages entre processus
- Vérifie l'exactitude des données transmises
- S'assure que la synchronisation entre processus fonctionne correctement
- Utilise MPI_Send, MPI_Recv et MPI_Barrier

### Test13_RobinBoundary
**Description**: Test des conditions aux limites de Robin
- Teste l'implémentation des conditions aux limites de Robin dans le produit matrice-vecteur
- Compare les résultats avec et sans conditions de Robin
- Vérifie que les conditions de Robin modifient effectivement les calculs

### Test14_FileIO
**Description**: Test d'entrée/sortie fichier
- Teste la fonction `sauvegarder_resultats` pour l'écriture des résultats
- Vérifie la création et l'accès aux fichiers
- Teste la lecture des données sauvegardées
- Gère les différents chemins d'accès possibles selon l'environnement d'exécution

### Test15_NumericalStability
**Description**: Test de stabilité numérique
- Vérifie la robustesse numérique des fonctions de calcul de termes source
- Teste les termes source avec le cas de la gaussienne
- S'assure que les calculs ne produisent ni NaN ni valeurs infinies
- Vérifie que les dimensions des résultats sont correctes

## Utilisation

Les tests peuvent être exécutés en parallèle avec MPI :
```bash
mpirun -np [nombre_de_processus] ./run_tests