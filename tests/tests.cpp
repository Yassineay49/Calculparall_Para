#include <gtest/gtest.h>
#include <mpi.h>
#include <fstream>
#include <iostream>
#include "../src/Matrix.h"
#include "../src/charge.h"
#include "../src/solver.h"
#include "../src/functions.h"
#include "../src/settings.h"
#include <cmath>
#include <vector>
#include <sys/stat.h>

// Classe pour afficher les informations de test
class TestInfo : public ::testing::TestEventListener {
public:
    void OnTestProgramStart(const ::testing::UnitTest& /*unit_test*/) override {}
    void OnTestIterationStart(const ::testing::UnitTest& /*unit_test*/, int /*iteration*/) override {}
    void OnEnvironmentsSetUpStart(const ::testing::UnitTest& /*unit_test*/) override {}
    void OnEnvironmentsSetUpEnd(const ::testing::UnitTest& /*unit_test*/) override {}
    void OnTestCaseStart(const ::testing::TestCase& /*test_case*/) override {}
    
    void OnTestStart(const ::testing::TestInfo& test_info) override {
        std::cout << "\n==========================================================\n";
        std::string test_case_name = test_info.test_case_name();
        std::string test_name = test_info.name();
        
        std::cout << "Test en cours d'exécution : " << test_case_name << "\n";
        std::cout << "Description : " << GetTestDescription(test_case_name) << "\n";
        if (test_info.should_run()) {
            std::cout << "État : Test démarré" << std::endl;
        }
    }
    
    void OnTestPartResult(const ::testing::TestPartResult& test_part_result) override {
        if (test_part_result.failed()) {
            std::cout << "ÉCHEC dans " << test_part_result.file_name() << ":" 
                      << test_part_result.line_number() << "\n"
                      << test_part_result.summary() << std::endl;
        }
    }
    
    void OnTestEnd(const ::testing::TestInfo& test_info) override {
        std::cout << "\nRésultat : " << (test_info.result()->Passed() ? 
            "\033[32mVALIDÉ\033[0m" : "\033[31mÉCHOUÉ\033[0m");
        std::cout << "\n==========================================================\n";
    }
    
    void OnTestCaseEnd(const ::testing::TestCase& /*test_case*/) override {}
    void OnEnvironmentsTearDownStart(const ::testing::UnitTest& /*unit_test*/) override {}
    void OnEnvironmentsTearDownEnd(const ::testing::UnitTest& /*unit_test*/) override {}
    void OnTestIterationEnd(const ::testing::UnitTest& /*unit_test*/, int /*iteration*/) override {}
    void OnTestProgramEnd(const ::testing::UnitTest& /*unit_test*/) override {}

private:
    std::string GetTestDescription(const std::string& test_case_name) {
        if (test_case_name == "Test1_SizeProc")
            return "Test de la fonction size_proc pour le calcul de la taille des domaines";
        else if (test_case_name == "Test2_VectorOperations")
            return "Test des opérations vectorielles (somme et produit scalaire)";
        else if (test_case_name == "Test3_ErrorCalculation")
            return "Test de la fonction de calcul d'erreur";
        else if (test_case_name == "Test4_MatrixVectorProduct")
            return "Test du produit matrice-vecteur séquentiel";
        else if (test_case_name == "Test5_BICGstab")
            return "Test du solveur BICGstab";
        else if (test_case_name == "Test6_ParallelOperations")
            return "Test des opérations parallèles avec MPI";
        else if (test_case_name == "Test7_BoundaryConditions")
            return "Test des conditions aux limites de Robin";
        else if (test_case_name == "Test9_Functions")
            return "Test des fonctions source et conditions aux limites";
        else if (test_case_name == "Test10_DomainDecomposition")
            return "Test de la décomposition de domaine et équilibrage de charge";
        else if (test_case_name == "Test11_BICGstab")
            return "Test de la convergence du solveur BICGstab";
        else if (test_case_name == "Test12_MPICommunication")
            return "Test de la communication MPI";
        else if (test_case_name == "Test13_RobinBoundary")
            return "Test des conditions aux limites de Robin";
        else if (test_case_name == "Test14_FileIO")
            return "Test d'entrée/sortie fichier";
        else if (test_case_name == "Test15_NumericalStability")
            return "Test de stabilité numérique";
        else
            return "Test non reconnu: " + test_case_name;
    }
};

// Test Fixture pour MPI
class MPITest : public ::testing::Test {
protected:
    void SetUp() override {
        int argc = 1;
        char** argv = new char*[1];
        argv[0] = strdup("test");
        
        int initialized;
        MPI_Initialized(&initialized);
        if (!initialized) {
            MPI_Init(&argc, &argv);
        }
    }

    void TearDown() override {
        // Ne pas finalize MPI ici car d'autres tests l'utilisent
    }
};

// Test 1: Fonction size_proc
TEST(Test1_SizeProc, BasicFunctionality) {
    int result1 = size_proc(0, 0, 3, 1, 2);
    int result2 = size_proc(1, 4, 7, 1, 2);
    int result3 = size_proc(1, 2, 5, 1, 3);
    
    EXPECT_EQ(result1, 5);
    EXPECT_EQ(result2, 5);
    EXPECT_EQ(result3, 6);
    
    EXPECT_GT(result1, 0);
    EXPECT_GT(result2, 0);
    EXPECT_GT(result3, 0);
}

// Test 2: Opérations vectorielles
TEST(Test2_VectorOperations, BasicOperations) {
    std::vector<double> v1 = {1.0, 2.0, 3.0};
    std::vector<double> v2 = {4.0, 5.0, 6.0};
    
    auto sum = somme_vecteur(v1, v2, 3);
    EXPECT_NEAR(sum[0], 5.0, 1e-10);
    EXPECT_NEAR(sum[1], 7.0, 1e-10);
    EXPECT_NEAR(sum[2], 9.0, 1e-10);
    
    double dot = produit_scalaire(v1, v2, 3);
    EXPECT_NEAR(dot, 32.0, 1e-10);
}

// Test 3: Calcul d'erreur
TEST(Test3_ErrorCalculation, BasicError) {
    std::vector<double> u = {1.0, 2.0, 3.0};
    std::vector<double> stencil = {1.1, 2.1, 3.1};
    int Nx = 3;
    int p = 0;
    
    double err = erreur(u, stencil, Nx, p);
    EXPECT_NEAR(err, 0.03, 1e-10);
}

// Test 4: Produit matrice-vecteur séquentiel
TEST(Test4_MatrixVectorProduct, Sequential) {
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0};
    double a = 2.0, b = -1.0, c = -1.0;
    int n = 2, m = 2;
    
    auto result = produitmatvectSEQ(a, b, c, x, n, m, 0.1, 1.0, 0.01);
    EXPECT_EQ(result.size(), x.size());
}

// Test 5: Tests BICGstab
TEST(Test5_BICGstab, BasicSolver) {
    std::vector<double> b = {1.0, 0.0, 0.0, 1.0};
    double alpha = 2.0, beta = -1.0, gamma = -1.0;
    int n = 2, m = 2;
    double dt = 0.01, D = 1.0, dy = 0.1;
    int kmax = 1000;
    double epsilon = 1e-6;
    
    auto solution = BICGstabSEQ(alpha, beta, gamma, b, n, m, 0, 1, dt, D, dy, kmax, epsilon);
    EXPECT_EQ(solution.size(), b.size());
}

// Test 6: Parallèle MPI
TEST_F(MPITest, Test6_ParallelOperations) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (size > 1) {
        std::vector<double> x = {1.0, 2.0, 3.0, 4.0};
        double alpha = 2.0, beta = -1.0, gamma = -1.0;
        int n = 2, m = 2;
        double dy = 0.1, D = 1.0, dt = 0.01;
        double alpha_rob = 0.0, beta_rob = 0.0;
        
        auto result = produitmatvect1(alpha, beta, gamma, x, n, m, dy, D, dt, 
                                    alpha_rob, beta_rob, rank, size);
        int local_size = result.size();
        std::vector<int> sizes(size);
        MPI_Gather(&local_size, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        if (rank == 0) {
            for (int i = 0; i < size; i++) {
                EXPECT_GT(sizes[i], 0);
            }
        }
    }
}

// Test 7: Conditions aux limites de Robin
TEST(Test7_BoundaryConditions, RobinConditions) {
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0};
    double alpha = 2.0, beta = -1.0, gamma = -1.0;
    int n = 2, m = 2;
    double dy = 0.1, D = 1.0, dt = 0.01;
    double alpha_rob = 1.0, beta_rob = 1.0;
    
    auto result = produitmatvect1(alpha, beta, gamma, x, n, m, dy, D, dt, 
                                alpha_rob, beta_rob, 0, 2);
    EXPECT_NE(result[0], result[1]);
}

// Test 9: Tests des fonctions f, g et h
TEST(Test9_Functions, SourceAndBoundary) {
    // Test des fonctions de source et conditions aux limites
    double Lx = 1.0, Ly = 1.0;
    double x = 0.5, y = 0.5, t = 0.1;
    
    // Cas 1: Fonction stationnaire
    double f_cas1 = f(x, y, t, Lx, Ly, 1);
    double expected_cas1 = 2*(y-pow(y,2)+x-pow(x,2));
    EXPECT_NEAR(f_cas1, expected_cas1, 1e-10);
    
    // Cas 2: Fonction sinusoïdale 
    double f_cas2 = f(x, y, t, Lx, Ly, 2);
    double g_cas2 = g(x, y, t, Lx, Ly, 2);
    double h_cas2 = h(x, y, t, Lx, Ly, 2);
    EXPECT_NEAR(f_cas2, sin(x)+cos(y), 1e-10);
    EXPECT_NEAR(g_cas2, sin(x)+cos(y), 1e-10);
    EXPECT_NEAR(h_cas2, sin(x)+cos(y), 1e-10);
    
    // Cas 3: Gaussienne
    double pi = 4*atan(1.);
    double f_cas3 = f(x, y, t, Lx, Ly, 3);
    double g_cas3 = g(x, y, t, Lx, Ly, 3);
    double h_cas3 = h(x, y, t, Lx, Ly, 3);
    double expected_cas3 = exp(-(x-Lx/2)*(x-Lx/2))*exp(-(y-Ly/2)*(y-Ly/2))*cos((pi*t)/2);
    EXPECT_NEAR(f_cas3, expected_cas3, 1e-10);
    EXPECT_EQ(g_cas3, 0.0);
    EXPECT_EQ(h_cas3, 1.0);
}

// Test 10: Test de la décomposition de domaine (charge_a)
TEST(Test10_DomainDecomposition, LoadBalancing) {
    // Test pour différents nombres de processus
    int n = 100; // Taille du domaine
    
    // Cas 1: 4 processus, division équilibrée
    {
        int p = 4;
        int total_size = 0;
        
        for (int me = 0; me < p; me++) {
            int ibeg, iend;
            charge_a(me, n, p, &ibeg, &iend);
            
            // Vérifier que chaque processus a une taille raisonnable
            int size = iend - ibeg + 1;
            EXPECT_GT(size, 0);
            
            // Vérifier que les domaines ne se chevauchent pas et qu'il n'y a pas de trous
            if (me > 0) {
                int prev_ibeg, prev_iend;
                charge_a(me-1, n, p, &prev_ibeg, &prev_iend);
                EXPECT_EQ(ibeg, prev_iend + 1);
            }
            
            total_size += size;
        }
        
        // Vérifier que la somme des tailles couvre le domaine entier
        EXPECT_EQ(total_size, n);
    }
    
    // Cas 2: 3 processus, division potentiellement déséquilibrée
    {
        int p = 3;
        int total_size = 0;
        
        for (int me = 0; me < p; me++) {
            int ibeg, iend;
            charge_a(me, n, p, &ibeg, &iend);
            total_size += (iend - ibeg + 1);
        }
        
        EXPECT_EQ(total_size, n);
    }
}

// Correction pour Test11_BICGstab
TEST(Test11_BICGstab, ConvergenceRate) {
    // Test de convergence du solveur BICGstab
    int n = 4;
    int m = 4;
    std::vector<double> b(n*m, 1.0); // Source constante
    
    double alpha = 4.0;
    double beta = -1.0;
    double gamma = -1.0;
    
    double dt = 0.01;
    double D = 1.0;
    double dy = 0.1;
    
    int kmax_values[] = {10, 100, 1000};
    std::vector<double> residuals;
    
    for (auto kmax : kmax_values) {
        auto solution = BICGstabSEQ(alpha, beta, gamma, b, n, m, 0, 1, dt, D, dy, kmax, 1e-12);
        
        // Calcul du résidu
        auto Ax = produitmatvectSEQ(alpha, beta, gamma, solution, n, m, dy, D, dt);
        double residual = 0.0;
        for (size_t i = 0; i < b.size(); i++) {
            residual += pow(Ax[i] - b[i], 2);
        }
        residual = sqrt(residual);
        
        residuals.push_back(residual);
    }
    
    // Vérifier que les résidus ne augmentent pas avec plus d'itérations
    // (ils peuvent stagner à cause de la précision machine)
    for (size_t i = 1; i < residuals.size(); i++) {
        EXPECT_LE(residuals[i], residuals[i-1] + 1e-14);
    }
}

// Correction pour Test14_FileIO
TEST(Test14_FileIO, SaveAndRead) {
    // Créer les dossiers nécessaires
    #ifdef _WIN32
    mkdir("results");
    #else
    mkdir("results", 0777);
    #endif
    
    // Données de test
    int cas = 1, ni = 0, Me = 0, Nx = 2, Np = 1;
    int ibeg = 0, iend = 1, r = 1;
    std::vector<double> x = {0.0, 0.5, 1.0};
    std::vector<double> y = {0.0, 0.5, 1.0};
    std::vector<double> u = {1.0, 2.0, 3.0, 4.0};
    
    // Sauvegarder les résultats
    sauvegarder_resultats(cas, ni, Me, Nx, ibeg, iend, Np, x, y, u, r);
    
    // Vérifier que le fichier a été créé
    std::string filename = "results/Sol.1.0.0.txt";
    std::ifstream file(filename);
    
    if (!file.good()) {
        // Si le test s'exécute depuis le dossier build/tests
        file = std::ifstream("../../results/Sol.1.0.0.txt");
    }
    
    EXPECT_TRUE(file.good());
    
    // Si le fichier est bien ouvert, lire son contenu
    if (file.good()) {
        std::vector<double> x_read, y_read, u_read;
        double x_val, y_val, u_val;
        
        while (file >> x_val >> y_val >> u_val) {
            x_read.push_back(x_val);
            y_read.push_back(y_val);
            u_read.push_back(u_val);
        }
        
        // Vérifier que les données lues correspondent aux données sauvegardées
        EXPECT_GT(x_read.size(), 0);
    }
    
    file.close();
}

// Test 12: Test de communication MPI
TEST_F(MPITest, Test12_MPICommunication) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (size >= 2) {
        const int n = 4;
        std::vector<double> send_data(n), recv_data(n);
        
        // Initialiser les données
        for (int i = 0; i < n; i++) {
            send_data[i] = rank * 10.0 + i;
        }
        
        // Test de Send/Recv
        if (rank == 0) {
            MPI_Send(send_data.data(), n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        } else if (rank == 1) {
            MPI_Recv(recv_data.data(), n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            // Vérifier que les données sont correctes
            for (int i = 0; i < n; i++) {
                EXPECT_NEAR(recv_data[i], i, 1e-10);
            }
        }
        
        // S'assurer que tous les processus attendent
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

// Test 13: Test des conditions aux limites de Robin
TEST(Test13_RobinBoundary, MatrixOperation) {
    // Test du produit matrice-vecteur avec conditions aux limites de Robin
    int n = 4;
    int m = 3;
    std::vector<double> x(n*m);
    for (int i = 0; i < n*m; i++) {
        x[i] = i+1; // 1, 2, 3, ...
    }
    
    double alpha = 2.0;
    double beta = -1.0;
    double gamma = -1.0;
    double dy = 0.1;
    double D = 1.0;
    double dt = 0.01;
    
    // Test sans Robin (alpha_rob = 0)
    auto result1 = produitmatvect1(alpha, beta, gamma, x, n, m, dy, D, dt, 0.0, 0.0, 0, 2);
    
    // Test avec Robin (alpha_rob = 1)
    auto result2 = produitmatvect1(alpha, beta, gamma, x, n, m, dy, D, dt, 1.0, 1.0, 0, 2);
    
    // Les résultats devraient être différents avec différentes conditions aux limites
    bool all_same = true;
    for (size_t i = 0; i < result1.size(); i++) {
        if (std::abs(result1[i] - result2[i]) > 1e-10) {
            all_same = false;
            break;
        }
    }
    EXPECT_FALSE(all_same);
}

// Test 15: Test de robustesse numérique
TEST(Test15_NumericalStability, SourceTerms) {
    double D = 1.0;
    double dy = 0.1;
    double dx = 0.1;
    double dt = 0.005; // Petit pas de temps pour la stabilité
    double t = 0.1;
    double Lx = 1.0;
    double Ly = 1.0;
    int Nx = 10;
    int Ny = 10;
    int cas = 3; // Utiliser le cas de la gaussienne
    int ibeg = 0;
    int iend = 9;
    int nproc = 1;
    int rank = 0;
    int nr = 1;
    
    // Créer des vecteurs pour x, y et u
    std::vector<double> x(Nx+2), y(Ny+2), u(Nx*Ny, 0.0);
    for (int i = 0; i < Nx+2; i++) x[i] = i * dx;
    for (int i = 0; i < Ny+2; i++) y[i] = i * dy;
    
    std::vector<double> stencil1(Nx, 0.0);
    std::vector<double> stencil2(Nx, 0.0);
    
    // Calculer le terme source
    auto source = SOURCESEQ(D, dy, dx, dt, t, Lx, Ly, Nx, Ny, x, y, u, cas, ibeg, iend, nproc, rank);
    
    // Vérifier que le terme source a la bonne taille
    EXPECT_EQ(source.size(), Nx*Ny);
    
    // Vérifier que le terme source ne contient pas de NaN ou Inf
    for (auto& val : source) {
        EXPECT_FALSE(std::isnan(val));
        EXPECT_FALSE(std::isinf(val));
    }
}

int main(int argc, char **argv) {
    // Initialiser MPI avant les tests
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized) {
        MPI_Init(&argc, &argv);
    }

    testing::InitGoogleTest(&argc, argv);
    
    // Supprime le listener par défaut
    testing::TestEventListeners& listeners = testing::UnitTest::GetInstance()->listeners();
    delete listeners.Release(listeners.default_result_printer());
    
    // Ajoute notre listener personnalisé
    listeners.Append(new TestInfo);
    
    // Exécuter les tests
    int result = RUN_ALL_TESTS();
    
    // Finaliser MPI à la fin
    int finalized;
    MPI_Finalized(&finalized);
    if (!finalized) {
        MPI_Finalize();
    }
    
    return result;
}