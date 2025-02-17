#include <gtest/gtest.h>
#include <mpi.h>
#include <fstream>
#include <iostream>
#include "Matrix.h"
#include "solver.h"
#include "functions.h"
#include "settings.h"
#include <cmath>
#include <vector>

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
        else if (test_case_name == "Test7_FileIO")
            return "Test des opérations de fichiers";
        else if (test_case_name == "Test8_BoundaryConditions")
            return "Test des conditions aux limites de Robin";
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
        MPI_Init(&argc, &argv);
    }

    void TearDown() override {
        MPI_Finalize();
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

// Test 7: Sauvegarde des résultats
TEST(Test7_FileIO, SaveResults) {
    std::vector<double> x = {0.0, 0.1, 0.2};
    std::vector<double> y = {0.0, 0.1, 0.2};
    std::vector<double> u = {1.0, 2.0, 3.0, 4.0};
    
    sauvegarder_resultats(1, 0, 0, 2, 1, 2, 1, x, y, u, 1);
    std::ifstream file("Sol.1.0.0.txt");
    EXPECT_TRUE(file.good());
    file.close();
}

// Test 8: Conditions aux limites de Robin
TEST(Test8_BoundaryConditions, RobinConditions) {
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0};
    double alpha = 2.0, beta = -1.0, gamma = -1.0;
    int n = 2, m = 2;
    double dy = 0.1, D = 1.0, dt = 0.01;
    double alpha_rob = 1.0, beta_rob = 1.0;
    
    auto result = produitmatvect1(alpha, beta, gamma, x, n, m, dy, D, dt, 
                                alpha_rob, beta_rob, 0, 2);
    EXPECT_NE(result[0], result[1]);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    
    // Supprime le listener par défaut
    testing::TestEventListeners& listeners = testing::UnitTest::GetInstance()->listeners();
    delete listeners.Release(listeners.default_result_printer());
    
    // Ajoute notre listener personnalisé
    listeners.Append(new TestInfo);
    
    return RUN_ALL_TESTS();
}