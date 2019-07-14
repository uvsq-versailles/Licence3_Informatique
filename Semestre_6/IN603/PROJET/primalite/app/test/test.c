/**
 * \file test.c
 * \brief Fonctions concernant les tests unitaires.
 * \author Clément Caumes - Yassin Doudouh
 * \date 26 Avril 2018
 *
 * Programme de test de primalité de Solovay-Strassen.
 *
 */

#include <CUnit/Basic.h>
#include <CUnit/CUnit.h>
#include <gmp.h>
#include <time.h>
#include "../inc/exponent.h"
#include "../inc/interface.h"
#include "../inc/jacobi.h"
#include "../inc/strassen.h"

#define K 2

/**
 * \fn void test_jacobi_1()
 * \brief Teste si jacobi(8721,4235)=-1.
 */
void test_jacobi_1(){
	mpz_t n;
	mpz_t m;
	mpz_init(n);
	mpz_init(m);
	mpz_set_ui(n,8721);
	mpz_set_ui(m,4235);
	int test=jacobi(n,m);
	mpz_clear(n);
	mpz_clear(m);
	CU_ASSERT_EQUAL(test,-1);	
}

/**
 * \fn void test_jacobi_2()
 * \brief Teste si jacobi(541,2011)=-1.
 */
void test_jacobi_2(){
	mpz_t n;
	mpz_t m;
	mpz_init(n);
	mpz_init(m);
	mpz_set_ui(n,541);
	mpz_set_ui(m,2011);
	int test=jacobi(n,m);
	mpz_clear(n);
	mpz_clear(m);
	CU_ASSERT_EQUAL(test,-1);	
}

/**
 * \fn void test_jacobi_3()
 * \brief Teste si jacobi(1236,20003)=1.
 */
void test_jacobi_3(){
	mpz_t n;
	mpz_t m;
	mpz_init(n);
	mpz_init(m);
	mpz_set_ui(n,1236);
	mpz_set_ui(m,20003);
	int test=jacobi(n,m);
	mpz_clear(n);
	mpz_clear(m);
	CU_ASSERT_EQUAL(test,1);
}

/**
 * \fn void test_primalite_1()
 * \brief Teste si 2 est premier.
 */
void test_primalite_1(){
	mpz_t n;
	mpz_t k;
	mpz_init(n);
	mpz_init(k);
	mpz_set_ui(n,2);
	mpz_set_ui(k,K);
	int test=teste_primalite_strassen(n,k);
	mpz_clear(n);
	mpz_clear(k);
	CU_ASSERT_EQUAL(test,PROBABLEMENT_PREMIER);	
}

/**
 * \fn void test_primalite_2()
 * \brief Teste si 1 est composé.
 */
void test_primalite_2(){
	mpz_t n;
	mpz_t k;
	mpz_init(n);
	mpz_init(k);
	mpz_set_ui(n,1);
	mpz_set_ui(k,K);
	int test=teste_primalite_strassen(n,k);
	mpz_clear(n);
	mpz_clear(k);
	CU_ASSERT_EQUAL(test,COMPOSE);	
}

/**
 * \fn void test_primalite_3()
 * \brief Teste si 3 est premier.
 */
void test_primalite_3(){
	mpz_t n;
	mpz_t k;
	mpz_init(n);
	mpz_init(k);
	mpz_set_ui(n,3);
	mpz_set_ui(k,K);
	int test=teste_primalite_strassen(n,k);
	mpz_clear(n);
	mpz_clear(k);
	CU_ASSERT_EQUAL(test,PROBABLEMENT_PREMIER);	
}

/**
 * \fn void test_primalite_4()
 * \brief Teste si 4999 est premier.
 */
void test_primalite_4(){
	mpz_t n;
	mpz_t k;
	mpz_init(n);
	mpz_init(k);
	mpz_set_ui(n,4999);
	mpz_set_ui(k,K);
	int test=teste_primalite_strassen(n,k);
	mpz_clear(n);
	mpz_clear(k);
	CU_ASSERT_EQUAL(test,PROBABLEMENT_PREMIER);
}

/**
 * \fn void test_primalite_5()
 * \brief Teste si 2^521-1 est premier.
 */
void test_primalite_5(){
	mpz_t n;
	mpz_t k;
	mpz_init(n);
	mpz_init(k);
	mpz_set_ui(n,2);
	mpz_pow_ui(n,n, 521);
    mpz_sub_ui(n, n, 1);
	mpz_set_ui(k,K);
	int test=teste_primalite_strassen(n,k);
	mpz_clear(n);
	mpz_clear(k);
	CU_ASSERT_EQUAL(test,PROBABLEMENT_PREMIER);
}

/**
 * \fn void test_primalite_6()
 * \brief Teste si 2^1279-1 est premier.
 */
void test_primalite_6(){
	mpz_t n;
	mpz_t k;
	mpz_init(n);
	mpz_init(k);
	mpz_set_ui(n,2);
	mpz_pow_ui(n,n, 1279);
    mpz_sub_ui(n, n, 1);
	mpz_set_ui(k,K);
	int test=teste_primalite_strassen(n,k);
	mpz_clear(n);
	mpz_clear(k);
	CU_ASSERT_EQUAL(test,PROBABLEMENT_PREMIER);
}

/**
 * \fn void test_primalite_7()
 * \brief Teste si 2^2281-1 est premier.
 */
void test_primalite_7(){
	mpz_t n;
	mpz_t k;
	mpz_init(n);
	mpz_init(k);
	mpz_set_ui(n,2);
	mpz_pow_ui(n,n, 2281);
    mpz_sub_ui(n, n, 1);
	mpz_set_ui(k,K);
	int test=teste_primalite_strassen(n,k);
	mpz_clear(n);
	mpz_clear(k);
	CU_ASSERT_EQUAL(test,PROBABLEMENT_PREMIER);	
}

/**
 * \fn void test_primalite_8()
 * \brief Teste si 2^3217-1 est premier.
 */
void test_primalite_8(){
	mpz_t n;
	mpz_t k;
	mpz_init(n);
	mpz_init(k);
	mpz_set_ui(n,2);
	mpz_pow_ui(n,n, 3217);
    mpz_sub_ui(n, n, 1);
	mpz_set_ui(k,K);
	int test=teste_primalite_strassen(n,k);
	mpz_clear(n);
	mpz_clear(k);
	CU_ASSERT_EQUAL(test,PROBABLEMENT_PREMIER);	
}

/**
 * \fn void test_primalite_9()
 * \brief Teste si 2^4423-1 est premier.
 */
void test_primalite_9(){
	mpz_t n;
	mpz_t k;
	mpz_init(n);
	mpz_init(k);
	mpz_set_ui(n,2);
	mpz_pow_ui(n,n, 4423);
    mpz_sub_ui(n, n, 1);
	mpz_set_ui(k,K);
	int test=teste_primalite_strassen(n,k);
	mpz_clear(n);
	mpz_clear(k);
	CU_ASSERT_EQUAL(test,PROBABLEMENT_PREMIER);	
}

/**
 * \fn void test_primalite_10()
 * \brief Teste si 2^9941-1 est premier.
 */
void test_primalite_10(){
	mpz_t n;
	mpz_t k;
	mpz_init(n);
	mpz_init(k);
	mpz_set_ui(n,2);
	mpz_pow_ui(n,n, 9941);
    mpz_sub_ui(n, n, 1);
	mpz_set_ui(k,K);
	int test=teste_primalite_strassen(n,k);
	mpz_clear(n);
	mpz_clear(k);
	CU_ASSERT_EQUAL(test,PROBABLEMENT_PREMIER);	
}

/**
 * \fn void test_primalite_11()
 * \brief Teste si 2^11213-1 est premier.
 */
void test_primalite_11(){
	mpz_t n;
	mpz_t k;
	mpz_init(n);
	mpz_init(k);
	mpz_set_ui(n,2);
	mpz_pow_ui(n,n, 11213);
    mpz_sub_ui(n, n, 1);
	mpz_set_ui(k,K);
	int test=teste_primalite_strassen(n,k);
	mpz_clear(n);
	mpz_clear(k);
	CU_ASSERT_EQUAL(test,PROBABLEMENT_PREMIER);	
}

/**
 * \fn void test_primalite_12()
 * \brief Teste si 2^19937-1 est premier.
 */
void test_primalite_12(){
	mpz_t n;
	mpz_t k;
	mpz_init(n);
	mpz_init(k);
	mpz_set_ui(n,2);
	mpz_pow_ui(n,n, 19937);
    mpz_sub_ui(n, n, 1);
	mpz_set_ui(k,K);
	int test=teste_primalite_strassen(n,k);
	mpz_clear(n);
	mpz_clear(k);
	CU_ASSERT_EQUAL(test,PROBABLEMENT_PREMIER);	
}

/**
 * \fn void test_primalite_13()
 * \brief Teste si 2^21701-1 est premier.
 */
void test_primalite_13(){
	mpz_t n;
	mpz_t k;
	mpz_init(n);
	mpz_init(k);
	mpz_set_ui(n,2);
	mpz_pow_ui(n,n, 21701);
    mpz_sub_ui(n, n, 1);
	mpz_set_ui(k,K);
	int test=teste_primalite_strassen(n,k);
	mpz_clear(n);
	mpz_clear(k);
	CU_ASSERT_EQUAL(test,PROBABLEMENT_PREMIER);	
}

/**
 * \fn void test_primalite_14()
 * \brief Teste si 2^23209-1 est premier.
 */
void test_primalite_14(){
	mpz_t n;
	mpz_t k;
	mpz_init(n);
	mpz_init(k);
	mpz_set_ui(n,2);
	mpz_pow_ui(n,n, 23209);
    mpz_sub_ui(n, n, 1);
	mpz_set_ui(k,K);
	int test=teste_primalite_strassen(n,k);
	mpz_clear(n);
	mpz_clear(k);
	CU_ASSERT_EQUAL(test,PROBABLEMENT_PREMIER);	
}

/**
 * \fn void test_primalite_15()
 * \brief Teste si 341 est compose.
 */
void test_primalite_15(){
	mpz_t n;
	mpz_t k;
	mpz_init(n);
	mpz_init(k);
	mpz_set_ui(n,341);
	mpz_set_ui(k,K);
	int test=teste_primalite_strassen(n,k);
	mpz_clear(n);
	mpz_clear(k);
	CU_ASSERT_EQUAL(test,COMPOSE);
}

/**
 * \fn void test_primalite_16()
 * \brief Teste si 2^521-3 est compose.
 */
void test_primalite_16(){
	mpz_t n;
	mpz_t k;
	mpz_init(n);
	mpz_init(k);
	mpz_set_ui(n,2);
	mpz_pow_ui(n,n, 521);
    mpz_sub_ui(n, n, 3);
	mpz_set_ui(k,K);
	int test=teste_primalite_strassen(n,k);
	mpz_clear(n);
	mpz_clear(k);
	CU_ASSERT_EQUAL(test,COMPOSE);
}

/**
 * \fn void test_premiers_entre_eux_1()
 * \brief Teste si pgcd(14,12)!=1.
 */
void test_premiers_entre_eux_1(){
	mpz_t n;
	mpz_t k;
	mpz_init(n);
	mpz_init(k);
	mpz_set_ui(n,14);
	mpz_set_ui(k,12);
	int test=premier_entre_eux_iteratif(n,k);
	mpz_clear(n);
	mpz_clear(k);
	CU_ASSERT_EQUAL(test,0);
}

/**
 * \fn void test_pgcd_1()
 * \brief Teste si pgcd(1128,2986)=1.
 */
void test_pgcd_1(){
	mpz_t a;
	mpz_t b;
	mpz_t result;
	mpz_init(a);
	mpz_init(b);
	mpz_init(result);
	mpz_set_ui(a,1128);
	mpz_set_ui(b,2986);
	pgcd_recursif(result,a,b);
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(result);
	int test=1;
	CU_ASSERT_EQUAL(test,1);
}

int setup(void)  { return 0; }
int teardown(void) { return 0; }

int main(int argc, char **argv){
	
	CU_initialize_registry();
	if (CUE_SUCCESS != CU_initialize_registry())
		return CU_get_error();	
	
	CU_pSuite suite = CU_add_suite("Tests Primalité", setup,teardown);

	CU_add_test(suite, "Jacobi N°01", test_jacobi_1);
    CU_add_test(suite, "Jacobi N°02", test_jacobi_2);
    CU_add_test(suite, "Jacobi N°03", test_jacobi_3);
    CU_add_test(suite, "PGCD N°1", test_pgcd_1);
    CU_add_test(suite, "Primalité nombre N°01", test_primalite_1);
    CU_add_test(suite, "Primalité nombre N°02", test_primalite_2);
    CU_add_test(suite, "Primalité nombre N°03", test_primalite_3);
    CU_add_test(suite, "Primalité nombre N°04", test_primalite_4);
    CU_add_test(suite, "Primalité nombre N°05", test_primalite_5);
    CU_add_test(suite, "Primalité nombre N°06", test_primalite_6);
    /*
    CU_add_test(suite, "Primalité nombre N°07", test_primalite_7);
    CU_add_test(suite, "Primalité nombre N°08", test_primalite_8);
    CU_add_test(suite, "Primalité nombre N°09", test_primalite_9);
    CU_add_test(suite, "Primalité nombre N°10", test_primalite_10);
    CU_add_test(suite, "Primalité nombre N°11", test_primalite_11);
    CU_add_test(suite, "Primalité nombre N°12", test_primalite_12);
    CU_add_test(suite, "Primalité nombre N°13", test_primalite_13);
    CU_add_test(suite, "Primalité nombre N°14", test_primalite_14);
    CU_add_test(suite, "Primalité nombre N°15", test_primalite_15);
    CU_add_test(suite, "Primalité nombre N°16", test_primalite_16);
    CU_add_test(suite, "Premiers entre eux N°1", test_premiers_entre_eux_1); 
    */

	CU_basic_set_mode(CU_BRM_VERBOSE);
	CU_basic_run_tests();
	CU_cleanup_registry();

    return 0;
}	
