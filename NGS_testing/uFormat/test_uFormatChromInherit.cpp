/**< Test Common inherited functions */
#include "NGS++.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <time.h>
#include <thread>
#include <random>
#include <string.h>
#include <time.h>
#include <gtest/gtest.h>



using namespace std;
using namespace NGS;



/*
 * Setters/Getters testing - start/end (positions)
 *	Valid cases:
 *		SETGETSTART
 *		SETGETEND
 *	Invalid cases:
 *		SETSTARTILLEGAL
 *		SETSENDILLEGAL
 */
TEST(uBasicNGSCHROM_GetChrom, NONAMECHROM){
    //CHROMOSOME sans nom ""
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHROM_GetChrom, NOCHROMTHROWEXC){
    //Demande de chrom existe pas, throw
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHROM_GetChrom, VALIDCHROM){
    //VALID CHROM
       ASSERT_TRUE(false);
 }
/**<  */
TEST(uBasicNGSCHROM_getChromP, NONAMECHROM){
    //Chromosome pointer no name and works
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHROM_getChromP, NOCHROMTHROWEXC){
    //Chromosome pointer not exist and throw
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHROM_getChromP, VALIDCHROM){
    // Chromosome pointer, valid and works
       ASSERT_TRUE(false);
 }
/**<  */
 TEST(uBasicNGSCHROM_getConstChromP, NONAMECHROM){
    //Chromosome pointer no name and works
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHROM_getConstChromP, NOCHROMTHROWEXC){
    //Chromosome pointer not exist and throw
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSCHROM_getConstChromP, VALIDCHROM){
    // Chromosome pointer, valid and works
       ASSERT_TRUE(false);
 }
/**<  */

TEST(uBasicNGSCHROM_getConstChromP, VALIDCHROM){
    // Chromosome pointer, valid and works
       ASSERT_TRUE(false);
 }


