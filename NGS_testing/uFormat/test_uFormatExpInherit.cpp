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
TEST(uBasicNGSEXP_GetChrom, NONAMECHROM){
    //CHROMOSOME sans nom ""
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_GetChrom, NOCHROMTHROWEXC){
    //Demande de chrom existe pas, throw
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_GetChrom, VALIDCHROM){
    //VALID CHROM
       ASSERT_TRUE(false);
 }
/**<  */
TEST(uBasicNGSEXP_getChromP, NONAMECHROM){
    //Chromosome pointer no name and works
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_getChromP, NOCHROMTHROWEXC){
    //Chromosome pointer not exist and throw
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_getChromP, VALIDCHROM){
    // Chromosome pointer, valid and works
       ASSERT_TRUE(false);
 }
/**<  */
 TEST(uBasicNGSEXP_getConstChromP, NONAMECHROM){
    //Chromosome pointer no name and works
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_getConstChromP, NOCHROMTHROWEXC){
    //Chromosome pointer not exist and throw
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_getConstChromP, VALIDCHROM){
    // Chromosome pointer, valid and works
       ASSERT_TRUE(false);
 }
/**<  */

/**< GetSite */
TEST(uBasicNGSEXP_getSite, VALID){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_getSite, OUTOFBOUND){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSEXP_getSite, BELOW0){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSEXP_getSite, VALIDITERRATOR){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSEXP_getSite, INVALIDITERATOR){
       ASSERT_TRUE(false);
 }
 /**< GetOverlapping */
 TEST(uBasicNGSEXP_getOverlapping, VALIDEXP){
       ASSERT_TRUE(false);
 }
  TEST(uBasicNGSEXP_getOverlapping, VALIDCHROM){
       ASSERT_TRUE(false);
 }
  TEST(uBasicNGSEXP_getOverlapping, VALIDPOS){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_getOverlapping, EMPTYEXP){
       ASSERT_TRUE(false);
 }
  TEST(uBasicNGSEXP_getOverlapping, EMPTYCHROM){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_getOverlapping, EMPTYTHISEXP){
       ASSERT_TRUE(false);
 }
  TEST(uBasicNGSEXP_getOverlapping, EMPTYTHISCHROM){
       ASSERT_TRUE(false);
 }
  TEST(uBasicNGSEXP_getOverlapping, EMPTYTHISANDPOS){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_getOverlapping, POLYMORPHICHROM){
       ASSERT_TRUE(false);
 }
TEST(uBasicNGSEXP_getOverlapping, POLYMORPHICEXP){
       ASSERT_TRUE(false);
 }

/**< set/get Chr Size */
TEST(uBasicNGSEXP_getChrSize, VALID){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSEXP_getChrSize, NOTVALID){
       ASSERT_TRUE(false);
 }
  TEST(uBasicNGSEXP_getChrSize, EMPTY){
       ASSERT_TRUE(false);
 }

 TEST(uBasicNGSEXP_setOverlapping, INVALIDVALUE){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSEXP_setOverlapping, INVALUECHR){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSEXP_setOverlapping, VALID){
       ASSERT_TRUE(false);
 }

/**< getSubsetCount */
 TEST(uBasicNGSEXP_getSubsetCount, POSITIONS){
       ASSERT_TRUE(false);
 }
  TEST(uBasicNGSEXP_getSubsetCount, ELEMENT){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSEXP_getSubsetCount, CHROMNOEXIST){
       ASSERT_TRUE(false);
 }
 TEST(uBasicNGSEXP_getSubsetCount, NONAMECHROM){
       ASSERT_TRUE(false);
 }
