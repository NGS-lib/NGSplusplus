#include "utility.h"
#include "gtest.h"

#include "uGeneException.h"

void FirstError();
void SecondLevel();
void ThirdLevel();

void FirstError(){

    std::cerr << "In first" <<std::endl;

        throw ugene_exception_base()<<test_string_error("happy") ;



}

void SecondLevel(){

    try {
        std::cerr << "In Second" <<std::endl;
    FirstError();


    }
    catch(boost::exception & e){
         std::cerr << "Throwing in second" <<std::endl;
        e << test_string_error(*boost::get_error_info<test_string_error>(e)+" happy2");
        throw;
    }

}


void ThirdLevel(){

    try {
         std::cerr << "In trird" <<std::endl;
        SecondLevel();

    }
        catch(ugene_exception_base & e){
            std::cerr << "Catching  in third" <<std::endl;

        if (std::string const * fn = boost::get_error_info<test_string_error>(e))
            std::cerr << "I am happy??? " << *fn<< std::endl;

        utility::pause_input();
    }


}


TEST(ERROR, ErrorTest)
{

    ThirdLevel();

}
