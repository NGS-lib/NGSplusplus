#ifndef UGENEERROR_H
#define UGENEERROR_H

#include <iostream>
#include <string>
#include <stdexcept>
#include <boost/exception/all.hpp>

class uRegion;
class uTags;
class uGenericNGS;

typedef boost::error_info<struct string_info,std::string> string_error;
typedef boost::error_info<struct region_info,uRegion> region_error;
typedef boost::error_info<struct tag_info,uTags> tag_error;
typedef boost::error_info<struct generic_info,uGenericNGS> generic_error;

typedef boost::error_info<struct skip_tags,std::vector<uTags>> skipped_tags;
typedef boost::error_info<struct skip_regions,std::vector<uRegion>> skipped_regions;

struct ugene_exception_base : virtual std::exception, virtual boost::exception
{};


struct elem_throw : virtual ugene_exception_base{};
struct skipped_elem_throw : virtual ugene_exception_base{};
struct construct_elem_throw : virtual elem_throw{};
struct format_parsing_error : virtual ugene_exception_base{};
#endif // UGENEERROR_H
