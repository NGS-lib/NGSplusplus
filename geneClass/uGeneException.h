#ifndef UGENEERROR_H
#define UGENEERROR_H

#include <iostream>
#include <string>
#include <stdexcept>
#include "boost/exception/all.hpp"
#include <vector>

namespace NGS {
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
struct skipped_elem_throw : virtual elem_throw{};
struct construct_elem_throw : virtual elem_throw{};
struct param_throw : virtual elem_throw{};
struct format_parsing_error : virtual ugene_exception_base{};
struct ugene_operation_throw : virtual ugene_exception_base{};
struct uChrom_operation_throw : virtual ugene_operation_throw{};
struct uExp_operation_throw : virtual uChrom_operation_throw{};

/**< File IO exceptions */



/**< uToken Exceptions */
struct uToken_exception_base : virtual ugene_exception_base{};
struct invalid_uToken_throw : virtual uToken_exception_base{};
struct invalid_token_param_throw : virtual uToken_exception_base{};
struct invalid_value_throw : virtual uToken_exception_base{};
struct param_not_found : virtual invalid_uToken_throw{};



 /**< uParser Exceptions */
struct uParser_exception_base : virtual ugene_exception_base{};
struct uParser_invalid_header : virtual uParser_exception_base{};
struct end_of_file_throw : virtual uParser_exception_base{};
struct uParser_missing_mandatory_values : virtual uParser_exception_base{};
struct uParser_missing_mandatory_header : virtual uParser_invalid_header{};
struct uParser_invalid_line : virtual uParser_exception_base{};
struct customParser_missing_mandatory_values : virtual uParser_missing_mandatory_values{};
/**< Sam Parser exception */
struct uParser_invalid_Sam_header : virtual uParser_invalid_header{};
struct uParser_invalid_Sam_line : virtual uParser_invalid_line{};

/**< uHeader Exceptions */

struct invalid_header_param_throw : virtual uParser_exception_base{};

/**< uWrite exception */
struct uWriter_exception_base : virtual ugene_exception_base{};
struct no_fields_names : virtual uWriter_exception_base{};

// Util functions
static inline void addStringError(ugene_exception_base & e, const std::string & err){
	std::string trace;
	if (std::string const * ste =boost::get_error_info<string_error>(e) )
        trace=*ste;
	e << string_error(trace+err+"\n");
}

static inline std::string fetchStringError(ugene_exception_base& e) {
	std::string trace;
	if (std::string const * ste =boost::get_error_info<string_error>(e))
		trace = *ste;
	return trace;
}
} // End of namespace NGS
#endif // UGENEERROR_H

