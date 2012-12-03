#ifndef UGENEERROR_H
#define UGENEERROR_H

#include <iostream>
#include <string>
#include <stdexcept>
#include "boost-include/boost/exception/all.hpp"
#include <vector>

namespace NGS {
class uRegion;
class uTags;
class uGenericNGS;
class uToken;
class uBasicNGS;
enum class token_param;


typedef boost::error_info<struct string_info,std::string> string_error;
typedef boost::error_info<struct region_info,uRegion> region_error;
typedef boost::error_info<struct basic_info,uBasicNGS> basic_error;
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



/**< To pass a token*/
typedef boost::error_info<struct token_info,uToken> token_error;
/**< To pass one or more param type */
typedef boost::error_info<struct token_param_info,token_param> token_param_error;


/**< uToken Exceptions */
struct uToken_exception_base : virtual ugene_exception_base{};
struct invalid_uToken_throw : virtual uToken_exception_base{};
struct invalid_token_param_throw : virtual uToken_exception_base{};
struct invalid_value_throw : virtual uToken_exception_base{};
struct param_not_found : virtual invalid_uToken_throw{};


 /**< uParser Exceptions */
struct uParser_exception_base : virtual ugene_exception_base{};
struct uParser_invalid_type_instance : virtual uParser_exception_base{};
struct uParser_invalid_header : virtual uParser_exception_base{};
struct end_of_file_throw : virtual uParser_exception_base{};
struct uParser_missing_mandatory_values : virtual uParser_exception_base{};
struct uParser_missing_mandatory_header : virtual uParser_invalid_header{};
struct uParser_invalid_line : virtual uParser_exception_base{};
struct customParser_missing_mandatory_values : virtual uParser_missing_mandatory_values{};
struct uParserBed_invalid_number_of_columns : virtual uParser_exception_base{};
/**< Sam Parser exception */
struct uParser_invalid_Sam_header : virtual uParser_invalid_header{};
struct uParser_invalid_Sam_line : virtual uParser_invalid_line{};
/**< GFF Parser exception */
struct uParser_invalid_GFF_line : virtual uParser_invalid_line{};
/**< GTF Parser exception */
struct uParser_invalid_GTF_line : virtual uParser_invalid_line{};
/**<  BedGraph exceptions*/
struct uParser_invalid_BedGraph_line : virtual uParser_invalid_line{};

/**< uHeader Exceptions */

struct invalid_header_param_throw : virtual uParser_exception_base{};


/**< uWrite exception */
struct uWriter_exception_base : virtual ugene_exception_base{};
struct uWriter_invalid_type_instance : virtual uWriter_exception_base{};
struct no_fields_names : virtual uWriter_exception_base{};
struct uWriter_missing_mandatory_param : virtual uWriter_exception_base{};
struct uWriter_missing_mandatory_header : virtual uWriter_exception_base{};

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

