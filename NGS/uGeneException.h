// ***************************************************************************
// uGeneException.h (c) 2013
// Alexei Nordell-Markovits : Sherbrooke University
// Charles Joly Beauparlant : Laval University
//
//       This file is part of the NGS++ library.
//
//    The NGS++ library is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU  Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public License
//    along with this program (lgpl-3.0.txt).  If not, see <http://www.gnu.org/licenses/>.
// ***************************************************************************


#ifndef UGENEERROR_H
#define UGENEERROR_H

#include <iostream>
#include <string>
#include <stdexcept>
#include <boost/exception/all.hpp>
#include <vector>

namespace NGS {
class uRegion;
class uTags;
class uToken;
class uBasicNGS;
class uGene;
enum class token_param;


typedef boost::error_info<struct string_info,std::string> string_error;
typedef boost::error_info<struct region_info,uRegion> region_error;
typedef boost::error_info<struct basic_info,uBasicNGS> basic_error;
typedef boost::error_info<struct tag_info,uTags> tag_error;
typedef boost::error_info<struct skip_tags,std::vector<uTags>> skipped_tags;
typedef boost::error_info<struct skip_genes,std::vector<uGene>> skipped_genes;
typedef boost::error_info<struct skip_regions,std::vector<uRegion>> skipped_regions;
typedef boost::error_info<struct skip_regions,std::vector<uBasicNGS>> skipped_Basic;

struct ugene_exception_base : virtual std::exception, virtual boost::exception
{};

struct unsorted_throw : virtual ugene_exception_base{};
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
/**< GenePred exceptions */
struct uParser_invalid_GenePred_line  : virtual uParser_invalid_line{};

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

