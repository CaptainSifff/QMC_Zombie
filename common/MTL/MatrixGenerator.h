/***************************************************************************
 *   Copyright (C) 2006 by Florian Goth   *
 *   CaptainSifff@gmx.de   *
 *                                                                         *
 *   Permission is hereby granted, free of charge, to any person obtaining *
 *   a copy of this software and associated documentation files (the       *
 *   "Software"), to deal in the Software without restriction, including   *
 *   without limitation the rights to use, copy, modify, merge, publish,   *
 *   distribute, sublicense, and/or sell copies of the Software, and to    *
 *   permit persons to whom the Software is furnished to do so, subject to *
 *   the following conditions:                                             *
 *                                                                         *
 *   The above copyright notice and this permission notice shall be        *
 *   included in all copies or substantial portions of the Software.       *
 *                                                                         *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       *
 *   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    *
 *   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. *
 *   IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR     *
 *   OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, *
 *   ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR *
 *   OTHER DEALINGS IN THE SOFTWARE.                                       *
 ***************************************************************************/
#ifndef MATRIX_GENERATOR_H
#define MATRIX_GENERATOR_H
#include "DSL.h"

/*enum {
do_all,
defaults_and_assemble,
assemble_components };*/

using namespace MatrixDSL;

/*template < class InputDSL = matrix<> , int WhatToDo = do_all >
class MATRIX_GENERATOR
{
// parse InputDSL
typedef SWITCH < WhatToDo
                ,CASE<assemble_components     ,           matrix<>
                ,CASE<defaults_and_assemble,                     matrix<>
           , DEFAULT<                                                  InputDSL>
 > > >::RET DSL_Description;

typedef typename MATRIX_DSL_PARSER < DSL_Description >::RET ParsedDSL__;
};*/

#endif
