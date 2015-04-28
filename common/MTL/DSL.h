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
//DSL-Specification:
#ifndef DSL_H
#define DSL_H
namespace MatrixDSL
{
//Declarations
//Matrix : matrix[Elementtype , BoundsChecking, Size]
/*template < class Elementtype , class BoundsChecking , class Size >
struct matrix;*/
//Elementtype : float | double | ...
//They're built-in types -> nothing to declare
//BoundsChecking : CheckBounds | NoBoundsCheck
template < class dummy > class CheckBounds;
template < class dummy > class NoBoundsCheck;
//Size : whatever | { n ,m }
template < class dummy >
class Size;

struct whatever
    {
enum {
whatever_id = -1,
check_bounds_id,
no_check_bounds_id };
}
;

//Implementation:
template < typename Elementtype = whatever , 
                 class BoundsChecking = whatever, 
                 class Size = whatever >
struct matrix
{
    typedef Elementtype elementtype;
    typedef BoundsChecking boundschecking;
    typedef Size size;
};

template < class dummy = whatever > struct CheckBounds : public whatever
{//It inherits the Declaration of check_bounds_id from whatever
enum { id = check_bounds_id };
};

template < class dummy = whatever > class NoBoundsCheck : public whatever
{
enum { id = no_check_bounds_id };
};

template < class dummy = whatever >
class Size
{
    enum { rows = 0 ,columns = 0};
};

};
#endif
