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
#ifndef META_H
#define META_H
/* A Selection of all used helper functions*/
//A template to select a type on a condition
template <bool flag, typename T, typename U>
struct Select//Select template
{
    typedef T Result;
};

template <typename T, typename U>
struct Select<false, T, U>
{
    typedef U Result;
};


const int DEFAULT = ~(~0u >> 1); //initialize with the smallest int

struct NilCase
{};

template <int tag_, class Type_, class Next_ = NilCase>
struct CASE
{ enum { tag = tag_ };
  typedef Type_ Type;
  typedef Next_ Next;
};


// SWITCH<>
template<int tag, class Case>
class SWITCH
{ typedef typename Case::Next NextCase;
  enum {
    caseTag = Case::tag,
    found   = (caseTag == tag || caseTag == DEFAULT)
  };
public:
  typedef typename Select<found,
             typename Case::Type,
             typename SWITCH<tag, NextCase>::RET
            >::Result RET;
};

template<int tag>
class SWITCH<tag, NilCase>
{public:
  typedef NilCase RET;
};

template <typename T>
class TypeTraits
{// A template to strip the constness from a type from "Modern C++ Design" by Andrei Alexandrescu
private:
    template <class U>
    struct UnConst
    {
        typedef U Result;
    };
    template <class U>
    struct UnConst<const U>
    {
        typedef U Result;
    };
public:
    typedef typename UnConst<T>::Result NonConstType;
};

template <int v>
struct Int2Type
{// A template to Conserve an Int in a type from "Modern C++ Design" by Andrei Alexandrescu
    enum { value = v };
};

template <typename T>
struct Type2Type
{
   typedef T OriginalType;
};

template<bool>
struct CompileTimeChecker
{
    CompileTimeChecker(...);
};
/*
template<>
struct CompileTimeChecker<false>
    { }
;//From Modern C++ Design
#define STATIC_CHECK(expr, msg) \
   {\
       class ERROR_##msg {}; \
       (void)sizeof(CompileTimeChecker<(expr) != 0>((ERROR_##msg())));\
   }
*/
template <class T, class U>
class Conversion
{//from "Modern C++ Design" by Andrei Alexandrescu
    typedef char Small;
    class Big
    {
        char dummy[2];
    };
    static Small Test(U);
    static Big Test(...);
    static T MakeT();
public:
    enum { exists =
               sizeof(Test(MakeT())) == sizeof(Small),
           sametype = false
         };
};

template <class T>
class Conversion<T,T>
{
public:
    enum { exists = 1, sametype = 1 };
};

#endif
