#ifndef TLIST_H
#define TLIST_H
struct NullType {};

template <class T, class U>
struct TList
{
    typedef T Head;
    typedef U Tail;
};

template <class T>
struct TListLength {};

template <class T, class U>
struct TListLength<TList<T,U> >
{
    enum
    {
        Ret = 1 + TListLength<U>::Ret
    };
};

template <>
struct TListLength<NullType>
{
    enum
    {
        Ret = 0
    };
};
#endif
