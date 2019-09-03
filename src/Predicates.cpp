#include "Predicates.h"

// statics

//template <> struct CtorConsts<float> {
//    /*  cword = 4223; */
//    static const int CWORD = 4210;
//};
//
//template <> struct CtorConsts<double> {
//    /*  cword = 4735; */
//    static const int CWORD = 4722;
//};

template<> const int CtorConsts<float>::CWORD = 4210;
template<> const int CtorConsts<double>::CWORD = 4722;

template<typename Precision> Precision PredicatesBase<Precision>::splitter;

template<typename Precision> Precision PredicatesBase<Precision>::epsilon = 99999;

template<typename Precision> Precision PredicatesBase<Precision>::resulterrbound;

template<typename Precision> Precision PredicatesBase<Precision>::ccwerrboundA;
template<typename Precision> Precision PredicatesBase<Precision>::ccwerrboundB;
template<typename Precision> Precision PredicatesBase<Precision>::ccwerrboundC;

template<typename Precision> Precision PredicatesBase<Precision>::o3derrboundA;
template<typename Precision> Precision PredicatesBase<Precision>::o3derrboundB;
template<typename Precision> Precision PredicatesBase<Precision>::o3derrboundC;

template<typename Precision> Precision PredicatesBase<Precision>::iccerrboundA;
template<typename Precision> Precision PredicatesBase<Precision>::iccerrboundB;
template<typename Precision> Precision PredicatesBase<Precision>::iccerrboundC;

template<typename Precision> Precision PredicatesBase<Precision>::isperrboundA;
template<typename Precision> Precision PredicatesBase<Precision>::isperrboundB;
template<typename Precision> Precision PredicatesBase<Precision>::isperrboundC;

template<typename Precision> typename PredicatesBase<Precision>::constructor PredicatesBase<Precision>::ctor;

// specilizations

template
class PredicatesBase<float>;

template
class PredicatesBase<double>;
