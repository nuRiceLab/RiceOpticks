//
// Created by Ilker Parmaksiz on 8/31/26.
// Implementation of G4OpWLS for RiceOpticks
//

#pragma once

#include <string>
#include "QUDARAP_API_EXPORT.hh"
#include "plog/Severity.h"

struct dim3 ;
struct NP ;
template <typename T> struct QTex ;
struct qoptwls ;

struct QUDARAP_API QOptWLS
{
    static const plog::Severity LEVEL ;
    static QTex<float>* MakeWLSQTex(const NP* src, unsigned hd_factor);

    const NP*      dsrc ;
    const NP*      src ;
    QTex<float>*    tex ;
    qoptwls*       optwls ;
    qoptwls*       d_optwls ;

    QOptWLS(const NP* icdf, unsigned hd_factor);

    void init();
    std::string desc() const ;

    void configureLaunch( dim3& numBlocks, dim3& threadsPerBlock, unsigned width, unsigned height );

    void check();
    NP*  lookup();

    void lookup( float* lookup, unsigned num_lookup, unsigned width, unsigned height );
    void dump(   float* lookup, unsigned num_lookup, unsigned edgeitems=10 );

};

thread_local QOptWLS qwls;
