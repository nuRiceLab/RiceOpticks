//
// Created by Ilker Parmaksiz on 8/31/26.
// Implementation of G4OpWLS for Opticks
//


#pragma once
/**
qoptwls.h
==================
**/

#if defined(__CUDACC__) || defined(__CUDABE__)
   #define QOPTWLS_METHOD __device__
#else
   #define QOPTWLS_METHOD
#endif


#include "qrng.h"
struct quad4 ;
struct quad6 ;
struct sphoton ;

#include "OpticksPhoton.h"
struct qoptwls
{
    cudaTextureObject_t qoptwls_tex ;
    quad4*              qoptwls_meta ; // HUH: not used ?
    unsigned            hd_factor ;

#if defined(__CUDACC__) || defined(__CUDABE__) || defined(MOCK_CURAND) || defined(MOCK_CUDA)
    QOPTWLS_METHOD void    generate( sphoton& p, RNG& rng, const quad6& gs, unsigned long long photon_id, int genstep_id ) const ;
    QOPTWLS_METHOD void    reemit(   sphoton& p, RNG& rng, float scintillationTime) const ;
    QOPTWLS_METHOD void    momw_polw(sphoton& p, RNG& rng) const ;
    // sets direction, polarization and wavelength as needed by both generate and reemit

    QOPTWLS_METHOD float   wavelength(     const float& u0) const ;
    QOPTWLS_METHOD float   wavelength_hd0( const float& u0) const ;
    QOPTWLS_METHOD float   wavelength_hd10(const float& u0) const ;
    QOPTWLS_METHOD float   wavelength_hd20(const float& u0) const ;

#endif

};


#if defined(__CUDACC__) || defined(__CUDABE__) || defined(MOCK_CURAND) || defined(MOCK_CUDA)

//#include "sscint.h"
#include "qoptwls.h"
/**
qoptwls::generate_photon
------------------------

**/

inline QOPTWLS_METHOD void qoptwls::generate(
    sphoton& p,
    RNG& rng,
    const quad6& _gs,
    unsigned long long photon_id,
    int genstep_id ) const
{
    momw_polw(p, rng );

    const sscint& gs = (const sscint&)_gs ;

    float fraction = gs.charge == 0.f  ? 1.f : curand_uniform(&rng) ;
    p.pos = gs.pos + fraction*gs.DeltaPosition ;

    float u4 = curand_uniform(&rng) ;
    float deltaTime = fraction*gs.step_length/gs.meanVelocity - gs.ScintillationTime*logf(u4) ;

    p.time = gs.time + deltaTime ;
    p.zero_flags();
    p.set_flag(SCINTILLATION) ;
    p.set_PID(gs.ParentId); // For LArSoft
}




inline QOPTWLS_METHOD float qoptwls::wavelength(const float& u0) const
{
    float wl ;
    switch(hd_factor)
    {
        case 0:  wl = wavelength_hd0(u0)  ; break ;
        case 10: wl = wavelength_hd10(u0) ; break ;
        case 20: wl = wavelength_hd20(u0) ; break ;
        default: wl = 0.f ;
    }
    //printf("//qoptwls::wavelength wl %10.4f hd %d \n", wl, hd_factor );
    return wl ;
}


inline QOPTWLS_METHOD float qoptwls::wavelength_hd0(const float& u0) const
{
    constexpr float y0 = 0.5f/3.f ;
    return tex2D<float>(scint_tex, u0, y0 );
}

/**
qoptwls::wavelength_hd10
--------------------------------------------------

Idea is to improve handling of extremes by throwing ten times the bins
at those regions, using simple and cheap linear mappings.

TODO: move hd "layers" into float4 payload so the 2d cerenkov and 1d scint
icdf texture can share some of teh implementation

**/

inline QOPTWLS_METHOD float qoptwls::wavelength_hd10(const float& u0) const
{
    float wl ;

    constexpr float y0 = 0.5f/3.f ;
    constexpr float y1 = 1.5f/3.f ;
    constexpr float y2 = 2.5f/3.f ;

    if( u0 < 0.1f )
    {
        wl = tex2D<float>(scint_tex, u0*10.f , y1 );
    }
    else if ( u0 > 0.9f )
    {
        wl = tex2D<float>(scint_tex, (u0 - 0.9f)*10.f , y2 );
    }
    else
    {
        wl = tex2D<float>(scint_tex, u0,  y0 );
    }
    return wl ;
}



inline QOPTWLS_METHOD float qoptwls::wavelength_hd20(const float& u0) const
{
    float wl ;

    constexpr float y0 = 0.5f/3.f ;
    constexpr float y1 = 1.5f/3.f ;
    constexpr float y2 = 2.5f/3.f ;

    if( u0 < 0.05f )
    {
        wl = tex2D<float>(scint_tex, u0*20.f , y1 );
    }
    else if ( u0 > 0.95f )
    {
        wl = tex2D<float>(scint_tex, (u0 - 0.95f)*20.f , y2 );
    }
    else
    {
        wl = tex2D<float>(scint_tex, u0,  y0 );
    }
    return wl ;
}


#endif



