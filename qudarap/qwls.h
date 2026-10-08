//
// Created by Ilker Parmaksiz on 8/31/26.
// Implementation of G4OpWLS for Opticks
//


#pragma once
//#include "QWLS.hh"
/**
wls.h
==================
**/
struct QWLS;
#if defined(__CUDACC__) || defined(__CUDABE__)
   #define WLS_METHOD __device__
#else
   #define WLS_METHOD
#endif


#include "qrng.h"
struct quad4 ;
struct quad6 ;
struct sphoton ;

#include "OpticksPhoton.h"
struct qwls
{
    cudaTextureObject_t qwls_tex ;
    quad4*              qwls_meta ; // HUH: not used ?  q0.f.x (WLS Time),q0.f.y (Mean Photons), q0.u.z (delta or exponential time)
    unsigned            hd_factor ;
    float              time_constant ;
    float              mean_number_photons ;


#if defined(__CUDACC__) || defined(__CUDABE__) || defined(MOCK_CURAND) || defined(MOCK_CUDA)
    WLS_METHOD void    wlsemit(   sphoton& p, RNG& rng) const ;
    // sets direction, polarization and wavelength as needed by both generate and reemit
    WLS_METHOD float   wavelength(     const float& u0) const ;
    WLS_METHOD float   wavelength_hd0( const float& u0) const ;
    WLS_METHOD float   wavelength_hd10(const float& u0) const ;
    WLS_METHOD float   wavelength_hd20(const float& u0) const ;


#endif

};


#if defined(__CUDACC__) || defined(__CUDABE__) || defined(MOCK_CURAND) || defined(MOCK_CUDA)


/**
qwls::wlsemit : dir,pol and wavelength do not depend on genstep param
--------------------------------------------------------------------------------

Translation of "jcv DsG4Scintillation"

**/

inline WLS_METHOD void qwls::wlsemit(sphoton& p, RNG& rng) const
{
    float u0 = curand_uniform(&rng);
    float u1 = curand_uniform(&rng);
    float u2 = curand_uniform(&rng);
    float u3 = curand_uniform(&rng);
    float u4 = curand_uniform(&rng);
    float cost = 1.f - 2.f*u0;
    float sint = sqrt((1.f-cost)*(1.f+cost));
    float phi = 2.f*M_PIf*u1;
    float sinp = sin(phi);
    float cosp = cos(phi);

    p.mom.x = sint*cosp;
    p.mom.y = sint*sinp;
    p.mom.z = cost ;
    p.orient_iindex = 0u ;

    // Determine polarization of new photon
    p.pol.x = cost*cosp ;
    p.pol.y = cost*sinp ;
    p.pol.z = -sint ;

    phi = 2.f*M_PIf*u2 ;
    sinp = sin(phi);
    cosp = cos(phi);

    p.pol = normalize( cosp*p.pol + sinp*cross(p.mom, p.pol) ) ;
    p.wavelength = wavelength(u3);

    // Add exponential decay time here if needed, otherwise delta time is used
    if(qwls_meta->q0.u.z) p.time += -logf(u4) * time_constant;// HUH: q0.u.z is delta or exponential time, not used in G4OpWLS

}


inline WLS_METHOD float qwls::wavelength(const float& u0) const
{
    float wl ;
    switch(hd_factor)
    {
        case 0:  wl = wavelength_hd0(u0)  ; break ;
        case 10: wl = wavelength_hd10(u0) ; break ;
        case 20: wl = wavelength_hd20(u0) ; break ;
        default: wl = 0.f ;
    }
    //printf("//wls::wavelength wl %10.4f hd %d \n", wl, hd_factor );
    return wl ;
}


inline WLS_METHOD float qwls::wavelength_hd0(const float& u0) const
{
    constexpr float y0 = 0.5f/3.f ;
    return tex2D<float>(qwls_tex, u0, y0 );
}

/**
wls::wavelength_hd10
--------------------------------------------------

Idea is to improve handling of extremes by throwing ten times the bins
at those regions, using simple and cheap linear mappings.

TODO: move hd "layers" into float4 payload so the 2d cerenkov and 1d scint
icdf texture can share some of teh implementation

**/

inline WLS_METHOD float qwls::wavelength_hd10(const float& u0) const
{
    float wl ;

    constexpr float y0 = 0.5f/3.f ;
    constexpr float y1 = 1.5f/3.f ;
    constexpr float y2 = 2.5f/3.f ;

    if( u0 < 0.1f )
    {
        wl = tex2D<float>(qwls_tex, u0*10.f , y1 );
    }
    else if ( u0 > 0.9f )
    {
        wl = tex2D<float>(qwls_tex, (u0 - 0.9f)*10.f , y2 );
    }
    else
    {
        wl = tex2D<float>(qwls_tex, u0,  y0 );
    }
    return wl ;
}



inline WLS_METHOD float qwls::wavelength_hd20(const float& u0) const
{
    float wl ;

    constexpr float y0 = 0.5f/3.f ;
    constexpr float y1 = 1.5f/3.f ;
    constexpr float y2 = 2.5f/3.f ;

    if( u0 < 0.05f )
    {
        wl = tex2D<float>(qwls_tex, u0*20.f , y1 );
    }
    else if ( u0 > 0.95f )
    {
        wl = tex2D<float>(qwls_tex, (u0 - 0.95f)*20.f , y2 );
    }
    else
    {
        wl = tex2D<float>(qwls_tex, u0,  y0 );
    }
    return wl ;
}


#endif



