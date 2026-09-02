#pragma once
/**
U4OptWLS
================
WLS : Wavelength Shifting

History:
Ilker Parmaksiz, 08/31/2026

**/

#include <string>
#include <iomanip>

#include "Randomize.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "ssys.h"
#include "NPFold.h"
#include "U4MaterialPropertyVector.h"

struct U4OptWLS
{
    static constexpr const bool VERBOSE = false ;
    static constexpr const char* PROPS = "WLSMEANNUMBERPHOTONS,WLSABSLENGTH,WLSCOMPONENT,WLSTIMECONSTANT" ;
    static U4OptWLS* Create(const NPFold* materials );

    const NPFold* optwls ;
    const char* name ;

    const NP* wlsabslength ;
    const NP* wlscomponent ;
    const NP* wlstimeconstant ;

    const double epsilon ;

    // Variables to collect from Geant4
    const G4MaterialPropertyVector* WLSABSVector ;
    const G4MaterialPropertyVector* WLSComponentVector ;
    const G4MaterialPropertyVector* WLSIntegral ;
    const G4double MeanNumberPhotons ;
    const G4double TimeConstant;
    const NP* icdf ;

    const int num_wlsamp ;
    const NP* wlsamp ;


    U4OptWLS(const NPFold* fold, const char* name);
    void init();

    std::string desc() const ;
    NPFold* make_fold() const ;
    void save(const char* base, const char* rel=nullptr ) const ;

    NP* createWavelengthSamples( int num_samples );
    NP* createGeant4InterpolatedInverseCDF(
        int num_bins=4096,
        int hd_factor=20,
        const char* material_name="LS",
        bool energy_not_wavelength=false );

    static G4MaterialPropertyVector* Integral( const G4MaterialPropertyVector* theFastLightVector ) ;
    static NP* CreateWavelengthSamples(
        const G4MaterialPropertyVector* WLSIntegral,
        int num_samples
        );
    static NP* CreateGeant4InterpolatedInverseCDF(
        const G4MaterialPropertyVector* WLSIntegral,
        int num_bins,
        int hd_factor,
        const char* name,
        bool energy_not_wavelength
        );

};

inline U4OptWLS* U4OptWLS::Create(const NPFold* materials ) // static
{
    std::vector<const NPFold*> subs ;
    std::vector<std::string> names ;
    materials->find_subfold_with_all_keys( subs, names, PROPS );

    int num_subs = subs.size();
    int num_names = names.size() ;
    assert( num_subs == num_names );

    const char* name = num_names > 0 ? names[0].c_str() : nullptr ;
    const NPFold* sub = num_subs > 0 ? subs[0] : nullptr ;
    bool with_qoptwls = name && sub ;

    return with_qoptwls ? new U4OptWLS(sub, name) : nullptr ;
}


inline U4OptWLS::U4OptWLS(const NPFold* optwls_, const char* name_)
    :
    optwls(optwls_),
    name(strdup(name_)),
    wlsabslength(optwls->get("WLSABSLENGTH")),
    wlscomponent(optwls->get("WLSCOMPONENT")),
    wlstimeconstant(optwls->get("WLSTIMECONSTANT")),
    MeanNumberPhotons(optwls->get("WLSMEANNUMBERPHOTONS")),
    epsilon(0.),
    WLSABSVector(U4MaterialPropertyVector::FromArray(wlsabslength)),
    WLSComponentVector(U4MaterialPropertyVector::FromArray(wlscomponent)),
    WLSIntegral(Integral(WLSComponentVector)),
    icdf(nullptr),
    num_wlsamp(ssys::getenvint("U4OptWLS__num_wlsamp", 0)),
    wlsamp(nullptr)
{
    init();
}

inline void U4OptWLS::init()
{
    int num_bins = 4096 ;
    int hd_factor = 20 ;
    icdf = createGeant4InterpolatedInverseCDF(num_bins, hd_factor, name) ;
    wlsamp = createWavelengthSamples(num_wlsamp) ;
}



inline std::string U4OptWLS::desc() const
{
    std::stringstream ss ;
    ss << "U4OptWLS::desc" << std::endl
       << " name " << name
       << " WLSLength " << ( wlsabslength ? wlsabslength->sstr() : "-" )
       << " WLSComponent " << ( wlscomponent ? wlscomponent->sstr() : "-" )
       << " WLSTimeConstant " << ( wlstimeconstant ? wlstimeconstant->sstr() : "-" )
       << " icdf " << ( icdf ? icdf->sstr() : "-" )
       << " wlsamp " << ( wlsamp ? wlsamp->sstr() : "-" )
       << std::endl
       ;

    std::string str = ss.str();
    return str ;
}

inline NPFold* U4OptWLS::make_fold() const
{
    NPFold* fold = new NPFold ;
    fold->add("wlsabslength", wlsabslength) ;
    fold->add("wlscomponent", wlscomponent) ;
    fold->add("wlstimeconstant", wlstimeconstant) ;
    fold->add("icdf", icdf) ;
    if(wlsamp) fold->add("wlsamp", wlsamp) ;
    return fold ;
}

inline void U4OptWLS::save(const char* base, const char* rel ) const
{
    NPFold* fold = make_fold();
    fold->save(base, rel);
}


inline NP* U4OptWLS::createWavelengthSamples( int num_samples )
{
    return CreateWavelengthSamples(WLSIntegral, num_samples );
}

inline NP* U4OptWLS::createGeant4InterpolatedInverseCDF(
    int num_bins,
    int hd_factor,
    const char* material_name,
    bool energy_not_wavelength
    )
{
    return CreateGeant4InterpolatedInverseCDF(
               WLSIntegral,
               num_bins,
               hd_factor,
               material_name,
               energy_not_wavelength  );
}


/**
U4OptWLS::Integral
---------------------------

Returns cumulative sum of the input property on the same energy domain,
with values starting at 0. and increasing monotonically.

The is using trapezoidal numerical integration.

**/

inline G4MaterialPropertyVector* U4OptWLS::Integral( const G4MaterialPropertyVector* theFastLightVector )
{
     G4MaterialPropertyVector* aMaterialPropertyVector = new G4MaterialPropertyVector();

          if (theFastLightVector) {

               G4double currentIN = (*theFastLightVector)[0];

                if (currentIN >= 0.0) {

                    // Create first (photon energy, Scintillation
                    // Integral pair

                    G4double currentPM = theFastLightVector->
                        Energy(0);

                    G4double currentCII = 0.0;

                    aMaterialPropertyVector->
                        InsertValues(currentPM , currentCII);

                    // Set previous values to current ones prior to loop

                    G4double prevPM  = currentPM;
                    G4double prevCII = currentCII;
                    G4double prevIN  = currentIN;

                    // loop over all (photon energy, intensity)
                    // pairs stored for this material

                    for(size_t ii = 1;
                              ii < theFastLightVector->GetVectorLength();
                              ++ii)
                    {
                        currentPM = theFastLightVector->Energy(ii);

                        currentIN= (*theFastLightVector)[ii];

                        currentCII = 0.5 * (prevIN + currentIN);

                        currentCII = prevCII +
                            (currentPM - prevPM) * currentCII;

                        aMaterialPropertyVector->
                            InsertValues(currentPM, currentCII);

                        prevPM  = currentPM;
                        prevCII = currentCII;
                        prevIN  = currentIN;
                    }
               }
            }

    return aMaterialPropertyVector ;
}


inline NP* U4OptWLS::CreateWavelengthSamples(
    const G4MaterialPropertyVector* WLSIntegral_,
    int num_samples )
{
    if(num_samples == 0) return nullptr ;

    G4MaterialPropertyVector* WLSIntegral = const_cast<G4MaterialPropertyVector*>(WLSIntegral_) ;

    double mx = WLSIntegral->GetMaxValue() ;
    std::cerr
        << "U4OptWLS::CreateWavelengthSamples"
        << " WLSIntegral.max*1e9 "
        << std::fixed << std::setw(10) << std::setprecision(4) << mx*1e9
        ;

    NP* wl = NP::Make<double>(num_samples);
    wl->fill<double>(0.);
    double* wl_v = wl->values<double>() ;

    for(int i=0 ; i < num_samples ; i++)
    {
        G4double u = G4UniformRand() ;
        G4double CIIvalue = u*mx;
        G4double sampledEnergy = WLSIntegral->GetEnergy(CIIvalue);  // from value to domain

        G4double sampledWavelength_nm = h_Planck*c_light/sampledEnergy/nm ;

        wl_v[i] = sampledWavelength_nm ;

        if( i < 10 ) std::cout
            << " sampledEnergy/eV "
            << std::fixed << std::setw(10) << std::setprecision(4) << sampledEnergy/eV
            << " sampledWavelength_nm "
            <<  std::fixed << std::setw(10) << std::setprecision(4) << sampledWavelength_nm
            << std::endl
            ;
    }
    return wl ;
}


/**
U4OptWLS::CreateGeant4InterpolatedInverseCDF
-----------------------------------------------------

Reproducing the results of Geant4 dynamic bin finding interpolation
using GPU texture lookups demands very high resolution textures for some
ICDF shapes. This function prepares a three item buffer that can be used
to create a 2D texture that effectively mimmicks variable bin sizing even
though GPU hardware does not support that without paying the cost of
high resolution across the entire range.

* item 0 : full range "standard" resolution
* item 1: left hand side high resolution
* item 2: right hand side high resolution

::

    hd_factor                LHS            RHS
    10          10x bins:    0.00->0.10     0.90->1.00
    20          20x bins:    0.00->0.05     0.95->1.00


The ICDF is formed using Geant4s "domain lookup from value" functionality
in the form of G4MaterialPropertyVector::GetEnergy

::

    g4-cls G4MaterialPropertyVector

    096 G4double G4PhysicsOrderedFreeVector::GetEnergy(G4double aValue)
     97 {
     98         G4double e;
     99         if (aValue <= GetMinValue()) {
    100           e = edgeMin;
    101         } else if (aValue >= GetMaxValue()) {
    102           e = edgeMax;
    103         } else {
    104           size_t closestBin = FindValueBinLocation(aValue);
    105           e = LinearInterpolationOfEnergy(aValue, closestBin);
    106     }
    107         return e;
    108 }

    118 G4double G4PhysicsOrderedFreeVector::LinearInterpolationOfEnergy(G4double aValue,
    119                                  size_t bin)
    120 {
    121         G4double res = binVector[bin];
    122         G4double del = dataVector[bin+1] - dataVector[bin];
    123         if(del > 0.0) {
    124           res += (aValue - dataVector[bin])*(binVector[bin+1] - res)/del;
    125         }
    126         return res;
    127 }




                                                        1  (x1,y1)     (  binVector[bin+1], dataVector[bin+1] )
                                                       /
                                                      /
                                                     *  ( xv,yv )       ( res, aValue )
                                                    /
                                                   /
                                                  0  (x0,y0)          (  binVector[bin], dataVector[bin] )


              Similar triangles::

                 xv - x0       x1 - x0
               ---------- =   -----------
                 yv - y0       y1 - y0




                  res - binVector[bin]             binVector[bin+1] - binVector[bin]
               ----------------------------  =     -----------------------------------
                 aValue - dataVector[bin]          dataVector[bin+1] - dataVector[bin]


                                                                              binVector[bin+1] - binVector[bin]
                   res  = binVector[bin] +  ( aValue - dataVector[bin] ) *  -------------------------------------
                                                                              dataVector[bin+1] - dataVector[bin]

                                                   x1 - x0
                   xv  =    x0  +   (yv - y0) *  --------------
                                                   y1 - y0

**/

inline NP* U4OptWLS::CreateGeant4InterpolatedInverseCDF(
       const G4MaterialPropertyVector* WLSIntegral_,
       int num_bins,
       int hd_factor,
       const char* material_name,
       bool energy_not_wavelength
)   // static
{

    assert( material_name );

    G4MaterialPropertyVector* WLSIntegral= const_cast<G4MaterialPropertyVector*>(WLSIntegral_) ;  // tut tut : G4 GetMaxValue() GetEnergy() non-const
    double mx = WLSIntegral->GetMaxValue() ;   // dataVector.back(); because its **ORDERED** to be increasing on Insert


    // hmm more extensible (eg for Cerenkov [BetaInverse,u,payload] icdf)
    // with the 3 for the different resolutions to be in the payload rather than as separate items ?
    // would of course use 4 to map to float4 after narrowing


    NP* icdf = NP::Make<double>(3, num_bins, 1);
    icdf->fill<double>(0.);

    int ni = icdf->shape[0];
    int nj = icdf->shape[1];
    int nk = icdf->shape[2];

    assert( ni == 3 );
    assert( nk == 1 );
    int k = 0 ;

    assert( hd_factor == 10 || hd_factor == 20 );
    double edge = 1./double(hd_factor) ;

    icdf->names.push_back(material_name) ;  // match X4/GGeo
    icdf->set_meta<std::string>("name", material_name );

    icdf->set_meta<std::string>("creator", "U4OptWLS::CreateGeant4InterpolatedInverseCDF" );
    icdf->set_meta<int>("hd_factor", hd_factor );
    icdf->set_meta<int>("num_bins", num_bins );
    icdf->set_meta<double>("edge", edge );


    if(VERBOSE) std::cerr
        << "U4OptWLS::CreateGeant4InterpolatedInverseCDF"
        << " num_bins " << num_bins
        << " hd_factor " << hd_factor
        << " mx " << std::fixed << std::setw(10) << std::setprecision(4) << mx
        << " mx*1e9 " << std::fixed << std::setw(10) << std::setprecision(4) << mx*1e9
        << " edge " << std::fixed << std::setw(10) << std::setprecision(4) << edge
        << " icdf " << icdf->sstr()
        << std::endl
        ;

    for(int j=0 ; j < nj ; j++)
    {
        double u_all = double(j)/double(nj) ;                         // 0 -> (nj-1)/nj = 1-1/nj
        double u_lhs = double(j)/double(hd_factor*nj) ;               // hd_factor=10(20) u_lhs=0.0->0.1  (0.0  ->~0.05)
        double u_rhs = 1. - edge + double(j)/double(hd_factor*nj) ;   // hd_factor=10(20) u_rhs=0.9->~1.0 (0.95 ->~1.0)
        // u_lhs and u_rhs mean there are hd_factor more lookups at the extremes

        double energy_all = WLSIntegral->GetEnergy( u_all*mx );
        double energy_lhs = WLSIntegral->GetEnergy( u_lhs*mx );
        double energy_rhs = WLSIntegral->GetEnergy( u_rhs*mx );

        double wavelength_all = h_Planck*c_light/energy_all/nm ;
        double wavelength_lhs = h_Planck*c_light/energy_lhs/nm ;
        double wavelength_rhs = h_Planck*c_light/energy_rhs/nm ;

        double v_all = energy_not_wavelength ? energy_all :  wavelength_all ;
        double v_lhs = energy_not_wavelength ? energy_lhs :  wavelength_lhs ;
        double v_rhs = energy_not_wavelength ? energy_rhs :  wavelength_rhs ;

        icdf->set<double>(v_all, 0, j, k );
        icdf->set<double>(v_lhs, 1, j, k );
        icdf->set<double>(v_rhs, 2, j, k );
    }
    return icdf ;
}

