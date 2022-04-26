//******************************************************************************************
//* Sakktábla
//* Készítette: Gaf
//* Web: www.raytracer.hu
//* 2008.05.18
//******************************************************************************************

#include "colors.inc"
#include "shapes.inc"

#declare CameraFrom  = <6.75, 5.25, -9>;        // Kamera honnan néz
#declare CameraAt    = <2, .2, -3>;             // Kamera hová néz
#declare CameraAngle = 49;                      // Kamera látószöge
#declare LightFrom = <-17, 10.4, 10.4>;         // Fényforrás helye
#declare NumberOfPhotons = 150000;              // Fotonok mennyisége


global_settings { 
 max_trace_level 10 
 photons {
 count NumberOfPhotons
  spacing 0.00003
  autostop 0
  jitter .4 }} 


camera { 
 location CameraFrom 
 look_at CameraAt 
 angle CameraAngle 
 right 4/3*x 
// aperture 0.0625          // Rekesznyílás
// blur_samples 32          // Elmosás minõsége (mintavételezési sûrûség)
// focal_point <1,0,-1>     // Ez az ami éles
// confidence 0.95          // Felsõ korlát mintavételezéshez
}


sky_sphere {
 pigment {
  gradient z
  color_map {
   [0.0 rgb <0.6,0.7,1.0>]
   [0.2 rgb <0.2,0.3,0.9>] }}}


light_source {
 LightFrom
 color rgb <1, 1, .9>
 adaptive 50
 spotlight
 radius 7
 falloff 22
 point_at <3, 0, -3>
 area_light <15, 0, 0>, <0, 0, 15>, 12, 12 
 photons{
  reflection on
  refraction on }}   


// *** FIGURÁK ANYAGA ***
// Sötét figurák
#declare DarkPlayerTexture = texture {  
 pigment { color rgb<1,0.8392,0> filter 0.98 }
 finish { ambient .2 diffuse .6 phong .15 phong_size 5 reflection 0.15 } };
#declare DarkPlayerInterior = interior { ior 1.51714 caustics 1.0 fade_power 12 };
// Világos figurák
#declare LightPlayerTexture = texture {  
 pigment { color rgb<1,1,1> filter 0.98 }
 finish { ambient .2 diffuse .6 phong .15 phong_size 5 reflection 0.15 } };
#declare LightPlayerInterior = interior { ior 1.51714 caustics 1.0 fade_power 12 };


// *** FIGURÁK ***
// Testek modellezése
#declare Figura_Paraszt = lathe { bezier_spline 20, <0, 9>, <3, 9>, <2.5, 5>, <1, 5>, <1, 5>, <1, 2>, <2, 0>, <3, 0>, <3, 0>, <3.5, 0>, <3.5, -0.5>, <3, -0.5>, <3, -0.5>, <3, -1>, <3.5, -2>, <3, -2>, <3, -2>, <1, -2>, <1, -2>, <0, -2> translate <0,2,0> };
#declare Figura_Futo = lathe { bezier_spline 24, <0, 12>, <3, 12>, <0.7, 11>, <0.7, 9>, <0.7, 9>, <3, 9>, <2.5, 5>, <1, 5>,  <1, 5>, <1, 2>, <2, 0>, <3, 0>, <3, 0>, <3.5, 0>, <3.5, -0.5>, <3, -0.5>, <3, -0.5>, <3, -1>, <3.5, -2>, <3, -2>, <3, -2>, <1, -2>, <1, -2>, <0, -2> translate <0,2,0> };
#declare Figura_Bastya = lathe { bezier_spline 28, <0, 9>, <2, 9>, <0, 9>, <2.5, 9>, <2.5, 9>, <3, 5.5>, <2.75, 7>, <3, 5.5>, <3, 5.5>, <1.25, 5>, <3, 5.6>, <1.25, 5>, <1.25, 5>, <1.25, 2>, <2, 0>, <3, 0>, <3, 0>, <3.5, 0>, <3.5, -0.5>, <3, -0.5>, <3, -0.5>, <3, -1>, <3.5, -2>, <3, -2>, <3, -2>, <1, -2>, <1, -2>, <0, -2> translate <0,2,0> }; 
#declare Figura_Queen = lathe { bezier_spline 28, <0, 14>, <2, 14>, <2.5, 12>, <0.5, 12>, <0.5, 12>, <5, 12>, <1.25,10>, <1.25, 6>, <1.25, 6>, <1.5, 6>, <3.5, 4>, <1.25, 4>, <1.25, 4>, <1.25, 2>, <2, 0>, <3, 0>, <3, 0>, <3.5, 0>, <3.5, -0.5>, <3, -0.5>, <3, -0.5>, <3, -1>, <3.5, -2>, <3, -2>, <3, -2>, <1, -2>, <1, -2>, <0, -2> translate <0,2,0> }; 
#declare Figura_King = merge { object { lathe { bezier_spline 28, <0.5, 12>, <0.6, 12>, <0, 12>, <0.5, 12>, <0.5, 12>, <5, 12>, <1.25,10>, <1.25, 6>, <1.25, 6>, <1.5, 6>, <3.5, 4>, <1.25, 4>, <1.25, 4>, <1.25, 2>, <2, 0>, <3, 0>, <3, 0>, <3.5, 0>, <3.5, -0.5>, <3, -0.5>, <3, -0.5>, <3, -1>, <3.5, -2>, <3, -2>, <3, -2>,  <1, -2>, <1, -2>, <0, -2> translate <0,2,0> } } cylinder { <0,12,0>, <0,16.5,0>, 0.5} cylinder { <-1,15.5,0>, <1,15.5,0>, 0.5}};
#declare Figura_Horse = difference {
 object { merge { object { lathe { bezier_spline 12, <3, -2>, <1, -2>, <1, -2>, <0, -2>, <3, -2>, <4, -2>, <4, -1>, <3, -1>, <3, -1>, <2, -1>, <1.5, 1>, <1.5, 3> translate <0,2,0> } } sphere { <0, 6,0>, 3 } sphere { <0, 11,0>, 3 } cylinder { <0, 6,0>, <0, 11,0>, 3 }} }
 cylinder { <21, 14,-4>, <21, 14, 4>, 20 }                                      // Jobb oldal leválasztása
 cylinder { <-21, 14,-4>, <-21, 14, 4>, 20 }                                    // Bal oldal leválasztása
 cylinder { <-4, 15,-3.7>, <4, 15, -3.7>, 3.4 }                                 // orrnyereg kifaragása
 cylinder { <-4, 14,-1>, <4, 14, -1>, 1 }                                       // Fülek faragása
 cylinder { <-4, 15.5, 4>, <4, 15.5, 4>, 4 }                                    // Tarkó kifaragása
 cylinder { <-4, 10,-0.6>, <4, 10, -0.6>, .5 }                                  // Az álla...
 box { <-4, 10.5, 2>, <4, 9.5, -2> rotate <-15,0,0> translate <0, -0.2, 0> }    // Mellkas 1
 box { <-4, 10.5, 2>, <4, 9.5, -2> rotate <-45,0,0> translate <0, 1.5, 5> }     // Mellkas 2
 cylinder { <-4, 9,-3>, <4, 9, -3>, 1 }                                         // Bugfix a mellkason
 cylinder { <0, 14,-4>, <0, 14, 4>, .5 }};                                      // Fülek szétválasztása


// Komplett testek definiálása (modell + anyag)
#declare Sotet_Paraszt = object { Figura_Paraszt texture { DarkPlayerTexture } hollow on interior { DarkPlayerInterior } photons { target reflection on refraction on } scale .1 };
#declare Sotet_Futo    = object { Figura_Futo texture { DarkPlayerTexture } hollow on interior { DarkPlayerInterior } photons { target reflection on refraction on } scale .1 };
#declare Sotet_Bastya  = object { Figura_Bastya texture { DarkPlayerTexture } hollow on interior { DarkPlayerInterior } photons { target reflection on refraction on } scale .1 };
#declare Sotet_Queen   = object { Figura_Queen texture { DarkPlayerTexture } hollow on interior { DarkPlayerInterior } photons { target reflection on refraction on } scale .1 }; 
#declare Sotet_King    = object { Figura_King texture { DarkPlayerTexture } hollow on interior { DarkPlayerInterior } photons { target reflection on refraction on } scale .1 }; 
#declare Sotet_Horse   = object { Figura_Horse texture { DarkPlayerTexture } hollow on interior { DarkPlayerInterior } photons { target reflection on refraction on } scale .1 }; 
#declare Vilagos_Paraszt = object { Figura_Paraszt texture { LightPlayerTexture } hollow on interior { LightPlayerInterior } photons { target reflection on refraction on } scale .1 }; 
#declare Vilagos_Futo    = object { Figura_Futo texture { LightPlayerTexture } hollow on interior { LightPlayerInterior } photons { target reflection on refraction on } scale .1 }; 
#declare Vilagos_Bastya  = object { Figura_Bastya texture { LightPlayerTexture } hollow on interior { LightPlayerInterior } photons { target reflection on refraction on } scale .1 }; 
#declare Vilagos_Queen   = object { Figura_Queen texture { LightPlayerTexture } hollow on interior { LightPlayerInterior } photons { target reflection on refraction on } scale .1 }; 
#declare Vilagos_King    = object { Figura_King texture { LightPlayerTexture } hollow on interior { LightPlayerInterior } photons { target reflection on refraction on } scale .1 }; 
#declare Vilagos_Horse   = object { Figura_Horse texture { LightPlayerTexture } hollow on interior { LightPlayerInterior } photons { target reflection on refraction on } scale .1 };


// *** SZÍNPAD ***
// Sötét bábúk
object { Sotet_Paraszt translate <3.5, 0, 1.5> }
object { Sotet_Paraszt translate <2.5, 0, 2.5> }
object { Sotet_Paraszt translate <1.5, 0, 2.5> }        
object { Sotet_Paraszt translate <0.5, 0, 1.5> }                        // E6
object { Sotet_Paraszt translate <-0.5, 0, 0.5> }                       // D5
object { Sotet_Paraszt translate <2.5, -0.3, -5.5> }                    // Kiütve       
object { Sotet_Paraszt translate <-2.5, 0, 2.5> }
object { Sotet_Paraszt translate <-3.5, 0, 2.5> }
object { Sotet_Horse rotate <0,-30,0> translate <-1.5, 0, 1.5> }        // C6
object { Sotet_Horse rotate <0,-15,0> translate <2.5, 0, 1.5> }         // G6 
object { Sotet_Bastya translate <-3.5, 0, 3.5> }
object { Sotet_Bastya translate <3.5, 0, 3.5> }
object { Sotet_Futo translate <-0.5, 0, 2.5> }                          // D7
object { Sotet_Futo translate <1, -0.3, -5.7> }                         // Kiütve
object { Sotet_Queen translate <-3.5, 0, 0.5> }                         // A5
object { Sotet_King translate <0.5, 0, 3.5> }
 
// Világos bábúk
object { Vilagos_Paraszt translate <3.5, 0, -1.5> }                     // H4
object { Vilagos_Paraszt translate <2.5, 0, -1.5> }                     // G4
object { Vilagos_Paraszt translate <1.5, 0, -2.5> }
object { Vilagos_Paraszt translate <0.5, 0, 0.5> }                      // E5
object { Vilagos_Paraszt translate <3.6, -0.3, -5.52> }                 // Kiütve
object { Vilagos_Paraszt translate <-1.5, 0, -1.5> }                    // C3
object { Vilagos_Paraszt translate <-2.5, 0, -2.5> }
object { Vilagos_Paraszt translate <-3.5, 0, -2.5> }
object { Vilagos_Horse rotate <0,110,0> translate <5.6, -0.3, -3.5>  }  // Kiütve
object { Vilagos_Horse rotate <0,195,0> translate <1.5, 0, -1.5> }      // F3 
object { Vilagos_Bastya translate <-3.5, 0, -3.5> }
object { Vilagos_Bastya translate <3.5, 0, -3.5> }
object { Vilagos_Futo translate <-0.5, 0, -1.5> }                       // D3
object { Vilagos_Futo translate <0.5, 0, -1.5> }                        // E3
object { Vilagos_Queen translate <-0.5, 0, -3.5> }      
object { Vilagos_King translate <0.5, 0, -3.5> }


// *** ASZTALLAP ***
cylinder { <0, -.3, 0>, <0, -2, 0>, 20
 pigment {
  image_map { tga "textura_asztallap.TGA"
  interpolate 2 }
  scale 10
  rotate <90,73,0> }
 normal {
  bump_map { tga "bumpmap_asztallap.TGA"
  interpolate 2 }
  scale 10
  rotate <90,73,0>
  bump_size 0.5 }
 finish { ambient .05 diffuse .5 reflection 0.05 }}

 
// *** HÁTTÉR ***
cylinder { <0, -25, 0>, <0, 25, 0>, 100
 pigment {
      image_map { tga "textura_hatter.TGA"
      interpolate 2
      map_type 2 }
      scale 40
      translate <0,-10,0> }
 hollow on
 finish { ambient .15 diffuse .5 }}

 
// *** SAKKTÁBLA ***
// Világos kockák felülete
#declare LightCubeTexture = texture {  
 pigment {
  image_map { tga "textura_light.TGA"
  interpolate 2 }
  scale 10
  rotate <90,0,0> }
 finish { ambient .2 diffuse .6 phong .85 phong_size 50 reflection 0.025 } };
// Sötét kockák felülete
#declare DarkCubeTexture = texture {  
 pigment {
  image_map { tga "textura_dark.TGA"
  interpolate 2 }
  scale 10
  rotate <90,0,0> }
 finish { ambient .2 diffuse .6 phong .85 phong_size 50 reflection 0.025 } };
// Keret felülete
#declare MidCubeTexture = texture {  
 pigment {
  image_map { tga "textura_mid.TGA"
  interpolate 2 }
  scale 10
  rotate <90,0,0> }
 finish { ambient .2 diffuse .6 phong .85 phong_size 50 reflection 0.025 } };

// Sakktábla kockái
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <3,0,4>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <1,0,4>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <-1,0,4>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <-3,0,4>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <4,0,4>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <2,0,4>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <0,0,4>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <-2,0,4>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <3,0,3>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <1,0,3>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <-1,0,3>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <-3,0,3>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <4,0,3>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <2,0,3>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <0,0,3>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <-2,0,3>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <3,0,2>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <1,0,2>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <-1,0,2>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <-3,0,2>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <4,0,2>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <2,0,2>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <0,0,2>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <-2,0,2>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <3,0,1>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <1,0,1>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <-1,0,1>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <-3,0,1>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <4,0,1>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <2,0,1>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <0,0,1>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <-2,0,1>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <3,0,0>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <1,0,0>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <-1,0,0>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <-3,0,0>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <4,0,0>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <2,0,0>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <0,0,0>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <-2,0,0>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <3,0,-1>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <1,0,-1>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <-1,0,-1>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <-3,0,-1>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <4,0,-1>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <2,0,-1>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <0,0,-1>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <-2,0,-1>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <3,0,-2>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <1,0,-2>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <-1,0,-2>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <-3,0,-2>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <4,0,-2>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <2,0,-2>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <0,0,-2>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <-2,0,-2>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <3,0,-3>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <1,0,-3>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <-1,0,-3>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { DarkCubeTexture } translate <-3,0,-3>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <4,0,-3>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <2,0,-3>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <0,0,-3>}
object { Round_Box_Union(<-1, -0.01, -1>, < 0,  -1,  0>, .01) texture { LightCubeTexture } translate <-2,0,-3>}
// Keret
object { Wire_Box_Union(<-4.02, -0.02, -4.02>, < 4.02,  -1,  4.02>, .01) texture { MidCubeTexture } translate <0,0,0>}
difference {
 box { <-4.5, -0.02, -4.5>, < 4.5,  -1,  4.5> }
 box { <-4.01, 1, -4.01>, < 4.01,  -2,  4.01> }
 texture { MidCubeTexture }}
object { Wire_Box_Union(<-4.6, -0.05, -4.6>, < 4.6,  -1,  4.6>, .1) texture { MidCubeTexture } translate <0,0,0>}
 box { <-4.8, -0.1, -4.8>, < 4.8,  -2,  4.8> texture { MidCubeTexture } }
object { Wire_Box_Union(<-4.9, -0.1, -4.9>, < 4.9,  -0.3,  4.9>, .1) texture { MidCubeTexture } translate <0,0,0>}
