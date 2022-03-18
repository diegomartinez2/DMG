// Insert menu illustration scene
// Created June-August 2001 by Christoph Hormann

// ----- Shapes submenu -----

// -w120 -h48 +a0.1 +am2 -j +r3

#version 3.5;

#include "colors.inc"

/*
#declare Typ=1;    // sphere
#declare Typ=2;    // box
#declare Typ=3;    // cylinder
#declare Typ=4;    // cone
#declare Typ=5;    // disc
#declare Typ=6;    // plane
#declare Typ=7;    // torus
#declare Typ=8;    // superellipsoid
#declare Typ=9;    // mesh
#declare Typ=10;   // fractal
#declare Typ=11;   // blob
#declare Typ=12;   // height_field
#declare Typ=13;   // height_field_img
#declare Typ=14;   // polygon
#declare Typ=15;   // lathe
#declare Typ=16;   // prism
#declare Typ=17;   // text
#declare Typ=18;   // sor
#declare Typ=19;   // bicubic_patch
#declare Typ=20;   // quartic
#declare Typ=21;   // poly
#declare Typ=22;   // quadric
#declare Typ=23;   // isosurface
#declare Typ=24;   // sphere_sweep
#declare Typ=25;   // triangle
#declare Typ=26;   // smooth_triangle
#declare Typ=27;   // cubic
#declare Typ=28;   // parametric
*/

global_settings {
  assumed_gamma 1
  max_trace_level 5
}

light_source {
  <1.5, 0.35, 1.0>*10000
  color rgb 1.0
}

camera {
  location    <10, 20, 5.5>
  direction   y
  sky         z
  up          z
  right       (120/48)*x
  look_at     <0, 0, -0.2>
  angle       28
}

sky_sphere {
  pigment {
    gradient y
    color_map {
      [0.0 rgb <0.6,0.7,1.0>]
      [1.0 rgb <0.2,0.2,0.8>]
    }
  }
}

// ----------------------------------------

plane
{
  z, -1
  texture
  {
    pigment {checker color rgb 1 color rgb <0,0,0.1> scale 1.1}
    finish {
      diffuse 0.7
      specular 0.4
      roughness 0.01

      reflection { 0.5 , 1.0
        fresnel on
        metallic 0.8
      }
      conserve_energy

    }
  }
}

//---------------------------------------------------

#declare Tex1=
  texture {
    pigment { color ForestGreen }
    finish { specular 0.4 roughness 0.02 }
  }

#declare Tex2=
  texture {
    pigment { color Coral }
    finish { specular 0.4 roughness 0.02 }
  }

#declare Tex3=
  texture {
    pigment { color SpicyPink }
    finish { specular 0.4 roughness 0.02 }
  }


// =============================================

#if(Typ=1)

sphere { <0, 0, 0>, 1.2 texture { Tex1 } }

sphere { <0, 0, 0>, 1.2 texture { Tex2 }
  scale <0.5, 0.5, 1.2>
  translate <-3.5, 1, 0.1>
}

sphere { <0, 0, 0>, 1.2 texture { Tex3 }
  scale <1, 1, 0.5>
  translate <3.8, 0, 0>
}

#end

// =============================================

#if(Typ=2)

box { -0.9, 0.9 texture { Tex1 } }

box { -0.9, 0.9 texture { Tex2 }
  scale <0.5, 0.5, 1.3>
  translate <-3.5, 1, 0>
}

box { -0.9, 0.9 texture { Tex3 }
  scale <1, 1, 0.5>
  translate <3.8, 0, 0>
}

#end

// =============================================

#if(Typ=3)

cylinder { -0.9*z, 0.8*z, 1.2 texture { Tex1 } }

cylinder { -1.0*z, 1.2*z, 0.8
  texture { Tex2 }
  translate <-3.5, 1, 0>
}

cylinder { -0.6*z, 0.0*z, 1.8
  texture { Tex3 }
  scale <0.3, 1.0, 1.0>
  translate <3.8, 0, 0>
}

#end

// =============================================

#if(Typ=4)

cone { -0.9*z, 1.2, 0.8*z, 0.0 texture { Tex1 } }

cone { -1.0*z, 0.4, 1.2*z, 1.5
  texture { Tex2 }
  translate <-3.5, 1, 0>
}

cone { -0.6*z, 1.8, 0.8*z, 0.5
  texture { Tex3 }
  scale <0.4, 1.0, 1.0>
  translate <4, 0, 0>
}

#end

// =============================================

#if(Typ=5)

disc { <0, 0, -0.3>, z, 2.4, 0.6 texture { Tex1 } }

disc { <0, 0, 0>, z, 1.8, 0.0
  texture { Tex2 }
  translate <-3.5, 1, 0.4>
}

disc { <0, 0, 0>, z, 1.4, 0.4
  texture { Tex3 }
  rotate 45*y
  translate <3.8, 0, 0.3>
}

#end

// =============================================

#if(Typ=6)

plane {
  z, 0
  texture { Tex1 }
  rotate 20*y
  translate -1*z
}

#end

// =============================================

#if(Typ=7)

torus {
  1.4, 0.4
  rotate 90*x
  translate -0.3*z
  texture { Tex1 }
}

torus { 1.0, 0.2
  texture { Tex2 }
  translate <-3.5, 1, 0.2>
}

torus { 1.2, 0.3
  texture { Tex3 }
  rotate 90*z
  scale <1, 2, 0.6>
  rotate -30*y
  translate <3.8, 0, 0>
}

#end

// =============================================

#if(Typ=8)

superellipsoid { <0.8, 0.2>
  scale 1.3
  translate -0.1*z
  texture { Tex1 }
}

superellipsoid { <0.1, 1.0>
  texture { Tex2 }
  translate <-3.5, 1, 0.2>
}

superellipsoid { <2.0, 3.0>
  scale 1.3
  texture { Tex3 }
  translate <3.8, 0, 0>
}

#end

// =============================================

#if(Typ=9)

mesh {
  #include "..\scenes\advanced\newltpot\teapot_tri.inc"
  scale 0.8
  translate -0.8*z
  texture { Tex2 }
  translate <-2.8, 1, 0>
}


#declare WHITE_KNIGHT_1 = texture { Tex1 }
#declare BLACK_ROOK_1 = texture { Tex3 }

union {
 #include "..\scenes\objects\chess.inc"
 rotate 90*x
 rotate 90*z
 scale 0.25
 translate <2.4, 0, 0.25>
}

#end

// =============================================

#if(Typ=10)

julia_fractal
{
  <-0.083,0.0,-0.83,-0.025>
  hypercomplex
  sqr
  max_iteration 8
  precision 120

  scale 1.3
  translate -0.1*z
  texture { Tex1 }
}

julia_fractal {
  <.55,.5,-.34,.1>

  quaternion
  max_iteration 7
  precision 80
  texture { Tex2 }
  translate <-3.5, 1, 0.2>
}

julia_fractal {
  <-.2,.4,-.34,.2>

  quaternion
  sqr
  max_iteration 8
  precision 120

  scale 1.3
  texture { Tex3 }
  translate <3.4, 0, 0>
}

#end

// =============================================

#if(Typ=11)

blob {
  threshold 0.6
  sphere { <0.75, 0, 0>, 1.0, 1.0 }
  sphere { <-0.375, 0.64952, 0>, 1.0, 1.0 }
  sphere { <-0.375, -0.64952, 0>, 1.0, 1.0 }

  scale 1.5
  translate -0.1*z
  texture { Tex1 }
}


blob {
  threshold 0.5
  sphere { <0.8, 1.8, -0.6>, 0.8, 1.0 }
  sphere { <0, 1, 0>, 2.2, 1.0 }
  sphere { <0, 1, 1.5>, 1.0, -2.0 }

  texture { Tex2 }
  translate <-3.3, 1, 0.2>
}

blob {
  threshold 0.5
  cylinder { <-7, 0, 0>, <7, 0, 0>, 4, 2 }
  cylinder { <0, 0, -7>, <0, 0, 7>, 4, 2 }
  sphere { <0, 3, 0>, 2.5, -4 }

  scale 0.13

  texture { Tex3 }
  translate <3.8, 0, 0.3>
}

#end

// =============================================

#if(Typ=12)

height_field
{
  function 200, 200 {
    pigment {
    marble
    sine_wave
    color_map { [0 rgb 1][1 rgb 0] }
    turbulence 0.7
    }
  }
  translate -0.5
  rotate 90*x
  scale <3, 3, 0.8>
  scale 2.2
  texture { Tex1 }
}

#end

// =============================================

#if(Typ=13)

height_field
{
  png "plasma3.png"
  texture
  {
    pigment
    {
      image_map { png "plasma3.png" map_type 0 interpolate 2 once }
      rotate x*90
    }
  }
  translate -0.5
  rotate 90*x
  scale <3, 3, 0.8>
  scale 2.2
}

#end

// =============================================

#if(Typ=14)

polygon
{
  14,
  <0, 0>, <0, 6>, <1, 6>, <1, 5>, <4, 5>, <4, 2>, <1, 2>, <1, 0>, <0, 0>,
  <1, 3>, <1, 4>, <3, 4>, <3, 3>, <1, 3>

  scale 0.8
  translate <-1, -1.6, -0.2>
  texture { Tex1 }
}

polygon {
  16,
  <0, 0>, <1, 0>,<0, 1>,<0, 0>,
  <.3, .7>,<.5, .7>,<.3, .9>,<.3, .7>,
  <.5, .5>,<.8, .5>,<.5, .8>,<.5, .5>,
  <.2, .2>,<.5, .2>,<.2, .5>,<.2, .2>

  scale 3
  rotate 90*y
  texture { Tex2 }
  translate <-3.5, 0, 1.2>
}

polygon {
  7,
  <0, 0>, <1, 0>,<1.5, 0.5>,<1, 1>,<0, 1>,<-0.5, 0.5>,<0, 0>

  scale 1.6
  rotate 90*x
  texture { Tex3 }
  translate <3, 0, -0.8>
}

#end

// =============================================

#if(Typ=15)

lathe {
  linear_spline
  12,
  <2, 1>, <2, -1>, <3, -1>, <3.4, -2>, <4, -1.1>, <3.6, -0.9>,
  <2.6, 0>, <3.6.9>, <4, 1.1>, <3.4, 2>, <3, 1>, <2, 1>

  rotate 90*x
  scale 0.4
  translate -0.1*z
  texture { Tex1 }
}


lathe{
  cubic_spline
  9,
  <0.870490,0.570719>,
  <1.325912,0.622602>,
  <1.579565,0.403538>,
  <1.412384,0.069178>,
  <0.709075,0.011530>,
  <0.593778,0.236358>,
  <0.853196,0.357420>,
  <1.124143,0.184475>,
  <1.274028,0.149886>

  rotate -90*z
  scale 0.8
  texture { Tex2 }
  translate <-3.6, 1, 0.4>
}


lathe{
  cubic_spline
  9,
  <0.326558,-0.139954>,
  <0.439854,-0.046651>,
  <0.736422,0.643120>,
  <0.973011,0.999668>,
  <0.353216,0.599801>,
  <0.156615,0.236588>,
  <0.236588,0.139954>,
  <0.356548,0.403200>,
  <0.513163,0.466512>

  rotate 90*x
  scale 1.2
  texture { Tex3 }
  translate <3.8, 0, -1>
}

#end

// =============================================

#if(Typ=16)

prism {
  conic_sweep
  linear_spline
  0.6,
  1.0,
  10,
  <0.872689,0.705541>,
  <1.024002,0.510241>,
  <0.895561,0.293829>,
  <0.566544,0.555987>,
  <0.279753,0.446901>,
  <0.661554,0.012316>,
  <0.029911,0.168907>,
  <0.367726,0.710819>,
  <0.800551,0.578860>,
  <0.872689,0.705541>

  rotate 90*x
  rotate 90*z
  scale 4.4
  translate <-0.8,0, -3.4>
  texture { Tex1 }
}


prism {
  cubic_spline
  linear_sweep
  0.5,
  1.0,
  13,

  <-0.2, -1.0>,
  < 0.2, -1.0>, < 0.2,  0.2>, < 1.0, -0.2>, < 1.0,  0.2>, < 0.2,  1.0>,
  <-0.2,  1.0>, <-1.0,  0.2>, <-1.0, -0.2>, <-0.2,  0.2>, <-0.2, -1.0>,
  < 0.2, -1.0>,
  < 0.2,  0.2>

  rotate 90*z
  scale 1.0
  texture { Tex2 }
  translate <-3.8, 1, 0.0>
}


prism{
  cubic_spline
  0,
  0.5,
  12,
  <1.061539,0.669709>,
  <0.841635,0.487788>,
  <0.921600,0.027988>,
  <0.351847,0.075967>,
  <0.771665,0.267884>,
  <0.281878,0.247892>,
  <0.081964,0.433812>,
  <0.743677,0.355846>,
  <0.689701,0.673708>,
  <1.061539,0.669709>,
  <0.841635,0.487788>,
  <0.921600,0.027988>

  scale 2
  texture { Tex3 }
  translate <3.0, 0, -1>
}

#end

// =============================================

#if(Typ=17)

text {
  ttf
  "crystal.ttf",
  "text",
  0.3,
  0

  rotate 90*x
  rotate 90*z
  scale 3.0
  translate <-1.6, 0, -0.5>
  texture { Tex1 }
}

text {
  ttf
  "crystal.ttf",
  "PovRay",
  0.2,
  0

  rotate 90*x
  rotate 180*z
  scale 2.0
  translate <5.5, 0, -0.5>
  texture { Tex2 }
}

#end

// =============================================

#if(Typ=18)

sor{
  7,
  <0.306379,-0.135782>,
  <0.905210,-0.027853>,
  <0.988768,0.341195>,
  <0.313342,0.546608>,
  <0.295934,0.793800>,
  <0.295934,0.999213>,
  <0.748539,1.215071>

  rotate 90*x
  scale 2
  translate -0.9*z
  texture { Tex1 }
}


sor {
  13,
  <0.148225, 0.000000>
  <0.189979, 0.000000>
  <0.154489, 0.096033>
  <0.075157, 0.123173>
  <0.070981, 0.164927>
  <0.129436, 0.223382>
  <0.070981, 0.281837>
  <0.108559, 0.354906>
  <0.075157, 0.436326>
  <0.175365, 0.536534>
  <0.240084, 0.672234>
  <0.327766, 1.000000>
  <0.290188, 1.000000>
  open

  rotate 90*x
  scale 2.5
  texture { Tex2 }
  translate <-3.5, 1, -0.9>
}

sor {
  12,
  <.517241379, -.132625995>
  <.249336870, 0.000000>
  <.068965517, .031830239>
  <.021220159, .050397878>
  <.058355438, .347480106>
  <.132625995, .381962865>
  <.196286472, .464190981>
  <.238726790, .602122016>
  <.249336870, .721485411>
  <.233421751, .864721485>
  <.167108753, 1.000000000>
  <.084880637, 1.055702918>
  open

  rotate 90*x
  scale 2.5
  translate -0.9*z
  translate <3.8, 0, 0>
  texture { Tex3 }
}

#end

// =============================================

#if(Typ=19)

bicubic_patch { type 1 flatness 0.0100 u_steps 3 v_steps 3,
  <-2.000000, -2.000000, 0.000000>, <-1.000000, -2.000000, 0.000000>, <0.000000, -2.000000, 0.000000>, <1.000000, -2.000000, 0.000000>,
  <-2.000000, -1.000000, 0.000000>, <-1.000000, -1.000000, 0.000000>, <0.078802, -1.000000, -5.358550>, <1.000000, -1.000000, 0.000000>,
  <-2.000000, 0.000000, 0.000000>, <-1.000000, 0.000000, 0.000000>, <0.000000, 0.000000, 0.000000>, <1.000000, 0.000000, 0.000000>,
  <-2.000000, 1.000000, 0.000000>, <-1.000000, 1.000000, 0.000000>, <0.000000, 1.000000, 0.000000>, <1.000000, 1.000000, 0.000000>


  scale <1, 1, -1>*1.2
  rotate 90*z
  translate <-0.6, 0, -0.5>
  texture { Tex1 }
}

bicubic_patch {
  type 1 flatness 0.1  u_steps 8  v_steps 8
  < 0.0, 0.0, 2.0>, < 1.0, 0.0, 0.0>, < 2.0, 0.0, 0.0>, < 3.0, 0.0, -2.0>,
  < 0.0, 1.0, 0.0>, < 1.0, 1.0, 0.0>, < 2.0, 1.0, 0.0>, < 3.0, 1.0,  0.0>,
  < 0.0, 2.0, 0.0>, < 1.0, 2.0, 0.0>, < 2.0, 2.0, 0.0>, < 3.0, 2.0,  0.0>,
  < 0.0, 3.0, 2.0>, < 1.0, 3.0, 0.0>, < 2.0, 3.0, 0.0>, < 3.0, 3.0, -2.0>

  scale 0.6
  translate 0.1*z
  texture { Tex2 }
  translate <-3.7, 1, 0>
}

bicubic_patch { type 1 flatness 0.0100 u_steps 3 v_steps 3,
  <-1.0, -2.0, -1.0>, <-1.0, -2.0, -0.0>, < 1.0, -2.0, -0.0>, < 1.0, -2.0, -1.0>,
  <-1.0,  1.0,  0.0>, <-1.0,  0.0,  0.0>, < 1.0,  0.0,  0.0>, < 1.0,  1.0,  0.0>,
  <-1.0,  1.0,  6.0>, <-1.0,  0.0,  6.0>, < 1.0,  0.0,  6.0>, < 1.0,  1.0,  6.0>,
  <-1.0, -2.0,  7.0>, <-1.0, -2.0,  6.0>, < 1.0, -2.0,  6.0>, < 1.0, -2.0,  7.0>


  scale 0.6
  rotate 90*x
  rotate 90*y
  texture { Tex3 }
  translate <4.2, 0.6, -0.2>
}

#end

// =============================================

#if(Typ=20)

quartic {
   < 0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.01,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.01, 0.0, -0.01 >

  scale 1.2
  translate 0.2*z
  texture { Tex1 }
}

quartic {
   < 4.0,  0.0,  0.0,  -4.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0,
     0.0,  0.0,  0.0,   0.0, 0.0,  0.0,  0.0,  0.0,  0.0, 0.0,
     0.0,  0.0,  0.0,   0.0, 0.0,  1.0,  0.0,  0.0,  0.0, 0.0,
     0.0,  0.0,  1.0,   0.0, 0.0 >

  scale 1.8
  rotate 90*y
  texture { Tex2 }
  translate <-3.7, 0.2, 0.6>
}

quartic {
  < 1.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0, -1.0,
    0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0, 0.0,  1.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  1.0,  0.0, 0.0 >

  scale 1.5
  texture { Tex3 }
  translate <3.6, 0, 0>
}

#end

// =============================================

#if(Typ=21)

 poly
  {6,
   < 4,   0,   0,   0,   0,   0,   0,   0,   0, -4,
     0,   0,   0,   0,   0,   0,   0,   0,   0,  0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,  0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,  0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,  0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,  0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,  0,
     1,   0,   0,   0,   0,   0,   0,   0,   0,  0,
     0,   1,   0,   0>

  scale 1.2
  translate 0.2*z
  texture { Tex1 }
}

 poly
  {6,
   < 4,   0,   0,   0,   0,   0,   0,   0,   0, -4,
     0,   1,   0,   0,   0,   0,   0,   0,   0,  0,
     0,   0,   0,   1,   0,   0,   0,   0,   0,  0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,  0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,  0,
     0,   0,   0,   0,   0,   0,   2,   0,   0,  0,
     0,   -2,   0,   0,   0,   -1,  0,   0,   0,  0,
     1,   0,   0,   0,   0,   0,   0,   0,   0,  0,
     0,   1,   0,   0>

  scale 1.2
  rotate 90*y
  rotate 90*z
  texture { Tex2 }
  translate <-3.7, 0.2, 0.0>
}


 poly
  {6,
   < 4,   0,   0,   0,   0,   0,   0,   0,   0, -4,
     0,   1,   0,   0,   0,   0,   0,   0,   0,  0,
     0,   0,   0,   1,   0,   0,   0,   0,   0,  0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,  0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,  0,
     0,   0,   0,   0,   0,   0,   2,   0,   0,  0,
     0,   -2,   0,   0,   0,   -1,  0,   0,   0,  0,
    -2,   0,   0,   0,   0,   0,   0,   0,   0,  0,
     0,   4,   0,   -1>

  scale 0.9
  rotate 90*x
  texture { Tex3 }
  translate <3.6, 0, 0.2>
}

#end

// =============================================

#if(Typ=22)

quadric {
   <1, 1, 0.4>, <0, 0.4, 0>, <1, 0.5, 0>, 0

  scale 1.4
  translate <1, 2, 0.4>
  texture { Tex1 }
}

quadric {
   <1, 0, 1>, <0, 0, 0>, <0, 0, 0>, -0.3

  scale 1.8
  rotate 90*y
  texture { Tex2 }
  translate <-3.7, 0.2, 0.6>
}

quadric {
   <1, 1, -0.1>, <0, -0.5, 1>, <0, 0, 0>, -0.1

  scale 1.5
  texture { Tex3 }
  translate <3.6, 0, 0>
}

#end

// =============================================

#if(Typ=23)

#include "functions.inc"

isosurface {
  function{ sqrt( x*x+y*y+z*z-1) - f_noise3d( x*3, y*3, z*3 )  }
  contained_by { box { -1.2, 1.2 } }
  accuracy 0.002
  max_gradient 4

  translate 0.1*z
  texture { Tex1 }
}

isosurface {
  function{
   (1+0.002)
   -pow(0.002,sqrt(x*x+y*y)-0.5)
   -pow(0.008,abs(x)+abs(z)-0.5)
   -pow(0.002,sqrt( y*y+z*z)-0.5)
  }
  contained_by { box { -1.2, 1.2 } }
  accuracy 0.003
  max_gradient 5

  texture { Tex2 }
  translate <-3.5, 1, 0.2>
}


isosurface {
  function{
   z-f_ridged_mf(x, y, z, 0.1, 2.8, 8, 0.7, 0.7, 1)*0.5
  }
  contained_by { box { -1.2, 1.2 } }
  accuracy 0.003
  max_gradient 5

  texture { Tex3 }
  translate <3.7, 0, 0>
}

#end

// =============================================

#if(Typ=24)

sphere_sweep {
  linear_spline
  4,
  <-5, -5, 1>, 1.5
  <-5, 5, 0>, 0.7
  < 5, 5, 0>, 0.7
  < 5, -5, 1>, 1.5

  scale 0.3
  translate -0.3*z
  texture { Tex1 }
}


sphere_sweep {
  cubic_spline
  8,
  <-5, -5, -4>, 1.6
  <-5,  5, -4>, 1.6
  < 5,  5,  0>, 1.6
  < 5, -5,  0>, 1.6
  <-5, -5,  4>, 1.6
  <-5,  5,  4>, 1.6
  < 5,  5,  8>, 1.6
  < 5, -5,  8>, 1.6

  rotate 90*y
  scale 0.18
  texture { Tex2 }
  translate <-3.5, 1, 0.3>
}


sphere_sweep {
  cubic_spline
  9,
  < 0,  1,  2>, 0.4
  < 1,  2,  0>, 0.4
  < 2,  0, -1>, 0.4
  < 0, -1, -2>, 0.4
  <-1, -2,  0>, 0.4
  <-2,  0,  1>, 0.4
  < 0,  1,  2>, 0.4
  < 1,  2,  0>, 0.4
  < 2,  0, -1>, 0.4

  scale 0.6
  rotate 85*z
  rotate -45*y
  texture { Tex3 }
  translate <3.8, 0, 0>
}

#end

// =============================================

#if(Typ=25)

triangle {
  <-2,  0, -1>,
  < 2,  0, -1>,
  < 0,  0,  2>

  rotate -60*z
  translate <0, 0.6, 0.1>
  texture { Tex1 }
}


triangle {
  < 1, -1, 0.9>,
  <-2,  0, 1.5>,
  < 2,  1, -1>

  texture { Tex2 }
  translate <-3.5, 1, 0.2>
}


triangle {
  < 1, -1, 0.5>,
  <-2,  3, 0.5>,
  < 2,  4, 0.5>

  texture { Tex3 }
  translate <3.8, 0, 0>
}

#end

// =============================================

#if(Typ=26)

smooth_triangle {
  <-2,  0.2, -1>, <-2, 1, 1>,
  < 2,  0.2, -1>, <2, 1,-1>,
  < 0,  0.4,  2>, <1, 2, 4>

  rotate -60*z
  translate <0, 0.6, 0.1>
  texture { Tex1 }
}


smooth_triangle {
  < 1, -1, 0.9>, <1, 0, 0>,
  <-2,  0, 1.5>, <0, 1, 0>,
  < 2,  1, -1>,  <0, 0, 1>

  texture { Tex2 }
  translate <-3.5, 1, 0.2>
}


smooth_triangle {
  < 1, -1, 0.5>, <-1.8,  -1, 1>,
  <-2,  3, 0.5>, <   0, 2.2, 1>,
  < 2,  4, 0.5>, < 2.2,   0, 1>

  texture { Tex3 }
  translate <3.8, 0, 0>
}

#end

// =============================================

#if(Typ=27)

cubic
{
  <
// x^3,      x^2y,     x^2z,     x^2,
   1,        0,        0,        0,
// xy^2,     xyz,      xy,       xz^2,
   0,        1,        0,        0,
// xz,       x,        y^3,      y^2z,
   0,        0,        1,        0,
// y^2,      yz^2,     yz,       y,
   0,        0,        0,        1,
// z^3,      z^2,      z,        C
   0,        0,        0,        1
  >
  sturm
  bounded_by {box {-1.2,1.2}}
   texture { Tex1 }
}

cubic
{
  <
// x^3,      x^2y,     x^2z,     x^2,
   0.5,        1,        0,        0,
// xy^2,     xyz,      xy,       xz^2,
   0.5,        0,        0.5,        1,
// xz,       x,        y^3,      y^2z,
   0,        1,        0,        0,
// y^2,      yz^2,     yz,       y,
   0,        1,        1,        0,
// z^3,      z^2,      z,        C
   1,        0,        1,        1
  >
  sturm
  bounded_by {box {-1.2,1.2}}
   texture { Tex2 }
 translate -3.8*x
}

cubic
{
  <
// x^3,      x^2y,     x^2z,     x^2,
   1,        0,        0,        0,
// xy^2,     xyz,      xy,       xz^2,
   0,        0,        0,        0,
// xz,       x,        y^3,      y^2z,
   0,        1,        0,        0,
// y^2,      yz^2,     yz,       y,
   0,        0,        0,        1,
// z^3,      z^2,      z,        C
   0,        0,        1,        1
  >
  sturm
  bounded_by {box {-1.2,1.2}}
   texture { Tex3 }
 translate 3.5*x
}

#end

// =============================================

#if(Typ=28)

// Parametric in a spiral or helical shape

parametric {
// (1) like Galileo's helio flying machine...?
  function { u*v*sin (15*v) },            // x-axis
  function { v },                         // y-axis
  function { u*v*cos (15*v) }             // z-axis
  <0,0>, <1,1>
 // contained_by { sphere { 0,1 }} // !texturing problem, use box instead to see!
  contained_by { box { <-1, -1, -1>, <1, 1, 1> } }
 // max_gradient 2
  accuracy 0.01 // 0.001 default, lower slower but better
  precompute 15 x,y,z // precompute normally gives faster rendering (<=20) but longer parsing

  bounded_by {box {-1.2,1.2}}
 scale 2.5 rotate -90*x translate 1.5*z
   texture { Tex1 }
}

// Parametric in a flower shape

parametric {
// (2) like a flower...?
  function { u*v*sin (10*v) },            // x-axis
  function { v*pow(2,u) },                // y-axis
  function { u*v*cos (10*v) }             // z-axis
  <0,0>, <5,5>
  contained_by { box { <-0.5,-1,-0.5>, <0.5,1,0.5> } }
 // max_gradient 2
  accuracy 0.01 // 0.001 default, lower slower but better
  precompute 15 x,y,z // precompute normally gives faster rendering (<=20) but longer parsing

  bounded_by {box {-1.2,1.2}}
 scale 2.5 rotate 90*x translate -1.5*z
   texture { Tex2 }
 translate -3.8*x
}

// Parametric in a kind of ornamental shape

parametric {
// (4) or maybe... ?? (can be very slow and inaccurate)
  function{ sin(v)*sin(u)*pow(cos(v), 3) },     // x-axis
  function{ cos(u)*pow(sin(v), 2) },            // y-axis
  function{ sin(v)*cos(u)*pow(cos(v), 3) }      // z-axis
  <0,0>, <5,5>
  contained_by { box { <-0.5,-1,-0.5>, <0.5,1,0.5> } }
 // max_gradient 2
  accuracy 0.01 // 0.001 default, lower slower but better
  precompute 15 x,y,z // precompute normally gives faster rendering (<=20) but longer parsing

  bounded_by {box {-1.2,1.2}}
 scale 2.5 rotate -90*x
   texture { Tex3 }
 translate 3.5*x
}

#end

// =============================================
