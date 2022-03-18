// Insert menu illustration scene
// Created June-August 2001 by Christoph Hormann
// Modified (added to) October 2001 by Bob Hughes

// ----- patterns submenu -----

// -w120 -h48 +a0.1 +am2 -j +r3

#version 3.5;

#include "colors.inc"

/*
#declare Typ=1;   // agate
#declare Typ=2;   // boxed
#declare Typ=3;   // bozo
#declare Typ=4;   // brick
#declare Typ=5;   // bumps
#declare Typ=6;   // cells
#declare Typ=7;   // checker
#declare Typ=8;   // crackle
#declare Typ=9;   // cylindrical
#declare Typ=10;  // density_file
#declare Typ=11;  // dents
#declare Typ=12;  // fractal
#declare Typ=13;  // function
#declare Typ=14;  // gradient
#declare Typ=15;  // granite
#declare Typ=16;  // hexagon
#declare Typ=17;  // image_pattern
#declare Typ=18;  // leopard
#declare Typ=19;  // marble
#declare Typ=20;  // onion
#declare Typ=21;  // planar
#declare Typ=22;  // quilted
#declare Typ=23;  // radial
#declare Typ=24;  // ripples
#declare Typ=25;  // spherical
#declare Typ=26;  // spiral1
#declare Typ=27;  // spiral2
#declare Typ=28;  // spotted
#declare Typ=29;  // waves
#declare Typ=30;  // wood
#declare Typ=31;  // wrinkles
#declare Typ=32;  // facets
#declare Typ=33;  // average
#declare Typ=34;  // object
#declare Typ=35;  // pigment_pattern
#declare Typ=36;  // slope
*/

#if (Typ=17)     // image_pattern
  #include "logo.inc"
#end

global_settings {
  assumed_gamma 1.0
}

camera {
  orthographic
  location <0,0,1>
  look_at  <0,0,0>
  right    -(120/48)*x
  up y
}

light_source {
  <-1000, 2500, 2000>
  color rgb 1.2
}

#declare P1=pigment {checker
         color rgb <0.85, 0.55, 0.45>, color rgb <0.00, 0.00, 0.42> scale 0.2}
#declare P2=pigment {hexagon
         color rgb <0.85, 0.55, 0.45>, color rgb <0.00, 0.00, 0.42>, color rgb <0.25, 0.40, 0.35> scale 0.3 rotate 90*x}

#if (Typ=36)
union {
  sphere { 0, 0.47 translate <0.5, 0, 0>}
  sphere { 0, 1.7 }

  cylinder { -z*0.5+0.3*(x+y), z*0.5-0.3*(x+y), 0.2 translate <-0.5, 0, 0>}
#else
plane {
  z, 0
#end
  texture {
    pigment {
      #switch (Typ)
        #case (1)
          agate
          #break
        #case (2)
          boxed
          #break
        #case (3)
          bozo
          #break
        #case (4)
          brick
          #break
        #case (5)
          bumps
          #break
        #case (6)
          cells
          #break
        #case (7)
          checker
          #break
        #case (8)
          crackle
          #break
        #case (9)
          cylindrical
          #break
        #case (10)
          density_file df3 "spiral.df3"
          #break
        #case (11)
          dents
          #break
        #case (12)
          mandel 50
            interior 2,.5
            exterior 1,.01
          #break
        #case (13)
          function { pow(min(min(abs(x) + y, (x*2)*(x*2)),y*5),2) }
          #break
        #case (14)
          gradient x
          #break
        #case (15)
          granite
          #break
        #case (16)
          hexagon
          #break
        #case (17)
          image_pattern {
            //png "mtmandj.png"
            function 50, 50 {
              pattern {
                object {
                  Povray_Logo_Prism
                  scale <1,-1,1>*0.5
                  translate 0.5*(x+y)
                }
              }
            }
            interpolate 2
          }
          #break
        #case (18)
          leopard
          #break
        #case (19)
          marble
          #break
        #case (20)
          onion
          #break
        #case (21)
          planar
          #break
        #case (22)
          quilted
          #break
        #case (23)
          radial frequency 6
          #break
        #case (24)
          ripples
          #break
        #case (25)
          spherical
          #break
        #case (26)
          spiral1 5
          #break
        #case (27)
          spiral2 5
          #break
        #case (28)
          spotted
          #break
        #case (29)
          waves
          #break
        #case (30)
          wood
          #break
        #case (31)
          wrinkles
          #break
        #case (32)
          // facets
          #break
        #case (33)
          average
          pigment_map {
                [1, P1]
                [1, P2]
                }
          #break
        #case (34)
          object {
            object {
              union {
                cone {-y,1,y/2,0}
                difference {
                  sphere {y,1}
                  cylinder {-z,z,0.25 translate y*2}
                  cylinder {-z,z,0.25 scale <1,0.5,1> translate <0.5,1.5,0>}
                }
                sphere {0,0.125 scale <2,1,1> translate <1,2,0> }
              }
              scale <0.3,0.6,0.3> rotate 90*z translate 0.25*x
            }
            color rgb <0.85, 0.55, 0.45>
            color rgb <0.00, 0.00, 0.42>
          }
          #break
        #case (35)
          pigment_pattern {
            checker
            pigment { gradient x color_map {[0 color rgb 0.0][1 color rgb 0.5]} scale 0.34 rotate 45*z },
            pigment { marble  color_map {[0 color rgb 0.3][1 color rgb 1.0]} scale 2 turbulence 1 }
            scale 0.5
          }
          #break
        #case (36)
          slope {<0,1,0>,0,1 altitude <0,0,0>,0,1}
          #break
      #end

     #if ((Typ=33) | (Typ=34))
     #else
      #if ((Typ=4) | (Typ=7) | (Typ=16))
        color rgb <0.85, 0.55, 0.45>,
        color rgb <0.00, 0.00, 0.42>,
        #if (Typ=16)
          color rgb <0.25, 0.40, 0.35>
        #end
      #else
        #if (Typ=32)
          color rgb <0.25, 0.40, 0.35>
        #else
        color_map {
          [0.0 color rgb <0.85, 0.20, 0.10> ]
          [0.1 color rgb <0.85, 0.55, 0.45> ]
          [0.5 color rgb <0.25, 0.40, 0.35> ]
          [0.9 color rgb <0.00, 0.00, 0.42> ]
          [1.0 color rgb <0.00, 0.00, 0.22> ]
        }
        #end
      #end

      scale 0.7

      #if (Typ=17)
        translate 1 scale 0.9  // image_pattern
      #end
      #if (Typ=10)
        translate -<0.4, 0.4, 0.4> scale 1.3 // density_file
      #end
      #if (Typ=4)
        scale 0.19  // brick
      #end
      #if ((Typ=6) | (Typ=7) | (Typ=16) |
           (Typ=18) | (Typ=21) | (Typ=24) | (Typ=29))
        scale 0.3   // cells / checker / hexagon / leopard / planar / ripples / waves
      #end
      #if (Typ=0)
        rotate 90*y // agate
      #end
      #if ((Typ=16) | (Typ=23) | (Typ=24))
        rotate 90*x // hexagon / radial / ripples
      #end
     #end
    } // close pigment
    #if (Typ=32)  // facets
      normal {
        facets 0.5
        coords 1.5
        scale 3
      }
    #end
    finish {
      #if (Typ=32)
        specular 0.4  // facets
      #else
        diffuse 0.0
        ambient 1.0
      #end
    }
  } // close texture
} // close plane
