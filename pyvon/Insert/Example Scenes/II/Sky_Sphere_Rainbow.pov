// Persistence of Vision Ray Tracer Scene Description File
// File: im_sky_sphere_rainbow.pov
// Vers: 3.5
// Desc: sky and rainbow previews
// Date: August 2001
// Auth: Bob Hughes
//       modified Aug. 25 by Christoph Hormann
// Note: pot o' gold was added, animatable
// [revised]
// +w120 +h40 +a0.2 +r2 +am2 +j0.5

#version 3.5;

/*
#declare Typ=1;     // sky_sphere
#declare Typ=2;     // rainbow
*/

#declare PotOGold=no;

#if (Typ=1)
  #declare Rainbow=no;
#else
  #declare Rainbow=yes;
#end

global_settings {
  assumed_gamma 1
  max_trace_level 15
}

// ----------------------------------------

camera { //fisheye
  location  <0, 0, -15>
  right     3*x
  rotate 15*y
  look_at   0//-0.25*y
  angle 89//36
}

#declare Clock=0; // clock

// ----------------------------------------

#if (PotOGold=yes)

union {
//sphere {0,1 clipped_by {plane {y,.25}}
lathe {
        cubic_spline
        13,
        <-0.1,0.1>,
        <0,0>,
        <0.2,0.1>,
        <0.5,0.3>,
        <0.7,0.7>,
        <0.65,0.9>,
        <0.7,0.95>,
        <0.75,1>,
        <0.6,0.9>,
        <0.7,0.7>
        <0.5,0.25>,
        <0.15,0.1>
        <-0.1,0.5>
       // pigment {rgb 1}
       translate -y/1.5 scale <1.5,2,1.5>
hollow

pigment {rgb .333}
normal {dents 2 scale .05}
finish {diffuse .5 specular .2 roughness .1 phong .5 phong_size 123}
}
//plane {y,.3 clipped_by {sphere {0,.99}}
height_field {
        pattern 200,200 {
                hf_gray_16
                crackle offset .125 metric 2
                color_map {
                        [.15 rgb 0]
                        [.2 rgb 1]
        } scale .25}
         clipped_by {sphere {y/.75,.9}}
pigment {rgb <1.3,1.1,.25>}
normal {crackle 2 metric 1 scale.125}
finish {ambient .125 diffuse .5 specular 2 roughness .075 phong .5 phong_size 10 brilliance 2
        reflection {.1,.6 metallic}
        }
translate -.5 scale 1.4 rotate 180*y rotate 15*z translate <.33,-.1,0>
}
height_field {
        pattern 200,200 {
                hf_gray_16
                crackle offset .125 metric 2
                color_map {
                        [.15 rgb 0]
                        [.2 rgb 1]
        } scale .25}
         clipped_by {sphere {y/.67,.75}}
pigment {rgb <1.3,1.1,.25>}
normal {crackle 2 metric 1 scale.125}
finish {ambient .125 diffuse .5 specular 2 roughness .075 phong .5 phong_size 10 brilliance 2
        reflection {.1,.6 metallic}
        }
translate -.5 scale <1.5,.75,1.5> rotate -90*y rotate 15*z translate <.5,-1,-.25>
}
 scale 10 rotate -15*x translate <33,-7.5,25>
}


#declare P1=
        pigment {spherical
                color_map {
                        [0 rgb 1][1 rgbf 1]
        }frequency -1 }
#declare P2=
        pigment {radial
                color_map {
                        [0 rgb 1][1 rgbf 1]
        } scallop_wave frequency 5 rotate 90*x}
#declare P3=
        pigment {rgbf 1}

disc {0,z,1 hollow
        pigment {onion
                pigment_map {
                        [0 P1][.5 P2][1 P3]
        }}
        finish {ambient 1 diffuse 0}
 scale 4*(.1+Clock) rotate z*60*Clock rotate <5,15,0> translate <25,-1,12.5>
  no_shadow
}

#end

//////////////////////////////////////////////////////////////////////////////////

/* sky */

sky_sphere {
  pigment {
    gradient y
    color_map {
      [0.0 rgb <0.7,0.8,1.1>]
      [0.1 rgb <0.2,0.5,1.2>]
      [0.6 rgb <0.1,0.4,0.8>]
    } turbulence .1 translate -y/7
  }
  pigment {
    bozo
    color_map {
      [0.0 rgbf <0.67,0.5,1.0,0>]
      [0.3 rgbf <0.9,0.9,0.9,.1>]
      [0.4 rgbf <1,1,1,1>]
    } turbulence 0.5 scale <.5,.125,.25>
  }
  pigment {
    gradient y
    color_map {
      [0.5 rgbf <1.0,1.0,1.0,1.0>]
      [0.5 rgb <0.8,0.6,0.2>]
      [0.9 rgb <0.2,0.6,0.3>]
      [1.0 rgb <0.6,0.9,0.1>]
    } cubic_wave turbulence .1 translate -y/7
  }
}

light_source {
  0,
  color rgb 1  // light's color
  translate <-333, 333, -333>
 // media_interaction on media_attenuation on
}


//////////////////////////////////////////////////////////////////////////////////

#if (Rainbow=on)

// use clock here if animating otherwise use static number
#declare C=1; // orientation/size
#declare CC=0.5; // color intensity change (higher is lighter)

/* bright single rainbow */

  rainbow {
    direction <.2,-.36,.8>
    angle 35
    width 10
    distance 5
    color_map {
  [0.0,0.12 color rgbt<.5,.5,.5,1> color rgbt<1.5,0,2.5,1.25>/(2-CC) ]
  [0.12,0.23 color rgbt<1.5,0,2.5,1.25>/(2-CC) color rgbt<0,1,2,.975>/(2-CC) ]
  [0.23,0.3 color rgbt<0,1,2,.975>/(2-CC) color rgbt<0,1.5,0.5,.925>/(2-CC) ]
  [0.3,0.45 color rgbt<0,1.5,0.5,.925>/(2-CC) color rgbt<1.5,2,0,.95>/(2-CC)]
  [0.45,0.6 color rgbt<1.5,2,0,.95>/(2-CC) color rgbt<2,1.5,0,.925>/(2-CC)]
  [0.6,0.75 color rgbt<2,1.5,0,.925>/(2-CC) color rgbt<2.5,1,0,.9125>/(2-CC)]
  [0.75,0.9 color rgbt<2.5,1,0,.9125>/(2-CC) color rgbt<3,0,0,1.25>/(2-CC)]
  [0.9,1.0 color rgbt<3,0,0,1.25>/(2-CC) color rgbt<.5,.5,.5,1> ]
	}
    jitter 0.02
    up <0,1.2-(2-2*C),-.8+(2-2*C)>  // 1*y is normal up (any x or z value will tilt)
    arc_angle (30+(150*C))-(180*Clock)
    falloff_angle (30+(30*C))-(60*Clock)
  }

#end

