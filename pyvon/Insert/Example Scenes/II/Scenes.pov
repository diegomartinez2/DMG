// Insert menu illustration scene
// Created June-August 2001 by Christoph Hormann

// ----- scene templates -----

// -w120 -h90 +a0.2 +am2 -j +r3

// Typ4:
// -w90 -h90 +a0.2 +am2 -j +r3

/*
#declare Typ=1;     // basic
#declare Typ=2;     // radiosity
#declare Typ=3;     // checkered floor
#declare Typ=4;     // orthographic
#declare Typ=5;     // TTF
#declare Typ=6;     // photons
#declare Typ=7;     // image map
*/

// =============================================

#if (Typ=1)

#version 3.5;

#include "colors.inc"

global_settings {
  assumed_gamma 1.0
}

// ----------------------------------------

camera {
  location  <0.0, 0.5, -4.0>
  direction 1.5*z
  right     4/3*x
  look_at   <0.0, 0.0,  0.0>
}

sky_sphere {
  pigment {
    gradient y
    color_map {
      [0.0 rgb <0.6,0.7,1.0>]
      [0.7 rgb <0.0,0.1,0.8>]
    }
  }
}

light_source {
  <0, 0, 0>            // light's position (translated below)
  color rgb <1, 1, 1>  // light's color
  translate <-30, 30, -30>
}

// ----------------------------------------

plane {
  y, -1
  pigment { color rgb <0.7,0.5,0.3> }
}

sphere {
  0.0, 1
  texture {
    pigment {
      radial
      frequency 8
      color_map {
        [0.00 color rgb <1.0,0.4,0.2> ]
        [0.33 color rgb <0.2,0.4,1.0> ]
        [0.66 color rgb <0.4,1.0,0.2> ]
        [1.00 color rgb <1.0,0.4,0.2> ]
      }
    }
    finish{
      specular 0.6
    }
  }
}

#end

// =============================================

#if (Typ=2)

#version 3.5;

#declare Radiosity=on;    // use +QR on commandline

global_settings {
  assumed_gamma 1.0
  //max_trace_level 25
  #if (Radiosity)
    radiosity {
      pretrace_start 0.08           // start pretrace at this size
      pretrace_end   0.04           // end pretrace at this size
      count 35                      // higher -> higher quality (1..1600) [35]
      nearest_count 5               // higher -> higher quality (1..10) [5]
      error_bound 1.8               // higher -> smoother, less accurate [1.8]
      recursion_limit 3             // how much interreflections are calculated (1..5+) [3]
      low_error_factor .5           // reduce error_bound during last pretrace step
      gray_threshold 0.0            // increase for weakening colors (0..1) [0]
      minimum_reuse 0.015           // reuse of old radiosity samples [0.015]
      brightness 1                  // brightness of radiosity effects (0..1) [1]

      adc_bailout 0.01/2
      //normal on                   // take surface normals into account [off]
      //media on                    // take media into account [off]
      //save_file "file_name"       // save radiosity data
      //load_file "file_name"       // load saved radiosity data
      //always_sample off           // turn sampling in final trace off [on]
      //max_sample 1.0              // maximum brightness of samples
    }
  #end
}

#default {
  texture {
    pigment {rgb 1}
    #if (Radiosity)
      finish {
        ambient 0.0
        diffuse 0.6
        specular 0.3
      }
    #else
      finish {
        ambient 0.1
        diffuse 0.6
        specular 0.3
      }
    #end
  }
}

// ----------------------------------------

camera {
  right x*image_width/image_height
  location  <0,1.5,-4>
  look_at   <0,1,0>
}

light_source {
  <500,500,-500>       // light's position
  color rgb <1, 1, 1>  // light's color
}

sky_sphere {
  pigment {
    gradient y
    color_map {
      [0.0 rgb <0.6,0.7,1.0>]
      [0.7 rgb <0.0,0.1,0.8>]
    }
  }
}

// ----------------------------------------

plane {
  y, 0
  texture {
    pigment {
      checker
      color rgb <1.0, 0.8, 0.6>
      color rgb <1.0, 0.0, 0.0>
      scale 0.5
    }
  }
}

sphere {
  <0,1,0>, 1
}

#end

// =============================================

#if (Typ=3)

#version 3.5;

#include "colors.inc"

global_settings {
  assumed_gamma 1.0
  max_trace_level 5
}

// ----------------------------------------

camera {
  location  <0.0, 0.5, -4.0>
  direction 1.5*z
  right     4/3*x
  look_at   <0.0, 0.0,  0.0>
}

sky_sphere {
  pigment {
    gradient y
    color_map {
      [0.0 rgb <0.6,0.7,1.0>]
      [0.7 rgb <0.0,0.1,0.8>]
    }
  }
}

light_source {
  <0, 0, 0>            // light's position (translated below)
  color rgb <1, 1, 1>  // light's color
  translate <-30, 30, -30>
}

// ----------------------------------------

plane {               // checkered floor
  y, -1
  texture
  {
    pigment {
      checker
      color rgb 1
      color blue 1
      scale 0.5
    }
    finish{
      diffuse 0.8
      ambient 0.1
    }
  }
}

sphere {              // reflective sphere
  0.0, 1
  texture {
    pigment {
      color rgb <0.8,0.8,1.0>
    }
    finish{
      diffuse 0.3
      ambient 0.0
      specular 0.6
      reflection {
        0.8
        metallic
      }
      conserve_energy
    }
  }
}

#end

// =============================================

#if (Typ=4)

#version 3.5;

global_settings {
  assumed_gamma 1.0
  //hf_gray_16 on
}

// ----------------------------------------

camera {
  orthographic
  location <0,0,1>     // position & direction of view
  look_at  <0,0,0>
  right 1*x            // horizontal size of view  \___ to be rendered at square size
  up 1*y               // vertical size of view    /
}

// ----------------------------------------

box {                  // this box fits exactly in view
  <-0.5, -0.5, 0>, <0.5, 0.5, 0>
  texture {
    pigment {
      agate
      color_map {
        [0.0 color rgb 0.0 ]
        [1.0 color rgb 1.0 ]
      }
    }
    finish {
      ambient 1.0
      diffuse 0.0
    }
  }
}

#end

// =============================================

#if (Typ=5)

#version 3.5;

global_settings {
  assumed_gamma 1.0
}

// ----------------------------------------

camera {
  location  <0.0, 2.0, -6.0>
  direction 1.5*z
  right     4/3*x
  look_at   <0.0, 0.0,  0.0>
}

sky_sphere {
  pigment {
    gradient y
    color_map {
      [0.0 color blue 0.6]
      [1.0 color rgb 1]
    }
  }
}

light_source {
  <0, 0, 0>            // light's position (translated below)
  color rgb <1, 1, 1>  // light's color
  translate <-30, 30, -30>
}

// ----------------------------------------

#declare Text_Tex = texture {
  pigment { granite scale 0.5 }
  finish { specular 0.7 }
}

text {
  ttf "crystal.ttf", "Hello",
  2, // depth
  0  // spacing
  texture { Text_Tex }
  rotate <0, -20, 0>
  translate <-1, 0, -3>
}

text {
  ttf "crystal.ttf", "Virtual World!",
  1, // depth
  0  // spacing
  scale <1, 2, 1> // stretch it taller
  texture { Text_Tex }
  rotate <0, -30, 0>
  translate <-3, 0, 3>
}


plane { y, 0 pigment { color rgb <0.7,0.5,0.3> } }

#end

// =============================================

#if (Typ=6)

#version 3.5;

#declare Photons=on;

global_settings {
  assumed_gamma 1.0
  max_trace_level 5
  #if (Photons)          // global photon block
    photons {
      spacing 0.02                 // specify the density of photons
      //count 100000               // alternatively use a total number of photons

      //gather min, max            // amount of photons gathered during render [20, 100]
      //media max_steps [,factor]  // media photons
      //jitter 1.0                 // jitter phor photon rays
      //max_trace_level 5          // optional separate max_trace_level
      //adc_bailout 1/255          // see global adc_bailout
      //save_file "filename"       // save photons to file
      //load_file "filename"       // load photons from file
      //autostop 0                 // photon autostop option
      //radius 10                  // manually specified search radius
      // (---Adaptive Search Radius---)
      //steps 1
      //expand_thresholds 0.2, 40
    }

  #end
}

// ----------------------------------------

camera {
  right x*image_width/image_height
  location  <0,1.6,-5>
  look_at   <0,0.75,0>
}

light_source {
  <500,500,150>       // light's position
  color rgb 1.3       // light's color
  photons {           // photon block for a light source
    refraction on
    reflection on
  }
}

sky_sphere {
  pigment {
    gradient y
    color_map {
      [0.0 rgb <0.6,0.7,1.0>]
      [0.7 rgb <0.0,0.1,0.8>]
    }
  }
}

// ----------------------------------------

plane {
  y, 0
  texture {
    pigment { color rgb <1.0, 0.8, 0.6> }
  }
}


#declare M_Glass=    // Glass material
material {
  texture {
    pigment {rgbt 1}
    finish {
      ambient 0.0
      diffuse 0.05
      specular 0.6
      roughness 0.005
      reflection {
        0.1, 1.0
        fresnel on
      }
      conserve_energy
    }
  }
  interior {
    ior 1.5
    fade_power 1001
    fade_distance 0.9
    fade_color <0.5,0.8,0.6>
  }
}


sphere {
  <0,1,0>, 1
  translate <1.0,0,-1.3>
  material { M_Glass }

  photons {  // photon block for an object
    target 1.0
    refraction on
    reflection on
  }
}

cylinder {
  <0,0.01,0>, <0,2.5,0>, 1
  translate <-1.2,0,0.8>
  material { M_Glass }

  photons {  // photon block for an object
    target 1.0
    refraction on
    reflection on
  }
}

#end

// =============================================

#if (Typ=7)

#version 3.5;

#include "colors.inc"

global_settings {
  assumed_gamma 1.0
}

// ----------------------------------------

camera {
  location  <0.0, 0.0, -4.0>
  direction 2*z
  right     4/3*x
  look_at   <0.0, 0.0,  0.0>
}

sky_sphere {
  pigment {
    gradient y
    color_map {
      [0.0 color blue 0.6]
      [1.0 color rgb 1]
    }
  }
}

light_source {
  <0, 0, 0>            // light's position (translated below)
  color rgb <1, 1, 1>  // light's color
  translate <-30, 30, -30>
}

// ----------------------------------------

plane {
  y, -1
  texture {
    pigment { checker color rgb 1 color blue 1 scale 0.5 }
    finish { reflection 0.2 }
  }
}

plane {
  z, -1
  texture {
    pigment {
      image_map {
        png "test.png"
        interpolate 2 // smooth it
        once   // don't tile image, just one copy
        filter 0 0.8  // make 1st color mostly transparent
        filter 1 0.8  // make 2nd color mostly transparent
      }
      // transform it to unit-size (-1 to +1)
      translate -0.5*(x+y) // center on the origin
      scale 2              // make it unit-sized
    }
    finish { ambient 0.3 }
  }
}

#end

// =============================================
