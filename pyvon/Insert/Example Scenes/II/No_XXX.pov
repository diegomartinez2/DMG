// Insert menu illustration scene
// Created June-August 2001 by Christoph Hormann

// ----- no_image / no_reflection / no_shadow -----

// -w120 -h48 +a0.1 +am2 -j +r3

#version 3.5;

#include "colors.inc"

/*
#declare Typ=1;     // no_image
#declare Typ=2;     // no_reflection
#declare Typ=3;     // no_shadow
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
  location    <10, 20, 7.5>
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
    pigment { color rgb 1 }
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
    pigment { color MediumAquamarine }
    finish { specular 0.4 roughness 0.02 }
  }

#declare Tex2=
  texture {
    pigment { color OrangeRed }
    finish { specular 0.4 roughness 0.02 }
  }

#declare Tex3=
  texture {
    pigment { color DarkOrchid }
    finish { specular 0.4 roughness 0.02 }
  }


// =============================================

union {
  sphere { <0.5, 0, 0>, 1.2 texture { Tex1 } }

  box { -0.9, 0.9 texture { Tex2 }
    translate <-3.8, 0.6, 0>
  }

  cylinder { -z, z, 0.8 texture { Tex3 }
    translate <4, -0.6, 0>
  }
  #if (Typ=1)
    no_image
  #end
  #if (Typ=2)
    no_reflection
  #end
  #if (Typ=3)
    no_shadow
  #end
}

// =============================================


