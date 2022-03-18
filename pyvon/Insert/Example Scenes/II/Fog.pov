// Insert menu illustration scene
// Created June-August 2001 by Christoph Hormann

// ----- fog -----

// -w120 -h48 +a0.1 +am2 -j +r3

#version 3.5;

#include "colors.inc"

/*
#declare Typ=1;     // ambient fog
#declare Typ=2;     // ground fog
*/

global_settings {
  assumed_gamma 1
  max_trace_level 5
}

light_source {
  <1.5, 2.5, 2.5>*10000
  color rgb 1.0
}


#if (Typ=1)

fog {
  color rgb <0.9, 0.8, 0.3>
  distance 360
}

#else

fog {
  fog_type 2
  fog_alt 0.5
  fog_offset 0
  color rgb <0.9, 0.8, 0.3>
  distance 60
  up z
}

#end

camera {
  location    <2.0, 40, 3.2>
  direction   y
  sky         z
  up          z
  right       (120/48)*x
  look_at     <0, 0, 2.7>
  angle       7
}

sky_sphere {
  pigment {
    gradient z
    color_map {
      [0.0 rgb <0.6,0.7,1.0>]
      [1.0 rgb <0.2,0.2,0.8>]
    }
  }
}

// ----------------------------------------
plane
{
  z, 0
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


// =============================================


#declare Cnt=0;

#while (Cnt < 80)

  object {
    cylinder { -z, z*2.2, 0.1  }
    texture {
      pigment { color rgb <1.0,0.2,0.2> }
      finish {
        diffuse 0.7
        specular 0.4
        roughness 0.01
      }
    }
    no_shadow
    translate y*(-Cnt*12)-3*x
  }
  #declare Cnt=Cnt+1;
#end

