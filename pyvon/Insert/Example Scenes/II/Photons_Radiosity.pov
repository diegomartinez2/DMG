// Persistence of Vision Ray Tracer Scene Description File
// File: im_photons_radiosity.pov
// Vers: 3.5
// Desc: CSG operations: cutaway_textures, difference, intersection,
// inverse, merge, union
// Also Material: interior, media
// Date: August 2001
// Auth: Bob Hughes
//       modified Aug. 25 by Christoph Hormann
// [revised]
// +w120 +h60 +a0.2 +r2 +am2 +j0.5
// +w120 +h60 +a0.2 +r2 +am2 +j0.5

#version 3.5;

/*
#declare Typ=1;  // 1=radiosity
#declare Typ=2;  // 2=photons
*/

#declare MenuPreview=Typ;

global_settings
 {
  #if (MenuPreview=1)
        ambient_light 1
        radiosity {
         // pretrace_start 1
         // pretrace_end 1
          error_bound 0.6
          count 50
         // brightness 1
         // gray_threshold 0
         // media on
         // normal on
          recursion_limit 2
        }
  #else
        max_trace_level 15
        photons {
                spacing 0.03
                // save_file "menupreview.ph"
                // load_file "menupreview.ph"
        }
  #end
}

// ----------------------------------------

#declare Count=1;

// radiosity object
#declare IterativeShape=
union
{
#while (Count<4)

difference
{
        box {0,2}
        sphere {1,1.5 inverse}
        sphere {1,1.25}
 #switch (Count)
 #case (1)
        pigment {color red 1.25}
 #break
 #case (2)
        pigment {color rgb 1.25}
 #break
 #case (3)
        pigment {color green 1.25}
 #break
 #end
        finish {
                diffuse 1.3333-(Count/3)

                ambient 0
                }
 translate -1 // center on world origin
 scale Count translate Count*y rotate Count*72*y
}

#declare Count=Count+1;

#end
}


#declare Count=0;

// photon object
#declare Crystal=

intersection
{

 #while (Count<9)

        box
        {
        -1,1
  rotate <45,0,45> rotate Count*40*y
        }
 #declare Count=Count+1;

 #end
}


//////////////////////////////////////////////////////////////////////////////////

#switch (MenuPreview)

#case (1) /* radiosity */

sky_sphere {
  pigment {
    color rgb <0.7,0.8,1.0>*1.15
  }
}

camera {
  location  <0, 7.5, -15>
  right     x*2 // 120X60, x*2
  rotate 12*y
  look_at   <-1.25,0.75,0>
  angle 33
}

union
{

plane // floor
{
  y, -1.001
  texture
  {
    pigment { color rgb 1 }
    finish {
      ambient 0.0
      diffuse 1
      specular 0.4
      roughness 0.01

    }
  }
}

object
{
        IterativeShape
}

}

#break

//////////////////////////////////////////////////////////////////////////////////

#case (2) /* photons */

sky_sphere {
  pigment {
    gradient y
    color_map {
      [0.0 rgb <0.6,0.8,1.0>]
      [0.3 rgb <0.3,0.6,0.9>]
      [1.0 rgb <0.1,0.4,0.8>]
    }
  }
}

camera {
  location  <0, 7.5, -15>
  right     x*2 // 120X60, x*2
  rotate 6*y
  look_at   <0.25,0,0>  
  angle 30  
}

union
{

light_source {
  0,
  color rgb 1
  translate <-200, 150, -200>
}


sphere // floor
{
        0,1
        clipped_by {plane{y,0}}
        pigment {rgb 0.9}
        finish {diffuse 0.9}
        scale <10,2,10>
        translate y
}

object  // transparent refractive crystal
{
        Crystal
        pigment {rgbft <1,1,1,0.9,0.1>}
        finish {specular 0.3 roughness 0.009
                reflection {0.1}
        }
        interior {ior 1.33 dispersion 1.2 dispersion_samples 12}
  translate <1.75,0.5,0>
        photons {target 0.9 refraction on reflection on}
}

object  // opaque reflective crystal
{
        Crystal
        pigment {rgb 0.1}
        finish {diffuse 0.2 specular 0.6 roughness 0.01
                reflection {0.9}
        }
  translate <-1.75,0.5,0>
        photons {target 0.9 refraction off reflection on}
}

}

#break

#end
