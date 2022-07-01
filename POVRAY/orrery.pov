#include "colors.inc"    // The include files contain
#include "shapes.inc"    // pre-defined scene elements
#include "textures.inc" 

/* ****************************************************************************
* Orrery in POV-Ray
* Created By: Jason Vroustouris
* http://www.jasonv.com/
******************************************************************************/
// 1) Global Settings

	global_settings {
		assumed_gamma 1
		radiosity {
			pretrace_start 0.08
			pretrace_end 0.02
			count 50
			error_bound 0.5
			recursion_limit 1
		}
	}
// 2) Lighting
	light_source { <7, 5, 0> color White area_light <0,0,1>, <1,0,0>, 10, 10 } 
	light_source { <0, 50, 0> color White } 
	//light_source { <4.5, 0, 0> color Red } 

// 3) Camera Location
	camera {location  <0, 6, -12> look_at   <1, 1.5,  0>} // Grandstand Cam
	//camera {location  <13, 10, -10> look_at   <13, 2.5,  0>} // Saturn Cam
	//camera {location  <0, 10, -10> look_at   <0, 0,0>} // Origin Cam
	//camera {location  <0, 10, -2> look_at   <0, 0,0>} // Down Origin Cam


// 4) Plane
	//plane { <1, 0, 0>,-18 pigment {checker color White color White} }
	//plane { <0, 1, 0>,-20 pigment {checker color White color White} }
	//plane { <0, 0, 1>,-18 pigment {checker color White color White} }

// 5) Planet Texture
#declare Planet_T = texture{Ruby_Glass}
//#declare Planet_T = texture{DMFWood6}
#declare Our_Metal = texture{Polished_Chrome}
//#declare Our_Metal = texture{Polished_Chrome}

/* ****************************************************************************
* Planet and Arm
**************************************************************************** */
#macro planet_and_arm(x1,y1,xlength,ylength,planet_size)
union {
	cylinder {<0,0,0>, <xlength,0,0>, .1     texture {Our_Metal}}
	cylinder {<0,-.05,0>, <0,.05,0>, .2      texture {Our_Metal}}	
	cylinder {<xlength,0,0>, <xlength,ylength-(planet_size/2),0>, .1
		texture {Our_Metal}
	}
	sphere {<xlength,0,0>, .1                texture {Our_Metal}}
	sphere {<xlength,ylength,0>, planet_size texture {Planet_T}}
	translate <x1,y1,0>
}
#end

#macro moon_and_arm(x1,y1,xlength,ylength,rot,planet_size)
union {
	cylinder {
		<0,0,0>, <xlength,0,0>, .05
		texture {Our_Metal}
	}
	cylinder {
		<0,-.05,0>, <0,.05,0>, .2
		texture {Our_Metal}
	}	
	cylinder {
		<xlength,0,0>, <xlength,ylength,0>, .05
		texture {Our_Metal}
	}
	sphere {
		<xlength,0,0>, .05
		texture {Our_Metal}
	}
	sphere {
		<xlength,ylength,0>, planet_size
		texture {Planet_T}
	}
   	rotate <0,rot*4,0>
	translate <x1,y1,0>
}
#end

#macro planet_system(x1,y1,xlength,ylength,rot,planet_size)
union {
	planet_and_arm(x1,y1,xlength,ylength,planet_size)
   	rotate <0,rot*4,0>
}
#end

#macro earth_system(x1,y1,xlength,ylength,rot,planet_size)
union {
	planet_and_arm(x1,y1,xlength,ylength,planet_size)
	moon_and_arm(xlength,y1+ylength-.5,.5,.5,0,.1)
   	rotate <0,rot*4,0>
}
#end

#macro jupitor_system(x1,y1,xlength,ylength,rot,planet_size)
union {
	planet_and_arm(x1,y1,xlength,ylength,planet_size)
	moon_and_arm(xlength,y1+ylength-0.5,0.5,0.5,00,.1)
	moon_and_arm(xlength,y1+ylength-1.0,1.0,1.0,20,.1)
	moon_and_arm(xlength,y1+ylength-1.5,1.5,1.5,40,.1)
	moon_and_arm(xlength,y1+ylength-2.0,2.0,2.0,80,.1)
   	rotate <0,rot*4,0>
}
#end

#macro saturn_system(x1,y1,xlength,ylength,rot,planet_size)
union {
	planet_and_arm(x1,y1,xlength,ylength,planet_size)
	moon_and_arm(xlength,y1+ylength-0.50,.25+0.50,0.50,0,.1)
	moon_and_arm(xlength,y1+ylength-0.75,.25+0.75,0.75,20,.1)
	moon_and_arm(xlength,y1+ylength-1.00,.25+1.00,1.00,40,.1)
	moon_and_arm(xlength,y1+ylength-1.25,.25+1.25,1.25,60,.1)
	moon_and_arm(xlength,y1+ylength-1.50,.25+1.50,1.50,80,.1)
	moon_and_arm(xlength,y1+ylength-1.75,.25+1.75,1.75,100,.1)
	moon_and_arm(xlength,y1+ylength-2.00,.25+2.00,2.00,120,.1)
	cylinder {
		<xlength-x1,ylength+y1-.01,0>,
		<xlength-x1,ylength+y1+.01,0>, .5
		texture {Planet_T}
	}
   	rotate <0,rot*4,0>
}
#end

/* ****************************************************************************
* Orrery
**************************************************************************** */
#macro orrery()
/*cylinder { <0,0,0>, <0,1.0,0>, .2*5 texture {Our_Metal} }
cylinder { <0,0,0>, <0,0.9,0>, .3*5 texture {Our_Metal} }
cylinder { <0,0,0>, <0,0.8,0>, .4*5 texture {Our_Metal} }
cylinder { <0,0,0>, <0,0.7,0>, .5*5 texture {Our_Metal} }
cylinder { <0,0,0>, <0,0.6,0>, .6*5 texture {Our_Metal} }
cylinder { <0,0,0>, <0,0.5,0>, .7*5 texture {Our_Metal} }
cylinder { <0,0,0>, <0,0.4,0>, .8*5 texture {Our_Metal} }
cylinder { <0,0,0>, <0,0.3,0>, .9*5 texture {Our_Metal} }
cylinder { <0,0,0>, <0,0.2,0>,1.0*5 texture {Our_Metal} }
*/
cylinder { <0,0,0>, <0,4.0,0>, .1  texture {Our_Metal} }
cylinder { <0,-4,0>, <0,8.0,0>, 100  texture {Our_Metal} }
lathe {
  quadratic_spline
  6,
  <0, 0>, <3, 0>, <3,.3>,<2,.5>,<0, 1>, <0, 0>
  texture {Our_Metal}
 }
lathe {
  quadratic_spline
  6,
  <0, 0>, <3, 0>, <3,.3>,<2,.5>,<0, 1>, <0, 0>
  texture {Our_Metal}
 }


/*sphere   { <0,4.5,0>,          .3 
	pigment {
	color White
	}
	finish { 
	diffuse 1 ambient 1 }  
}*/


sphere   { <0,4.5,0>,          .5 
texture{NBglass}
}

planet_system (0,  4.0,  2,  0.5,  0, .3)// Mercury
planet_system (0,  3.5,  3,  1,    10, .3)// Venus
earth_system  (0,  3.0,  4,  1.5,  20, .3)// Earth
planet_system (0,  2.5,  5,  2.0,  30, .3)// Mars
jupitor_system(0,  2.0,  8,  2.5,  60, .3)// Jupitor
saturn_system (0,  1.5,  14, 3.0,  80, .3)// Saturn
#end

/* ****************
* Wood Case
**************************************************************************** */
#macro glasscase()
lathe {
  linear_spline
  9,
  <0,0>,<17,0>,<17,-.2>,<17.2,-.2>,<17.2,0>,<17.5,0>,<17.5,-.3>,<0,-.3>,<0,0>
  texture {DMFWood6}
}
/*lathe {
  linear_spline
  9,
  <0,5.5>,<17.5,5.5>,<17.5,5.4>,<17.4,5.4>,<17.4,0>,<17.3,0>,<17.3,5.4>,<0,5.4>,<0,5.5>
  texture {NBglass}
}*/


#end
/* ****************************************************************************
* Main
**************************************************************************** */
orrery()
glasscase()

