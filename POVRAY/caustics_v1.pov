/////////////////////////////////////////////
//
//     ~~ [ Caustics ] ~~
//        version 1 out of 1
//
//  by Michael Scharrer
//  https://mscharrer.net
//
// CC-BY-4.0 license
// https://creativecommons.org/licenses/by/4.0/
//
/////////////////////////////////////////////


camera{
	right x*image_width/image_height
 location <2,7,-10>
 look_at <0,0,0>
}

global_settings {
 max_trace_level 7
 photons {
  count 200000
  media 500
  autostop 1000
  jitter .2
 }
}

#declare r = seed(69);
#declare i=0;
#while (i<15)
 light_source{
  <30-60*rand(r),10,30-60*rand(r)>
  color rgb <0.5*rand(r),0.5*rand(r),0.5*rand(r)>
 }
 #declare i = i+1;
#end


plane { 
 <0, 1, 0>, -2
 hollow
 pigment{
  checker
  pigment {
   granite
   turbulence 3
   color_map {
    [0.0 color <0.05,0.05,0.05>]
    [0.2 color <0.3,0.2,0>]
    [0.4 color <0.05,0.05,0.05>]
    [0.6 color <0.2,0.2,0.2>]
    [1.0 color <0.05,0.05,0.05>]
    }
   }
  pigment {
   marble
   turbulence 2
   color_map {
    [0.0 color <0.95,0.95,0.9>]
    [0.8 color <0.6,0.6,0.6>]
    [1.0 color <0.3,0.2,0.2>]
    }
   }
  turbulence 0.01
 }
 finish{ reflection 0.4 }
 /*photons {
  target
  refraction on
  reflection on
 }*/
}

plane { 
 <0, 1, 0>, 20
 hollow
 pigment{
  checker
  color <0.6,0.6,0.6>
  color <0.4,0.4,0.4>
  turbulence 0.02
 }
}

sphere{
 <0,0,12> 2
 hollow
 pigment { color <1,1,1,1> }
 interior{
  media{
   emission 0.15
  }
 }
 photons {
  target
  refraction on
  reflection on
 }
 photons {
  target
  refraction on
  reflection on
 }
}

sphere{
 <-4,0,5> 2
 pigment { color <0.95,1,0.95,0.9> }
 finish {
  ior 1.5
  reflection 0.1
 }
 photons {
  target
  refraction on
  reflection on
 }
}

sphere{
 <4,0,5> 2
 pigment { color <0.3,0.3,0.4> }
 finish {
  reflection 0.8
 }
 photons {
  target
  refraction on
  reflection on
 }
}


sphere{
 <5,0,-1> 2
 pigment { color <0.95,1,0.95,0.9> }
 finish {
  ior 1.4
  reflection 0.1
 }
 normal{
  bumps 1/50
  scale 1/5
 }
 photons {
  target
  refraction on
  reflection on
 }
}

sphere{
 <-5,0,-1> 2
 pigment { color <0.3,0.3,0.4> }
 finish {
  reflection 0.8
 }
 normal{
  bumps 1/15
  scale 1/5.5
 }
 photons {
  target
  refraction on
  reflection on
 }
}
