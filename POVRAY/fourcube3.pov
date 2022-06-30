//EXAMPLE OF PROJECTION OF FOUR-CUBE

//Files with predefined colors and textures
#include "colors.inc"
#include "glass.inc"
#include "golds.inc"
#include "metals.inc"
#include "stones.inc"
#include "woods.inc"

//Place the camera
camera {
   sky <0,0,1>          //Don't change this
   direction <-1,0,0>   //Don't change this
   right <-4/3,0,0>     //Don't change this
   location  <10,20,30>  //Camera location
   look_at   <0,0,0>    //Where camera is pointing
   angle 10       //Angle of the view
}

//Ambient light to "brighten up" darker pictures
global_settings { ambient_light White }
global_settings { max_trace_level 10 }


//Place a light
light_source {
   <10,20,30>
   color White*2
}

//Set a background color
background { color White }

//List the vertices of the hypercube
#declare p =
  array[16][4]
  {
   {-1,-1,-1,-1}, //0
   {-1,-1,-1,1},  //1
   {-1,-1,1,-1},  //2
   {-1,-1,1,1},   //3
   {-1,1,-1,-1},  //4
   {-1,1,-1,1},   //5
   {-1,1,1,-1},   //6
   {-1,1,1,1},    //7
   {1,-1,-1,-1},  //8
   {1,-1,-1,1},   //9
   {1,-1,1,-1},   //10
   {1,-1,1,1},    //11
   {1,1,-1,-1},   //12
   {1,1,-1,1},    //13
   {1,1,1,-1},    //14
   {1,1,1,1}      //15
  };

//Create the array to collect the projected points
#declare q = array[16];

//Specify the direction of projection
#declare a=1;
#declare b=1;
#declare c=1;
#declare d=1;

//Project the points
#declare l=sqrt(a*a+b*b+c*c+d*d);
#declare i=0;
#while(i<16)
  #declare q[i]=
    < (d*p[i][0]+c*p[i][1]-b*p[i][2]-a*p[i][3])/l,
      (-c*p[i][0]+d*p[i][1]+a*p[i][2]-b*p[i][3])/l,
      (b*p[i][0]-a*p[i][1]+d*p[i][2]-c*p[i][3])/l >;
  #declare i=i+1;
#end

//Create the polygons of the resulting polytope
#declare f = array[24];

#declare f[0] = polygon { 5, q[0], q[1], q[3], q[2], q[0] };
#declare f[1] = polygon { 5, q[4], q[5], q[7], q[6], q[4] };
#declare f[2] = polygon { 5, q[8], q[9], q[11], q[10], q[8] };
#declare f[3] = polygon { 5, q[12], q[13], q[15], q[14], q[12] };

#declare f[4] = polygon { 5, q[0], q[1], q[5], q[4], q[0] };
#declare f[5] = polygon { 5, q[2], q[3], q[7], q[6], q[2] };
#declare f[6] = polygon { 5, q[8], q[9], q[13], q[12], q[8] };
#declare f[7] = polygon { 5, q[10], q[11], q[15], q[14], q[10] };

#declare f[8] = polygon { 5, q[0], q[2], q[6], q[4], q[0] };
#declare f[9] = polygon { 5, q[1], q[3], q[7], q[5], q[1] };
#declare f[10] = polygon { 5, q[8], q[10], q[14], q[12], q[8] };
#declare f[11] = polygon { 5, q[9], q[11], q[15], q[13], q[9] };

#declare f[12] = polygon { 5, q[0], q[1], q[9], q[8], q[0] };
#declare f[13] = polygon { 5, q[2], q[3], q[11], q[10], q[2] };
#declare f[14] = polygon { 5, q[4], q[5], q[13], q[12], q[4] };
#declare f[15] = polygon { 5, q[6], q[7], q[15], q[14], q[6] };

#declare f[16] = polygon { 5. q[0], q[2], q[10], q[8], q[0] };
#declare f[17] = polygon { 5, q[1], q[3], q[11], q[9], q[1] };
#declare f[18] = polygon { 5, q[4], q[6], q[14], q[12], q[4] };
#declare f[19] = polygon { 5, q[5], q[7], q[15], q[13], q[5] };

#declare f[20] = polygon { 5, q[0], q[4], q[12], q[8], q[0] };
#declare f[21] = polygon { 5, q[1], q[5], q[13], q[9], q[1] };
#declare f[22] = polygon { 5, q[2], q[6], q[14], q[10], q[2] };
#declare f[23] = polygon { 5, q[3], q[7], q[15], q[11], q[3] };

#declare hypercube = object { union {
   object{f[0]}
   object{f[1]}
   object{f[2]}
   object{f[3]}
   object{f[4]}
   object{f[5]}
   object{f[6]}
   object{f[7]}
   object{f[8]}
   object{f[9]}
   object{f[10]}
   object{f[11]}
   object{f[12]}
   object{f[13]}
   object{f[14]}
   object{f[15]}
   object{f[16]}
   object{f[17]}
   object{f[18]}
   object{f[19]}
   object{f[20]}
   object{f[21]}
   object{f[22]}
   object{f[23]}
   }
   //texture { pigment { color rgbf <1,0,0,.5> }}
   texture {T_Gold_1A}
   };

//display it
hypercube





