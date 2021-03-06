#include "colors.inc"

camera {
    location <-15, 0, 0>
    direction <0, 0, 5>
    look_at <0, 0, 0>
    aperture 0.4
    blur_samples 100
    focal_point <-1, 0, 0>
}

light_source { <-9, 7, -6> color White }
light_source { <9, -7, 6> color White }

#declare DiceColor = color red 1 green .95 blue .65;
#declare DotColor = color red .1 green .1 blue .1;


#declare DiceBody = intersection {
    box { <-1, -1, -1>, <1, 1, 1> scale 0.5 }
    superellipsoid { <0.7, 0.7>  scale 0.63 }
}

#declare Middle = sphere { <0, 0.6, 0>, 0.13 }

#declare Corners1 = union {
    sphere { <-.25, .6, -.25>, 0.13 }
    sphere { <.25, .6, .25>, 0.13 }
}

#declare Corners2 = union {
    sphere { <-.25, .6, .25>, 0.13 }
    sphere { <.25, .6, -.25>, 0.13 }
}

#declare Middles = union {
    sphere { <-.25, .6, 0>, 0.13 }
    sphere { <.25, .6, 0>, 0.13 }
}

#declare One = Middle

#declare Two = Corners1

#declare Three = union {
    object { Middle }
    object { Corners1 }
}

#declare Four = union {
    object { Corners1 }
    object { Corners2 }
}

#declare Five = union {
    object { Four }
    object { One }
}

#declare Six = union {
    object { Corners1 }
    object { Corners2 }
    object { Middles }
}

#declare DiceInterior = interior { ior 1.5 }
#declare DiceFinish = finish { phong 0.1 specular 0.5 ambient 0.4 }

#macro Dice(Color)
difference {
    object {
        DiceBody
        pigment { color Color filter 0.4 transmit 0.3 }
        interior { DiceInterior }
        finish { DiceFinish }
    }
    union {
        object { One rotate -90*z }
        object { Two }
        object { Three rotate -90*x }
        object { Four rotate 90*x }
        object { Five rotate 180*x }
        object { Six rotate 90*z }
        pigment { White }
        finish { ambient 0.5 roughness 0.5 }

    }
    bounded_by { box { <-0.52, -0.52, -0.52>, <0.52, 0.52, 0.52> } }
}
#end

object { Dice(color rgb <0.7, 0, 0>)  rotate <195, -30, 10> }
object { Dice(color rgb <0, 0, 0.7>)  rotate <30, 40, 50> translate <3.5, 1, 1> }
object { Dice(color rgb <0, 0.5, 0>)  rotate <-40, 20, -120> translate <4.5, 1, -1> }
object { Dice(color rgb <0.5, 0.5, 0>)  rotate <-10, 290, -30> translate <5.5 ,-0.8, 0> }

