// Persistence of Vision Ray Tracer Scene Description File
// File: images.pov
// Vers: 3.5
// Desc: wrapping file for rendering all Insert menu pictures
// Date: August 2001
// Auth: Christoph Hormann

// Modified numbering for easier customizing. August 30th, Bob Hughes
// When adding to the menu previews simply append to the existing scene files
// and use next number in series.  Or add new file and insert a new range.
// Be sure to also add entries into Insert Menu\image.ini

// aspect ratios are not standardized,
// for correct options see file 'im.ini'

#switch (Switch)
  #range (1,28)      // --- shapes ---

    #declare Typ=Switch;
    #include "Shapes.pov"

    #break
  #range (101,136)     // --- patterns ---

    #declare Typ=Switch-100;
    #include "Patterns.pov"

    #break
  #range (201,208)     // --- pattern texture attributes ---

    #declare Typ=Switch-200;
    #include "Attributes.pov"

    #break
  #range (301,304)     // --- transformations ---

    #declare Typ=Switch-300;
    #include "Transform.pov"

    #break
  #range (401,402)     // --- sky_sphere and rainbow ---

    #declare Typ=Switch-400;
    #include "Sky_Sphere_Rainbow.pov"

    #break
  #range (501,502)     // --- photons and radiosity ---

    #declare Typ=Switch-500;
    #include "Photons_Radiosity.pov"

    #break
  #range (601,602)     // --- fog ---

    #declare Typ=Switch-600;
    #include "Fog.pov"

    #break
  #range (701,702)     // --- material ---

    #declare Typ=Switch-700;
    #include "Material.pov"

    #break
  #range (801,803)     // --- no_image / no_reflection / no_shadow ---

    #declare Typ=Switch-800;
    #include "No_XXX.pov"

    #break
  #range (901,908)     // --- CSG / material ---

    #declare Typ=Switch-900;
    #include "CSG_Material.pov"

    #break
  #range (1001,1009)     // --- lights ---

    #declare Typ=Switch-1000;
    #include "Lights.pov"

    #break
  #range (1101,1106)     // --- include ---

    #declare Typ=Switch-1100;
    #include "Include.pov"

    #break
  #range (1201,1207)     // --- scenes ---

    #declare Typ=Switch-1200;
    #include "Scenes.pov"

    #break
  #range (1301,1311)     // --- camera ---

    #declare Typ=Switch-1300;
    #include "Camera.pov"

#end
