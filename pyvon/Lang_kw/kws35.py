color_tags=['keyword','keywordbold','pattern','num','slcomment','mlcomment','str','squig','math','dirc','highl','pattern','brace']

ini_list=    [["Location : ","",11],	# 0
              ["Help file : ","",12],	# 1
              ["Help viewer : ","",14],	# 2
              ["Parameters : ","",13],	# 3
              ["Directory : ","",12],	# 4
              ["Font : ","",7],		# 5
              ["Font size : ","",12],	# 6
              ["Auto save : ","",12],	# 7
              ["Enable coloring : ","",18],# 8
              ["Language : ","",11],	# 9
              ["keywords : ","",11],	# 10
	      ["bold keywords : ","",16], # 11
	      ["numbers & math : ","",17],#12
	      ["comments : ","",11],	# 13
	      ["strings : ","",10],	# 14
	      ["highlighting : ","",15],#15
	      ["patterns : ","",11],	#16
	      ["directives : ","",13],	#17
	      ["braces : ","",9],	#18
	      ["editor background : ","",20], #19
	      ["identifiers : ","",14], #20
	      ["tab length : ","",13], #21
	      ["show help bubbles : ","",20], # 22
	      ["auto expand keywords : ","",23]] # 23
 



kws=['aa_level','aa_threshold','abs','absorption','accuracy','acos','acosh','adaptive','adc_bailout','agate',
     'agate_turb','all','all_intersections','alpha','altitude','always_sample','ambient','ambient_light','angle',
     'aperture','append','arc_angle','area_light','array','asc','ascii','asin','asinh','assumed_gamma','atan','atan2',
     'atanh','autostop','average','b_spline','background','bezier_spline','bicubic_patch','black_hole','blob','blue',
     'blur_samples','bounded_by','box','boxed','bozo','break','brick','brick_size','brightness','brilliance','bump_map',
     'bump_size','bumps','camera','case','catmull_rom_spline','caustics','ceil','cells','charset','checker','chr','circular',
     'clipped_by','clock','clock_delta','clock_on','collect','color','color_map','colour','colour_map','component','composite',
     'concat','cone','confidence','conic_sweep','conserve_energy','contained_by','control0','control1','coords','cos','cosh',
     'count','crackle','crand','cube','cubic','cubic_spline','cubic_wave','cutaway_textures','cylinder','cylindrical',
     'debug','declare','default','defined','degrees','density','density_file','density_map','dents','df3','difference',
     'diffuse','dimension_size','dimensions','direction','disc','dispersion','dispersion_samples','dist_exp','distance',
     'div','double_illuminate','eccentricity','else','emission','end','error','error_bound','evaluate','exp',
     'expand_thresholds','exponent','exterior','extinction','face_indices','facets','fade_color','fade_colour','fade_distance',
     'fade_power','falloff','falloff_angle','false','fclose','file_rexists','filter','final_clock','final_frame','finish','fisheye',
     'flatness','flip','floor','focal_point','fog','fog_alt','fog_offset','fog_type','fopen','form','frame_number','frequency',
     'fresnel','function','gather','gif','global_lights','global_settings','gradient','granite','gray','gray_threshold',
     'green','h_angle','height_field','hexagon','hf_gray_16','hierarchy','hollow','hypercomplex','if','ifdef','iff','ifndef'
     ,'image_height','image_map','image_pattern','image_width','include','initial_clock','initial_frame','inside','int',
     'interior','interior_texture','internal','interpolate','intersection','intervals','inverse','ior','irid',
     'irid_wavelength','isosurface','jitter','jpeg','julia','julia_fractal','lambda','lathe','leopard','light_group',
     'light_source','linear_spline','linear_sweep','ln','load_file','local','location','log','look_at','looks_like',
     'low_error_factor','macro','magnet','major_radius','mandel','map_type','marble','material','material_map','matrix',
     'max','max_extent','max_gradient','max_intersections','max_iteration','max_sample','max_trace','max_trace_level',
     'media','media_attenuation','media_interaction','merge','mesh','mesh2','metallic','method','metric','min',
     'min_extent','minimum_reuse','mod','mortar','nearest_count','no','no_bump_scale','no_image','no_reflection',
     'no_shadow','noise_generator','normal','normal_indices','normal_map','normal_vectors','number_of_waves',
     'object','octaves','off','offset','omega','omnimax','on','once','onion','open','orient','orientation','orthographic'
     ,'panoramic','parallel','parametric','pass_through','pattern','perspective','pgm','phase','phong',
     'phong_size','photons','pi','pigment','pigment_map','pigment_pattern','planar','plane','png','point_at',
     'poly','poly_wave','polygon','pot','pow','ppm','precision','precompute','pretrace_end','pretrace_start',
     'prism','projected_through','pwr','quadratic_spline','quadric','quartic','quaternion','quick_color','quick_colour',
     'quilted','radial','radians','radiosity','radius','rainbow','ramp_wave','rand','range','range_divider','ratio',
     'read','reciprocal','recursion_limit','red','reflection','reflection_exponent','refraction','render','repeat',
     'rgb','rgbf','rgbft','rgbt','right','ripples','rotate','roughness','samples','save_file','scale','scallop_wave',
     'scattering','seed','select','shadowless','sin','sine_wave','sinh','size','sky','sky_sphere','slice','slope',
     'slope_map','smooth','smooth_triangle','solid','sor','spacing','specular','sphere','sphere_sweep','spherical',
     'spiral1','spiral2','spline','split_union','spotlight','spotted','sqr','sqrt','statistics','str','strcmp','strength',
     'strlen','strlwr','strupr','sturm','substr','superellipsoid','switch','sys','t','tan','tanh','target','text','texture',
     'texture_list','texture_map','tga','thickness','threshold','tiff','tightness','tile2','tiles','tolerance','toroidal',
     'torus','trace','transform','translate','transmit','triangle','triangle_wave','true','ttf','turb_depth','turbulence',
     'type','u','u_steps','ultra_wide_angle','undef','union','up','use_alpha','use_color','use_colour','use_index','utf8',
     'uv_indices','uv_mapping','uv_vectors','v','v_angle','v_steps','val','variance','vaxis_rotate','vcross','vdot',
     'version','vertex_vectors','vlength','vnormalize','vrotate','vstr','vturbulence','warning','warp','water_level',
     'waves','while','width','wood','wrinkles','write','x','y','yes','z']



keywords={
    "obj"   :[
             'aa_level','aa_threshold','abs','absorption','accuracy','acos','acosh','adaptive','adc_bailout',
             'agate_turb','all','all_intersections','alpha','altitude','always_sample','ambient','ambient_light','angle',
             'aperture','append','arc_angle','area_light','array','asc','ascii','asin','asinh','assumed_gamma','atan','atan2',
             'atanh','autostop','b_spline','bezier_spline','black_hole','blue',
             'blur_samples','bounded_by','break','brick_size','brightness','brilliance','bump_map',
             'bump_size','case','catmull_rom_spline','caustics','ceil','charset','chr','circular',
             'clipped_by','clock','clock_delta','clock_on','collect','color','color_map','colour','colour_map','component','composite',
             'concat','confidence','conic_sweep','conserve_energy','contained_by','control0','control1','coords','cos','cosh',
             'count','crand','cube','cubic','cubic_spline','cubic_wave','cutaway_textures',
             'debug','declare','default','defined','degrees','density','density_map','df3',
             'diffuse','dimension_size','dimensions','direction','dispersion','dispersion_samples','dist_exp','distance',
             'div','double_illuminate','eccentricity','else','emission','end','error','error_bound','evaluate','exp',
             'expand_thresholds','exponent','exterior','extinction','face_indices','fade_color','fade_colour','fade_distance',
             'fade_power','falloff','falloff_angle','false','fclose','file_exists','filter','final_clock','final_frame','finish','fisheye',
             'flatness','flip','floor','focal_point','fog_alt','fog_offset','fog_type','fopen','frame_number','frequency',
             'fresnel','function','gather','gif','global_lights','gray','gray_threshold',
             'green','h_angle','hf_gray_16','hierarchy','hollow','hypercomplex','if','ifdef','iff','ifndef',
             'image_height','image_map','image_width','include','initial_clock','initial_frame','inside','int',
             'interior','interior_texture','internal','interpolate','intervals','inverse','ior','irid',
             'irid_wavelength','jitter','jpeg','lambda','light_group',
             'linear_spline','linear_sweep','ln','load_file','local','location','log','look_at','looks_like',
             'low_error_factor','macro','major_radius','map_type','material','material_map','matrix',
             'max','max_extent','max_gradient','max_intersections','max_iteration','max_sample','max_trace','max_trace_level',
             'media','media_attenuation','media_interaction','metallic','method','min',
             'min_extent','minimum_reuse','mod','mortar','nearest_count','no','no_bump_scale','no_image','no_reflection',
             'no_shadow','noise_generator','normal','normal_indices','normal_map','normal_vectors','number_of_waves',
             'octaves','off','omega','omnimax','on','once','open','orient','orientation','orthographic',
             'panoramic','parallel','parametric','pass_through','perspective','pgm','phase','phong',
             'phong_size','photons','pi','pigment','pigment_map','png','point_at',
             'poly_wave','pot','pow','ppm','precision','precompute','pretrace_end','pretrace_start',
             'projected_through','pwr','quadratic_spline','quaternion','quick_color','quick_colour',
             'radians','radius','ramp_wave','rand','range','range_divider','ratio',
             'read','reciprocal','recursion_limit','red','reflection','reflection_exponent','refraction','render','repeat',
             'rgb','rgbf','rgbft','rgbt','right','rotate','roughness','samples','save_file','scale','scallop_wave',
             'scattering','seed','select','shadowless','sin','sine_wave','sinh','size','sky','slice',
             'slope_map','smooth','sor','spacing','specular',
             'spline','split_union','spotlight','sqr','sqrt','statistics','str','strcmp','strength',
             'strlen','strlwr','strupr','sturm','substr','switch','sys','t','tan','tanh','target','texture',
             'texture_list','texture_map','tga','thickness','threshold','tiff','tightness','tolerance','toroidal',
             'trace','transform','translate','transmit','triangle_wave','true','ttf','turb_depth','turbulence',
             'type','u','u_steps','ultra_wide_angle','undef','up','use_alpha','use_color','use_colour','use_index','utf8',
             'uv_indices','uv_mapping','uv_vectors','v','v_angle','v_steps','val','variance','vaxis_rotate','vcross','vdot',
             'version','vertex_vectors','vlength','vnormalize','vrotate','vstr','vturbulence','warning','warp','water_level',
             'while','width','write','x','y','yes','z'],




    "objbold"  :[
            'bicubic_patch', 'box','blob','camera','cone','cylinder','difference','disc','height_field','julia_fractal','lathe',
            'prism','sphere','sphere_sweep','superellipsoid','text','torus','mesh','mesh2','polygon','triangle',
            'smooth_triangle','plane','poly','quadric','quartic','isosurface','sor',
            'union','merge','intersection','object',
            'global_settings','background','radiosity','light_source',
            'fog','sky_sphere','rainbow'],
    


   "pattern" :['agate','average','boxed','bozo','brick','bumps','cells',
              'checker','crackle','form','metric','offset','solid',
              'cylindrical','density_file','dents','facets','mandel','julia',
              'magnet','gradient','granite','hexagon','image_pattern','leopard',
              'marble','onion','pattern','pigment_pattern','planar','quilted',
              'radial','ripples','slope','spherical','spiral1','spiral2',
              'spotted','waves','wood','wrinkles','tile2','tiles']

    }

directives={
    "dirc"   :['declare','default','include','local','macro','switch','undef',
              'version','break','case','else','end','if','ifdef','ifndef','range','while',
              'fclose','fopen','read','write','debug','error','render','statistics','warning']
    }
