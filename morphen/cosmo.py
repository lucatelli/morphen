"""
 ____  _     _
|  _ \(_)___| |_ __ _ _ __   ___ ___  ___
| | | | / __| __/ _` | '_ \ / __/ _ \/ __|
| |_| | \__ \ || (_| | | | | (_|  __/\__ \
|____/|_|___/\__\__,_|_| |_|\___\___||___/
  ____                                     _                         
 / ___|   ___    ___   _ __ ___     ___   | |   ___     __ _   _   _ 
| |      / _ \  / __| | '_ ` _ \   / _ \  | |  / _ \   / _` | | | | |
| |___  | (_) | \__ \ | | | | | | | (_) | | | | (_) | | (_| | | |_| |
 \____|  \___/  |___/ |_| |_| |_|  \___/  |_|  \___/   \__, |  \__, |
                                                       |___/   |___/ 

"""



def luminosity_distance_cosmo(z,Om0=0.308):
    h = 67.8  # * (h1 + h2) / 2
    cosmo = FlatLambdaCDM(H0=h, Om0=Om0)
    D_L = cosmo.luminosity_distance(z=z).value
    # print('D_l = ', D_L)  # 946.9318492873492 Mpc
    return(D_L)

def comoving_distance_cosmo(z,Om0=0.308):
    h = 67.8  # * (h1 + h2) / 2
    cosmo = FlatLambdaCDM(H0=h, Om0=Om0)
    D_C = cosmo.comoving_distance(z=z).value
    # print('D_l = ', D_L)  # 946.9318492873492 Mpc
    return(D_C)

def angular_distance_cosmo(z, Om0=0.308):
    h = 67.8# * (h1 + h2) / 2
    cosmo = FlatLambdaCDM(H0=h, Om0=Om0)
    d_A = cosmo.angular_diameter_distance(z=z)
    # print('D_a = ', d_A)  # 946.9318492873492 Mpc
    return(d_A)

def arcsec_to_pc(z, cell_size=1, Om0=0.308):
    h = 67.8# * (h1 + h2) / 2
    cosmo = FlatLambdaCDM(H0=h, Om0=Om0)
    d_A = cosmo.angular_diameter_distance(z=z)
    # print('D_a = ', d_A)  # 946.9318492873492 Mpc
    theta = 1 * u.arcsec
    distance_pc = (theta * d_A).to(u.pc, u.dimensionless_angles())
    # unit is Mpc only now
    # print('Linear Distance = ', distance_pc)  # 3.384745689510495 Mpc
    return (distance_pc)

    
    

def pc_to_arcsec(parsecs, redshift, h=67.8, Om0=0.308):
    """
    Convert a distance in parsecs to an angular size in arcseconds,
    based on the angular diameter distance to an object at redshift z.
    
    Parameters:
    parsecs (float): Distance in parsecs
    redshift (float): Redshift of the object
    h (float): Hubble constant in km/s/Mpc
    Om0 (float): Matter density parameter
    
    Returns:
    float: Angular size in arcseconds
    """
    # Create the cosmology object
    cosmo = FlatLambdaCDM(H0=h, Om0=Om0)
    
    # Calculate the angular diameter distance in Mpc
    d_A = cosmo.angular_diameter_distance(z=redshift)
    
    # Convert parsecs to the same units as d_A (Mpc)
    distance_mpc = parsecs * u.pc.to(u.Mpc)
    
    # Calculate the angle in radians: theta = physical_size / d_A
    theta_rad = distance_mpc / d_A.value
    
    # Convert radians to arcseconds
    theta_arcsec = theta_rad * u.rad.to(u.arcsec)
    
    return theta_arcsec

def pixsize_to_pc(z, cell_size, Om0=0.308):
    h = 67.8  # * (h1 + h2) / 2
    cosmo = FlatLambdaCDM(H0=h, Om0=Om0)
    d_A = cosmo.angular_diameter_distance(z=z)
    # print('D_a = ', d_A)  # 946.9318492873492 Mpc
    theta = cell_size * u.arcsec
    distance_pc = (theta * d_A).to(u.pc, u.dimensionless_angles())  # unit is Mpc only now

    # print('Linear Distance = ', distance_pc)  # 3.384745689510495 Mpc
    return (distance_pc.value)



def cosmo_stats(imagename,z,results=None):
    """
    Get beam shape info in physical units.
    """
    if results is None:
        results = {}
        results['#imagename'] = os.path.basename(imagename)
    pc_scale, bmaj, bmin, BA_pc = beam_physical_area(imagename, z=z)
    results['arcsec_to_pc'] = pc_scale.value
    results['bmaj_pc'] = bmaj.value
    results['bmin_pc'] = bmin.value
    results['BA_pc'] = BA_pc.value
    return(results)

# def find_z_NED(source_name):
#     """
#     Find the redshift of a source (by name) from NED.

#     Parameters
#     ----------
#     source_name : str
#         Source name.
#             Example: 'VV705'

#     Returns
#     -------
#     redshift_NED : float, None
#     """
#     from astroquery.ipac.ned import Ned
#     result_table = Ned.query_object(source_name)
#     redshift_NED = result_table['Redshift'].data.data
#     if redshift_NED.shape[0] == 0:
#         return None
#     else:
#         return redshift_NED[0]

def find_z_NED(source_name, return_luminosity_distance=False):
    """
    Find the redshift of a source (by name) from NED.

    Parameters
    ----------
    source_name : str
        Source name.
            Example: 'VV705'
    return_luminosity_distance : bool, optional
        If True, also calculate and return the luminosity distance.
        Default: False

    Returns
    -------
    redshift_NED : float, None
        The redshift value from NED
    luminosity_distance : astropy.units.Quantity, optional
        The luminosity distance (only returned if return_luminosity_distance=True)
    """
    from astroquery.ipac.ned import Ned
    
    result_table = Ned.query_object(source_name)
    redshift_NED = result_table['Redshift'].data.data
    
    if redshift_NED.shape[0] == 0:
        return None
    else:
        z = redshift_NED[0]
        
        if return_luminosity_distance:
            from astropy.cosmology import Planck18
            
            # Calculate luminosity distance
            lum_dist = Planck18.luminosity_distance(z)
            
            # Print the luminosity distance
            print(f"Luminosity distance for {source_name} (z={z:.4f}): {lum_dist:.2f}")
            
            return z, lum_dist
        else:
            return z
