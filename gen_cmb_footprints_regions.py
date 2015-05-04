import numpy as np
import pylab as pl
import healpy as H

h = 360.0 / 24.0
m = h / 60.0
s = m / 60.0

d = 1.0
am = 1.0 / 60.0

def gen_bicep2(fp):
    #fn_b2 = 'maps/BICEP2/B2_3yr_373sqdeg_field_20140509.txt'
    #vtx_b2 = np.loadtxt(fn_b2,delimiter=',')

    vtx_ra = np.array([-17.53,-21.81,-25.53,-27.73,-27.73,-25.53,-21.81,-17.53,17.53,21.81,25.53,27.73,27.73,25.53,
                        21.81, 17.53])
    vtx_dec = np.array([-50.88,-51.12,-51.96,-53.2,-61.56,-62.8,-63.64,-63.88,-63.88,-63.64,-62.8,-61.56,-53.2,
                        -51.96,-51.12,-50.88])
    
    vtx_b2_b = np.transpose([vtx_ra,vtx_dec])
    
    idx = np.array([True,False,False,True,True,False,False,True,True,False,False,True,True,False,False,True])

    pl.figure()
    pl.scatter(vtx_b2_b[idx,0],vtx_b2_b[idx,1])
    pl.show()

    map_b2 = fp.superimpose_bound_vtx(vtx_b2_b[idx,:],'BICEP2',color='cyan',return_map=True)
    map_b2 = H.ud_grade(map_b2,256)
    H.write_map('maps/BICEP2_CMB_hpx.fits',map_b2)

def gen_pb(fp):
    ra23_cen = (23*h + 1*m + 48*s, -32.0 - 48 / 60.0)
    ra12_cen = (11*h + 53*m + 0*s, -30.0 / 60.0)
    ra4p5_cen = (4*h + 40*m + 12*s, -45.0)

    ra23_size = (1.75,1.75)
    ra12_size = (1.75,1.75)
    ra4p5_size = (1.75,1.75)
   
    map_ra23 = fp.superimpose_bound_cen(ra23_cen,ra23_size,'RA23',color='blue',return_map=True)
    map_ra12 = fp.superimpose_bound_cen(ra12_cen,ra12_size,'RA12',color='red',return_map=True)
    map_ra4p5 = fp.superimpose_bound_cen(ra4p5_cen,ra4p5_size,'RA4p5',color='green',return_map=True)
  
    map_pb1 = np.zeros_like(map_ra23)

    idx = np.isfinite(map_ra23)
    map_pb1[idx] += map_ra23[idx]
    print np.sum(map_ra23[idx]) / len(map_ra23)*100.0, 8.8 / 41253.0 * 100.0
    idx = np.isfinite(map_ra12)
    map_pb1[idx] += map_ra12[idx]
    print np.sum(map_ra12[idx]) / len(map_ra12)*100.0, 8.7 / 41253.0 * 100.0
    idx = np.isfinite(map_ra4p5)
    map_pb1[idx] += map_ra4p5[idx]
    print np.sum(map_ra4p5[idx]) / len(map_ra4p5)*100.0, 7.0 / 41253.0 * 100.0

    map_pb1 = H.ud_grade(map_pb1,256)

    H.write_map('maps/PB1_S1_hpx.fits',map_pb1)
    
def gen_quiet(fp):
    quiet_rad = 8.9

    quiet_cmb1 = (12*h + 4*m, -39*d)
    quiet_cmb2 = (5*h + 12*m, -39*d)
    quiet_cmb3 = (0*h + 48*m, -48*d)
    quiet_cmb4 = (22*h + 44*m, -36*d)
    quiet_g1 = (16*h + 0*m, -53*d)
    quiet_g2 = (17*h + 46*m, -28*d - 56*am)

    map_quiet_cmb1 = fp.superimpose_bound_circ(quiet_cmb1,quiet_rad,'CMB1',color='blue',return_map=True)
    map_quiet_cmb2 = fp.superimpose_bound_circ(quiet_cmb2,quiet_rad,'CMB2',color='red',return_map=True)
    map_quiet_cmb3 = fp.superimpose_bound_circ(quiet_cmb3,quiet_rad,'CMB3',color='green',return_map=True)
    map_quiet_cmb4 = fp.superimpose_bound_circ(quiet_cmb4,quiet_rad,'CMB4',color='magenta',return_map=True)
    map_quiet_g1 = fp.superimpose_bound_circ(quiet_g1,quiet_rad,'G1',color='yellow',return_map=True)
    map_quiet_g2 = fp.superimpose_bound_circ(quiet_g2,quiet_rad,'G2',color='orange',return_map=True)
    
    map_quiet = np.zeros_like(map_quiet_cmb1)
    
    idx = np.isfinite(map_quiet_cmb1)
    map_quiet[idx] += map_quiet_cmb1[idx]
    print np.sum(map_quiet_cmb1[idx]) / len(map_quiet_cmb1)*100.0, 250.0 / 41253.0 * 100.0
    idx = np.isfinite(map_quiet_cmb2)
    map_quiet[idx] += map_quiet_cmb2[idx]
    idx = np.isfinite(map_quiet_cmb3)
    map_quiet[idx] += map_quiet_cmb3[idx]
    idx = np.isfinite(map_quiet_cmb4)
    map_quiet[idx] += map_quiet_cmb4[idx]
    idx = np.isfinite(map_quiet_g1)
    map_quiet[idx] += map_quiet_g1[idx]
    idx = np.isfinite(map_quiet_g2)
    map_quiet[idx] += map_quiet_g2[idx]
    
    map_quiet = H.ud_grade(map_quiet,256)
    
    H.write_map('maps/QUIET_hpx.fits',map_quiet)

def gen_B03(fp):
    B03_shallow = ( (61*d,-31*d), (104*d,-31*d), (104*d,-58*d), (61*d, -58*d) )
    B03_deep = ( (76*d, -40*d), (89*d, -40*d), (89*d, -50*d), (76*d, -50*d) )

    map_B03_shallow = fp.superimpose_bound_vtx(B03_shallow,'Shallow',color='blue',return_map=True)
    map_B03_deep = fp.superimpose_bound_vtx(B03_deep,'Deep',color='red',return_map=True)

    map_B03 = np.zeros_like(map_B03_shallow)
    
    idx = np.isfinite(map_B03_shallow)
    map_B03[idx] += map_B03_shallow[idx]
    idx = np.isfinite(map_B03_deep)
    map_B03[idx] += map_B03_deep[idx]
    
    map_B03 = H.ud_grade(map_B03,256)
    
    H.write_map('maps/BOOMERANG_hpx.fits',map_B03)

if __name__ == '__main__':
    fn_background = 'maps/Planck/COM_CompMap_DustPol-commander_1024_R2.00.fits'

    data = H.read_map(fn_background,field=(0,1))
    background_map = np.sqrt(data[0]**2 + data[1]**2)

    fp = footprint.SurveyStack(background_map, fignum=1, projection=H.mollview, coord_plot='C', rot=[0,0]) 
