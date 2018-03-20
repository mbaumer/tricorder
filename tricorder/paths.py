config_dir = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/'

corr_out_dir = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/test2/'
ang_out_dir = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/test3_angular/'

rm_y1 = [
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/a/buzzard_1.6-6a_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/b/buzzard_1.6-6b_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/c/buzzard_1.6-6c_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/d/buzzard_1.6-6d_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/e/buzzard_1.6-6e_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/f/buzzard_1.6-6f_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',

    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-1/a/buzzard-1_1.6-6a_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/rdmagic/y1/buzzard/flock/buzzard-1/b/buzzard-1_1.6-6b_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-1/c/buzzard-1_1.6-6c_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-1/d/buzzard-1_1.6-6d_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-1/e/buzzard-1_1.6-6e_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-1/f/buzzard-1_1.6-6f_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',

    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/a/buzzard2_1.6-6a_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/b/buzzard2_1.6-6b_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/c/buzzard2_1.6-6c_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/d/buzzard2_1.6-6d_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/e/buzzard2_1.6-6e_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/f/buzzard2_1.6-6f_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',

    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-3/a/buzzard-3_1.6-6a_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-3/b/buzzard-3_1.6-6b_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-3/c/buzzard-3_1.6-6c_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-3/d/buzzard-3_1.6-6d_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-3/e/buzzard-3_1.6-6e_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-3/f/buzzard-3_1.6-6f_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',

    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/a/buzzard5_1.6-6a_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/b/buzzard5_1.6-6b_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/c/buzzard5_1.6-6c_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/d/buzzard5_1.6-6d_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/e/buzzard5_1.6-6e_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/f/buzzard5_1.6-6f_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',

    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-21/a/buzzard21_1.6-6a_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-21/b/buzzard21_1.6-6b_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-21/c/buzzard21_1.6-6c_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-21/d/buzzard21_1.6-6d_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-21/e/buzzard21_1.6-6e_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-21/f/buzzard21_1.6-6f_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
]

lss_y1 = [
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-0/a/Buzzard_v1.6_Y1a_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-0/b/Buzzard_v1.6_Y1b_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-0/c/Buzzard_v1.6_Y1c_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-0/d/Buzzard_v1.6_Y1d_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-0/e/Buzzard_v1.6_Y1e_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-0/f/Buzzard_v1.6_Y1f_gold.fits',

    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-1/a/Buzzard_v1.6_Y1a_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-1/b/Buzzard_v1.6_Y1b_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-1/c/Buzzard_v1.6_Y1c_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-1/d/Buzzard_v1.6_Y1d_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-1/e/Buzzard_v1.6_Y1e_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-1/f/Buzzard_v1.6_Y1f_gold.fits',

    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-2/a/Buzzard_v1.6_Y1a_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-2/b/Buzzard_v1.6_Y1b_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-2/c/Buzzard_v1.6_Y1c_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-2/d/Buzzard_v1.6_Y1d_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-2/e/Buzzard_v1.6_Y1e_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-2/f/Buzzard_v1.6_Y1f_gold.fits',

    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-3/a/Buzzard_v1.6_Y1a_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-3/b/Buzzard_v1.6_Y1b_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-3/c/Buzzard_v1.6_Y1c_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-3/d/Buzzard_v1.6_Y1d_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-3/e/Buzzard_v1.6_Y1e_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-3/f/Buzzard_v1.6_Y1f_gold.fits',

    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-5/a/Buzzard_v1.6_Y1a_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-5/b/Buzzard_v1.6_Y1b_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-5/c/Buzzard_v1.6_Y1c_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-5/d/Buzzard_v1.6_Y1d_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-5/e/Buzzard_v1.6_Y1e_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-5/f/Buzzard_v1.6_Y1f_gold.fits',

    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-21/a/Buzzard_v1.6_Y1a_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-21/b/Buzzard_v1.6_Y1b_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-21/c/Buzzard_v1.6_Y1c_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-21/d/Buzzard_v1.6_Y1d_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-21/e/Buzzard_v1.6_Y1e_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-21/f/Buzzard_v1.6_Y1f_gold.fits',
]

dm_y1 = [
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-0/a/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-0/b/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-0/c/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-0/d/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-0/e/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-0/f/downsampled_particles.fits.downsample',

    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-1/a/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-1/b/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-1/c/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-1/d/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-1/e/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-1/f/downsampled_particles.fits.downsample',

    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-2/a/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-2/b/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-2/c/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-2/d/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-2/e/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-2/f/downsampled_particles.fits.downsample',

    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-5/a/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-5/b/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-5/c/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-5/d/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-5/e/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y1/buzzard/flock/buzzard-5/f/downsampled_particles.fits.downsample',
]
