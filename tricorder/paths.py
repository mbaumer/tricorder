
# CAUTION: Only use NFS paths!! Otherwise you'll have disk access problems...

config_dir = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/'

corr_out_dir = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/paper_validation'
ang_out_dir = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new3/tolerance_syst'

rm_y1_randoms = '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3a/buzzard-3_1.9.2+1-6a_run_redmapper_v6.4.22_redmagic_highdens_0.5-10_randoms.fit'
rm_y1_HL_randoms = '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3a/buzzard-3_1.9.2+1-6a_run_redmapper_v6.4.22_redmagic_highlum_1.0-04_randoms.fit'
rm_y1_HHL_randoms = '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3a/buzzard-3_1.9.2+1-6a_run_redmapper_v6.4.22_redmagic_higherlum_1.5-01_randoms.fit'

#dm_y1_randoms = '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/a/buzzard_1.6-6a_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01_randoms.fit'
#dm_y1_randoms = '/nfs/slac/des/fs1/g/sims/erykoff/clusters/mocks/Buzzard/buzzard-1.6/des-y1a1/redmapper_v6.4.18/redmagic_a/buzzard_1.6-6a_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01_randoms.fit'
dm_y1_randoms = '/nfs/slac/g/ki/ki19/des/mbaumer/dark_matter_joe/buzzard_1.6-6a_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01_randoms.fit'

rm_y1 = [
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3a/buzzard-3_1.9.2+1-6a_run_redmapper_v6.4.22_redmagic_highdens_0.5-10.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3b/buzzard-3_1.9.2+1-6b_run_redmapper_v6.4.22_redmagic_highdens_0.5-10.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3c/buzzard-3_1.9.2+1-6c_run_redmapper_v6.4.22_redmagic_highdens_0.5-10.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3d/buzzard-3_1.9.2+1-6d_run_redmapper_v6.4.22_redmagic_highdens_0.5-10.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3e/buzzard-3_1.9.2+1-6e_run_redmapper_v6.4.22_redmagic_highdens_0.5-10.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3f/buzzard-3_1.9.2+1-6f_run_redmapper_v6.4.22_redmagic_highdens_0.5-10.fit',

    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_4a/buzzard-4_1.9.2+1-6a_run_redmapper_v6.4.22_redmagic_highdens_0.5-10.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_4b/buzzard-4_1.9.2+1-6b_run_redmapper_v6.4.22_redmagic_highdens_0.5-10.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_4c/buzzard-4_1.9.2+1-6c_run_redmapper_v6.4.22_redmagic_highdens_0.5-10.fit',
    # this one is bad b/c of joe healpix error!
    #'/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_4d/buzzard-4_1.9.2+1-6d_run_redmapper_v6.4.22_redmagic_highdens_0.5-10.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_4e/buzzard-4_1.9.2+1-6e_run_redmapper_v6.4.22_redmagic_highdens_0.5-10.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_4f/buzzard-4_1.9.2+1-6f_run_redmapper_v6.4.22_redmagic_highdens_0.5-10.fit',
]

rm_y1_HL = [
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3a/buzzard-3_1.9.2+1-6a_run_redmapper_v6.4.22_redmagic_highlum_1.0-04.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3b/buzzard-3_1.9.2+1-6b_run_redmapper_v6.4.22_redmagic_highlum_1.0-04.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3c/buzzard-3_1.9.2+1-6c_run_redmapper_v6.4.22_redmagic_highlum_1.0-04.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3d/buzzard-3_1.9.2+1-6d_run_redmapper_v6.4.22_redmagic_highlum_1.0-04.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3e/buzzard-3_1.9.2+1-6e_run_redmapper_v6.4.22_redmagic_highlum_1.0-04.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3f/buzzard-3_1.9.2+1-6f_run_redmapper_v6.4.22_redmagic_highlum_1.0-04.fit',

    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_4a/buzzard-4_1.9.2+1-6a_run_redmapper_v6.4.22_redmagic_highlum_1.0-04.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_4b/buzzard-4_1.9.2+1-6b_run_redmapper_v6.4.22_redmagic_highlum_1.0-04.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_4c/buzzard-4_1.9.2+1-6c_run_redmapper_v6.4.22_redmagic_highlum_1.0-04.fit',
    #'/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_4d/buzzard-4_1.9.2+1-6d_run_redmapper_v6.4.22_redmagic_highlum_1.0-04.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_4e/buzzard-4_1.9.2+1-6e_run_redmapper_v6.4.22_redmagic_highlum_1.0-04.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_4f/buzzard-4_1.9.2+1-6f_run_redmapper_v6.4.22_redmagic_highlum_1.0-04.fit',
]

rm_y1_HHL = [
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3a/buzzard-3_1.9.2+1-6a_run_redmapper_v6.4.22_redmagic_higherlum_1.5-01.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3b/buzzard-3_1.9.2+1-6b_run_redmapper_v6.4.22_redmagic_higherlum_1.5-01.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3c/buzzard-3_1.9.2+1-6c_run_redmapper_v6.4.22_redmagic_higherlum_1.5-01.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3d/buzzard-3_1.9.2+1-6d_run_redmapper_v6.4.22_redmagic_higherlum_1.5-01.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3e/buzzard-3_1.9.2+1-6e_run_redmapper_v6.4.22_redmagic_higherlum_1.5-01.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_3f/buzzard-3_1.9.2+1-6f_run_redmapper_v6.4.22_redmagic_higherlum_1.5-01.fit',

    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_4a/buzzard-4_1.9.2+1-6a_run_redmapper_v6.4.22_redmagic_higherlum_1.5-01.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_4b/buzzard-4_1.9.2+1-6b_run_redmapper_v6.4.22_redmagic_higherlum_1.5-01.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_4c/buzzard-4_1.9.2+1-6c_run_redmapper_v6.4.22_redmagic_higherlum_1.5-01.fit',
    #'/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_4d/buzzard-4_1.9.2+1-6d_run_redmapper_v6.4.22_redmagic_higherlum_1.5-01.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_4e/buzzard-4_1.9.2+1-6e_run_redmapper_v6.4.22_redmagic_higherlum_1.5-01.fit',
    '/nfs/slac/g/ki/ki19/des/erykoff/clusters/mocks/Buzzard/buzzard-1.9.2+1/des-y1a1/redmapper_v6.4.22/redmagic_4f/buzzard-4_1.9.2+1-6f_run_redmapper_v6.4.22_redmagic_higherlum_1.5-01.fit',
]

dm_y1 = [
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-0/downsample_particles/y1a/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-0/downsample_particles/y1b/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-0/downsample_particles/y1c/downsampled_particles.fits',
    #'/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-0/downsample_particles/y1d/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-0/downsample_particles/y1e/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-0/downsample_particles/y1f/downsampled_particles.fits',

    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-1/downsample_particles/y1a/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-1/downsample_particles/y1b/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-1/downsample_particles/y1c/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-1/downsample_particles/y1d/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-1/downsample_particles/y1e/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-1/downsample_particles/y1f/downsampled_particles.fits',

    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-2/downsample_particles/y1a/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-2/downsample_particles/y1b/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-2/downsample_particles/y1c/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-2/downsample_particles/y1d/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-2/downsample_particles/y1e/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-2/downsample_particles/y1f/downsampled_particles.fits',

    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-5/downsample_particles/y1a/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-5/downsample_particles/y1b/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-5/downsample_particles/y1c/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-5/downsample_particles/y1d/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-5/downsample_particles/y1e/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-5/downsample_particles/y1f/downsampled_particles.fits',

]
