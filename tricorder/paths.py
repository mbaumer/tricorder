config_dir = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/configs/'

corr_out_dir = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/test5/'
ang_out_dir = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/paper3/'

rm_y1_randoms = '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/a/buzzard_1.6-6a_run_redmapper_v6.4.18_redmagic_highdens_0.5-10_randoms.fit'
rm_y1_HL_randoms = '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/a/buzzard_1.6-6a_run_redmapper_v6.4.18_redmagic_highlum_1.0-04_randoms.fit'
rm_y1_HHL_randoms = '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/a/buzzard_1.6-6a_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01_randoms.fit'

rm_y3_randoms = '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-0/a/buzzard-0_1.6_y3_run_redmapper_v6.4.20_redmagic_highdens_0.5-10_randoms.fit'
rm_y3_HL_randoms = '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-0/a/buzzard-0_1.6_y3_run_redmapper_v6.4.20_redmagic_highlum_1.0-04_randoms.fit'
rm_y3_HHL_randoms = '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-0/a/buzzard-0_1.6_y3_run_redmapper_v6.4.20_redmagic_higherlum_1.5-01_randoms.fit'

dm_y1_randoms = '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/a/buzzard_1.6-6a_run_redmapper_v6.4.18_redmagic_highdens_0.5-10_randoms.fit'
dm_y3_randoms = '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-0/a/buzzard-0_1.6_y3_run_redmapper_v6.4.20_redmagic_highdens_0.5-10_randoms.fit'
lss_y1_randoms = '/nfs/slac/des/fs1/g/sims/mbaumer/3pt_sims/new2/lss_y1_randoms.fits'

rm_mice_y1 = [
    '/nfs/slac/g/ki/ki23/des/jderose/des/MICE/redmagic/mice2_des_run_redmapper_v6.4.16_redmagic_highdens_0.5-10.fit']

rm_y1 = [
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/a/buzzard_1.6-6a_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/b/buzzard_1.6-6b_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/c/buzzard_1.6-6c_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/d/buzzard_1.6-6d_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/e/buzzard_1.6-6e_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/f/buzzard_1.6-6f_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',

    # 2pt bias visibly low at small angles for these runs
    #'/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-1/a/buzzard-1_1.6-6a_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    #'/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-1/b/buzzard-1_1.6-6b_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    #'/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-1/c/buzzard-1_1.6-6c_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    #'/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-1/d/buzzard-1_1.6-6d_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    #'/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-1/e/buzzard-1_1.6-6e_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    #'/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-1/f/buzzard-1_1.6-6f_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',

    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/a/buzzard2_1.6-6a_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/b/buzzard2_1.6-6b_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/c/buzzard2_1.6-6c_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/d/buzzard2_1.6-6d_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/e/buzzard2_1.6-6e_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/f/buzzard2_1.6-6f_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',

    # 2pt bias also visibly low for these:
    #'/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-3/a/buzzard-3_1.6-6a_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    #'/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-3/b/buzzard-3_1.6-6b_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    #'/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-3/c/buzzard-3_1.6-6c_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    #'/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-3/d/buzzard-3_1.6-6d_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    #'/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-3/e/buzzard-3_1.6-6e_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    #'/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-3/f/buzzard-3_1.6-6f_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',

    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/a/buzzard5_1.6-6a_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/b/buzzard5_1.6-6b_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/c/buzzard5_1.6-6c_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/d/buzzard5_1.6-6d_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/e/buzzard5_1.6-6e_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/f/buzzard5_1.6-6f_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',

    # this one not visibly bad, but since there are no DM catalogs for it;
    # '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-21/a/buzzard21_1.6-6a_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    # '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-21/b/buzzard21_1.6-6b_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    # '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-21/c/buzzard21_1.6-6c_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    # '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-21/d/buzzard21_1.6-6d_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    # '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-21/e/buzzard21_1.6-6e_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
    # '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-21/f/buzzard21_1.6-6f_run_redmapper_v6.4.18_redmagic_highdens_0.5-10.fit',
]

rm_y1_HL = [
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/a/buzzard_1.6-6a_run_redmapper_v6.4.18_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/b/buzzard_1.6-6b_run_redmapper_v6.4.18_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/c/buzzard_1.6-6c_run_redmapper_v6.4.18_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/d/buzzard_1.6-6d_run_redmapper_v6.4.18_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/e/buzzard_1.6-6e_run_redmapper_v6.4.18_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/f/buzzard_1.6-6f_run_redmapper_v6.4.18_redmagic_highlum_1.0-04.fit',

    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/a/buzzard2_1.6-6a_run_redmapper_v6.4.18_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/b/buzzard2_1.6-6b_run_redmapper_v6.4.18_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/c/buzzard2_1.6-6c_run_redmapper_v6.4.18_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/d/buzzard2_1.6-6d_run_redmapper_v6.4.18_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/e/buzzard2_1.6-6e_run_redmapper_v6.4.18_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/f/buzzard2_1.6-6f_run_redmapper_v6.4.18_redmagic_highlum_1.0-04.fit',

    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/a/buzzard5_1.6-6a_run_redmapper_v6.4.18_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/b/buzzard5_1.6-6b_run_redmapper_v6.4.18_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/c/buzzard5_1.6-6c_run_redmapper_v6.4.18_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/d/buzzard5_1.6-6d_run_redmapper_v6.4.18_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/e/buzzard5_1.6-6e_run_redmapper_v6.4.18_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/f/buzzard5_1.6-6f_run_redmapper_v6.4.18_redmagic_highlum_1.0-04.fit',
]

rm_y1_HHL = [
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/a/buzzard_1.6-6a_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/b/buzzard_1.6-6b_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/c/buzzard_1.6-6c_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/d/buzzard_1.6-6d_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/e/buzzard_1.6-6e_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-0/f/buzzard_1.6-6f_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01.fit',

    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/a/buzzard2_1.6-6a_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/b/buzzard2_1.6-6b_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/c/buzzard2_1.6-6c_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/d/buzzard2_1.6-6d_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/e/buzzard2_1.6-6e_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-2/f/buzzard2_1.6-6f_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01.fit',

    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/a/buzzard5_1.6-6a_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/b/buzzard5_1.6-6b_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/c/buzzard5_1.6-6c_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/d/buzzard5_1.6-6d_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/e/buzzard5_1.6-6e_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y1/buzzard/flock/buzzard-5/f/buzzard5_1.6-6f_run_redmapper_v6.4.18_redmagic_higherlum_1.5-01.fit',

]

rm_y3 = [
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-0/a/buzzard-0_1.6_y3_run_redmapper_v6.4.20_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-1/a/buzzard-1_1.6_y3_run_redmapper_v6.4.20_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-2/a/buzzard-2_1.6_y3_run_redmapper_v6.4.20_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-3/a/buzzard-3_1.6_y3_run_redmapper_v6.4.20_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-5/a/buzzard-5_1.6_y3_run_redmapper_v6.4.20_redmagic_highdens_0.5-10.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-21/a/buzzard-21_1.6_y3_run_redmapper_v6.4.20_redmagic_highdens_0.5-10.fit',
]

rm_y3_HL = [
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-0/a/buzzard-0_1.6_y3_run_redmapper_v6.4.20_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-1/a/buzzard-1_1.6_y3_run_redmapper_v6.4.20_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-2/a/buzzard-2_1.6_y3_run_redmapper_v6.4.20_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-3/a/buzzard-3_1.6_y3_run_redmapper_v6.4.20_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-5/a/buzzard-5_1.6_y3_run_redmapper_v6.4.20_redmagic_highlum_1.0-04.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-21/a/buzzard-21_1.6_y3_run_redmapper_v6.4.20_redmagic_highlum_1.0-04.fit'
]

rm_y3_HHL = [
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-0/a/buzzard-0_1.6_y3_run_redmapper_v6.4.20_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-1/a/buzzard-1_1.6_y3_run_redmapper_v6.4.20_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-2/a/buzzard-2_1.6_y3_run_redmapper_v6.4.20_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-3/a/buzzard-3_1.6_y3_run_redmapper_v6.4.20_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-5/a/buzzard-5_1.6_y3_run_redmapper_v6.4.20_redmagic_higherlum_1.5-01.fit',
    '/u/ki/jderose/public_html/bcc/catalog/redmagic/y3/buzzard/flock/buzzard-21/a/buzzard-21_1.6_y3_run_redmapper_v6.4.20_redmagic_higherlum_1.5-01.fit'
]

lss_y1 = [
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-0/a/Buzzard_v1.6_Y1a_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-0/b/Buzzard_v1.6_Y1b_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-0/c/Buzzard_v1.6_Y1c_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-0/d/Buzzard_v1.6_Y1d_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-0/e/Buzzard_v1.6_Y1e_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-0/f/Buzzard_v1.6_Y1f_gold.fits',

    # '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-1/a/Buzzard_v1.6_Y1a_gold.fits',
    # '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-1/b/Buzzard_v1.6_Y1b_gold.fits',
    # '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-1/c/Buzzard_v1.6_Y1c_gold.fits',
    # '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-1/d/Buzzard_v1.6_Y1d_gold.fits',
    # '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-1/e/Buzzard_v1.6_Y1e_gold.fits',
    # '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-1/f/Buzzard_v1.6_Y1f_gold.fits',

    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-2/a/Buzzard_v1.6_Y1a_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-2/b/Buzzard_v1.6_Y1b_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-2/c/Buzzard_v1.6_Y1c_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-2/d/Buzzard_v1.6_Y1d_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-2/e/Buzzard_v1.6_Y1e_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-2/f/Buzzard_v1.6_Y1f_gold.fits',

    # '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-3/a/Buzzard_v1.6_Y1a_gold.fits',
    # '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-3/b/Buzzard_v1.6_Y1b_gold.fits',
    # '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-3/c/Buzzard_v1.6_Y1c_gold.fits',
    # '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-3/d/Buzzard_v1.6_Y1d_gold.fits',
    # '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-3/e/Buzzard_v1.6_Y1e_gold.fits',
    # '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-3/f/Buzzard_v1.6_Y1f_gold.fits',

    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-5/a/Buzzard_v1.6_Y1a_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-5/b/Buzzard_v1.6_Y1b_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-5/c/Buzzard_v1.6_Y1c_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-5/d/Buzzard_v1.6_Y1d_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-5/e/Buzzard_v1.6_Y1e_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-5/f/Buzzard_v1.6_Y1f_gold.fits',


    # this one not visibly bad, but since there are no DM catalogs for it;
    # '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-21/a/Buzzard_v1.6_Y1a_gold.fits',
    # '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-21/b/Buzzard_v1.6_Y1b_gold.fits',
    # '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-21/c/Buzzard_v1.6_Y1c_gold.fits',
    # '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-21/d/Buzzard_v1.6_Y1d_gold.fits',
    # '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-21/e/Buzzard_v1.6_Y1e_gold.fits',
    # '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y1/buzzard/flock/buzzard-21/f/Buzzard_v1.6_Y1f_gold.fits',
]

dm_y1 = [
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-0/downsample_particles/y1a/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-0/downsample_particles/y1b/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-0/downsample_particles/y1c/downsampled_particles.fits',
    '/nfs/slac/des/fs1/g/sims/jderose/BCC/Chinchilla/Herd//Chinchilla-0/downsample_particles/y1d/downsampled_particles.fits',
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

dm_y3 = [
    '/u/ki/jderose/public_html/bcc/catalog/particles/y3/buzzard/flock/buzzard-0/a/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y3/buzzard/flock/buzzard-1/a/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y3/buzzard/flock/buzzard-2/a/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y3/buzzard/flock/buzzard-3/a/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y3/buzzard/flock/buzzard-5/a/downsampled_particles.fits.downsample',
    '/u/ki/jderose/public_html/bcc/catalog/particles/y3/buzzard/flock/buzzard-21/a/downsampled_particles.fits.downsample',
]

lss_y3 = [
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y3/buzzard/flock/buzzard-21/a/Buzzard_v1.6_Y3_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y3/buzzard/flock/buzzard-5/a/Buzzard_v1.6_Y3_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y3/buzzard/flock/buzzard-3/a/Buzzard_v1.6_Y3_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y3/buzzard/flock/buzzard-2/a/Buzzard_v1.6_Y3_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y3/buzzard/flock/buzzard-1/a/Buzzard_v1.6_Y3_gold.fits',
    '/u/ki/jderose/public_html/bcc/catalog/mergedcats/y3/buzzard/flock/buzzard-0/a/Buzzard_v1.6_Y3_gold.fits',
]
