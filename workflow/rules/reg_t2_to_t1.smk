


rule n4biasfield:
    input: 
        img = lambda wildcards: config['in_raw'][wildcards.modality],
    output:
        img = bids(root='results/preproc',suffix='{modality}.nii.gz',desc='n4',subject='{subject}'),
    threads: 1
    container: config['singularity']['prepdwi']
    group: 'preproc'
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'N4BiasFieldCorrection -d 3 -i {input.img} -o {output}'


rule reg_aladin_t2w_to_t1:
    input: 
        ref = bids(root='results/preproc',suffix='T1w.nii.gz',desc='n4',subject='{subject}'),
        flo = bids(root='results/preproc',suffix='T2w.nii.gz',desc='n4',subject='{subject}'),
    output: 
        warped = bids(root='results/preproc',suffix='T2w.nii.gz',space='T1w',desc='n4',subject='{subject}'),
        xfm_ras = bids(root='results/preproc',suffix='xfm.txt',from_='T2w',to='T1w',type_='ras',subject='{subject}'),
    container: config['singularity']['prepdwi']
    group: 'preproc'
    shell:
        'reg_aladin -rigOnly -flo {input.flo} -ref {input.ref} -res {output.warped} -aff {output.xfm_ras}'

