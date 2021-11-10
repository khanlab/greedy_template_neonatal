
# these rules are for registering the newly-generated cohort templates back 
# to the initial (MNI152) template, and composing the warps from template-building


rule reg_cohort_template_to_std:
    input:
        std_template = [ config['std_template'][channel]  for channel in channels ],
        cohort_template =  expand('results/cohort-{{cohort}}/iter_{iteration}/template_{channel}.nii.gz',iteration=config['max_iters'],channel=channels)
    params:
        input_fixed_moving = lambda wildcards, input: [f'-i {fixed} {moving}' for fixed,moving in zip(input.std_template, input.cohort_template) ],
        input_moving_warped = lambda wildcards, input, output: [f'-rm {moving} {warped}' for moving,warped in zip(input.cohort_template,output.warped) ],
    output:
        warp = 'results/cohort-{cohort}/reg_to_{std_template}/cohort-{cohort}_to-{std_template}_1Warp.nii.gz',
        invwarp = 'results/cohort-{cohort}/reg_to_{std_template}/cohort-{cohort}_to-{std_template}_1InverseWarp.nii.gz',
        affine_xfm_ras = 'results/cohort-{cohort}/reg_to_{std_template}/cohort-{cohort}_to-{std_template}_affine_ras.txt',
        warped = expand('results/cohort-{cohort}/reg_to_{std_template}/cohort-{cohort}_to-{std_template}_WarpedToTemplate_{channel}.nii.gz',channel=channels,allow_missing=True)
        
    log: 'logs/reg_cohort_template_to_std/cohort-{cohort}_{std_template}.log'
    threads: 8
    container: config['singularity']['itksnap']
    resources:
        # this is assuming 1mm
        mem_mb = 32000,
        time = 60
    group: 'reg_to_std'
    shell: 
        #affine first
        'greedy -d 3 -threads {threads} -a -m NCC 2x2x2 {params.input_fixed_moving} -o {output.affine_xfm_ras} -ia-image-centers -n 100x50x10 &> {log} && '
        #then deformable:
        'greedy -d 3 -threads {threads} -m NCC 2x2x2 {params.input_fixed_moving} -it {output.affine_xfm_ras} -o {output.warp} -oinv {output.invwarp} -n 100x50x10 &>> {log} && '

        #and finally warp the moving image
        'greedy -d 3 -threads {threads} -rf {input.std_template[0]} {params.input_moving_warped} -r {output.warp} {output.affine_xfm_ras} &>> {log}'



def get_inputs_composite_subj_to_std (wildcards):
    """ Function for setting all the inputs
        Needed since cohort isn't in the output filename, and is determined
        by looking at the input lists
    """
    for c in cohorts:
        if wildcards.subject in subjects[c]:
            cohort = c
    std_template = wildcards.std_template
    subject = wildcards.subject
    iteration=config['max_iters']

    return {
        'cohort2std_warp': f'results/cohort-{cohort}/reg_to_{std_template}/cohort-{cohort}_to-{std_template}_1Warp.nii.gz',
        'cohort2std_affine_xfm_ras': f'results/cohort-{cohort}/reg_to_{std_template}/cohort-{cohort}_to-{std_template}_affine_ras.txt',
        'subj2cohort_warp': f'results/cohort-{cohort}/iter_{iteration}/sub-{subject}_1Warp.nii.gz',
        'subj2cohort_affine_xfm_ras': f'results/cohort-{cohort}/iter_{iteration}/sub-{subject}_affine_ras.txt',
        'ref_std': config['init_template'][channels[0]] }
      

rule create_composite_subj_to_std:
    """ This concatenates the subject to cohort to mni warps/affines to get a single warp from subject to mni """
    input: unpack(get_inputs_composite_subj_to_std)
    output:
        subj2std_warp = 'results/composite/sub-{subject}_to-{std_template}_via-cohort_CompositeWarp.nii.gz'
    group: 'composite'
    threads: 8
    shell: 'greedy -d 3 -threads {threads} -rf {input.ref_std} '
          ' -r {input.cohort2std_warp} {input.cohort2std_affine_xfm_ras} '
          '  {input.subj2cohort_warp} {input.subj2cohort_affine_xfm_ras} '
          ' -rc {output.subj2std_warp}'



def get_inputs_composite_subj_to_std_inverse (wildcards):
    """ Function for setting all the inputs
        Needed since cohort isn't in the output filename, and is determined
        by looking at the input lists
    """
    for c in cohorts:
        if wildcards.subject in subjects[c]:
            cohort = c
    std_template = wildcards.std_template
    subject = wildcards.subject
    iteration=config['max_iters']

    return {
        'cohort2std_invwarp': f'results/cohort-{cohort}/reg_to_{std_template}/cohort-{cohort}_to-{std_template}_1InverseWarp.nii.gz',
        'cohort2std_affine_xfm_ras': f'results/cohort-{cohort}/reg_to_{std_template}/cohort-{cohort}_to-{std_template}_affine_ras.txt',
        'subj2cohort_invwarp': f'results/cohort-{cohort}/iter_{iteration}/sub-{subject}_1InverseWarp.nii.gz',
        'subj2cohort_affine_xfm_ras': f'results/cohort-{cohort}/iter_{iteration}/sub-{subject}_affine_ras.txt',
        'ref_subj':   config['in_images'][channels[0]] }
      



rule create_composite_subj_to_std_inverse:
    """ This concatenates the subject to cohort to mni warps/affines to get a single warp from subject to std_template"""
    input: unpack(get_inputs_composite_subj_to_std_inverse)
    output:
        subj2std_invwarp = 'results/composite/sub-{subject}_to-{std_template}_via-cohort_CompositeInverseWarp.nii.gz'
    group: 'composite'
    threads: 8
    shell: 'greedy -d 3  -threads {threads} -rf {input.ref_subj} -r '
          ' {input.subj2cohort_affine_xfm_ras},-1 '
          ' {input.subj2cohort_invwarp}'
          ' {input.cohort2std_affine_xfm_ras},-1 '
          ' {input.cohort2std_invwarp} '
          ' -rc {output.subj2std_invwarp}'



