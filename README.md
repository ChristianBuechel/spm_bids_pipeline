SPM based data import, preprocessing and analyses for data in BIDS format.  
Includes a function to download the data from the ISN DICOM data base to the project directory  

> **VERSION 3 (OVERHAUL) - WORK IN PROGRESS, NOT YET FOR GENERAL USE**  
> This branch contains a major rewrite of the import, the file naming and most preprocessing modules.  
> `sbp_analyses.m` has only been patched for the new file names, its full refactoring (and the one of the QC and brainstem modules) is still to be done.  
> Until this is finished and merged into `main`, please use the `main` branch (tag `Version-2.1`).  
> Scripts, `get_study_specs.m` and data written by Version 2 are **not compatible** with this version, see [Changes in Version 3](#changes-in-version-3).

## Organisation  

**The data are organized like this:**  

c:\lcpa - the top data directory  
inside you will get a directory \derivatives  
The data import routine (sbp_import_data.m) will put all images in /derivatives/spm_preprocessing but will also make a  
gzipped copy in the main study directory  
  
You can clone this repository to any directory, but should keep it as is  

**It is best not to edit files inside spm_bids_pipeline**  
  
It is advisable to have a "scripts" directory somewhere where your wrapper functions are located  
In this wrapper function make sure to start with an   
**addpath(genpath('your_path_2_spm_bids_pipeline'))**

inside derivatives you have spm_preprocessing where all the data goes using \sub-01 etc  
inside sub-01 you have \ses-01 etc. and in there \anat and \func (and \fmap if specified)  

Before you start you need to create a get_study_specs.m in your "scripts" folder. There is an example inside the repository "get_study_specs.txt"  
Copy it to "scripts" edit it and rename it to "get_study_specs.m"  

Then you should create a wrapper file such as "do_em_all.m" that will sequentially call the subfunctions in the repository  
There is a commented example in the repo (with ext .txt), "do_em_all.txt"

The toolbox uses a nonlinear coregistration of the EPI to the T1w (as in CAT12) or, alternatively, a correction with fieldmaps using VDM files.  
Both workflows are described in **do_em_all.txt**

This toolbox also includes quality control routines found in \qc by Alexandra Tinnermann  
See full_diagn_qc.txt for how to call these functions. 

## File names (light weight BIDS)

All names are produced by the import (`bids_name`), later steps only prepend fixed prefixes:

| data | file name (in sub-XX/ses-YY/...) |
|---|---|
| functional magnitude | `func/sub-01_ses-01_task-<task>_run-01_part-mag_bold.nii` |
| functional phase | `func/sub-01_ses-01_task-<task>_run-01_part-phase_bold.nii` |
| T1w | `anat/sub-01_ses-01_run-01_T1w.nii` (always with `run-01`) |
| fieldmap | `fmap/sub-01_ses-01_run-01_magnitude1.nii`, `_magnitude2.nii`, `_phasediff.nii` (+ json) |

- Magnitude and phase of the same acquisition get the same run number (`part-mag`, `part-phase`). Data without phase images have no `part-` tag.
- `cbv` / `bold` series of VASO acquisitions are of type `cbv` / `bold`. Short reference runs get their own `task` (e.g. `vasoref`).
- Brain and spinal EPIs can be distinguished by an `acq-` tag (`vars.acq_type`).
- Runs that are shorter than expected can be left out at the import (`import.data(n).cond`, or in participants.tsv, see below).

## Preprocessing steps and image prefixes

The steps that are run on the EPIs are listed **in the order they are run** in `get_study_specs.m`:

```matlab
vars.preproc.steps = {'nordic','slicetiming','realign','smooth'};   % any of nordic, slicetiming, realign, smooth
```

`nordic` and `slicetiming` come before `realign`, `smooth` after it (`realign` can be left out, it is then put after the other steps). Each step has a fixed output prefix and SPM prepends it:

| step | function | prefix |
|---|---|---|
| NORDIC | `sbp_nordic` | n |
| slice timing | `sbp_slice_timing` | a |
| realignment (reslicing) | `sbp_realign_2_2` | r |
| smoothing | `sbp_smooth` | s |

So `{'nordic','slicetiming','realign'}` gives `n<epi>`, `an<epi>`, resliced `ran<epi>` (these are the images of the analysis), and `{'nordic','realign','smooth'}` gives `n<epi>`, `rn<epi>`, `srn<epi>`.  
Everything derived from the realignment is named after the prefix `P` of the images that were realigned: `meanP<epi>`, `rp_P<epi>.txt`, `cmeanP<epi>`, `y_meanP<epi>`, `wmeanP<epi>`, `noise_wm_csf_rP<epi>.mat` ...  
All names are built with the helpers `epi_file`, `epi_frames`, `epi_prefix`, `epi_steps` (in `utils`), so a new study only needs the list above.  
Helpers to find files: `get_func_files`, `get_ref_epi`, `get_run_numbers`, `get_anat`.

## Workflows

Typical order (see `do_em_all.txt`):

| | step | what it does |
|---|---|---|
| 1 | `sbp_import_data` | download from the DICOM server, convert, name, make 4D EPIs |
| 2 | `sbp_average_T1` | optional: average several T1w of a session, or of all sessions (`vars.average_T1.across_sessions`) |
| 3 | `sbp_align_images` | coarse rigid alignment of all images of a session to the template (separately for each session) |
| 4 | `sbp_nordic` | optional: NORDIC denoising of magnitude + phase, only the denoised magnitude is written |
| 5 | `sbp_slice_timing`, `sbp_realign_1_2` | slice timing, first realignment (estimate, mean image) |
| 6 | `sbp_coreg_nlin` | nonlinear coregistration EPI --> T1 (also segments the T1) |
| 7 | `sbp_skullstrip`, `sbp_create_dartel`, `sbp_create_trans`, `sbp_create_mask` | skull strip, Dartel T1 --> template, deformations EPI --> T1 --> template, brain mask |
| 8 | `sbp_realign_2_2` | second realignment with the brain mask, reslices all images |
| 9 | `sbp_warp`, `sbp_create_noise` | mean EPI and skull strip in template space, WM / CSF / ventricle noise regressors |
| 10 | `sbp_smooth` | optional: smooth the resliced images |
| 11 | `sbp_create_means`, `sbp_analyses` | group means, 1st / 2nd level analyses |

Fieldmap (VDM) workflow: `sbp_create_vdm_fmap`, `sbp_coreg_lin`, `sbp_seg_T1`, `sbp_create_trans_vdm`, `sbp_warp_vdm`, `sbp_realign_2_2_vdm`, `sbp_apply_vdm` instead of steps 6-9.  
`sbp_shift_images` moves images by a given 4x4 matrix if the automatic alignment fails.

### Import (`sbp_import_data(all_sub_ids, dry_run)`)

- Several subjects are downloaded and converted in parallel (SPM batch system, `import.max_network_procs` at a time). With `dry_run = 1` the plan (which series, which names and run numbers) is printed and only the first files are fetched. A protocol table and a check of the task labels are written.
- `import.data(n).seq` and `.cond` select the series, `.task` sets the task name of reference runs, `vars.acq_type` says whether brain and/or spinal data were acquired.
- participants.tsv with a column `scan_id` (see below) decides which subjects and sessions are imported.

### NORDIC (`sbp_nordic`)

Magnitude and phase are used as input, only the denoised magnitude `n<epi>.nii` is written (the phase is not needed later). Noise volumes at the end of a series (no excitation) are detected automatically (`detect_noise_volumes`), split off (`noise_<name>.nii/.json`) and used for the noise level. Runs without noise volumes use the noise level of the reference run of the same session if there is one, otherwise NORDIC assumes unit noise. **Check the tSNR maps of your data**, NORDIC can produce regions with lower tSNR than the raw data in SMS acquisitions (under investigation).

### VASO data (slab, 3D EPI, cbv + bold)

- `sbp_vaso_T1` preprocesses the short reference runs (cbv and bold, magnitude + phase) and builds T1-like images in the anat folder (`acq-vasomedian`, `acq-vasomeanstd`, `acq-vasosigned`, `acq-vasosignedpos`). Run it before `sbp_nordic` (it provides the noise level).
- `sbp_coreg_nlin_vaso_T1` replaces `sbp_coreg_nlin`: the clipped signed T1-like image is registered to the T1w and the deformation is transferred to the mean EPI of the main runs.
- `sbp_coreg_nlin` itself uses fixed settings that were tested on whole brain 1.5 mm SMS EPI: T1 segments smoothed with 1 mm, sampling distance and smoothing 2 mm, bias correction on, regularisation x1.

### Analysis

`sbp_analyses(all_sub_ids, idx)` runs analysis number `idx` of the `analysis` struct array. The images that are analysed are those after all steps (e.g. `rn` or `srn`, or `analysis.prefix`). The analysis folder is called `<noise_corr>_<name_ana>_<steps>_<shift>_<cvi>`, where `<steps>` lists the steps of the images in the order they were run (e.g. `nr`, `nrs`, `ar`). **`sbp_analyses` was only patched for the new file names, a refactoring is planned.**  
Noise regressors (`sbp_create_noise`): principal components of the WM and CSF voxels (and of two ventricle ROIs; ROIs outside the field of view of a slab are empty and give a message). Optional `vars.noise.smooth_fwhm` and `vars.noise.n_erode`.

## Changes in Version 3

Breaking changes compared to Version 2:

1) **New importer** (several subjects in parallel, plan / dry run, protocol, mag / phase pairing, spinal / brain, `task` for reference runs, T1w always `run-01`, fieldmap magnitude1 / magnitude2). The old importer is in the git history (tag `Version-2.1`).
2) **File names and prefixes** are set by `vars.preproc.steps` (see above) and no longer hard coded (`a`, `ra`, ...). All modules use `epi_file` / `epi_frames` / `get_anat` etc.
3) **New modules:** `sbp_align_images`, `sbp_average_T1`, `sbp_nordic`, `sbp_smooth`, `sbp_vaso_T1`, `sbp_coreg_nlin_vaso_T1`.
4) **New / renamed fields in get_study_specs:** `vars.templateStripID` (was `vars.templateStrip`), `vars.nT1`, `vars.average_T1.across_sessions`, `vars.acq_type`, `vars.preproc.steps`, `vars.smooth.fwhm`, `vars.ref_task`, `vars.nordic.vasoref_noisevols`, `vars.nordic.func_noisevols`, `vars.nordic.ref_noise`, `vars.noise.smooth_fwhm`, `vars.noise.n_erode`, `import.max_network_procs`, `import.data(n).task`, `vars.input_prefix.<stage>` (override). See `get_study_specs.txt`.
5) **Changed results:** `sbp_coreg_nlin` uses new settings (see above), the analysis folder names contain the tag of the steps, `sbp_create_noise` can skip empty ROIs.
6) Runners (`run_spm_parallel` etc.) ignore empty jobs, `spm_BIDS` indexes fieldmap magnitude images.

Known issues / to do: refactoring of `sbp_analyses`, the QC modules (`qc/`) and the brainstem modules (`bs/`); cleaning up dead code; NORDIC on SMS BOLD data; validation of the coregistration settings on more subjects.

## New features

### 10/26 (Version 3, branch v3-overhaul)

See [Changes in Version 3](#changes-in-version-3).

### 10/25

1) Multiple Analysis Definition: ``analysis`` is now a struct array instead of single struct, allowing multiple analyses to be run sequentially. This includes a major overhaul of the get_study_specs structure. All changes are backwards compatible

2) Overhauled PPI analysis: now includes reuse of existing base models, improved consistency checks across subjects, and better user interaction for managing corrupted analyses

3) Critical bug fix to physio regressor creation pipeline. 

### 06/25

1) The new version of the toolbox now also incorporates correction with fieldmaps using VDM files.

2) The new version uses an extension to participants.tsv using a column labeled 'scan_id' (the old import still works)
These contain the individual PRISMA numbers. If there are more sessions, seperate the PRISMA numbers by spaces NOT tabs !
In addition you can specify series that should NOT be imported by adding them in () with NO space in between e.g.

participant_id	sex	age	scan_id
sub-09	F	20	25879 25878(13)

imports sub-09, with sessions 25879 and 25878. Series 13 in 25878 will be omitted

3) all major functions were now renamed to start with sbp_
4) Can now do PPI analyses (Thanks to Lars Keuter for implementing this)
5) Can now incorporate physio (retroicor type) correction with the TAPAS toolbox (Thanks to Leo Asan for implementing this)
    

## Pipeline for physio regressor creation 


IMPORTANT: 
- you need to define path.physioDir, which is baseDir/physio.

- in PhysioDir, you need to have a folder 'mat' or a folder 'smr', depending on which format your file is.
- Your smr files need to be stored as such:
	path.physioDir/physio/smr/sub-XX_ses-YY_physio.smr
- Your mat files need to be stored as such:
	path.physioDir/physio/mat/sub-XX_ses-YY_physio.mat
	
- if you have session-wise files and wish to exclude certain runs, you need to add a json sidecar, which has the following structure:
{
  "includedRuns": [1, NaN, 2, 3, 4]
  }
 where 1-4 are the runs that are included, and NaN is the run you want to exclude.
 save this json sidecar with the same name and in the same location as your session file, with a .json extension.

= OPERATIONS (functions found in utils) =:
1. convert_and_clean_physio: flexibly reads in mat- or smr- session- or run-files, cleans them and converts them to BIDS-compatible TSV files and 
2. create_tapas: Uses tapas toolbox to create txt file of 18 physiological noise regressors in func dirs of subjects.


= OUTPUTS =:
-BIDS-compatible run-wise tsv files of physiological data with jsons in func dir of each sub
-txt file of 18 physiological noise regressors in func dirs of each sub.
-info_physio mat file containing some info coming from tapas toolbox, in func dir of each sub


= REQUIREMENTS =:

1. TAPAS PhysIO Toolbox
	Download:
    **HAS TO BE DOWNLOADED AND INITIALIZED**

    https://github.com/translationalneuromodeling/tapas/releases
    go to tapas_master and type tapas_init() in command line
    you need to do this only once. It will be set up for good in your spm12 after


## PPI Analysis

Psychophysiological interaction (PPI) analysis examines how functional connectivity between brain regions changes depending on psychological context.

### Requirements and Setup

- Every PPI analysis requires a base first-level model defined in `get_study_specs.m`. When first running the PPI analysis, this base model will be estimated. For subsequent PPIs using the same base model, it won't be estimated again.
  - Multiple PPIs (using different ROIs or conditions) can be run sequentially by defining multiple elements in the analysis struct array. These PPIs will use the same base model
  - **Important**: If some subjects don't have the base model yet, you must first estimate the base analysis for these subjects before running the PPI

- PPI Configuration: PPI settings (ROI definition, condition selection etc.) are defined in `get_study_specs.m`

### Best practise

- Verify that condition indices in your PPI specification reference valid conditions in your base model
- For two-pass ROI optimization, make sure the contrast vector matches your model dimensions
