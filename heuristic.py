import os, glob, re, json, pandas as pd
from frozendict import frozendict
import nibabel.nicom.dicomwrappers as nb_dw
from collections import OrderedDict
import string, random


def create_key(subdir, file_suffix , outtype = ("nii.gz"), annotation_classes = None, prefix = ""):

    '''
    Creates the bids key

    Args:
        subdir
        file_suffix
        outtype
        annotation_classes
    
    Return:
        template
        outtype
        annotation_classes
    '''
    
    if not subdir:
        raise ValueError("subdir must be a valid format string")
    
    template = os.path.join(
        prefix,
        "{bids_subject_session_dir}",
        subdir,
        "{bids_subject_session_prefix}_%s" % file_suffix,
    )

   
    return template, outtype, annotation_classes


def custom_seqinfo(wrapper, series_files):

    '''
    Loads info from dicoms

    Args:
        wrapper:
        series_files:
    
    Return:
        custom_info: Dictionary containing values pulled from dicoms
    '''
    slice_orient = wrapper.dcm_data.get([0x0051,0x100e])
    study = wrapper.dcm_data.get([0x0020, 0x000d])
    receive_coil = wrapper.dcm_data.get((0x0051,0x100f))
    contrast_bol = wrapper.dcm_data.get((0x0018, 0x0010))
    contrast_vol = wrapper.dcm_data.get((0x0018,0x1041))
    options = wrapper.dcm_data.get(("ScanOptions"))
    strength = wrapper.dcm_data.get(("MagneticFieldStrength"))
    
    if strength:
        if strength == 15000.0:
            strength = 1.5
        strength = 1.5 if strength <= 1.5 else 3.0
    device = wrapper.dcm_data.get(("DeviceSerialNumber"))
    pixSpacing = wrapper.dcm_data.get(("PixelSpacing"))
    vsd3 = wrapper.dcm_data.get(("SliceThickness"))
    body_part = wrapper.dcm_data.get("BodyPartExamined")

    custom_info = frozendict({
        'study_uid': str(study.value) if study else None,
        'contrast_bol': str(contrast_bol.value) if contrast_bol else None,
        'contrast_vol': str(contrast_vol.value) if contrast_vol else None,
        "strength":  float(round(strength,1)) if strength else None,
        'options': str(options) if options else None,
        'device': str(device) if device else None,
        'voxel': round(float(pixSpacing[0]),3) if pixSpacing else None,
        'vsd3': round(float(vsd3),1) if vsd3 else None,
        'slice_orient': str(slice_orient.value) if slice_orient else None,
        'echo_number': str(wrapper.dcm_data.get("EchoNumber", None)),
        'receive_coil': str(receive_coil.value) if receive_coil else None,
        'body_part': str(body_part).lower() if body_part else None,
    })

    return custom_info


def getDictInfo(dtype):

    '''
    Loads dictionaries

    Args:
        dtype: Dictionary type
    
    Return:
        info_dict: Dictionary 
    '''
    
    if dtype == "acq":
        info = "/app/scripts/acq_dict.json"
    elif dtype == "mod":
        info = "/app/scripts/mod_dict.json"
    elif dtype == "label":
        info = "/app/scripts/label_dict.json"
    
    with open(info) as json_file:
        info_dict = json.load(json_file)

    return info_dict


def get_seq_bids_info(s):

    '''
    Get BIDS info from a series of dicoms

    Args:
        s: Dictionary with series' dicom information
    
    Return:
        seq: Dictionary of information about a series of dicoms 
                    to be used to create a BIDS descriptive filename with BIDS key/value pairs as each entry
        seq_extra: Dictionary of information about a series of dicoms 
                    to be used to create a BIDS descriptive filename with BIDS key/value pairs as each entry
    '''
    
    sdescrip = s.series_description.lower()
    protocol = s.protocol_name.lower()

    # remove "pat2" and "bandwidth" from descriptions and protocol names since if will cause misclassification
    
    sdescrip = sdescrip.replace("pat2","").replace("bandwidth","")
    protocol = protocol.replace("pat2","").replace("bandwidth","")
    
    # Possible mods: anat, dwi, func, fmap, perf, swi, mra
    mod_dict = getDictInfo("mod")
    
    label_dict = getDictInfo("label")
    
    
    mod = None if not [k for k in mod_dict if(k in sdescrip or k in protocol)] else mod_dict[[k for k in mod_dict if(k in sdescrip or k in protocol)][0]]
        

    seq = {
        "type": mod,  
        "label": None,
    }

    seq_extra = {}

    scan_options = s.custom['options']
    fs = False
    if scan_options and ("FS" in scan_options or "FS" == scan_options):
        fs = True
    print(scan_options, fs)


    # recon pair
    if all(x in s.image_type for x in ['NORM', 'DIS2D']):
        seq_extra["rec"] = "DIS2D"
    elif all(x in s.image_type for x in ['NORM', 'DIS3D']):
        seq_extra["rec"] = "DIS3D"
    elif all(x in s.image_type for x in ['NORM', 'DIS2D',"MFSPLIT"]):
        seq_extra["rec"] = "MFSPLIT"
    elif all(x in s.image_type for x in ['NORM', 'FM3_2', 'FIL']):
        seq_extra["rec"] = "FIL"
    elif all(x in s.image_type for x in ['NORM']):
        seq_extra["rec"] = "NORM"
    elif all(x in s.image_type for x in ['PROPELLER']):
        seq_extra["rec"] = "PROPELLER"
    
    # ce pair
    if s.custom["contrast_bol"] or s.custom["contrast_vol"]:
        seq_extra["ce"] = "gad"
    
    # part pair
    part_dict = {
        "P": "phase",
        "M": "mag",
        "R": "real",
        "I": "imag"
    }
    seq_extra["part"] = None if not [x for x in s.image_type if (x in part_dict and (seq["type"] not in ["anat","dwi"]))] else part_dict[[x for x in s.image_type if (x in part_dict and (seq["type"] not in ["anat","dwi"]))][0]]

    # dir pair
    try:
        pedir = s.custom['pe_dir']
        if "COL" in pedir:
            pedir = "AP"
        else:
            pedir = "LR"
        pedir_pos = bool(
            s.custom['pe_dir_pos']
        )

        seq["dir"] = pedir if pedir_pos else pedir[::-1]
    except:
        pass

    # sample pair
    ## we can use this to label fetal scans?
    seq_extra["sample"] = "fetal" if [x for x in ["fetal","utero","gestation"] if (s.study_description and x in s.study_description.lower())] else None
    
    # acq pair 
    ## Current convention: ProtocolFieldStrengthDeviceDlicesVoxel
    acq_dict = getDictInfo("acq")

    ## Identify the CHOP standard MPRAGE
    if s.TR == 2.05 and s.custom["vsd3"] == 0.9:
        if sdescrip == "T1 MPR SAG 0.9 MM":
            acq_parts = [
                "MPRAGEStandardized",
                None if not s.custom["strength"] else "%sT" % str(s.custom["strength"]),
                None if not s.custom["device"] else "ScannerId%s" % s.custom["device"],
                None if not s.dim3 else "Slices%d" % s.dim3,
                None if not s.custom["voxel"] else "Voxel%smm" % str(s.custom["voxel"]),
                ]
        else: 
            acq_parts = [
                "MPRAGEStandardizedVariant",
                None if not s.custom["strength"] else "%sT" % str(s.custom["strength"]),
                None if not s.custom["device"] else "ScannerId%s" % s.custom["device"],
                None if not s.dim3 else "Slices%d" % s.dim3,
                None if not s.custom["voxel"] else "Voxel%smm" % s.custom["voxel"],
                ]
        
    else:

        acq_parts = [
                    None if not [k for k in acq_dict if(k in sdescrip or k in protocol)] else "%s" % acq_dict[[k for k in acq_dict if(k in sdescrip or k in protocol)][0]],
                    ]
    ## filter those which are None, and join 
    seq["acq"] = "".join(filter(bool, acq_parts)).replace(".","p").replace("3p0T","3T")
    seq["acq"] = seq["acq"]+"FatSat" if fs else seq["acq"]
    
    # mt pair
    if scan_options and "MT" in scan_options:
        seq["mt"] = "on"


    # inv pair
    if "INV1" in sdescrip or "INV1" in protocol:
        seq["inv"] = 1
    elif "INV2" in sdescrip or "INV2" in protocol:
        seq["inv"] = 2
    elif "UNI" in s.image_type:
        # seq['acq'] = 'UNI'
        seq["label"] = "UNIT1"


    # Get scan label/suffix
    seq["label"] = "unknown" if not [k for k in label_dict if(k in sdescrip or k in protocol)] else label_dict[[k for k in label_dict if(k in sdescrip or k in protocol)][0]]
    
    # Get modality
    if not mod:
        # Have all unknown data in "anat" for now, will move them to the "other" directory after CuBIDS
        if seq["label"] == "unknown":
            seq["type"] = "anat"
        else:
            # Check if the label can classify the data
            seq["type"] = "anat" if not [k for k in mod_dict if(k == seq["label"])] else label_dict[[k for k in mod_dict if(k == seq["label"])][0]]

        

    return seq, seq_extra


def id_generator(s):
    '''
    Generate filename for text file to keep track of BIDS filenames 
    to appropriately iterate run number to prevent file overwriting 

    Args:
        s:
    
    Return:
        filename
    '''
    
    # val =  ''.join(random.choice(chars) for _ in range(size))
    pat = s.patient_id.replace("HM9","HM") 
    val = pat + "_" + s.custom["study_uid"].rsplit(".",1)[1] + ".txt"
    print(f"file name {val}")
    return val


def getRun(bids_key, s, fname):
    
    '''
    Get run number for scans

    Args:
        bids_key: BIDS filename to add run to
        s:
        fname: name of text file keeping track of BIDS filenames 

    Return:
        run number
    '''

    path = f'/mnt/isilon/bgdlab_processing/Data/tmp_heudiconv/{fname}'

    # First append to file
    if not os.path.exists(path):
        os.system(f"touch {path}")
    with open(path,"a") as f:
        f.write(bids_key)
        f.write("\n")

    # Next, open the file and count how many instances match the current bids_key
    with open(path, 'r') as file:
        matches = [line for line in file if(bids_key in line)]
        
    # The count will be returned and used as the run number
    return str(len(matches)).zfill(3)


def generate_bids_key(s, seq_type, prefix, bids_info, 
                      rand_id, show_part=False,
                      show_dir=False, 
                      outtype=("nii.gz",), **bids_extra):
    
    '''
    Generate the bids key

    Args:
        s:
        seq_type: Type of sequence (T1w, T2w, etc)
        prefix:
        bids_info: Dictionary containing BIDS key/value pairs for a sequence
        show_part: True/False if the BIDS key 'part' should be used in the filename
        show_dir: True/False if the BIDS key 'dir' should be used in the filename
        outtype: Filetypes to output
        bids_extra:
    
    Return:
        BIDS filename
    '''
    
    # sub-<label>[_ses-<label>][sample-<label>][_task-<label>][_acq-<label>][_ce-<label>][_rec-<label>][_inv-<index>][_run-<index>][_echo-<index>][_flip-<index>][_part-<mag|phase|real|imag>][_mt-<label>]_<suffix>.nii[.gz]
    bids_info.update(bids_extra)
    # First without run
    suffix_parts = [
        None if not bids_info.get("sample") else "sample-%s" % bids_info["sample"],
        # None if not bids_info.get("task") else "task-%s" % bids_info["task"],
        None if not bids_info.get("acq") else "acq-%s" % bids_info["acq"],
        None if not bids_info.get("ce") else "ce-%s" % bids_info["ce"],
        None
        if not (bids_info.get("dir") and show_dir)
        else "dir-%s" % bids_info["dir"],
        None if not bids_info.get("rec") else "rec-%s" % bids_info["rec"],
        None if not bids_info.get("inv") else "inv-%d" % bids_info["inv"],
        None if not bids_info.get("echo") else "echo-%d" % int(bids_info["echo"]),
        None if not bids_info.get("flip") else "flip-%d" % int(bids_info["flip"]),
        None if not bids_info.get("mt") else "mt-%s" % bids_info["mt"],
        None 
        if not (bids_info.get("part") and show_part)
        else "part-%s" % bids_info["part"],
        "%s" % bids_info["label"],
    ]
     # filter those which are None, and join with _
    suffix = "_".join(filter(bool, suffix_parts))
    print("bids key w/o run:", suffix)

    suffix_parts = [
        None if not bids_info.get("sample") else "sample-%s" % bids_info["sample"],
        # None if not bids_info.get("task") else "task-%s" % bids_info["task"],
        None if not bids_info.get("acq") else "acq-%s" % bids_info["acq"],
        None if not bids_info.get("ce") else "ce-%s" % bids_info["ce"],
        None
        if not (bids_info.get("dir") and show_dir)
        else "dir-%s" % bids_info["dir"],
        None if not bids_info.get("rec") else "rec-%s" % bids_info["rec"],
        None if not bids_info.get("inv") else "inv-%d" % bids_info["inv"],
        "run-%s" % getRun(suffix, s, rand_id),
        None if not bids_info.get("echo") else "echo-%d" % int(bids_info["echo"]),
        None if not bids_info.get("flip") else "flip-%d" % int(bids_info["flip"]),
        None if not bids_info.get("mt") else "mt-%s" % bids_info["mt"],
        None 
        if not (bids_info.get("part") and show_part)
        else "part-%s" % bids_info["part"],
        "%s" % bids_info["label"],
    ]
    # filter those which are None, and join with _
    suffix = "_".join(filter(bool, suffix_parts))
    print("bids key with run:", suffix)
    return create_key(seq_type, suffix, prefix=prefix, outtype=outtype)


# Heudiconv active run
def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where

    allowed template fields - follow python string module:

    item: index within category
    subject: participant id
    seqitem: run number during scanning
    subindex: sub index within group
    """

    # Create the template output filename for each modality
    
    # Set up the dictionary of modality types
    info = OrderedDict()
    all_bids_infos = {}
    
    prefix = ""
    
    
    for s in seqinfo:
        """
        The namedtuple `s` contains the following fields:

        * total_files_till_now
        * example_dcm_file
        * series_id
        * dcm_dir_name
        * unspecified2
        * unspecified3
        * dim1
        * dim2
        * dim3
        * dim4
        * TR
        * TE
        * protocol_name
        * is_motion_corrected
        * is_derived
        * patient_id
        * study_description
        * referring_physician_name
        * series_description
        * image_type
        * custom
        """

        # Create string to create runs
        rand_id = id_generator(s)

        protocol = s.protocol_name.lower()
        description = s.series_description.lower()
        print("-------------------------------------------")
        print("Description:",description)
        print()


        # Get BIDS info
        bids_info, bids_extra = get_seq_bids_info(s)
        all_bids_infos[s.series_id] = (bids_info, bids_extra)

        # Classify as anat, dwi, swi
        seq_type = bids_info["type"]
        seq_label = bids_info["label"]

        # Let's start filtering out scans earlier
        isAScan = True if s.is_derived == False else False
        # https://github.com/vistalab/scitran-data/blob/master/scitran/data/medimg/dcm/mr/siemens.py
        # suggests there are a few types of "image_type" that should be excluded for siemens images
        typesToExclude = [
            # "DERIVED",
            "PROJECTION IMAGE",
            "MIP",
            "CSAPARALLEL",
            "MR3DNCOMP",
            "CSA 3D",
            "CSAFUSED",
            "CSA MPR",
            "CSA PARALLEL",
            "CSA MPR THICK",
            "ADC",
            # "MOCO",
            "PROJECTION IMAGE",
            "FA",
            # "DIS3D",
            # "DIS2D",
            # "MFSPLIT",
            "SPEC",
            "PERFUSION"
        ]
        excludeTypes = [i for i in typesToExclude if i in s.image_type]
        namesToExclude = ['local', 'scout', 'jj/fairest', "mrv", "hemo", 
                          "localizer",'corticospinal', 'spine', 
                          ' 2D ', "crusher", "wip", "trace",
                          #"t2_tse3dvfl_tra_st2_p2_WIP",
                          "ax dwi all b-1000 hyperband",
                          "cervical", "lumbar", "cardio", "aorta",
                          "svs", "edc",   # single voxel spectroscopy
                          "neck", "loc", "hip", # hippocampus focused
                          "documents", "csf flow", "csfflow", "venc", # CSF Flow
                          "mip", "carotids", "iac", "ratio (a-b)/(c-d)", "brain/scout",
                          "routine_brain/diffusion_fd", "3plane", "3 plane",
                          "mocoseries", "csi_se", "tof"] 
        excludeProtocols = [i for i in namesToExclude if i in protocol]
        excludeDescriptions = [i for i in namesToExclude if i in description]
        print("excludeTypes", len(excludeTypes),*excludeTypes)
        print("namesToExclude", (len(excludeProtocols)+len(excludeDescriptions)),*excludeProtocols,*excludeDescriptions)
        print("isDerived",s.is_derived)
        print("BodyPartExamined:",s.custom["body_part"])
        print("Dimensions:", s.dim1, s.dim2, s.dim3, s.dim4)
        
        if (
            len(excludeProtocols) > 0
            or len(excludeDescriptions) > 0
            or len(excludeTypes) > 0
            or s.custom["body_part"] == "CHEST"
        ):
            isAScan = False
        
    

        # Known issue with derived Siemens diffusion scans
        if (
            "DERIVED" in s.image_type 
            and "DIFFUSION" in s.image_type
        ):
            isAScan = False

        # Only show dir_<label> for dwi and field map data
        show_dir = seq_type in ["fmap", "dwi"]
        # Only show part_<label> for fmap, swi and other data
        show_part = seq_type in ["fmap", "swi","other"]

        # Only concerned with anatomical and diffusion
        print(f"type {seq_type} label {seq_label}")
        if seq_type in ["anat","dwi","swi","other"] and isAScan:
            
            # All images should be 3D or 4D
            if s.dim1 > 1 and s.dim2 > 1 and s.dim3 > 1:
                if s.dim4 > 1:

                    # DWI or fMRI
                    if seq_type == "dwi"  and s.series_files >= 7:
                        print("Found 4D sequence",seq_type, seq_label)
                        
                        template = generate_bids_key(s, seq_type, prefix, bids_info, rand_id, show_part, show_dir, outtype=("nii.gz",),**bids_extra)
                        print(template)
                        if template not in info:
                            info[template] = []
                        info[template].append(s.series_id)
                        
                        # Convert to dict since outside functionality depends on it being a basic dict
                        info = dict(info)  
                    
                    else:
                        print("Other 4D image found:")
                        print("ImageType:", s.image_type)
                        print(protocol)
                        print(description)
                        print(s.dim4)

                elif s.dim1 > 10 and s.dim2 > 10 and s.dim3 > 10 :  # Only 3 dimensions, larger than 10x10x10
                    print("Found",seq_type, seq_label)
                    template = generate_bids_key(s, seq_type, prefix, bids_info, rand_id, show_part, show_dir, outtype=("nii.gz",),**bids_extra)
                    print(template)    

                    if template not in info:
                        info[template] = []
                    info[template].append(s.series_id)
                    
                    # Convert to dict since outside functionality depends on it being a basic dict
                    info = dict(info)  


    return info
                
    