# This is a sample Python script.


from typing import (Tuple, Optional)

import pyradise.data as ps_data
import pyradise.fileio as ps_io
import pyradise.process as ps_proc
from pyradise.data import Modality
import pydicom
import os.path as osp
import os


def modify_operator_name(path: str):
    """
    modify operator names to make the rtss files compatible with pyradise.
    we use the Eclipse output nametags to differentiate the rtss files.
    :param path:
    :param name:
    an absolute path for the files
    :return:
    """
    files = os.listdir(path)
    for f in files:
        if osp.isdir(osp.join(path, f)):
            modify_operator_name(osp.join(path, f))
        elif f.startswith('RS') and f.endswith('dcm'):
            elements = f.split(".")[:]
            name = elements[-2]
            mod_f = ".".join(elements)
            with pydicom.dcmread(osp.join(path, f)) as ds:
                ds.OperatorsName = name
                ds.save_as(osp.join(path, mod_f))
            #os.remove(osp.join(path, f))


class ModalityExtractor(ps_io.ModalityExtractor):

    def extract_from_dicom(self, path: str) -> Optional[ps_data.Modality]:
        tags = (ps_io.Tag((0x0008, 0x0060)),
                ps_io.Tag((0x0008, 0x103e)),)
        dataset_dict = self._load_dicom_attributes(tags, path)

        # Idenfity the modality rule-based
        modality = dataset_dict.get("Modality", {}).get("value", None)
        series_desc = dataset_dict.get("Series Description", {}).get("value", "")
        if modality == "MR":
            if "t1" in series_desc.lower():
                return ps_data.Modality("T1")
            elif "t2" in series_desc.lower():
                return ps_data.Modality("T2")
            else:
                return None
        elif modality == "CT":
            return ps_data.Modality("CT")
        elif modality == "RTDOSE":
            return ps_data.Modality("RTDOSE")
        else:
            return None

    def extract_from_path(self, path: str) -> Optional[Modality]:
        # We can skip the implementation of this method, because we work
        # exclusively with DICOM files
        return None


def get_pipeline(output_size,
                 output_spacing,
                 reference_modality: str = 'CT'
                 ) -> ps_proc.FilterPipeline:
    # Create an empty filter pipeline
    pipeline = ps_proc.FilterPipeline()

    # Add an orientation filter to the pipeline
    #orientation_params = ps_proc.OrientationFilterParams(output_orientation='RAS')
    #pipeline.add_filter(ps_proc.OrientationFilter(), orientation_params)

    # Add a resampling filter to the pipeline
    resample_params = ps_proc.ResampleFilterParams(output_size,
                                                   output_spacing,
                                                   reference_modality=reference_modality,
                                                   centering_method='reference')
    pipeline.add_filter(ps_proc.ResampleFilter(), resample_params)

    return pipeline


def get_ct_size_spacing(series_info):
    series = [x for x in series_info if x.dicom_modality == "CT"]
    if len(series) > 1:
        raise ValueError("Multiple CT series found, currently not supported")
    ct_im = pydicom.dcmread(series[0].path[0])
    size = tuple(map(int, [ct_im.Columns, ct_im.Rows, len(series[0].path)]))
    spacing = [float(x) for x in ct_im.PixelSpacing]
    spacing.append(float(ct_im.SliceThickness))
    return size, tuple(spacing)


def convert_dicom_to_nifti_with_modality_extractor(input_path: str, output_path: str) -> None:
    loader = ps_io.SubjectLoader()
    writer = ps_io.SubjectWriter()
    expected_modalities = ("CT", "RTDOSE")
    modality_selection = ps_io.ModalityInfoSelector(expected_modalities)
    crawler = ps_io.DatasetDicomCrawler(input_path,
                                        modality_extractor=ModalityExtractor(),
                                        modality_config_file_name="doesnotexist.json")
    for series_info in crawler:  # Pyradise does not parse the RT related files
        # Keep just the selected modalities for loading
        series_info = modality_selection.execute(series_info)
        # Load the subject from the series info
        size, spacing = get_ct_size_spacing(series_info)
        pipeline = get_pipeline(output_size=size, output_spacing=spacing)
        subject = loader.load(series_info)
        print(f'Processing subject {subject.get_name()}...')
        subject = pipeline.execute(subject)
        writer.write_to_subject_folder(output_path, subject, write_transforms=False)


# Press the green button in the gutter to run the script.

input_path = "./examples/dose_image_read/patient_seg_data/"
output_path = "./examples/dose_image_read/nifti_data/"
modify_operator_name(input_path)
convert_dicom_to_nifti_with_modality_extractor(input_path, output_path)
