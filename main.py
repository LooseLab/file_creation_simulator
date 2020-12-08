import pickle
import time
from pathlib import Path
# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import os
from pathlib import Path
from collections import namedtuple, deque
import shutil
import pyfastx
import datetime


def get_files_in_origin_dir(origin_dir_path):
    """
    Iterare a given directory path and return a set of Named tuples containing file path and file creation time
    Parameters
    ----------
    origin_dir_path: pathlib.PosixPath
        The filepath to the directory
    Returns
    -------
    list
        List of all files a pathlib paths

    """
    FileInfo = namedtuple("FileInfo", ["file_name", "file_creation_time"])
    file_list = []
    for root, directories, files in os.walk(origin_dir_path):
        for filey in files:
            file_path = Path(filey)
            if ".fastq" in file_path.suffixes and file_path.suffix != ".fxi":
                full_path = Path(root).joinpath(file_path)
                fq = pyfastx.Fastq(str(full_path))
                times = []
                # todo add on the read length divided by 450 to get rough approximation of when read truly finished
                for read in fq:
                    times.append(parse_fastq_description(read.description))
                times.sort()
                print(times[-1])
                file_list.append(FileInfo(full_path, times[-1]))
    file_list.sort(reverse=True, key=lambda x: x.file_creation_time)
    return file_list


def parse_fastq_description(description):
    """
    Parse the description found in a fastq reads header

    Parameters
    ----------
    description: str
        A string of the fastq reads description header

    Returns
    -------
    description_dict: dict
        A dictionary containing the keys and values found in the fastq read headers
    """
    description_dict = {}
    descriptors = description.split(" ")
    # Delete symbol for header
    for item in descriptors:
        if "=" in item:
            bits = item.split("=")
            description_dict[bits[0]] = bits[1]
    try:
        return description_dict["start_time"]
    except KeyError:
        return "-1"


def move_files_in_time(files, destination_dir):
    """
    Pop first two files from deque, return
    Parameters
    ----------
    files: collections.deque
        deque of named tuples with file path and file creation time values
    destination_dir: pathlib.PosixPath
        The destination directory to move to

    Returns
    -------

    """
    base_time = 0
    while files:
        file_1 = files.popleft()
        print(file_1.file_name)
        (destination_dir / file_1.file_name.parts[-2]).mkdir(exist_ok=True, parents=True)
        print(f"Waiting time is {base_time}")
        time.sleep(base_time)
        shutil.copy(file_1.file_name, destination_dir/file_1.file_name.parts[-2])
        if files:
            base_time = (datetime.datetime.strptime(
                        file_1.file_creation_time, "%Y-%m-%dT%H:%M:%S%z"
                ) - datetime.datetime.strptime(files[0].file_creation_time, "%Y-%m-%dT%H:%M:%S%z")).seconds
        print(f"Moved file {file_1.file_name.name}")


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    if not Path("/media/rory/Bioinformatics/Baseline_NB_Single/sorted_files.pickle").exists():
        files = deque(get_files_in_origin_dir(Path("/media/rory/Bioinformatics/Baseline_NB_Single")))
        with open("/media/rory/Bioinformatics/Baseline_NB_Single/sorted_files.pickle", "wb") as fh:
            pickle.dump(files, fh)
    else:
        with open("/media/rory/Bioinformatics/Baseline_NB_Single/sorted_files.pickle", "rb") as fh:
            files = pickle.load(fh)
    move_files_in_time(files, Path("/media/rory/Bioinformatics/fake_data_dir"))