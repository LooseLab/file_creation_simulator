import gzip
import pickle
import time
import yaml
import os
from pathlib import Path
from collections import namedtuple, deque
import shutil
import configargparse
import datetime
from dateutil import parser as date_parser
import logging.config
import random

with open("logging.yaml", "r") as fh:
    configDict = yaml.load(fh, Loader=yaml.FullLoader)
logging.config.dictConfig(configDict)
logger = logging.getLogger("simulator")


def readfq(fp):  # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in ">@":  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        desc, name, seqs, last = last[1:], last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in "@+>":
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != "+":  # this is a fasta record
            yield desc, name, "".join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = "".join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield desc, name, seq, "".join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield desc, name, seq, None  # yield a fasta record instead
                break


def get_files_in_src_dir(origin_dir_path):
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
    file_list = []
    for root, directories, files in os.walk(origin_dir_path):
        for filey in files:
            file_path = Path(filey)
            if not {".fastq", ".fq", ".fna", ".fa", ".fsa", ".fasta"}.intersection(
                set(file_path.suffixes)
            ):
                continue
            handler = (
                gzip.open
                if {".fastq", ".fq"}.intersection(set(file_path.suffixes))
                else open
            )
            if file_path.suffix != ".fxi":
                full_path = Path(root).joinpath(file_path)
                times = []
                with handler(str(full_path), "rt") as fh:
                    for desc, name, seq, qual in readfq(fh):
                        # todo add on the read length divided by 450 to get rough approximation of when read truly finished
                        read_gen_time = len(seq) / 450
                        times.append(parse_fastq_description(desc, int(read_gen_time), read_name=name, file_name=filey))
                times.sort(key=lambda x: x.isoformat())
                file_list.append(FileInfo(full_path, times[-1]))
    file_list.sort(key=lambda x: x.file_creation_time.isoformat())
    return file_list


def parse_fastq_description(description, read_gen_time, **kwargs):
    """
    Parse the description found in a fastq reads header

    Parameters
    ----------
    description: str
        A string of the fastq reads description header
    read_gen_time: int
        Number of seconds taken for the read to actually generate - sequence length divided by 450
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
        t = date_parser.isoparse(description_dict["start_time"]) + datetime.timedelta(
            seconds=read_gen_time
        )
        return t
    except KeyError:
        logger.warning(
            f"No start time for read {kwargs['read_name']} in file {kwargs['file_name']}, randomly assigning time...."
        )
        return datetime.datetime.now() - datetime.timedelta(
            seconds=random.randint(0, 60)
        )


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
        logger.info(file_1.file_name)
        (destination_dir / file_1.file_name.parts[-2]).mkdir(
            exist_ok=True, parents=True
        )
        logger.info(f"Waiting time is {base_time}")
        time.sleep(base_time)
        shutil.copy(file_1.file_name, destination_dir / file_1.file_name.parts[-2])
        if files:
            base_time = abs((
                file_1.file_creation_time - files[0].file_creation_time
            ).total_seconds())
            base_time = base_time * 3 if base_time < 5 else base_time
        logger.info(f"Moved file {file_1.file_name.name}")


def config_parser(parser):
    """
    Configure the arg parse space for this command
    Parameters
    ----------
    parser: argparse.ArgumentParser
        Arguments from the command line
    Returns
    -------

    """
    parser.add_argument(
        "-sd",
        "--src_dir",
        type=str,
        required=True,
        help="Absolute or relative path to the directory containing the FASTQ or FASTA files"
        " to be moved as if they were being written out.",
    )
    parser.add_argument(
        "-dd",
        "--dest_dir",
        type=str,
        required=True,
        help="Directory to move the files from the src directory to,"
        " using last read start time as a proxy for file creation time.",
    )
    parser.add_argument(
        "-c",
        "--config-file",
        required=False,
        is_config_file=True,
        help="config file path",
    )
    return parser


def validate_args(args, parser):
    """
    Validate the arguments given to the parser.
    Parameters
    ----------
    args: argparse.NameSpace
        The chosen command line arguments
    parser: argparse.ArgumentParser
        Parser class instance
    Returns
    -------

    """
    if not Path(args.dest_dir).exists():
        parser.error(
            "Specified destination directory does not exist. Please double check given path."
        )
    if not Path(args.src_dir).exists():
        parser.error(
            "Specified source directory does not exist. Please double check given path."
        )


def get_pickle_if_run_before(src_dir):
    """
    If we have run this before we have created a pickled list of files sorted by their latest read start time,
     it will be fetched if it exists
    Parameters
    ----------
    src_dir: str
        file path to src directory

    Returns
    -------
    list of namedtuple
        List of FileInfo named Tuples, sorted by files latest read start time
    """
    if (Path(src_dir) / "sorted_files.pickle").exists():
        with open(f"{src_dir}/sorted_files.pickle", "rb") as fh:
            try:
                return pickle.load(fh)
            except EOFError as e:
                logger.error("Failed to load pickle file.")
                return None
    else:
        return None


def write_pickle(src_dir, files):
    """
    Write a pickle of sorted files to the src directory
    Parameters
    ----------
    src_dir: str
        Path to the source directory
    files: list of namedtuple
        The fiels to be written to

    Returns
    -------

    """
    with open(f"{src_dir}/sorted_files.pickle", "wb") as fh:
        pickle.dump(files, fh)


# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    logger.info("Beginning")
    FileInfo = namedtuple("FileInfo", ["file_name", "file_creation_time"])
    parser = configargparse.ArgParser(
        description="Move files from given source directory to destination directory "
        "using the the last read start time added to the file as a file creation time proxy."
    )
    parser = config_parser(parser)
    args = parser.parse_args()
    validate_args(args, parser)
    files = get_pickle_if_run_before(args.src_dir)
    if not files:
        files = deque(get_files_in_src_dir(Path(args.src_dir)))
        write_pickle(args.src_dir, files)
    move_files_in_time(files, Path(args.dest_dir))
