# File creation simulator

Command line client for copying FASTA or FASTQ files from their current **_source directory_** to a given **_destination directory_** using their latest read start time as a proxy for the file creation time. The files are then copied at the correct intervals as per their latest read start times.
Config.yaml currently lists the example directories as source and destination directories. 

### Installation

```bash
git clone https://github.com/LooseLab/file_creation_simulator
python -m venv file_creation_simulator
pip install -r requirements.txt
```

### Example commands

- Show help text and exit
    ```bash 
    python main.py -h
    ```

- Run with config file
    ```bash
    python main.py -c config.yaml  
    ```

- Run without config file
    ```bash
    python main.py -sd <path/to/src/directory> -dd <path/to/destination/directory
    ```

