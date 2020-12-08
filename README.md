# File creation simulator

Set the origin directory and destination directory in the script, and the fastq files in the origin directory will be iterated,
gathering all the read times for the file, and setting the last one as the file creation time.
Files are then moved one at a time to the destination directory, at interval between file last read times.

Requires pyfastx in the environment.

Run like

```bash 
python main.py
```